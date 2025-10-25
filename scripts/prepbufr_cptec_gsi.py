#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepbufr_cptec_gsi.py
=====================

Conversor de **BUFR** → **PREPBUFR** compatível com o **GSI** (NCEP).

Este utilitário lê um ou mais arquivos **BUFR** de entrada, extrai cabeçalho e
variáveis observadas, e escreve um arquivo **PREPBUFR** no formato e convenções
esperadas pelo GSI, incluindo as **pilhas de eventos** (Table D) por variável:

- `TEVN`  = [TOB, TQM, TPC, TRC]   (temperatura)
- `WEVN`  = [UOB, VOB, WQM, WPC, WRC]   (vento)
- `PEVN`  = [POB, PQM, PPC, PRC]   (pressão em nível; ADPUPA)
- `ZEVN`  = [ZOB, ZQM, ZPC, ZRC]   (altura/geopotencial; ADPUPA)
- `QEVN`  = [QOB, QQM, QPC, QRC]   (umidade específica)

O **GSI** usa o **último evento** de cada sequência como o estado “final” para
assimilação. Este script grava **um evento** por variável presente (efeito de
*pass-through*), preservando `QM/PC/RC` se existirem no BUFR de entrada, ou
usando `0` (neutro) quando ausentes.

Além disso, o cabeçalho mínimo por subset é escrito na forma:

`SID (opcional) XOB YOB DHR TYP ELV SAID T29`

> **Observação:** O exemplo de escrita do `SID` (IA5) está comentado, pois
> depende de funções de string expostas na sua build do `nceplibs-bufr`.
> O GSI costuma funcionar sem `SID`, mas recomendável incluí-lo quando possível.

-------------------------------------------------------------------------------
Funcionalidades principais
-------------------------------------------------------------------------------
- Lê 1+ BUFRs de entrada (mensagens e subsets).
- Usa a **DX do PREPBUFR** (`prepobs_prep.bufrtable`) para abrir a saída
  com `openbf(..., "OUT", ...)`, embutindo a tabela DX no arquivo.
- Calcula `DHR` (hora relativa ao ciclo) se ausente, a partir de
  `YEAR/MNTH/DAYS/HOUR/MINU/SECO`, em relação a `--idate`.
- Monta as sequências de evento por variável, com 1 evento final cada.
- Escreve campos não-evento lidos pelo GSI quando presentes (`PRSS`, `PWO`, `CAT`).
- Permite **forçar** `TYP/T29` via CLI (`--typ`, `--t29`) caso estejam ausentes
  ou inconsistentes no BUFR de entrada.
- Suporta “tipo lógico” do dado (`--type adpsfc|adpupa`) para decidir quais
  sequências adicionar (ex.: `PEVN/ZEVN` tipicamente em ADPUPA).

-------------------------------------------------------------------------------
Requisitos
-------------------------------------------------------------------------------
- Python 3.8+
- `nceplibs-bufr` (conda-forge) — módulo Python `ncepbufr`
- `numpy`

Instalação (sugestão):

    conda install -c conda-forge nceplibs-bufr numpy

-------------------------------------------------------------------------------
Uso
-------------------------------------------------------------------------------
Exemplo ADPSFC (copiando `TYP/T29` do input quando existirem):

    python prepbufr_cptec_gsi.py --type adpsfc --dx prepobs_prep.bufrtable \
        --idate 2025091100 --out saida.prepbufr in1.bufr in2.bufr

Exemplo ADPUPA forçando `TYP/T29`:

    python prepbufr_cptec_gsi.py --type adpupa --typ 120 --t29 0 \
        --dx prepobs_prep.bufrtable --idate 2025091100 --out saida.prepbufr in.bufr

-------------------------------------------------------------------------------
Notas e boas práticas
-------------------------------------------------------------------------------
- **DX correta**: use a `prepobs_prep.bufrtable` adequada ao seu pipeline.
- **DHR relativo**: `DHR = (t_obs - idate) [horas]`; se ausente, o script calcula.
- **Missing**: este script usa o missing da BUFRLIB (`1.0e10`). Evite sentinelas
  como `-9999`.
- **Unidades** (conforme DX): graus (XOB/YOB), K (TOB), m/s (U/V), hPa (POB/PRSS),
  m (ZOB), kg/kg (QOB).
- **TYP/T29**: mantenha coerentes; use `--typ/--t29` se precisar forçar.
- **SID**: opcional aqui; recomendo ativar (ver TODO no código) se a sua build
  do `nceplibs-bufr` expuser funções de string.

-------------------------------------------------------------------------------
Limitações
-------------------------------------------------------------------------------
- Não implementa QC físico/estatístico — apenas transporta/define `QM/PC/RC`.
- `SID` (IA5) comentado por depender de funções específicas da BUFRLIB.
- Não escreve erros de observação (`TOE/WOE/...`) por padrão (bloco comentado).

-------------------------------------------------------------------------------
Saída / retorno
-------------------------------------------------------------------------------
- Escreve o arquivo PREPBUFR indicado em `--out` e imprime a contagem de subsets
  gravados. Retorno de processo: `0` (sucesso) ou exceção em caso de erro.

"""

from __future__ import annotations
import argparse
import datetime as dt
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import ncepbufr

# Expondo a API Fortran da BUFRLIB via bindings do pacote nceplibs-bufr
BLIB = ncepbufr.bufrlib  # BUFRLIB (Fortran) exposta no Python


# ---------------- util ---------------- #

def parse_idate(s: str) -> dt.datetime:
    """Converte `YYYYMMDDHH` para `datetime`.

    Parameters
    ----------
    s : str
        String no formato `YYYYMMDDHH` (ex.: "2025091100").

    Returns
    -------
    datetime.datetime
        Data-centro do ciclo.

    Raises
    ------
    ValueError
        Se o formato não for `YYYYMMDDHH`.
    """
    return dt.datetime.strptime(s, "%Y%m%d%H")


def fnum(a) -> Optional[float]:
    """Converte qualquer valor numérico/array para `float` (primeiro elemento).

    Parameters
    ----------
    a : Any
        Escalar ou elemento indexável.

    Returns
    -------
    float | None
        Valor convertido ou `None` em caso de falha.
    """
    try:
        return float(a)
    except Exception:
        return None


def read1(bufr, mnem: str) -> Optional[float]:
    """Lê um mnemônico escalar do subset corrente.

    Usa `bufr.read_subset(mnem)` e retorna o **primeiro elemento** como `float`.

    Parameters
    ----------
    bufr : ncepbufr.Bufr
        Objeto de leitura do nceplibs-bufr posicionado no subset atual.
    mnem : str
        Nome do mnemônico (ex.: "TOB", "XOB").

    Returns
    -------
    float | None
        Valor como `float`, ou `None` se ausente ou não numérico.
    """
    arr = bufr.read_subset(mnem)
    if arr.size == 0:
        return None
    return fnum(arr[0][0])


def compute_dhr(bufr, idate: dt.datetime) -> Optional[float]:
    """Computa `DHR` (hora relativa ao ciclo) se `DHR` não existir no subset.

    A partir de `YEAR/MNTH/DAYS/HOUR/MINU/SECO`, calcula:
    `DHR = (t_obs - idate).total_seconds() / 3600.0`

    Parameters
    ----------
    bufr : ncepbufr.Bufr
        Leitor posicionado no subset corrente.
    idate : datetime.datetime
        Data-centro do ciclo.

    Returns
    -------
    float | None
        `DHR` em horas, ou `None` se campos de data/hora não estiverem completos.
    """
    y = read1(bufr, "YEAR"); m = read1(bufr, "MNTH"); d = read1(bufr, "DAYS")
    h = read1(bufr, "HOUR"); mi = read1(bufr, "MINU") or 0.0; se = read1(bufr, "SECO") or 0.0
    if None in (y, m, d, h):
        return None
    try:
        t = dt.datetime(int(y), int(m), int(d), int(h), int(mi), int(se))
    except Exception:
        return None
    return (t - idate).total_seconds() / 3600.0


def f_default(x: Optional[float], default: float) -> float:
    """Retorna `x` como `float` se finito; caso contrário retorna `default`."""
    return default if (x is None or not np.isfinite(x)) else float(x)


def f_missing(x: Optional[float], missing: float) -> float:
    """Retorna `x` se finito; caso contrário retorna o valor de `missing` (BUFR)."""
    return missing if (x is None or not np.isfinite(x)) else float(x)


# -------------- BUFRLIB wrappers -------------- #

def open_prepbufr(out_path: Path, dx_table: Path):
    """Abre arquivo PREPBUFR para escrita, embutindo a DX (tabela) na saída.

    Fluxo BUFRLIB (modo OUT):
    1) `datelen(10)`  → datas `YYYYMMDDHH`
    2) `setbmiss(1.0e10)`  → missing padrão
    3) `open(dx_table, L_TAB)`  → abre tabela DX
    4) `open(out_file, L_OUT, mode='wb')`  → cria arquivo binário
    5) `openbf(L_OUT, "OUT", L_TAB)`  → conecta e embute a DX na saída

    Parameters
    ----------
    out_path : pathlib.Path
        Caminho do arquivo PREPBUFR de saída.
    dx_table : pathlib.Path
        Caminho para a `prepobs_prep.bufrtable`.

    Returns
    -------
    tuple[int, int]
        (L_OUT, L_TAB) — unidades lógicas usadas pela BUFRLIB.
    """
    BLIB.datelen(10)           # YYYYMMDDHH
    BLIB.setbmiss(1.0e10)      # missing BUFRLIB
    L_OUT, L_TAB = 11, 12
    BLIB.open(str(dx_table).encode(), L_TAB)
    BLIB.open(str(out_path).encode(), L_OUT, mode='wb')
    BLIB.openbf(L_OUT, b"OUT", L_TAB)  # embute a DX no arquivo (estilo PREPBUFR)
    return L_OUT, L_TAB


def close_prepbufr(L_OUT: int, L_TAB: int):
    """Fecha mensagem corrente (se aberta), arquivo BUFR e unidade da tabela."""
    try:
        BLIB.closmg(L_OUT)
    except Exception:
        pass
    BLIB.closbf(L_OUT)
    try:
        BLIB.close(L_TAB)
    except Exception:
        pass


def ensure_open_message(curr: Optional[bytes], msgtype: bytes, L_OUT: int, idate_int: int) -> bytes:
    """Garante que a mensagem de saída do tipo `msgtype` esteja aberta.

    Fecha a mensagem anterior (se houver) e abre uma mensagem de tipo `msgtype`
    com `openmb(L_OUT, msgtype, idate_int)` quando necessário.

    Parameters
    ----------
    curr : bytes | None
        Tipo de mensagem atualmente aberta (ou `None` se nenhuma).
    msgtype : bytes
        Tipo de mensagem BUFR (ex.: b'ADPSFC', b'ADPUPA', etc.).
    L_OUT : int
        Unidade lógica do arquivo de saída.
    idate_int : int
        Data-centro (YYYYMMDDHH) como inteiro.

    Returns
    -------
    bytes
        O tipo de mensagem agora ativo.
    """
    if curr != msgtype:
        if curr is not None:
            BLIB.closmg(L_OUT)
        BLIB.openmb(L_OUT, msgtype, idate_int)
        return msgtype
    return curr


def ufbint_mat(L_OUT: int, names: bytes, mat: np.ndarray):
    """Escreve vetor/matriz (`ncampos x nreps`) usando `UFBINT`.

    Parameters
    ----------
    L_OUT : int
        Unidade lógica do arquivo de saída.
    names : bytes
        Lista de mnemonics separados por espaço (ex.: b"XOB YOB DHR").
    mat : numpy.ndarray
        Matriz 2D com shape `(ncampos, nreps)`.
    """
    BLIB.ufbint(L_OUT, mat, mat.shape[0], mat.shape[1], names)


def ufbseq_mat(L_OUT: int, seqname: bytes, mat: np.ndarray):
    """Escreve sequência Table D (`seqname`) com matriz `(ncampos x neventos)` via `UFBSEQ`.

    Parameters
    ----------
    L_OUT : int
        Unidade lógica do arquivo de saída.
    seqname : bytes
        Nome da sequência Table D (ex.: b"TEVN", b"WEVN").
    mat : numpy.ndarray
        Matriz 2D `(ncampos, neventos)` com os campos da sequência.
    """
    BLIB.ufbseq(L_OUT, mat, mat.shape[0], mat.shape[1], seqname)


# -------------- escrita de um subset -------------- #

# Cabeçalho mínimo (numérico). `SID` (IA5) é opcional (ver TODO logo abaixo).
HDR_NUM = b"XOB YOB DHR TYP ELV SAID T29"


def write_subset(
    bufr,
    L_OUT: int,
    idate_dt: dt.datetime,
    force_typ: Optional[float],
    force_t29: Optional[float],
    kind: str,
) -> bool:
    """Escreve **um subset** PREPBUFR a partir do subset BUFR atual (se possível).

    Fluxo:
      1) Monta e escreve o **cabeçalho** com `UFBINT`: `XOB YOB DHR TYP ELV SAID T29`
         (opcionalmente `SID`, comentado).
      2) Para cada variável presente, escreve **uma** sequência de **evento final**:
         `TEVN`, `WEVN`, `QEVN`, e (se `kind == "adpupa"`) `PEVN`, `ZEVN`.
         - `QM/PC/RC` são copiados do BUFR de entrada se existirem; caso contrário `0`.
      3) Escreve campos **não-evento** lidos pelo GSI quando presentes: `PRSS`, `PWO`, `CAT`.
      4) (Opcional) Exemplos comentados para erros de observação `TOE/WOE/...`.

    Parameters
    ----------
    bufr : ncepbufr.Bufr
        Objeto leitor posicionado no subset de entrada.
    L_OUT : int
        Unidade lógica do arquivo PREPBUFR de saída.
    idate_dt : datetime.datetime
        Data-centro do ciclo; usada para calcular `DHR` quando ausente.
    force_typ : float | None
        Valor para forçar `TYP` (quando `None`, usa `TYP` do input ou `0.0`).
    force_t29 : float | None
        Valor para forçar `T29` (quando `None`, usa `T29` do input ou `0.0`).
    kind : str
        Tipo lógico do dado de entrada (`"adpsfc"` ou `"adpupa"`).

    Returns
    -------
    bool
        `True` se o subset foi escrito; `False` se faltaram campos essenciais
        de cabeçalho (ex.: XOB/YOB/DHR).
    """
    # Header básico
    xob = read1(bufr, "XOB")
    yob = read1(bufr, "YOB")
    dhr = read1(bufr, "DHR") or compute_dhr(bufr, idate_dt)
    if xob is None or yob is None or dhr is None:
        # Sem lon/lat/tempo, não gravamos este subset
        return False

    typ_in = read1(bufr, "TYP")
    t29_in = read1(bufr, "T29")
    typ = f_default(force_typ if force_typ is not None else typ_in, 0.0)
    t29 = f_default(force_t29 if force_t29 is not None else t29_in, 0.0)
    elv = f_default(read1(bufr, "ELV"), 0.0)
    said = f_default(read1(bufr, "SAID"), 0.0)

    hdr = np.array([[xob, yob, dhr, typ, elv, said, t29]], dtype=np.float64).T
    ufbint_mat(L_OUT, HDR_NUM, hdr)

    # TODO (opcional): escrever SID (IA5).
    # Dependendo da sua build, pode existir BLIB.ufbstr(...) ou rotina similar.
    # Exemplo (quando disponível):
    # sid_arr = bufr.read_subset("SID")
    # if sid_arr.size:
    #     sid_str = sid_arr[0][0].tobytes().decode(errors="ignore").strip()
    #     BLIB.ufbstr(L_OUT, sid_str.encode(), 1, 1, b"SID")

    # Eventos finais (um evento por variável disponível)
    MISS = 1.0e10

    def ev4(obs_m, qm_m, pc_m, rc_m, seq_name: bytes):
        """Auxiliar para escrever uma sequência de 4 campos (ex.: TEVN/PEVN/ZEVN/QEVN)."""
        v = read1(bufr, obs_m)
        if v is None:
            return
        qm = f_default(read1(bufr, qm_m), 0.0)
        pc = f_default(read1(bufr, pc_m), 0.0)
        rc = f_default(read1(bufr, rc_m), 0.0)
        mat = np.array([[v, qm, pc, rc]], dtype=np.float64).T
        ufbseq_mat(L_OUT, seq_name, mat)

    # Temperatura
    ev4("TOB", "TQM", "TPC", "TRC", b"TEVN")

    # Vento (U/V podem faltar isoladamente)
    u, v = read1(bufr, "UOB"), read1(bufr, "VOB")
    if (u is not None) or (v is not None):
        wqm = f_default(read1(bufr, "WQM"), 0.0)
        wpc = f_default(read1(bufr, "WPC"), 0.0)
        wrc = f_default(read1(bufr, "WRC"), 0.0)
        mat = np.array([[f_missing(u, MISS), f_missing(v, MISS), wqm, wpc, wrc]], dtype=np.float64).T
        ufbseq_mat(L_OUT, b"WEVN", mat)

    # Pressão “em nível” (POB) e altura/geopotencial (ZOB) — típicos em ADPUPA
    if kind == "adpupa":
        ev4("POB", "PQM", "PPC", "PRC", b"PEVN")
        ev4("ZOB", "ZQM", "ZPC", "ZRC", b"ZEVN")

    # Umidade específica
    ev4("QOB", "QQM", "QPC", "QRC", b"QEVN")

    # Campos não-evento que o GSI também lê (quando existirem): PRSS/PWO/CAT
    for nm in (b"PRSS", b"PWO", b"CAT"):
        val = read1(bufr, nm.decode())
        if val is not None:
            arr = np.array([[val]], dtype=np.float64).T
            ufbint_mat(L_OUT, nm, arr)

    # (Opcional) erros de observação se estiverem no BUFR de entrada
    # for nm in (b"TOE", b"WOE", b"POE", b"ZOE", b"QOE"):
    #     v = read1(bufr, nm.decode())
    #     if v is not None:
    #         ufbint_mat(L_OUT, nm, np.array([[v]], dtype=np.float64).T)

    BLIB.writsb(L_OUT)
    return True


# -------------- driver -------------- #

def process(
    inputs: Sequence[Path],
    out_path: Path,
    dx_table: Path,
    idate: dt.datetime,
    kind: str,
    force_typ: Optional[float],
    force_t29: Optional[float],
):
    """Executa a conversão BUFR → PREPBUFR.

    Para cada arquivo BUFR de entrada:
      - percorre mensagens com `advance()`
      - ajusta o `msg_type` de saída com `ensure_open_message(...)`
      - itera subsets com `load_subset()`
      - escreve cada subset com `write_subset(...)`

    Parameters
    ----------
    inputs : Sequence[pathlib.Path]
        Lista de arquivos BUFR de entrada.
    out_path : pathlib.Path
        Caminho do arquivo PREPBUFR de saída.
    dx_table : pathlib.Path
        Caminho para a tabela DX (`prepobs_prep.bufrtable`).
    idate : datetime.datetime
        Data-centro do ciclo (para `DHR` e `OPENMB`).
    kind : {"adpsfc", "adpupa"}
        Tipo lógico do dado de entrada (controla sequências e extras).
    force_typ : float | None
        Força `TYP` no cabeçalho (se `None`, tenta usar o input).
    force_t29 : float | None
        Força `T29` no cabeçalho (se `None`, tenta usar o input).
    """
    L_OUT, L_TAB = open_prepbufr(out_path, dx_table)
    curr: Optional[bytes] = None
    idate_int = int(idate.strftime("%Y%m%d%H"))
    nout = 0
    try:
        for f in inputs:
            bufr = ncepbufr.open(str(f))
            while bufr.advance() == 0:
                msgtype: bytes = bufr.msg_type  # ex.: b'ADPSFC', b'ADPUPA', ...
                curr = ensure_open_message(msgtype, L_OUT, idate_int, curr)
                while bufr.load_subset() == 0:
                    if write_subset(bufr, L_OUT, idate, force_typ, force_t29, kind):
                        nout += 1
            bufr.close()
    finally:
        close_prepbufr(L_OUT, L_TAB)
    print(f"[OK] Escrevi {nout} subset(s) em {out_path}")


def main():
    """Interface de linha de comando.

    Argumentos:
      --type  : "adpsfc" | "adpupa" (requerido)
      --dx    : caminho para `prepobs_prep.bufrtable` (requerido)
      --idate : `YYYYMMDDHH` (requerido)
      --out   : arquivo PREPBUFR de saída (requerido)
      --typ   : força `TYP` (float)
      --t29   : força `T29` (float)
      inputs  : 1+ BUFRs de entrada
    """
    ap = argparse.ArgumentParser(description="BUFR → PREPBUFR (compatível GSI)")
    ap.add_argument("--type", choices=["adpsfc", "adpupa"], required=True,
                    help="Tipo lógico do dado de entrada (controla PEVN/ZEVN e extras).")
    ap.add_argument("--dx", required=True,
                    help="Caminho para a DX do PREPBUFR (prepobs_prep.bufrtable).")
    ap.add_argument("--idate", required=True,
                    help="Data-centro YYYYMMDDHH (para DHR e OPENMB).")
    ap.add_argument("--out", required=True,
                    help="Arquivo PREPBUFR de saída.")
    ap.add_argument("--typ", type=float, default=None,
                    help="Forçar TYP (se ausente/incorreto no input).")
    ap.add_argument("--t29", type=float, default=None,
                    help="Forçar T29 (se ausente/incorreto no input).")
    ap.add_argument("inputs", nargs="+", help="BUFR(s) de entrada.")
    args = ap.parse_args()

    process(
        inputs=[Path(x) for x in args.inputs],
        out_path=Path(args.out),
        dx_table=Path(args.dx),
        idate=parse_idate(args.idate),
        kind=args.type,
        force_typ=args.typ,
        force_t29=args.t29,
    )


if __name__ == "__main__":
    main()


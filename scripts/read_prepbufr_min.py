#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read_prepbufr_min.py
====================

Leitor **mínimo** (mas útil) de arquivos **PREPBUFR** para inspeção rápida.

Este utilitário percorre mensagens e subsets (relatórios) de um PREPBUFR e
imprime, por linha, informações essenciais do subset, incluindo cabeçalho e o
**último evento** das sequências de temperatura (`TEVN`) e vento (`WEVN`),
que é o que o **GSI** considera no momento da assimilação.

O script foi pensado para:
- validar rapidamente que campos-chave existem e possuem valores sensatos;
- inspecionar o *payload* de eventos finais (temperatura e vento);
- servir de exemplo mínimo de leitura com o pacote `ncepbufr`.

-------------------------------------------------------------------------------
Requisitos
-------------------------------------------------------------------------------
- Python 3.8+
- `nceplibs-bufr` (conda-forge) → módulo `ncepbufr`
- `numpy`

Instalação rápida:

    conda install -c conda-forge nceplibs-bufr numpy

-------------------------------------------------------------------------------
Conceitos rápidos
-------------------------------------------------------------------------------
- **DX embutida:** PREPBUFR costuma carregar a tabela DX embutida; portanto,
  `ncepbufr.open()` é suficiente na maioria dos casos.
- **Mensagem / subset:** uma *mensagem* contém N *subsets* (relatórios). Use
  `advance()` para ir à próxima mensagem e `load_subset()` para iterar subsets.
- **Header típico:** `SID XOB YOB DHR TYP ELV SAID T29`
- **Pilhas de eventos (Table D):**
    - `TEVN` → linhas = `[TOB, TQM, TPC, TRC]`, colunas = eventos (última coluna é o *final*)
    - `WEVN` → linhas = `[UOB, VOB, WQM, WPC, WRC]` (às vezes 5 linhas)
- **Último evento:** o GSI usa o **último** evento da sequência (coluna final).

-------------------------------------------------------------------------------
Uso
-------------------------------------------------------------------------------
Exemplo simples:

    python read_prepbufr_min.py input.prepbufr

Formatação de saída (por subset):

    [000123] SID=<texto> lon=<graus> lat=<graus> dhr=<horas> \
    T_last=<K> TQM=<flag> U_last=<m/s> V_last=<m/s> WQM=<flag>

Observações:
- Se um mnemônico não existir no subset, `read_subset(mnem)` retorna array vazio.
- Os `*_last` são os valores do **último** evento de `TEVN/WEVN` quando estas
  sequências existem; caso contrário, imprime `nan`.

-------------------------------------------------------------------------------
Limitações
-------------------------------------------------------------------------------
- Este leitor **não** executa QC físico/estatístico nem valida unidades.
- Apenas `TEVN`/`WEVN` são lidos como exemplo (adapte para `PEVN/ZEVN/QEVN`).
- `SID` é impresso se presente e convertido de IA5.

"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Tuple

import ncepbufr
import numpy as np


# ------------------------------- helpers ------------------------------------ #

def read_scalar(bufr, mnem: str) -> Optional[float]:
    """Lê um mnemônico escalar do subset atual e devolve `float` (ou `None`).

    Parameters
    ----------
    bufr : ncepbufr.Bufr
        Leitor posicionado no subset atual.
    mnem : str
        Nome do mnemônico (ex.: "XOB", "DHR").

    Returns
    -------
    float | None
        Valor convertido, ou `None` se ausente/não numérico.
    """
    arr = bufr.read_subset(mnem)
    if arr.size == 0:
        return None
    try:
        return float(arr[0][0])
    except Exception:
        return None


def read_sid(bufr) -> str:
    """Retorna o `SID` (IA5) do subset como `str`, ou string vazia se ausente.

    O `SID` pode vir como bytes IA5 em um array de 2D. Este helper trata a
    conversão de maneira resiliente.
    """
    sid = bufr.read_subset("SID")
    if not getattr(sid, "size", 0):
        return ""
    try:
        return sid[0][0].tobytes().decode(errors="ignore").strip()
    except Exception:
        # fallback: melhor esforço
        try:
            return str(sid[0][0]).strip()
        except Exception:
            return ""


def last_event_from_seq(seq: np.ndarray, *, expect_rows: Optional[int] = None) -> np.ndarray:
    """Retorna a **última coluna** de uma sequência Table D como vetor 1D.

    Parameters
    ----------
    seq : numpy.ndarray
        Matriz retornada por `read_subset("TEVN")`, `read_subset("WEVN")`, etc.
        Espera-se shape `(n_linhas, n_eventos)`.
    expect_rows : int | None
        Se informado, valida o mínimo de linhas. Ex.: `WEVN` pode ter 5 linhas.

    Returns
    -------
    numpy.ndarray
        Vetor 1D (shape `(n_linhas,)`) correspondente à **última coluna**.
        Se `seq` não for 2D ou não tiver colunas, levanta `ValueError`.
    """
    a = np.asarray(seq)
    if a.ndim != 2 or a.shape[1] == 0:
        raise ValueError("Sequência vazia ou com shape inesperado.")
    if expect_rows is not None and a.shape[0] < expect_rows:
        # Não é erro fatal; apenas segue — alguns DXs reduzem campos
        pass
    return a[:, -1]


# --------------------------------- main -------------------------------------- #

def main(argv=None) -> int:
    """Ponto de entrada da CLI.

    Parameters
    ----------
    argv : list[str] | None
        Argumentos de linha de comando (default: `sys.argv[1:]`).

    Returns
    -------
    int
        0 em caso de sucesso; diferente de 0 em falhas de I/O/formato.
    """
    ap = argparse.ArgumentParser(
        prog="read_prepbufr_min.py",
        description="Leitor mínimo de PREPBUFR: imprime cabeçalho e últimos eventos de TEVN/WEVN.",
    )
    ap.add_argument("file", help="Caminho do arquivo PREPBUFR (DX embutida).")
    args = ap.parse_args(argv)

    fn = args.file

    try:
        bufr = ncepbufr.open(fn)
    except Exception as e:
        print(f"[ERRO] Não foi possível abrir '{fn}': {e}", file=sys.stderr)
        return 2

    n = 0
    try:
        # Loop de mensagens → subsets
        while bufr.advance() == 0:
            while bufr.load_subset() == 0:
                n += 1

                # ---- Cabeçalho típico (GSI)
                # Nota: se algum mnemônico não existir, read_scalar retorna None
                sid = read_sid(bufr)
                xob = read_scalar(bufr, "XOB")
                yob = read_scalar(bufr, "YOB")
                dhr = read_scalar(bufr, "DHR")

                # ---- Sequência TEVN (temperatura)
                last_TOB = np.nan
                last_TQM = np.nan
                tev = bufr.read_subset("TEVN")
                if getattr(tev, "size", 0):
                    try:
                        v = last_event_from_seq(np.asarray(tev), expect_rows=4)
                        # TEVN linhas usuais: [TOB, TQM, TPC, TRC]
                        last_TOB = float(v[0]) if np.isfinite(v[0]) else np.nan
                        last_TQM = float(v[1]) if np.isfinite(v[1]) else np.nan
                    except Exception:
                        pass  # mantém nan

                # ---- Sequência WEVN (vento)
                last_U = np.nan
                last_V = np.nan
                last_WQM = np.nan
                wev = bufr.read_subset("WEVN")
                if getattr(wev, "size", 0):
                    try:
                        v = last_event_from_seq(np.asarray(wev))  # às vezes 5 linhas
                        # WEVN comum: [UOB, VOB, WQM, WPC, WRC]
                        # Tratamos defensivamente tamanho variável
                        if v.size >= 1 and np.isfinite(v[0]):
                            last_U = float(v[0])
                        if v.size >= 2 and np.isfinite(v[1]):
                            last_V = float(v[1])
                        if v.size >= 3 and np.isfinite(v[2]):
                            last_WQM = float(v[2])
                    except Exception:
                        pass  # mantém nan

                # ---- Impressão formatada (com tolerância a None)
                def ffmt(x, width=8, prec=3):
                    return f"{(x if x is not None else float('nan')):{width}.{prec}f}"

                print(
                    f"[{n:06d}] SID={sid:>10s} "
                    f"lon={ffmt(xob)} lat={ffmt(yob)} dhr={ffmt(dhr, width=6, prec=2)} "
                    f"T_last={(last_TOB if np.isfinite(last_TOB) else float('nan')):7.2f} "
                    f"TQM={(last_TQM if np.isfinite(last_TQM) else float('nan'))} "
                    f"U_last={(last_U if np.isfinite(last_U) else float('nan')):7.2f} "
                    f"V_last={(last_V if np.isfinite(last_V) else float('nan')):7.2f} "
                    f"WQM={(last_WQM if np.isfinite(last_WQM) else float('nan'))}"
                )
    finally:
        bufr.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


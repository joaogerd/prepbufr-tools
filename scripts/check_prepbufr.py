#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_prepbufr.py
=================

Validador “smoke test” para arquivos **PREPBUFR** visando compatibilidade com o GSI.

Este utilitário realiza checagens **estruturais** e de **plausibilidade** sobre cada
subset (relatório) em um arquivo PREPBUFR, ajudando a diagnosticar problemas comuns
em pipelines que geram PREPBUFR (ex.: conversões BUFR→PREPBUFR). Ele **não** executa
controle de qualidade físico/estatístico — seu objetivo é verificar se o arquivo
está **bem formado** e carrega os **metadados** e **sequências de evento** esperados
pelo GSI.

Principais verificações
-----------------------
1) **Cabeçalho mínimo** por subset:
   - Mnemonics requeridos: ``XOB YOB DHR TYP ELV SAID T29``.
   - Ausência de qualquer um conta como *header incompleto*.

2) **Janela de tempo**:
   - Verifica se ``|DHR| ≤ twind_h`` (hora relativa ao ciclo).
   - Se ``DHR`` estiver ausente ou fora da janela, é contabilizado.

3) **Sequências de eventos** (Table D) com *evento final*:
   - Procura ``TEVN`` (checando ``TQM``) ou ``WEVN`` (checando ``WQM``).
   - Se não encontrar, tenta ``PEVN/ZEVN/QEVN``.
   - Ausência total conta como *sem eventos finais*.

4) **Validação de unidades/valores plausíveis** (habilitado por padrão):
   - Checa ranges plausíveis para: ``XOB/YOB`` (graus),
     ``TOB`` (K), ``POB`` (hPa), ``ZOB`` (m), ``UOB/VOB`` (m s⁻¹), ``QOB`` (kg kg⁻¹).
   - Se o campo existir e estiver fora do intervalo esperado, o subset falha
     na checagem de unidades.
   - **Pressão**: aceite ``--pressure-unit {hpa,cb,pa,auto}`` (default: ``hpa``);
     em ``auto`` é feita uma **autodetecção** pela mediana de ``POB/PRSS``.

5) **Erros de observação** (opcional via ``--check-oberrs``):
   - Quando a variável existe, exige o respectivo erro:
     ``TOB→TOE``, ``U/V→WOE``, ``POB→POE``, ``ZOB→ZOE``, ``QOB→QOE``.

6) **Checagens específicas de ADPSFC** (opcional via ``--kind adpsfc``):
   - Requer presença de ``PRSS``, ``PWO`` e ``CAT`` (quando o subset é
     “significativo”: header ok + eventos).
   - Valida ranges plausíveis de ``PRSS`` e ``PWO``.

7) **Relatórios**:
   - ``--csv``: CSV **por msg_type** (total, aprovados, % aprovado, contagens de falha).
   - ``--where N``: imprime até **N exemplos** de “onde” ocorreram falhas (lon/lat/DHR).
   - ``--report-csv``: CSV **detalhado**, 1 linha por falha (tipo, msg_type, subset, var).
   - ``--vars-csv``: CSV **agregado** com contagem de falhas de **unidade** por variável/msg_type.

Critério de "aprovado"
----------------------
Um subset é considerado **aprovado** se **todas** as checagens aplicáveis a ele
passarem: header OK, dentro da janela, possui evento final, unidades plausíveis
e, se solicitado, erros de observação presentes e campos ADPSFC válidos.

Saída e código de retorno
-------------------------
- Mostra um resumo geral e inventário por ``msg_type``.
- Retorna **0** se o status geral for **PASSOU** (nenhuma falha nos contadores
  principais habilitados); caso contrário, retorna **1**.

Requisitos
----------
- Python 3.8+
- ``nceplibs-bufr`` (conda-forge) para o módulo ``ncepbufr``.
- Opcional: ``tqdm`` para a barra de progresso (``--progress bar``).

Exemplos
--------
Básico:

    python check_prepbufr.py seu.prepbufr

Com barra de progresso e janela de ±6 h:

    python check_prepbufr.py seu.prepbufr --progress bar --twind 6

Exigindo erros de observação, checando ADPSFC, autodetectando pressão e
emitindo relatórios:

    python check_prepbufr.py seu.prepbufr \
      --check-oberrs --kind adpsfc --pressure-unit auto \
      --where 5 --csv por_tipo.csv \
      --report-csv falhas_detalhe.csv --vars-csv falhas_unidades.csv
"""

from __future__ import annotations

import sys
import csv
import argparse
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from statistics import median
from collections import defaultdict
import ncepbufr  # conda-forge: nceplibs-bufr

# --------------------------------- Config -------------------------------------

#: Mnemonics de header que devem existir em cada subset.
REQ_HDR: tuple[str, ...] = tuple("XOB YOB DHR TYP ELV SAID T29".split())

#: Sequências Table D normalmente presentes em PREPBUFR (usadas na detecção de “evento final”).
REQ_SEQS: tuple[str, ...] = ("TEVN", "WEVN", "PEVN", "ZEVN", "QEVN")

#: Meia-largura padrão da janela de tempo ``|DHR|`` (em horas).
DEFAULT_TWIND_H: float = 3.0

#: Periodicidade de log do contador simples (quando ``--progress counter``).
DEFAULT_PROGRESS_EVERY: int = 10_000

#: Intervalos de **plausibilidade** (validação de unidades por variável).
#: Obs.: são limites operacionais razoáveis (não normas rígidas WMO).
RANGE = {
    "XOB": (-180.0, 180.0),    # graus (longitude)
    "YOB": (-90.0, 90.0),      # graus (latitude)
    "TOB": (180.0, 340.0),     # Kelvin (temperatura)
    "POB": (1.0, 1100.0),      # hPa (pressão em nível)
    "ZOB": (-500.0, 40000.0),  # m (altura/geopotencial) (teto mais realista p/ obs; 40 km)
    "UOB": (-120.0, 120.0),    # m/s (vento U)
    "VOB": (-120.0, 120.0),    # m/s (vento V)
    "QOB": (0.0, 0.05),        # kg/kg (umidade específica) (50 g/kg já é extremo tropical)
    # Extras típicos de ADPSFC:
    "PRSS": (500.0, 1085.0),   # hPa (pressão à superfície) (superfície realista; evita unidade errada)
    "PWO": (0.0, 100.0),       # mm (ou kg/m^2) — apenas plausibilidade
    "DHR":  (-6.0, 6.0),       # h (janela típica de ±6h)
    "HOVI": (0.0, 100000.0),   # m (visibilidade)
    "TDO":  (173.0, 333.0),    # K (ponto de orvalho; -100 a +60 °C ~ plausível)
    "MSST": (271.0, 313.0),    # K (SST de -2 a 40 °C)
}

# Unidade de pressão (para interpretar POB/PRSS do arquivo).
# A validação usa ranges base em hPa; aplicamos um fator de escala para
# converter o valor lido para hPa antes de validar.
PRESSURE_UNITS = ("hpa", "cb", "pa", "auto")
DEFAULT_PRESSURE_UNIT = "hpa"
PRESSURE_AUTODETECT_SAMPLES = 500  # nº máximo de amostras para auto
PRESSURE_AUTODETECT_THRESH = ( (50,1100), (5,110), (5000,110000) )  # hPa, cb, Pa janelas plausíveis

# Tipos de falha reconhecidos para o relatório detalhado
FAIL_HDR      = "bad_header"
FAIL_EVENTS   = "no_final_event"
FAIL_TIME     = "time_out_of_window"
FAIL_UNITS    = "bad_units"
FAIL_OBERRS   = "missing_oberrs"
FAIL_ADPSFC   = "missing_adpsfc"

@dataclass
class FailRecord:
    msg_type: str
    subset_idx: int
    fail_type: str
    var: str | None

# --------------------------------- Model --------------------------------------

@dataclass
class Summary:
    """Resumo final das checagens executadas sobre o arquivo PREPBUFR.

    Attributes
    ----------
    nmsg : int
        Número de mensagens BUFR no arquivo.
    nsub : int
        Número total de subsets (relatórios) processados.
    bad_hdr : int
        Subsets com header incompleto (faltando algum de REQ_HDR).
    no_events : int
        Subsets sem **qualquer** evento final detectado (TEVN/WEVN/PEVN/ZEVN/QEVN).
    time_oow : int
        Subsets cujo ``DHR`` está fora da janela especificada (|DHR| > twind_h) ou ausente.
    bad_units : int
        Subsets que apresentaram valores fora dos intervalos plausíveis (unidades).
    miss_oberrs : int
        Subsets em que faltam erros de observação (apenas se ``--check-oberrs``).
    miss_adpsfc : int
        Subsets ADPSFC faltando ``PRSS/PWO/CAT`` ou com ``PRSS/PWO`` fora de faixa (se ``--kind adpsfc``).
    types : Dict[object, int]
        Inventário bruto de subsets por ``msg_type`` (chave pode ser bytes).
    per_type_total : Dict[object, int]
        Total de subsets por ``msg_type``.
    per_type_approved : Dict[object, int]
        Subsets aprovados por ``msg_type`` (conforme critério descrito no docstring).
    per_type_bad_hdr : Dict[object, int]
        Contagem de header incompleto por ``msg_type``.
    per_type_no_events : Dict[object, int]
        Contagem de “sem eventos finais” por ``msg_type``.
    per_type_time_oow : Dict[object, int]
        Contagem de “fora da janela de tempo” por ``msg_type``.
    per_type_bad_units : Dict[object, int]
        Contagem de “unidades/valores implausíveis” por ``msg_type``.
    per_type_miss_oberrs : Dict[object, int]
        Contagem de “faltando erros de observação” por ``msg_type`` (se habilitado).
    per_type_miss_adpsfc : Dict[object, int]
        Contagem de “faltando/fora de faixa PRSS/PWO/CAT” por ``msg_type`` (ADPSFC).
    status : str
        "PASSOU" se todos os contadores aplicáveis forem zero; caso contrário "ATENÇÃO".
    pressure_unit_used : str
        Unidade de pressão efetivamente usada (após `auto` ou default).
    where_* : list[str]
        Exemplos "onde" (linhas já formatadas) para cada tipo de falha.
    fail_records : list[FailRecord]
        Registros detalhados de falhas (para os CSVs detalhados).
    """
    nmsg: int
    nsub: int
    bad_hdr: int
    no_events: int
    time_oow: int
    bad_units: int
    miss_oberrs: int
    miss_adpsfc: int
    types: Dict[object, int]
    per_type_total: Dict[object, int]
    per_type_approved: Dict[object, int]
    per_type_bad_hdr: Dict[object, int]
    per_type_no_events: Dict[object, int]
    per_type_time_oow: Dict[object, int]
    per_type_bad_units: Dict[object, int]
    per_type_miss_oberrs: Dict[object, int]
    per_type_miss_adpsfc: Dict[object, int]
    status: str
    pressure_unit_used: str
    # Amostras "onde" (até N por tipo de falha)
    where_hdr: list[str]
    where_events: list[str]
    where_time: list[str]
    where_units: list[str]
    # Registros detalhados para CSV
    fail_records: list[FailRecord]

# ------------------------------ Progress utils --------------------------------

def count_subsets(path: str) -> int:
    """Conta rapidamente o número total de subsets, somando ``nsubsets`` por mensagem.

    Faz um *pré-passo* sem carregar subsets individualmente, útil para configurar a
    barra de progresso com total conhecido.
    """
    b = ncepbufr.open(path)
    total = 0
    try:
        adv = b.advance
        while adv() == 0:
            total += getattr(b, "nsubsets", 0)
    finally:
        b.close()
    return total

def _maybe_get_tqdm():
    """Tenta importar ``tqdm`` e retorna o construtor da barra se disponível."""
    try:
        from tqdm import tqdm  # type: ignore
        return tqdm
    except Exception:
        return None

# ------------------------------ Low-level helpers -----------------------------

def _read1(read, mnem: str) -> Optional[np.ndarray]:
    """Lê um mnemônico com ``read`` e retorna um ``np.ndarray`` (ou ``None``)."""
    arr = read(mnem)
    if getattr(arr, "size", 0):
        return np.asarray(arr)
    return None

def _read_scalar(read, mnem: str) -> Optional[float]:
    """Lê um mnemônico escalar (primeiro elemento) e retorna ``float`` (ou ``None``)."""
    a = _read1(read, mnem)
    if a is None:
        return None
    try:
        return float(a.flat[0])
    except Exception:
        return None

def _has_header_fast(read, required: tuple[str, ...]) -> bool:
    """Retorna ``True`` se **todos** os mnemonics de ``required`` estiverem presentes."""
    for m in required:
        arr = read(m)
        if not getattr(arr, "size", 0):
            return False
    return True

def _last_event(seq: Optional[np.ndarray]) -> tuple[Optional[np.ndarray], int]:
    """Normaliza a sequência para 2D e retorna (matriz, n_eventos); (None, 0) se inválida."""
    if seq is None:
        return None, 0
    a = np.asarray(seq)
    if a.ndim == 2:
        return a, a.shape[1]
    return None, 0

def _in_range(val: Optional[float], lo: float, hi: float) -> bool:
    """Checa se ``val`` é finito e está no intervalo fechado ``[lo, hi]``."""
    if val is None or not np.isfinite(val):
        return False
    return (val >= lo) and (val <= hi)

# --------------------------------- Core ---------------------------------------

def _b2s(x) -> str:
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

def _fmt_where(mtyp: object, idx: int, xob: Optional[float], yob: Optional[float], dhr: Optional[float], extra: str="") -> str:
    xs = "nan" if xob is None else f"{xob:.3f}"
    ys = "nan" if yob is None else f"{yob:.3f}"
    hs = "nan" if dhr is None else f"{dhr:.2f}"
    return f"  - {_b2s(mtyp)} subset#{idx} @ lon={xs} lat={ys} dhr={hs}{extra}"

def _pressure_scale(unit: str) -> float:
    """Retorna fator p/ converter valores lidos → hPa.
       hpa: 1.0 ; cb: 10.0 (1 cb = 10 hPa) ; pa: 1/100 (100 Pa = 1 hPa)
    """
    unit = (unit or "").lower()
    if unit == "hpa":
        return 1.0
    if unit == "cb":
        return 10.0
    if unit == "pa":
        return 1.0/100.0
    return 1.0

def _autodetect_pressure_unit(path: str, max_samples: int = PRESSURE_AUTODETECT_SAMPLES) -> str:
    """Lê rapidamente `POB` e `PRSS` para estimar unidade (hpa/cb/pa).
       Heurística baseada na mediana de valores positivos.
    """
    vals = []
    b = ncepbufr.open(path)
    try:
        cnt = 0
        while b.advance() == 0 and cnt < max_samples:
            while b.load_subset() == 0 and cnt < max_samples:
                for nm in ("POB", "PRSS"):
                    arr = b.read_subset(nm)
                    if getattr(arr, "size", 0):
                        try:
                            v = float(arr[0][0])
                            if np.isfinite(v) and v > 0:
                                vals.append(v)
                                cnt += 1
                                if cnt >= max_samples:
                                    break
                        except Exception:
                            pass
    finally:
        b.close()
    if not vals:
        return "hpa"  # fallback razoável
    med = median(vals)
    hpa_lo, hpa_hi = PRESSURE_AUTODETECT_THRESH[0]
    cb_lo, cb_hi   = PRESSURE_AUTODETECT_THRESH[1]
    pa_lo, pa_hi   = PRESSURE_AUTODETECT_THRESH[2]
    if hpa_lo <= med <= hpa_hi: return "hpa"
    if cb_lo  <= med <= cb_hi:  return "cb"
    if pa_lo  <= med <= pa_hi:  return "pa"
    return "hpa"  # fallback

def check_file(
    path: str,
    twind_h: float = DEFAULT_TWIND_H,
    *,
    progress: str = "counter",
    progress_every: int = DEFAULT_PROGRESS_EVERY,
    quiet: bool = False,
    check_oberrs: bool = False,
    pressure_unit: str = DEFAULT_PRESSURE_UNIT,  # 'hpa' | 'cb' | 'pa' | 'auto'
    where_max: int = 0,                  # imprime até N exemplos "onde" por tipo
    report_csv: Optional[str] = None,    # CSV detalhado linha-a-linha das falhas
    vars_csv: Optional[str] = None,      # CSV agregado: contagem por variável/msg_type
    kind: Optional[str] = None,
    csv_out: Optional[str] = None,
) -> Summary:
    """Executa todas as checagens sobre um arquivo PREPBUFR."""
    # Configura barra de progresso se solicitado
    total_subsets: Optional[int] = None
    pbar = None
    tqdm = None

    if not quiet and progress == "bar":
        tqdm = _maybe_get_tqdm()
        if tqdm is None:
            print("[WARN] tqdm não encontrado; usando progresso 'counter'. "
                  "Instale com: pip install tqdm", file=sys.stderr)
            progress = "counter"
        else:
            # pré-passo para contar subsets
            total_subsets = count_subsets(path)
            pbar = tqdm(total=total_subsets, desc="Processando subsets", unit="subset")

    # Unidade de pressão: resolve 'auto' antes do processamento principal
    pressure_unit_used = pressure_unit
    if pressure_unit_used not in PRESSURE_UNITS:
        pressure_unit_used = DEFAULT_PRESSURE_UNIT
    if pressure_unit_used == "auto":
        pressure_unit_used = _autodetect_pressure_unit(path)
        if not quiet:
            print(f"[INFO] Unidade de pressão autodetectada: {pressure_unit_used}")

    b = ncepbufr.open(path)

    # buffers para "onde"
    where_hdr: list[str] = []
    where_events: list[str] = []
    where_time: list[str] = []
    where_units: list[str] = []
    fail_records: list[FailRecord] = []
    units_by_var_and_type = defaultdict(int)  # (msg_type,mnem) -> count

    # Contadores globais
    nmsg = nsub = 0
    bad_hdr = 0
    no_events = 0
    time_oow = 0
    bad_units = 0
    miss_oberrs = 0
    miss_adpsfc = 0

    # Inventários
    types: Dict[object, int] = {}
    per_type_total: Dict[object, int] = {}
    per_type_approved: Dict[object, int] = {}
    per_type_bad_hdr: Dict[object, int] = {}
    per_type_no_events: Dict[object, int] = {}
    per_type_time_oow: Dict[object, int] = {}
    per_type_bad_units: Dict[object, int] = {}
    per_type_miss_oberrs: Dict[object, int] = {}
    per_type_miss_adpsfc: Dict[object, int] = {}

    def inc(d: Dict[object, int], k: object, by: int = 1):
        """Incrementa d[k] em ``by`` (cria a chave se necessário)."""
        d[k] = d.get(k, 0) + by

    pscale = _pressure_scale(pressure_unit_used)  # escala para converter → hPa

    try:
        advance = b.advance
        load_subset = b.load_subset
        read = b.read_subset

        # Loop principal: mensagens → subsets
        while advance() == 0:
            nmsg += 1
            mtyp = b.msg_type
            if mtyp not in types:
                types[mtyp] = 0
            if mtyp not in per_type_total:
                per_type_total[mtyp] = 0
                per_type_approved[mtyp] = 0
                per_type_bad_hdr[mtyp] = 0
                per_type_no_events[mtyp] = 0
                per_type_time_oow[mtyp] = 0
                per_type_bad_units[mtyp] = 0
                per_type_miss_oberrs[mtyp] = 0
                per_type_miss_adpsfc[mtyp] = 0

            while load_subset() == 0:
                nsub += 1
                types[mtyp] += 1
                # Leitura básica para relatórios "onde"
                xob = _read_scalar(read, "XOB")
                yob = _read_scalar(read, "YOB")
                dhr_dbg = _read_scalar(read, "DHR")
                inc(per_type_total, mtyp)

                # --- Progresso
                if not quiet:
                    if pbar is not None:
                        pbar.update(1)
                    elif progress == "counter" and progress_every > 0 and (nsub % progress_every == 0):
                        print(f"[INFO] {nsub:,} subsets processados...", flush=True)

                # ---- Header
                subset_bad_hdr = False
                if not _has_header_fast(read, REQ_HDR):
                    bad_hdr += 1
                    inc(per_type_bad_hdr, mtyp)
                    subset_bad_hdr = True
                    if where_max and len(where_hdr) < where_max:
                        where_hdr.append(_fmt_where(mtyp, nsub, xob, yob, dhr_dbg))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_HDR, None))

                # ---- Janela de tempo
                subset_time_oow = False
                if dhr_dbg is None or not (-twind_h <= dhr_dbg <= twind_h):
                    time_oow += 1
                    inc(per_type_time_oow, mtyp)
                    subset_time_oow = True
                    if where_max and len(where_time) < where_max:
                        where_time.append(_fmt_where(mtyp, nsub, xob, yob, dhr_dbg))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_TIME, None))

                # ---- Eventos
                subset_no_events = False
                seq_ok = False

                # Tenta TEVN/WEVN com validação de QM
                for seq_name in ("TEVN", "WEVN"):
                    seq = _read1(read, seq_name)
                    a, nev = _last_event(seq)
                    if nev >= 1 and a is not None:
                        if seq_name == "TEVN":
                            # a[1, -1] = TQM
                            if a.shape[0] >= 2 and np.isfinite(a[1, -1]):
                                seq_ok = True
                                break
                        else:  # WEVN
                            # a[2, -1] = WQM
                            if a.shape[0] >= 3 and np.isfinite(a[2, -1]):
                                seq_ok = True
                                break

                # Se ainda não ok, aceita PEVN/ZEVN/QEVN (só presença de evento)
                if not seq_ok:
                    for seq_name in ("PEVN", "ZEVN", "QEVN"):
                        seq = _read1(read, seq_name)
                        _a, nev = _last_event(seq)
                        if nev >= 1:
                            seq_ok = True
                            break

                if not seq_ok:
                    no_events += 1
                    inc(per_type_no_events, mtyp)
                    subset_no_events = True
                    if where_max and len(where_events) < where_max:
                        where_events.append(_fmt_where(mtyp, nsub, xob, yob, dhr_dbg))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_EVENTS, None))

                # ---- Unidades/valores plausíveis
                subset_bad_units = False

                def check_range(name: str) -> bool:
                    """Valida um mnemônico escalar existente contra RANGE[name]."""
                    a = _read_scalar(read, name)
                    lo, hi = RANGE[name]
                    # Pressões precisam de conversão para hPa antes de validar
                    if name in ("POB", "PRSS") and a is not None:
                        a = a * pscale
                    return _in_range(a, lo, hi)

                # testa cada variável presente; registra **todas** que falharem
                for nm in ("XOB", "YOB", "TOB", "POB", "ZOB", "UOB", "VOB", "QOB"):
                    arr = _read1(read, nm)
                    if arr is not None:
                        if not check_range(nm):
                            subset_bad_units = True
                            units_by_var_and_type[(_b2s(mtyp), nm)] += 1
                            if where_max and len(where_units) < where_max:
                                where_units.append(_fmt_where(mtyp, nsub, xob, yob, dhr_dbg, extra=f" var={nm}"))
                            fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_UNITS, nm))
                            # (não faz break; queremos registrar todas as variáveis fora de faixa)

                if subset_bad_units:
                    bad_units += 1
                    inc(per_type_bad_units, mtyp)

                # ---- Erros de observação (opcional)
                subset_miss_oberrs = False
                if check_oberrs:
                    # Para cada variável presente, exige o erro correspondente
                    pairs = [
                        ("TOB", "TOE"),
                        (("UOB", "VOB"), "WOE"),
                        ("POB", "POE"),
                        ("ZOB", "ZOE"),
                        ("QOB", "QOE"),
                    ]
                    for obs, err in pairs:
                        if isinstance(obs, tuple):
                            present = any(_read1(read, o) is not None for o in obs)
                        else:
                            present = _read1(read, obs) is not None
                        if present and _read1(read, err) is None:
                            subset_miss_oberrs = True
                            break

                    if subset_miss_oberrs:
                        miss_oberrs += 1
                        inc(per_type_miss_oberrs, mtyp)
                        fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_OBERRS, None))

                # ---- ADPSFC: PRSS/PWO/CAT (opcional por --kind)
                subset_miss_adpsfc = False
                if kind == "adpsfc":
                    # Requer PRSS/PWO/CAT quando o subset é “significativo”
                    need_check = not (subset_bad_hdr or subset_no_events)
                    if need_check:
                        prss = _read1(read, "PRSS")
                        pwo = _read1(read, "PWO")
                        cat = _read1(read, "CAT")
                        if (prss is None) or (pwo is None) or (cat is None):
                            subset_miss_adpsfc = True
                        else:
                            # valida ranges se presentes
                            # PRSS pode não estar em hPa → converte
                            try:
                                prss_hpa = float(prss.flat[0]) * pscale
                            except Exception:
                                prss_hpa = None
                            if not _in_range(prss_hpa, *RANGE["PRSS"]):
                                subset_miss_adpsfc = True
                            if not _in_range(float(pwo.flat[0]), *RANGE["PWO"]):
                                subset_miss_adpsfc = True

                        if subset_miss_adpsfc:
                            miss_adpsfc += 1
                            inc(per_type_miss_adpsfc, mtyp)
                            fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_ADPSFC, None))

                # ---- Aprovação do subset
                approved = not any([
                    subset_bad_hdr,
                    subset_no_events,
                    subset_time_oow,
                    subset_bad_units,
                    subset_miss_oberrs,
                    subset_miss_adpsfc,
                ])
                if approved:
                    inc(per_type_approved, mtyp)

    finally:
        try:
            b.close()
        finally:
            if pbar is not None:
                pbar.close()

    # Determina o status geral
    status = "PASSOU" if all(x == 0 for x in (
        bad_hdr, no_events, time_oow, bad_units,
        (miss_oberrs if check_oberrs else 0),
        (miss_adpsfc if (kind == "adpsfc") else 0),
    )) else "ATENÇÃO"

    summary = Summary(
        nmsg=nmsg,
        nsub=nsub,
        bad_hdr=bad_hdr,
        no_events=no_events,
        time_oow=time_oow,
        bad_units=bad_units,
        miss_oberrs=miss_oberrs,
        miss_adpsfc=miss_adpsfc,
        types=types,
        per_type_total=per_type_total,
        per_type_approved=per_type_approved,
        per_type_bad_hdr=per_type_bad_hdr,
        per_type_no_events=per_type_no_events,
        per_type_time_oow=per_type_time_oow,
        per_type_bad_units=per_type_bad_units,
        per_type_miss_oberrs=per_type_miss_oberrs,
        per_type_miss_adpsfc=per_type_miss_adpsfc,
        status=status,
        pressure_unit_used=pressure_unit_used,
        where_hdr=where_hdr,
        where_events=where_events,
        where_time=where_time,
        where_units=where_units,
        fail_records=fail_records,
    )

    # CSVs (opcionais)
    if csv_out:
        _write_csv(csv_out, summary)
    if report_csv:
        _write_report_csv(report_csv, path, pressure_unit_used, summary)
    if vars_csv:
        _write_vars_csv(vars_csv, summary, units_by_var_and_type)

    return summary

# --------------------------------- CSV ----------------------------------------

def _write_csv(path: str, s: Summary) -> None:
    """Escreve um CSV por ``msg_type`` com totais, aprovados, % aprovado e contadores de falhas."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "msg_type", "total", "approved", "approved_pct",
            "bad_hdr", "no_events", "time_oow", "bad_units",
            "missing_oberrs", "missing_adpsfc_fields"
        ])
        # A ordenação abaixo mantém bytes primeiro, depois strings, de forma estável
        for k in sorted(s.per_type_total.keys(), key=lambda x: (isinstance(x, bytes), x)):
            mt = k.decode() if isinstance(k, (bytes, bytearray)) else str(k)
            tot = s.per_type_total.get(k, 0)
            app = s.per_type_approved.get(k, 0)
            pct = (100.0 * app / tot) if tot > 0 else 0.0
            w.writerow([
                mt, tot, app, f"{pct:.2f}",
                s.per_type_bad_hdr.get(k, 0),
                s.per_type_no_events.get(k, 0),
                s.per_type_time_oow.get(k, 0),
                s.per_type_bad_units.get(k, 0),
                s.per_type_miss_oberrs.get(k, 0),
                s.per_type_miss_adpsfc.get(k, 0),
            ])

def _write_report_csv(path: str, src_path: str, press_unit: str, s: Summary) -> None:
    """CSV detalhado com **todas** as falhas detectadas.
       Colunas: src_file, pressure_unit, msg_type, subset_idx, fail_type, var
    """
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["src_file", "pressure_unit", "msg_type", "subset_idx", "fail_type", "var"])
        for rec in s.fail_records:
            w.writerow([src_path, press_unit, rec.msg_type, rec.subset_idx, rec.fail_type, (rec.var or "")])

def _write_vars_csv(path: str, s: Summary, counter: dict) -> None:
    """CSV agregado: contagem de falhas de unidade por `msg_type` e variável."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["msg_type", "var", "bad_units_count"])
        # counter: (msg_type, var) -> count
        for (mt, var), cnt in sorted(counter.items()):
            w.writerow([mt, var, cnt])

# --------------------------------- CLI ----------------------------------------

def _fmt_types(types: Dict[object, int]) -> str:
    """Formata inventário por ``msg_type`` em linhas legíveis."""
    lines = []
    for k, v in sorted(types.items(), key=lambda kv: kv[0]):
        k_disp = k.decode() if isinstance(k, (bytes, bytearray)) else str(k)
        lines.append(f" - {k_disp}: {v} subsets")
    return "\n".join(lines)

def _build_parser() -> argparse.ArgumentParser:
    """Constrói o parser de argumentos da linha de comando."""
    p = argparse.ArgumentParser(
        prog="check_prepbufr.py",
        description=(
            "Checagens estruturais e de plausibilidade em PREPBUFR: "
            "header, eventos finais, janela de tempo, faixas/unidades, "
            "erros de observação e campos ADPSFC. Pode gerar CSV por msg_type, "
            "CSV detalhado de falhas e CSV agregado por variável."
        ),
    )
    p.add_argument("file", help="Caminho do PREPBUFR")
    p.add_argument(
        "--twind",
        type=float,
        default=DEFAULT_TWIND_H,
        metavar="HOURS",
        help=f"Meia-largura da janela |DHR| (padrão: ±{DEFAULT_TWIND_H} h).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Oculta saída humana; usa apenas o código de saída.",
    )
    p.add_argument(
        "--progress",
        choices=("off", "counter", "bar"),
        default="counter",
        help="Tipo de progresso: off | counter | bar (precisa do pacote tqdm).",
    )
    p.add_argument(
        "--progress-every",
        type=int,
        default=DEFAULT_PROGRESS_EVERY,
        metavar="N",
        help="(counter) imprime a cada N subsets (0 desliga).",
    )
    # Novos parâmetros
    p.add_argument(
        "--check-oberrs",
        action="store_true",
        help="Exige erros de observação (TOE/WOE/POE/ZOE/QOE) quando a variável está presente.",
    )
    p.add_argument(
        "--kind",
        choices=("adpsfc", "adpupa"),
        help="Checagens específicas por tipo lógico. Em 'adpsfc' confere PRSS/PWO/CAT.",
    )
    p.add_argument(
        "--pressure-unit",
        choices=PRESSURE_UNITS,
        default=DEFAULT_PRESSURE_UNIT,
        help=(
            "Unidade de pressão no arquivo: hpa (padrão), cb, pa ou auto. "
            "Usada para validar POB/PRSS com ranges base em hPa."
        ),
    )
    p.add_argument(
        "--where",
        type=int,
        default=0,
        metavar="N",
        help="Imprime até N exemplos 'onde' por tipo de falha (header/eventos/tempo/unidades).",
    )
    p.add_argument(
        "--report-csv",
        help="CSV detalhado com todas as falhas (msg_type, subset_idx, tipo, var).",
    )
    p.add_argument(
        "--vars-csv",
        help="CSV agregado com contagem de falhas de unidade por variável e msg_type.",
    )
    p.add_argument(
        "--csv",
        help="Escreve relatório CSV por msg_type (caminho do arquivo).",
    )
    return p

def main(argv: Optional[list[str]] = None) -> int:
    """Ponto de entrada do script (CLI).

    Lê os argumentos, executa as checagens e imprime o resumo.
    O código de retorno é 0 se o status geral for *PASSOU*; caso contrário 1.
    """
    args = _build_parser().parse_args(argv)
    # Se quiet, suprime progresso
    progress = "off" if args.quiet else args.progress

    res = check_file(
        args.file,
        twind_h=args.twind,
        progress=progress,
        progress_every=args.progress_every,
        quiet=args.quiet,
        check_oberrs=args.check_oberrs,
        pressure_unit=args.pressure_unit,
        where_max=args.where,
        report_csv=args.report_csv,
        vars_csv=args.vars_csv,
        kind=args.kind,
        csv_out=args.csv,
    )

    if not args.quiet:
        print(f"Arquivo: {args.file}")
        print(f"Unidade de pressão adotada: {res.pressure_unit_used}")
        print(f"Mensagens: {res.nmsg}  |  Subsets: {res.nsub}")
        print(_fmt_types(res.types))
        print("\nResumo:")
        print(f"  Header(s) incompleto(s): {res.bad_hdr}")
        print(f"  Subsets sem eventos finais (TEVN/WEVN/etc.): {res.no_events}")
        print(f"  Subsets fora da janela de tempo (±{args.twind}h): {res.time_oow}")
        print(f"  Subsets com unidades/valores implausíveis: {res.bad_units}")
        if args.check_oberrs:
            print(f"  Subsets faltando erros de observação: {res.miss_oberrs}")
        if args.kind == "adpsfc":
            print(f"  Subsets ADPSFC sem PRSS/PWO/CAT (ou fora de faixa): {res.miss_adpsfc}")

        # Onde (exemplos)
        if args.where:
            if res.where_hdr:
                print("\nOnde (header incompleto):")
                for line in res.where_hdr: print(line)
            if res.where_events:
                print("\nOnde (sem eventos finais):")
                for line in res.where_events: print(line)
            if res.where_time:
                print("\nOnde (fora da janela de tempo):")
                for line in res.where_time: print(line)
            if res.where_units:
                print("\nOnde (unidades/valores implausíveis):")
                for line in res.where_units: print(line)

        if args.report_csv:
            print(f"\nCSV detalhado de falhas: {args.report_csv}")
        if args.vars_csv:
            print(f"CSV agregado (falhas de unidade por variável): {args.vars_csv}")
        if args.csv:
            print(f"CSV (por msg_type) escrito em: {args.csv}")
        print(f"\nStatus: {res.status}")

    return 0 if res.status == "PASSOU" else 1

if __name__ == "__main__":
    sys.exit(main())


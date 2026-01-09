#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_prepbufr.py
=================

Validador unificado para arquivos PREPBUFR, com dois níveis de checagem:

1) **Estrutural (mínimo GSI)** → o que realmente importa para saber se o GSI
   consegue ler e usar o arquivo sem quebrar:

   - header essencial presente: XOB, YOB, DHR, TYP;
   - lon/lat grosseiramente válidos;
   - existe pelo menos uma variável "útil":
       TOB, QOB, UOB, VOB, POB, ZOB ou PRSS.

   Esse nível define:
   - `Status estrutural (GSI): OK / PROBLEMAS_ESTRUTURAIS`;
   - código de saída do programa (0 se OK, 1 se PROBLEMAS_ESTRUTURAIS);
   - campo `approved` e `% approved` no CSV (`--csv`).

2) **Diagnóstico** → checagens que **não impedem** o GSI de rodar, mas ajudam
   a entender se os dados estão "bonitos" e coerentes:

   - presença (ou não) de eventos finais (TEVN/WEVN/PEVN/ZEVN/QEVN);
   - janela de tempo: |DHR| <= twind (padrão ±3h);
   - faixas/unidades plausíveis de XOB/YOB/TOB/POB/ZOB/UOB/VOB/QOB/PRSS/PWO/etc;
   - (opcional via --check-oberrs) presença de erros de observação:
       TOB→TOE, U/V→WOE, POB→POE, ZOB→ZOE, QOB→QOE;
   - (opcional via --kind adpsfc) campos de superfície (PRSS/PWO/CAT).

   Esse nível define:
   - `Status diagnóstico: OK / ATENÇÃO`;
   - contadores em `Resumo diagnóstico`;
   - colunas extras no CSV e detalhes em `--report-csv` / `--vars-csv`.

A ideia é:

- **GSI** se preocupa com: "o arquivo é estruturalmente utilizável?"  
- Você, como produtor de PREPBUFR, pode ainda se preocupar com:
  "minhas variáveis/erros/eventos/unidades fazem sentido?".

Requisitos
----------
- Python 3.8+
- nceplibs-bufr (ncepbufr) via conda-forge
- numpy
- (opcional) tqdm para barra de progresso (--progress bar)
"""

from __future__ import annotations

import sys
import csv
import argparse
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List

import numpy as np
from statistics import median
from collections import defaultdict

import ncepbufr  # conda-forge: nceplibs-bufr

# -----------------------------------------------------------------------------
# Configuração geral
# -----------------------------------------------------------------------------

# Header mínimo para o GSI
# GSI REAL exige apenas posição e tempo
REQ_HDR_GSI: Tuple[str, ...] = ("XOB", "YOB", "DHR")

# Sequências Table D para eventos (diagnóstico)
REQ_SEQS: Tuple[str, ...] = ("TEVN", "WEVN", "PEVN", "ZEVN", "QEVN")

# Janela de tempo típica
DEFAULT_TWIND_H: float = 3.0

# Progresso
DEFAULT_PROGRESS_EVERY: int = 10_000

# Faixas de plausibilidade baseadas estritamente na tabela PREPBUFR

RANGE = {
    # Geolocalização
    "XOB": (-180.0, 180.0),          # DEG E
    "YOB": (-90.0, 90.0),            # DEG N

    # Tempo
    "DHR": (-24.0, 24.0),            # HOURS (tabela permite até ±24h)

    # Temperatura (DEG C)
    "TOB": (-90.0, 60.0),            # Atmosfera (°C)
    "TDO": (-100.0, 40.0),           # Ponto de orvalho (°C)

    # Temperatura de superfície (K)
    "TMSK": (200.0, 350.0),          # K (compatível com tabela)

    # Pressão
    "POB": (1.0, 1100.0),             # mb
    "PRSS": (50000.0, 108500.0),      # Pa (ATENÇÃO: Pascals!)

    # Altura
    "ZOB": (-500.0, 40000.0),         # m (ref = -1000 m)

    # Vento
    "UOB": (-150.0, 150.0),           # m/s
    "VOB": (-150.0, 150.0),           # m/s

    # Umidade específica
    "QOB": (0.0, 30000.0),            # mg/kg (0–30 g/kg)

    # Água precipitável
    "PWO": (0.0, 120.0),              # mm (valores extremos tropicais)
}


# Unidade de pressão (para interpretar POB/PRSS) e autodetecção
PRESSURE_UNITS = ("hpa", "cb", "pa", "auto")
DEFAULT_PRESSURE_UNIT = "hpa"
PRESSURE_AUTODETECT_SAMPLES = 500
PRESSURE_AUTODETECT_THRESH = (
    (50, 1100),     # hPa
    (5, 110),       # cb
    (5000, 110000)  # Pa
)

# Limites estruturais para lon/lat (mínimo GSI)
# GSI aceita lon 0–360 e normaliza internamente
LON_RANGE_STRUCT = (-180.0, 360.0)
LAT_RANGE_STRUCT = (-90.0, 90.0)

# Variáveis assimiláveis por tipo de mensagem (Tabela A / comportamento GSI)
USEFUL_BY_TYPE = {
    "ADPUPA": ("TOB", "QOB", "UOB", "VOB", "ZOB"),
    "AIRCAR": ("TOB", "UOB", "VOB"),
    "AIRCFT": ("TOB", "UOB", "VOB"),
    "SATWND": ("UOB", "VOB"),
    "PROFLR": ("UOB", "VOB"),
    "VADWND": ("UOB", "VOB"),
    "ADPSFC": ("PRSS", "TOB", "QOB", "UOB", "VOB"),
    "SFCSHP": ("PRSS", "TOB", "QOB", "UOB", "VOB"),
    # BUOY: GSI aceita só vento, só pressão ou só temperatura
    "SFCSHP": ("PRSS", "TOB", "QOB", "UOB", "VOB", "POB"),
    "SATEMP": ("TMBR",),
    "GOESND": ("TMBR",),
    "SPSSMI": ("PWO", "UOB", "VOB"),
    "QKSWND": ("UOB", "VOB"),
}


# Missing típico da BUFRLIB
try:
    BUFR_MISSING = float(ncepbufr.bufrlib.getbmiss())
except Exception:
    BUFR_MISSING = 1.0e10

# Tipos de falha (para report-csv)
FAIL_GSI_HDR = "gsi_bad_header"
FAIL_GSI_LONLAT = "gsi_bad_lonlat"
FAIL_GSI_NO_USEFUL = "gsi_no_useful"

FAIL_EVENTS = "no_final_event"
FAIL_TIME = "time_out_of_window"
FAIL_UNITS = "bad_units"
FAIL_OBERRS = "missing_oberrs"
FAIL_ADPSFC = "missing_adpsfc"


@dataclass
class FailRecord:
    msg_type: str
    subset_idx: int
    fail_type: str
    var: Optional[str]


@dataclass
class Summary:
    # Totais
    nmsg: int
    nsub: int

    # Estruturais (mínimo GSI)
    gsi_bad_hdr: int
    gsi_bad_lonlat: int
    gsi_no_useful: int
    per_type_total: Dict[object, int]
    per_type_gsi_bad_hdr: Dict[object, int]
    per_type_gsi_bad_lonlat: Dict[object, int]
    per_type_gsi_no_useful: Dict[object, int]
    per_type_gsi_approved: Dict[object, int]

    # Diagnósticos
    diag_no_events: int
    diag_time_oow: int
    diag_bad_units: int
    diag_miss_oberrs: int
    diag_miss_adpsfc: int
    per_type_no_events: Dict[object, int]
    per_type_time_oow: Dict[object, int]
    per_type_bad_units: Dict[object, int]
    per_type_miss_oberrs: Dict[object, int]
    per_type_miss_adpsfc: Dict[object, int]

    # Status
    status_struct: str
    status_diag: str
    pressure_unit_used: str

    # Amostras "onde"
    where_gsi_hdr: List[str]
    where_gsi_lonlat: List[str]
    where_gsi_no_useful: List[str]
    where_events: List[str]
    where_time: List[str]
    where_units: List[str]

    # Registros detalhados de falhas (para report-csv)
    fail_records: List[FailRecord]


# -----------------------------------------------------------------------------
# Helpers de leitura e missing
# -----------------------------------------------------------------------------

def _read1(read, mnem: str) -> Optional[np.ndarray]:
    arr = read(mnem)
    if getattr(arr, "size", 0):
        return np.asarray(arr)
    return None


def _read_scalar(read, mnem: str) -> Optional[float]:
    a = _read1(read, mnem)
    if a is None:
        return None
    try:
        return float(a.flat[0])
    except Exception:
        return None


def _is_missing(val: Optional[float]) -> bool:
    if val is None:
        return True
    if not np.isfinite(val):
        return True
    if abs(val) >= BUFR_MISSING * 0.9:
        return True
    return False


def _in_range(val: Optional[float], lo: float, hi: float) -> bool:
    if _is_missing(val):
        return False
    return lo <= val <= hi


def _last_event(seq: Optional[np.ndarray]) -> Tuple[Optional[np.ndarray], int]:
    if seq is None:
        return None, 0
    a = np.asarray(seq)
    if a.ndim == 2:
        return a, a.shape[1]
    return None, 0


def _has_gsi_header(read) -> bool:
    """Header mínimo estrutural para o GSI."""
    for m in REQ_HDR_GSI:
        arr = read(m)
        if not getattr(arr, "size", 0):
            return False
    return True


def _b2s(x) -> str:
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def _fmt_where(
    mtyp: object,
    idx: int,
    xob: Optional[float],
    yob: Optional[float],
    dhr: Optional[float],
    extra: str = "",
) -> str:
    xs = "nan" if xob is None else f"{xob:.3f}"
    ys = "nan" if yob is None else f"{yob:.3f}"
    hs = "nan" if dhr is None else f"{dhr:.2f}"
    return f"  - {_b2s(mtyp)} subset#{idx} @ lon={xs} lat={ys} dhr={hs}{extra}"


def _pressure_scale(unit: str) -> float:
    unit = (unit or "").lower()
    if unit == "hpa":
        return 1.0
    if unit == "cb":
        return 10.0       # 1 cb = 10 hPa
    if unit == "pa":
        return 1.0 / 100.0  # 100 Pa = 1 hPa
    return 1.0


# -----------------------------------------------------------------------------
# Progresso
# -----------------------------------------------------------------------------

def _maybe_get_tqdm():
    try:
        from tqdm import tqdm  # type: ignore
        return tqdm
    except Exception:
        return None


def _count_subsets(path: str) -> int:
    b = ncepbufr.open(path)
    total = 0
    try:
        while b.advance() == 0:
            total += getattr(b, "nsubsets", 0)
    finally:
        b.close()
    return total


# -----------------------------------------------------------------------------
# Autodeteção de unidade de pressão
# -----------------------------------------------------------------------------

def _autodetect_pressure_unit(path: str, max_samples: int = PRESSURE_AUTODETECT_SAMPLES) -> str:
    vals: List[float] = []
    b = ncepbufr.open(path)
    try:
        cnt = 0
        while b.advance() == 0 and cnt < max_samples:
            while b.load_subset() == 0 and cnt < max_samples:
                for nm in ("POB",):
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
        return "hpa"

    med = median(vals)
    hpa_lo, hpa_hi = PRESSURE_AUTODETECT_THRESH[0]
    cb_lo, cb_hi = PRESSURE_AUTODETECT_THRESH[1]
    pa_lo, pa_hi = PRESSURE_AUTODETECT_THRESH[2]

    if hpa_lo <= med <= hpa_hi:
        return "hpa"
    if cb_lo <= med <= cb_hi:
        return "cb"
    if pa_lo <= med <= pa_hi:
        return "pa"
    return "hpa"


# -----------------------------------------------------------------------------
# Núcleo: check_file
# -----------------------------------------------------------------------------

def check_file(
    path: str,
    twind_h: float = DEFAULT_TWIND_H,
    *,
    progress: str = "counter",
    progress_every: int = DEFAULT_PROGRESS_EVERY,
    quiet: bool = False,
    check_oberrs: bool = False,
    pressure_unit: str = DEFAULT_PRESSURE_UNIT,
    where_max: int = 0,
    report_csv: Optional[str] = None,
    vars_csv: Optional[str] = None,
    kind: Optional[str] = None,
    csv_out: Optional[str] = None,
) -> Summary:
    # Progresso
    tqdm = None
    pbar = None
    total_subsets: Optional[int] = None

    if not quiet and progress == "bar":
        tqdm = _maybe_get_tqdm()
        if tqdm is None:
            print("[WARN] tqdm não encontrado; usando progresso 'counter'.", file=sys.stderr)
            progress = "counter"
        else:
            total_subsets = _count_subsets(path)
            pbar = tqdm(total=total_subsets, desc="Processando subsets", unit="subset")

    # Unidade de pressão
    pressure_unit_used = pressure_unit if pressure_unit in PRESSURE_UNITS else DEFAULT_PRESSURE_UNIT
    if pressure_unit_used == "auto":
        pressure_unit_used = _autodetect_pressure_unit(path)
        if not quiet:
            print(f"[INFO] Unidade de pressão autodetectada: {pressure_unit_used}")

    pscale = _pressure_scale(pressure_unit_used)

    b = ncepbufr.open(path)

    # Contadores e estruturas
    nmsg = 0
    nsub = 0

    # Estrutural GSI
    gsi_bad_hdr = 0
    gsi_bad_lonlat = 0
    gsi_no_useful = 0
    per_type_total: Dict[object, int] = {}
    per_type_gsi_bad_hdr: Dict[object, int] = {}
    per_type_gsi_bad_lonlat: Dict[object, int] = {}
    per_type_gsi_no_useful: Dict[object, int] = {}
    per_type_gsi_approved: Dict[object, int] = {}

    # Diagnósticos
    diag_no_events = 0
    diag_time_oow = 0
    diag_bad_units = 0
    diag_miss_oberrs = 0
    diag_miss_adpsfc = 0
    per_type_no_events: Dict[object, int] = {}
    per_type_time_oow: Dict[object, int] = {}
    per_type_bad_units: Dict[object, int] = {}
    per_type_miss_oberrs: Dict[object, int] = {}
    per_type_miss_adpsfc: Dict[object, int] = {}

    # "Onde" e detalhes
    where_gsi_hdr: List[str] = []
    where_gsi_lonlat: List[str] = []
    where_gsi_no_useful: List[str] = []
    where_events: List[str] = []
    where_time: List[str] = []
    where_units: List[str] = []
    fail_records: List[FailRecord] = []

    units_by_var_and_type = defaultdict(int)  # (msg_type, var) -> count

    def inc(d: Dict[object, int], k: object, by: int = 1):
        d[k] = d.get(k, 0) + by

    try:
        advance = b.advance
        load_subset = b.load_subset
        read = b.read_subset

        while advance() == 0:
            nmsg += 1
            mtyp = b.msg_type
            if mtyp not in per_type_total:
                per_type_total[mtyp] = 0
                per_type_gsi_bad_hdr[mtyp] = 0
                per_type_gsi_bad_lonlat[mtyp] = 0
                per_type_gsi_no_useful[mtyp] = 0
                per_type_gsi_approved[mtyp] = 0
                per_type_no_events[mtyp] = 0
                per_type_time_oow[mtyp] = 0
                per_type_bad_units[mtyp] = 0
                per_type_miss_oberrs[mtyp] = 0
                per_type_miss_adpsfc[mtyp] = 0

            while load_subset() == 0:
                nsub += 1
                inc(per_type_total, mtyp)

                # progresso
                if not quiet:
                    if pbar is not None:
                        pbar.update(1)
                    elif progress == "counter" and progress_every > 0 and (nsub % progress_every == 0):
                        print(f"[INFO] {nsub:,} subsets processados...", flush=True)

                # leituras básicas para relatórios
                xob = _read_scalar(read, "XOB")
                yob = _read_scalar(read, "YOB")
                dhr = _read_scalar(read, "DHR")

                # -------------------------------
                # 1) Checagem estrutural (mínimo GSI)
                # -------------------------------
                subset_gsi_hdr_bad = False
                subset_gsi_lonlat_bad = False
                subset_gsi_no_useful = False

                # Header essencial
                if not _has_gsi_header(read):
                    subset_gsi_hdr_bad = True
                    gsi_bad_hdr += 1
                    inc(per_type_gsi_bad_hdr, mtyp)
                    if where_max and len(where_gsi_hdr) < where_max:
                        where_gsi_hdr.append(_fmt_where(mtyp, nsub, xob, yob, dhr))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_GSI_HDR, None))

                # lon/lat grosseiramente válidos (se header não está completamente quebrado)
                if not subset_gsi_hdr_bad:
                    if (not _in_range(xob, *LON_RANGE_STRUCT)) or (not _in_range(yob, *LAT_RANGE_STRUCT)):
                        subset_gsi_lonlat_bad = True
                        gsi_bad_lonlat += 1
                        inc(per_type_gsi_bad_lonlat, mtyp)
                        if where_max and len(where_gsi_lonlat) < where_max:
                            where_gsi_lonlat.append(_fmt_where(mtyp, nsub, xob, yob, dhr))
                        fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_GSI_LONLAT, None))

                # Pelo menos uma variável "útil" (dirigido por tipo)
                useful_mnems = USEFUL_BY_TYPE.get(_b2s(mtyp), ())
                has_useful = False
                for nm in useful_mnems:
                    arr = _read1(read, nm)
                    if arr is not None:
                        try:
                            v = float(arr.flat[0])
                            if not _is_missing(v):
                                has_useful = True
                                break
                        except Exception:
                            continue
                if not has_useful:
                    subset_gsi_no_useful = True
                    gsi_no_useful += 1
                    inc(per_type_gsi_no_useful, mtyp)
                    if where_max and len(where_gsi_no_useful) < where_max:
                        where_gsi_no_useful.append(_fmt_where(mtyp, nsub, xob, yob, dhr))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_GSI_NO_USEFUL, None))

                # aprovado estrutural GSI?
                gsi_ok = not (subset_gsi_hdr_bad or subset_gsi_lonlat_bad or subset_gsi_no_useful)
                if gsi_ok:
                    inc(per_type_gsi_approved, mtyp)

                # -------------------------------
                # 2) Diagnósticos
                # -------------------------------

                # 2.1. Eventos
                # NOTA: Eventos EVN são opcionais no GSI real → ignorados
                seq_ok = False
                subset_no_events = False

                # Tenta TEVN/WEVN com QM
                for seq_name in ("TEVN", "WEVN"):
                    seq = _read1(read, seq_name)
                    a, nev = _last_event(seq)
                    if nev >= 1 and a is not None:
                        if seq_name == "TEVN":
                            if a.shape[0] >= 2 and np.isfinite(a[1, -1]):
                                seq_ok = True
                                break
                        else:  # WEVN
                            if a.shape[0] >= 3 and np.isfinite(a[2, -1]):
                                seq_ok = True
                                break
                # fallback: PEVN/ZEVN/QEVN
                if not seq_ok:
                    for seq_name in ("PEVN", "ZEVN", "QEVN"):
                        seq = _read1(read, seq_name)
                        _a, nev = _last_event(seq)
                        if nev >= 1:
                            seq_ok = True
                            break

                if not seq_ok:
                    subset_no_events = True
                    diag_no_events += 1
                    inc(per_type_no_events, mtyp)
                    if where_max and len(where_events) < where_max:
                        where_events.append(_fmt_where(mtyp, nsub, xob, yob, dhr))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_EVENTS, None))

                # 2.2. Janela de tempo
                subset_time_oow = False
                if _is_missing(dhr) or not (-twind_h <= dhr <= twind_h):
                    subset_time_oow = True
                    diag_time_oow += 1
                    inc(per_type_time_oow, mtyp)
                    if where_max and len(where_time) < where_max:
                        where_time.append(_fmt_where(mtyp, nsub, xob, yob, dhr))
                    fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_TIME, None))

                # 2.3. Unidades / valores plausíveis
                subset_bad_units = False

                def check_range(name: str) -> bool:
                    a = _read_scalar(read, name)
                    
                    # Missing não é erro de unidade
                    if _is_missing(a):
                        return True
                    if name in ("POB", "PRSS") and a is not None:
                        a = a * pscale
                    lo, hi = RANGE[name]
                    return _in_range(a, lo, hi)

                for nm in ("XOB", "YOB", "TOB", "POB", "ZOB", "UOB", "VOB", "QOB"):
                    arr = _read1(read, nm)
                    if arr is not None and nm in RANGE:
                        if not check_range(nm):
                            subset_bad_units = True
                            units_by_var_and_type[(_b2s(mtyp), nm)] += 1
                            if where_max and len(where_units) < where_max:
                                where_units.append(_fmt_where(mtyp, nsub, xob, yob, dhr, extra=f" var={nm}"))
                            fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_UNITS, nm))

                if subset_bad_units:
                    diag_bad_units += 1
                    inc(per_type_bad_units, mtyp)

                # 2.4. Erros de observação (diagnóstico)
                subset_miss_oberrs = False
                if check_oberrs:
                    pairs = [
                        # *_OE são opcionais no GSI (erro default)
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
                        diag_miss_oberrs += 1
                        inc(per_type_miss_oberrs, mtyp)
                        fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_OBERRS, None))

                # 2.5. ADPSFC (diagnóstico, se --kind adpsfc)
                subset_miss_adpsfc = False
                # Diagnóstico ADPSFC só se msg_type == ADPSFC
                if kind == "adpsfc" and _b2s(mtyp) == "ADPSFC":
                    if not subset_gsi_hdr_bad:  # só olha ADPSFC se header básico ok
                        prss = _read1(read, "PRSS")
                        pwo = _read1(read, "PWO")
                        if (prss is None) or (pwo is None):
                            subset_miss_adpsfc = True
                        else:
                            try:
                                prss_hpa = float(prss.flat[0]) * pscale
                            except Exception:
                                prss_hpa = None
                            if not _in_range(prss_hpa, *RANGE["PRSS"]):
                                subset_miss_adpsfc = True
                            try:
                                pwo_val = float(pwo.flat[0])
                            except Exception:
                                pwo_val = None
                            if not _in_range(pwo_val, *RANGE["PWO"]):
                                subset_miss_adpsfc = True

                        if subset_miss_adpsfc:
                            diag_miss_adpsfc += 1
                            inc(per_type_miss_adpsfc, mtyp)
                            fail_records.append(FailRecord(_b2s(mtyp), nsub, FAIL_ADPSFC, None))

    finally:
        try:
            b.close()
        finally:
            if pbar is not None:
                pbar.close()

    # -------------------------------------------------------------------------
    # Status estrutural (mínimo GSI)
    #
    # Objetivo:
    #   Determinar se o arquivo PREPBUFR é ESTRUTURALMENTE UTILIZÁVEL pelo GSI,
    #   isto é, se o GSI consegue LER, INTERPRETAR e TENTAR ASSIMILAR os dados.
    #
    # Critérios adotados (estritamente alinhados ao GSI real):
    #
    #   - gsi_bad_hdr:
    #       Subsets sem header essencial (XOB, YOB, DHR, TYP).
    #       -> Sem essas chaves, o GSI NÃO consegue sequer posicionar a observação.
    #
    #   - gsi_bad_lonlat:
    #       Longitude ou latitude claramente inválidas.
    #       -> O GSI rejeita observações fora de limites geográficos básicos.
    #
    #   - gsi_no_useful:
    #       Subsets sem NENHUMA variável potencialmente assimilável
    #       (TOB, QOB, UOB, VOB, POB, ZOB ou PRSS).
    #       -> Observação sem conteúdo útil não entra no sistema.
    #
    # Regra:
    #   - status_struct = "OK"
    #       TODOS os critérios acima são satisfeitos para TODOS os subsets.
    #
    #   - status_struct = "PROBLEMAS_ESTRUTURAIS"
    #       Pelo menos um subset falhou em algum critério estrutural.
    #
    # Observação importante:
    #   Este status define o código de saída do programa:
    #     - 0 → OK (arquivo utilizável pelo GSI)
    #     - 1 → PROBLEMAS_ESTRUTURAIS (arquivo inadequado para o GSI)
    # -------------------------------------------------------------------------
    if gsi_bad_hdr == 0 and gsi_bad_lonlat == 0 and gsi_no_useful == 0:
        status_struct = "OK"
    else:
        status_struct = "PROBLEMAS_ESTRUTURAIS"

    # -------------------------------------------------------------------------
    # Status diagnóstico
    #
    # Objetivo:
    #   Elevar status_diag para "ATENÇÃO" apenas quando houver condições
    #   que possam causar REJEIÇÃO da observação ou IMPACTO OPERACIONAL REAL
    #   no GSI.
    #
    # Princípios adotados (alinhados ao comportamento do GSI):
    #
    #   - diag_no_events:
    #       -> EVN ausente (TEVN/WEVN/PEVN/ZEVN/QEVN)
    #       -> NÃO eleva status.
    #       -> O GSI não exige eventos para leitura ou assimilação.
    #
    #   - diag_miss_oberrs:
    #       -> Erros de observação ausentes (*_OE)
    #       -> NÃO eleva status.
    #       -> O GSI utiliza erros default definidos em convinfo.
    #
    #   - diag_time_oow:
    #       -> DHR fora da janela temporal
    #       -> ELEVA status (ATENÇÃO).
    #       -> Observações fora da janela podem ser descartadas pelo GSI.
    #
    #   - diag_bad_units:
    #       -> Valores ou unidades fisicamente implausíveis
    #       -> ELEVA status (ATENÇÃO).
    #       -> Podem causar rejeição posterior ou contaminar a análise.
    #
    #   - diag_miss_adpsfc:
    #       -> Inconsistência em campos críticos de superfície (PRSS/PWO/CAT)
    #       -> ELEVA status (ATENÇÃO).
    #       -> PRSS é essencial para observações de superfície no GSI.
    #
    # Resultado:
    #   - status_diag = "OK"
    #       Nenhuma condição com impacto real detectada (diagnóstico informativo).
    #   - status_diag = "ATENÇÃO"
    #       Existe pelo menos uma condição com potencial impacto na assimilação.
    # -------------------------------------------------------------------------
    if any([
        diag_time_oow,       # DHR fora da janela temporal
        diag_bad_units,      # valores/unidades fisicamente implausíveis
        diag_miss_adpsfc,    # campos críticos ausentes para ADPSFC
    ]):
        status_diag = "ATENÇÃO"
    else:
        status_diag = "OK"

    # -------------------------------------------------------------------------
    # Consolidação do resumo final da validação
    #
    # O objeto Summary concentra:
    #   - métricas estruturais mínimas (relevantes ao GSI)
    #   - métricas diagnósticas (informativas e/ou de ATENÇÃO)
    #   - status final (estrutural e diagnóstico)
    #   - amostras "onde" ocorreram os problemas
    #   - registros detalhados de falhas (para CSV detalhado)
    # -------------------------------------------------------------------------
    summary = Summary(
        # ---------------------------------------------------------------------
        # Totais globais
        # ---------------------------------------------------------------------
        nmsg=nmsg,                  # número total de mensagens BUFR no arquivo
        nsub=nsub,                  # número total de subsets processados
    
        # ---------------------------------------------------------------------
        # Métricas estruturais (mínimo GSI)
        # Estas métricas indicam se o GSI CONSEGUE ler e usar o PREPBUFR.
        # Qualquer uma delas diferente de zero implica status_struct != OK.
        # ---------------------------------------------------------------------
        gsi_bad_hdr=gsi_bad_hdr,    # subsets sem header essencial (XOB/YOB/DHR/TYP)
        gsi_bad_lonlat=gsi_bad_lonlat,  # lon/lat claramente inválidos (checagem grosseira)
        gsi_no_useful=gsi_no_useful,    # subsets sem NENHUMA variável assimilável
    
        # ---------------------------------------------------------------------
        # Contadores por tipo de mensagem (msg_type)
        # Permitem avaliar qualidade e cobertura por tipo (ADPSFC, SFCSHP, etc.)
        # ---------------------------------------------------------------------
        per_type_total=per_type_total,                      # total de subsets por msg_type
        per_type_gsi_bad_hdr=per_type_gsi_bad_hdr,          # header inválido por msg_type
        per_type_gsi_bad_lonlat=per_type_gsi_bad_lonlat,    # lon/lat inválido por msg_type
        per_type_gsi_no_useful=per_type_gsi_no_useful,      # sem variáveis úteis por msg_type
        per_type_gsi_approved=per_type_gsi_approved,        # subsets estruturalmente OK (GSI)
    
        # ---------------------------------------------------------------------
        # Métricas diagnósticas (não impedem leitura pelo GSI)
        # Servem para QA/QC, debug e avaliação da cadeia de produção.
        # ---------------------------------------------------------------------
        diag_no_events=diag_no_events,           # EVN ausente (TEVN/WEVN/PEVN/ZEVN/QEVN)
        diag_time_oow=diag_time_oow,             # DHR fora da janela temporal
        diag_bad_units=diag_bad_units,           # valores/unidades fisicamente implausíveis
        diag_miss_oberrs=diag_miss_oberrs,       # erros de observação (*_OE) ausentes
        diag_miss_adpsfc=diag_miss_adpsfc,       # inconsistência em PRSS/PWO/CAT (superfície)
    
        # ---------------------------------------------------------------------
        # Métricas diagnósticas por tipo de mensagem
        # Úteis para identificar problemas sistemáticos por origem/tipo.
        # ---------------------------------------------------------------------
        per_type_no_events=per_type_no_events,
        per_type_time_oow=per_type_time_oow,
        per_type_bad_units=per_type_bad_units,
        per_type_miss_oberrs=per_type_miss_oberrs,
        per_type_miss_adpsfc=per_type_miss_adpsfc,
    
        # ---------------------------------------------------------------------
        # Status finais
        # ---------------------------------------------------------------------
        status_struct=status_struct,              # OK / PROBLEMAS_ESTRUTURAIS
        status_diag=status_diag,                  # OK / ATENÇÃO (conforme regras definidas)
        pressure_unit_used=pressure_unit_used,    # unidade de pressão adotada (hpa/cb/pa)
    
        # ---------------------------------------------------------------------
        # Exemplos ("onde") para facilitar inspeção manual
        # Limitados por --where
        # ---------------------------------------------------------------------
        where_gsi_hdr=where_gsi_hdr,              # exemplos de header estrutural inválido
        where_gsi_lonlat=where_gsi_lonlat,        # exemplos de lon/lat inválidos
        where_gsi_no_useful=where_gsi_no_useful,  # exemplos sem variáveis úteis
        where_events=where_events,                # exemplos sem EVN
        where_time=where_time,                    # exemplos fora da janela temporal
        where_units=where_units,                  # exemplos com unidades implausíveis
    
        # ---------------------------------------------------------------------
        # Registro detalhado de falhas
        # Cada entrada representa uma falha específica em um subset.
        # Base para --report-csv
        # ---------------------------------------------------------------------
        fail_records=fail_records,
    )

    if csv_out:
        _write_csv(csv_out, path, summary)
    if report_csv:
        _write_report_csv(report_csv, path, summary)
    if vars_csv:
        _write_vars_csv(vars_csv, summary, units_by_var_and_type)

    return summary


# -----------------------------------------------------------------------------
# CSVs
# -----------------------------------------------------------------------------

def _write_csv(path: str, src: str, s: Summary) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "src_file",
            "msg_type",
            "total",
            "gsi_approved",
            "gsi_approved_pct",
            "gsi_bad_hdr",
            "gsi_bad_lonlat",
            "gsi_no_useful",
            "diag_no_events",
            "diag_time_oow",
            "diag_bad_units",
            "diag_missing_oberrs",
            "diag_missing_adpsfc_fields",
        ])
        for k in sorted(s.per_type_total.keys(), key=lambda x: (isinstance(x, bytes), x)):
            mt = _b2s(k)
            tot = s.per_type_total.get(k, 0)
            app = s.per_type_gsi_approved.get(k, 0)
            pct = (100.0 * app / tot) if tot > 0 else 0.0
            w.writerow([
                src,
                mt,
                tot,
                app,
                f"{pct:.2f}",
                s.per_type_gsi_bad_hdr.get(k, 0),
                s.per_type_gsi_bad_lonlat.get(k, 0),
                s.per_type_gsi_no_useful.get(k, 0),
                s.per_type_no_events.get(k, 0),
                s.per_type_time_oow.get(k, 0),
                s.per_type_bad_units.get(k, 0),
                s.per_type_miss_oberrs.get(k, 0),
                s.per_type_miss_adpsfc.get(k, 0),
            ])


def _write_report_csv(path: str, src: str, s: Summary) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["src_file", "msg_type", "subset_idx", "fail_type", "var"])
        for rec in s.fail_records:
            w.writerow([src, rec.msg_type, rec.subset_idx, rec.fail_type, rec.var or ""])


def _write_vars_csv(path: str, s: Summary, counter: dict) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["msg_type", "var", "bad_units_count"])
        for (mt, var), cnt in sorted(counter.items()):
            w.writerow([mt, var, cnt])


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="check_prepbufr.py",
        description=(
            "Validação estrutural mínima (GSI) + diagnósticos de eventos, tempo, "
            "unidades, erros de observação e ADPSFC em arquivos PREPBUFR."
        ),
    )
    p.add_argument("file", help="Caminho do arquivo PREPBUFR")
    p.add_argument(
        "--twind",
        type=float,
        default=DEFAULT_TWIND_H,
        metavar="HOURS",
        help=f"Meia-largura da janela |DHR| (diagnóstico, padrão: ±{DEFAULT_TWIND_H} h).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Oculta saída legível; usa apenas o código de saída.",
    )
    p.add_argument(
        "--progress",
        choices=("off", "counter", "bar"),
        default="counter",
        help="Tipo de progresso: off | counter | bar (precisa de tqdm para 'bar').",
    )
    p.add_argument(
        "--progress-every",
        type=int,
        default=DEFAULT_PROGRESS_EVERY,
        metavar="N",
        help="(counter) imprime a cada N subsets (0 desliga).",
    )
    p.add_argument(
        "--check-oberrs",
        action="store_true",
        help=(
            "Diagnóstico: verifica se existem erros de observação (TOE/WOE/POE/ZOE/QOE) "
            "quando as variáveis correspondentes estão presentes."
        ),
    )
    p.add_argument(
        "--pressure-unit",
        choices=PRESSURE_UNITS,
        default=DEFAULT_PRESSURE_UNIT,
        help="Unidade de pressão para POB/PRSS: hpa | cb | pa | auto.",
    )
    p.add_argument(
        "--where",
        type=int,
        default=0,
        metavar="N",
        help="Imprime até N exemplos por tipo de problema estrutural/diagnóstico. (0 desliga)",
    )
    p.add_argument(
        "--csv",
        dest="csv_out",
        metavar="FILE",
        help="Escreve CSV agregado por msg_type com totais, aprovados (GSI) e contadores de falhas.",
    )
    p.add_argument(
        "--report-csv",
        metavar="FILE",
        help="Escreve CSV detalhado com 1 linha por falha (estrutural + diagnóstica).",
    )
    p.add_argument(
        "--vars-csv",
        metavar="FILE",
        help="Escreve CSV com contagem de falhas de unidade (bad_units) por msg_type/var.",
    )
    p.add_argument(
        "--kind",
        choices=("adpsfc",),
        help="Ativa diagnósticos específicos para ADPSFC (PRSS/PWO/CAT).",
    )
    return p


def _fmt_types(per_type_total: Dict[object, int]) -> str:
    lines = []
    for k in sorted(per_type_total.keys(), key=lambda x: (isinstance(x, bytes), x)):
        mt = _b2s(k)
        lines.append(f" - {mt}: {per_type_total[k]} subsets")
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    summary = check_file(
        args.file,
        twind_h=args.twind,
        progress=args.progress,
        progress_every=args.progress_every,
        quiet=args.quiet,
        check_oberrs=args.check_oberrs,
        pressure_unit=args.pressure_unit,
        where_max=args.where,
        report_csv=args.report_csv,
        vars_csv=args.vars_csv,
        kind=args.kind,
        csv_out=args.csv_out,
    )

    if not args.quiet:
        print(f"Arquivo: {args.file}")
        print(f"Unidade de pressão adotada: {summary.pressure_unit_used}")
        print(f"Mensagens: {summary.nmsg:,}  |  Subsets: {summary.nsub:,}")
        if summary.per_type_total:
            print(_fmt_types(summary.per_type_total))

        # Resumo estrutural
        print("\nResumo estrutural (mínimo GSI):")
        print(f"  Subsets com header essencial faltando (XOB/YOB/DHR/TYP): {summary.gsi_bad_hdr:,}")
        print(f"  Subsets com lon/lat claramente inválidos: {summary.gsi_bad_lonlat:,}")
        print(f"  Subsets sem NENHUMA variável 'útil' (TOB/QOB/U/V/POB/ZOB/PRSS): {summary.gsi_no_useful:,}")

        # Resumo diagnóstico
        print("\nResumo diagnóstico (não afetam leitura pelo GSI):")
        print(f"  Subsets sem eventos finais (TEVN/WEVN/PEVN/ZEVN/QEVN): {summary.diag_no_events:,}")
        print(f"  Subsets fora da janela de tempo (±{args.twind:.1f}h): {summary.diag_time_oow:,}")
        print(f"  Subsets com unidades/valores implausíveis: {summary.diag_bad_units:,}")
        if args.check_oberrs:
            print(f"  Subsets faltando erros de observação: {summary.diag_miss_oberrs:,}")
        if args.kind == "adpsfc":
            print(f"  Subsets ADPSFC com problemas em PRSS/PWO/CAT: {summary.diag_miss_adpsfc:,}")

        # Onde (estrutural)
        if args.where and summary.where_gsi_hdr:
            print("\nOnde (header essencial faltando):")
            print("\n".join(summary.where_gsi_hdr))
        if args.where and summary.where_gsi_lonlat:
            print("\nOnde (lon/lat claramente inválidos):")
            print("\n".join(summary.where_gsi_lonlat))
        if args.where and summary.where_gsi_no_useful:
            print("\nOnde (sem variáveis úteis TOB/QOB/U/V/POB/ZOB/PRSS):")
            print("\n".join(summary.where_gsi_no_useful))

        # Onde (diagnóstico)
        if args.where and summary.where_events:
            print("\nOnde (sem eventos finais):")
            print("\n".join(summary.where_events))
        if args.where and summary.where_time:
            print("\nOnde (fora da janela de tempo):")
            print("\n".join(summary.where_time))
        if args.where and summary.where_units:
            print("\nOnde (unidades/valores implausíveis):")
            print("\n".join(summary.where_units))

        if args.report_csv:
            print(f"\nCSV detalhado de falhas: {args.report_csv}")
        if args.vars_csv:
            print(f"CSV agregado (falhas de unidade por variável): {args.vars_csv}")
        if args.csv_out:
            print(f"CSV (por msg_type) escrito em: {args.csv_out}")

        print(f"\nStatus estrutural (GSI): {summary.status_struct}")
        print(f"Status diagnóstico: {summary.status_diag}")

    # código de saída = 0 se estrutural OK (mínimo GSI), 1 se há problemas estruturais
    return 0 if summary.status_struct == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())


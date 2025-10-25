#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
write_prepbufr_min.py
=====================

Gerador **mínimo** de um arquivo **PREPBUFR** com 1 (um) subset, usando a
camada `bufrlib` exposta pelo pacote `nceplibs-bufr`.

O arquivo de saída é aberto no modo **OUT** com a **DX (tabela)** do PREPBUFR
embutida, exatamente como se espera para consumo por decodificadores do GSI.

Este exemplo:
- cria uma **mensagem** do tipo `ADPUPA` em `2025-01-01 00Z`;
- escreve um **subset** com o **cabeçalho mínimo**: `XOB YOB DHR TYP ELV SAID T29`;
- grava **1 evento de temperatura** (`TEVN = [TOB, TQM, TPC, TRC]`);
- finaliza (`WRITSB` subset, `CLOSMG` mensagem, `CLOsBF` arquivo).

> **Importante**: é um exemplo didático. Em produção, ajuste:
> - **DX** (use a `prepobs_prep.bufrtable` correta do seu ambiente);
> - `TYP/T29` coerentes com a sua rede/tipo de plataforma;
> - `QM/PC/RC` de acordo com sua lógica de QC;
> - conversões de unidade conforme sua DX (graus, K, m/s, hPa, etc.).

Requisitos
----------
- Python 3.8+
- `nceplibs-bufr` (conda-forge) — provê `import ncepbufr` e `ncepbufr.bufrlib`
- `numpy`

Instalação rápida:
    conda install -c conda-forge nceplibs-bufr numpy

Execução:
    python write_prepbufr_min.py

Saída esperada:
    - Cria `toy.prepbufr` contendo 1 subset ADPUPA com um evento TEVN.

Limitações
----------
- Não escreve `SID` (IA5) para manter o exemplo mínimo.
- Não grava eventos de vento (`WEVN`) nem outros (`PEVN/ZEVN/QEVN`).
- Não escreve erros de observação (`TOE/WOE/...`).

"""

from __future__ import annotations

import ncepbufr
import numpy as np


def open_prepbufr(out_fn: str, dx_table: str, lunit_out: int = 11, lunit_tab: int = 12):
    """Abre o arquivo PREPBUFR para **escrita** e embute a **DX** do PREPBUFR.

    Passos executados (BUFRLIB):
    1) `datelen(10)` → datas no formato `YYYYMMDDHH`
    2) `open(dx_table, LUNIT_TAB)` → abre a tabela DX
    3) `open(out_fn, LUNIT_OUT, mode='wb')` → cria o arquivo BUFR
    4) `openbf(LUNIT_OUT, "OUT", LUNIT_TAB)` → conecta e **embute** a DX

    Parameters
    ----------
    out_fn : str
        Caminho do arquivo PREPBUFR a ser criado (ex.: "toy.prepbufr").
    dx_table : str
        Caminho para a tabela DX do PREPBUFR (ex.: "prepobs_prep.bufrtable").
    lunit_out : int, optional
        Unidade lógica para o arquivo de saída (default 11).
    lunit_tab : int, optional
        Unidade lógica para a tabela DX (default 12).

    Returns
    -------
    tuple[ncepbufr.bufrlib, int, int]
        A referência à `bufrlib` e as unidades lógicas `(LUNIT_OUT, LUNIT_TAB)`.
    """
    blib = ncepbufr.bufrlib
    blib.datelen(10)  # YYYYMMDDHH
    blib.open(dx_table.encode(), lunit_tab)
    blib.open(out_fn.encode(), lunit_out, mode="wb")
    blib.openbf(lunit_out, b"OUT", lunit_tab)  # embute a DX no arquivo (estilo PREPBUFR)
    return blib, lunit_out, lunit_tab


def write_header_min(blib, lunit_out: int):
    """Escreve o **cabeçalho mínimo** de um subset PREPBUFR via `UFBINT`.

    Campos numéricos gravados:
    - `XOB` (lon, graus), `YOB` (lat, graus)
    - `DHR` (horas relativas ao ciclo)
    - `TYP` (código de tipo de relatório)
    - `ELV` (elevação em metros)
    - `SAID` (sensor/platform id, numérico)
    - `T29` (indicador interno; use 0 se não souber)

    Retorna o número de campos gravados (valor de retorno de `ufbint`).

    Notes
    -----
    - Exemplo didático: usa valores fixos para uma estação fictícia.
    - Ajuste conforme seu caso de uso (inclusive `TYP/T29`).
    """
    hdr_names = b"XOB YOB DHR TYP ELV SAID T29"
    hdr = np.zeros((7, 1), dtype=np.float64)
    hdr[0, 0] = -46.625  # XOB  (lon)
    hdr[1, 0] = -23.550  # YOB  (lat)
    hdr[2, 0] = 0.00     # DHR  (no centro do ciclo)
    hdr[3, 0] = 120.0    # TYP  (ex.: radiossonda; ajuste ao seu código)
    hdr[4, 0] = 10.0     # ELV  (m)
    hdr[5, 0] = 1.0      # SAID (id do sensor/plataforma)
    hdr[6, 0] = 0.0      # T29  (se não usar, mantenha 0)
    return blib.ufbint(lunit_out, hdr, hdr.shape[0], hdr.shape[1], hdr_names)


def write_tevent_min(blib, lunit_out: int):
    """Escreve **um evento de temperatura** (`TEVN`) via `UFBSEQ`.

    Sequência `TEVN` típica (4 linhas, 1 evento/coluna):
      - linha 1 → `TOB` (temperatura observada, K)
      - linha 2 → `TQM` (quality mark)
      - linha 3 → `TPC` (program code)
      - linha 4 → `TRC` (reason code)

    Para simplicidade, este exemplo usa:
      - `TOB = 290.15 K`
      - `TQM = 0` (neutro)
      - `TPC = 0`, `TRC = 0`

    Returns
    -------
    int
        Valor de retorno de `ufbseq`.
    """
    te_seq = b"TEVN"
    te = np.zeros((4, 1), dtype=np.float64)
    te[0, 0] = 290.15  # TOB (K)
    te[1, 0] = 0.0     # TQM
    te[2, 0] = 0.0     # TPC
    te[3, 0] = 0.0     # TRC
    return blib.ufbseq(lunit_out, te, te.shape[0], te.shape[1], te_seq)


def main() -> None:
    """Gera um PREPBUFR mínimo com 1 mensagem (ADPUPA) e 1 subset (TEVN)."""
    # ---------------------------------------------------------------------
    # Configuração do exemplo
    out_fn = "toy.prepbufr"
    dx_table = "prepobs_prep.bufrtable"  # ajuste para o caminho válido na sua máquina
    msgtype = b"ADPUPA"
    idate = 2025010100  # YYYYMMDDHH
    # ---------------------------------------------------------------------

    # 1) Abre saída e DX
    blib, LUNIT_OUT, LUNIT_TAB = open_prepbufr(out_fn, dx_table)

    # 2) Abre mensagem PREPBUFR (tipo/tempo)
    blib.openmb(LUNIT_OUT, msgtype, idate)

    # 3) Cabeçalho mínimo (numérico). Evito SID (IA5) neste exemplo.
    _ = write_header_min(blib, LUNIT_OUT)

    # 4) Evento TEVN (temperatura)
    _ = write_tevent_min(blib, LUNIT_OUT)

    # 5) Grava o subset e fecha mensagem/arquivo
    blib.writsb(LUNIT_OUT)   # escreve este subset
    blib.closmg(LUNIT_OUT)   # fecha a mensagem
    blib.closbf(LUNIT_OUT)   # fecha o arquivo BUFR

    print(f"Escrevi {out_fn} com 1 subset {msgtype.decode()}.")


if __name__ == "__main__":
    main()


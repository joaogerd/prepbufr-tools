# PREPBUFR Tools (CPTEC)

Ferramentas em Python para **ler**, **gerar** e **validar** arquivos **PREPBUFR** compatíveis com o **GSI** (NCEP/NCEPLIBS). Feitas para pipelines de assimilação no CPTEC, mas úteis para qualquer fluxo BUFR→PREPBUFR.

## Conteúdo

* `prepbufr_cptec_gsi.py` — **Gerador** de PREPBUFR a partir de 1+ BUFRs de entrada, montando **event stacks** (`TEVN/WEVN/PEVN/ZEVN/QEVN`) e cabeçalhos no formato que o GSI espera.
* `check_prepbufr.py` — **Validador** (smoke test) estrutural + plausibilidade: header mínimo, janela de tempo, presença de eventos, **unidades/intervalos**, **erros de observação**, e, para **ADPSFC**, checa `PRSS/PWO/CAT`. Gera **CSV** por `msg_type` com **% de subsets aprovados**.
* (Opcional) `read_prepbufr_min.py` — leitor mínimo para inspeções rápidas.
* (Opcional) `write_prepbufr_min.py` — escritor mínimo de PREPBUFR (toy).

> **Importante:** Estas ferramentas **não** executam QC físico/estatístico. O foco é **estrutura BUFR/PREPBUFR** e metadados/flags esperados pelo GSI.
---

## Instalação rápida

```bash
# Clone e entre no diretório
git clone https://github.com/joaogerd/prepbufr-tools.git
cd prepbufr-tools

# Crie o ambiente
mamba env create -f environment.yml

# Ative o ambiente
conda activate prepbufr-tools

# Execute qualquer ferramenta
python scripts/check_prepbufr.py arquivo.prepbufr

```

---

## Instalação

### Requisitos

* Python **3.8+**
* [conda-forge] `nceplibs-bufr` (fornece o módulo `ncepbufr`)
* (opcional) `tqdm` para barra de progresso

```bash
# Ambiente recomendado (conda ou mamba)
conda create -n prepbufr python=3.11 -y
conda activate prepbufr
conda install -c conda-forge nceplibs-bufr numpy -y
pip install tqdm  # opcional
```

---

## Tabelas DX

Para **gerar** PREPBUFR, use a DX do PREPBUFR (tipicamente `prepobs_prep.bufrtable`).
Coloque a tabela no repo (ou referencie seu caminho local) e passe com `--dx`.

---

## Uso Rápido

### 1) Gerar PREPBUFR (GSI-ready)

```bash
python prepbufr_cptec_gsi.py \
  --type adpsfc \
  --dx /caminho/prepobs_prep.bufrtable \
  --idate 2025091100 \
  --out saida.prepbufr \
  in1.bufr in2.bufr
```

* `--type`: `adpsfc` ou `adpupa` (ativa alguns campos/eventos por tipo)
* `--idate`: data-centro (YYYYMMDDHH) usada para calcular `DHR` relativo
* (opcional) `--typ 120 --t29 0` para forçar `TYP/T29` se necessário

### 2) Validar PREPBUFR

```bash
python check_prepbufr.py saida.prepbufr --progress bar --twind 6 \
  --check-oberrs --kind adpsfc --csv relatorio.csv
```

Saída típica:

```
Arquivo: saida.prepbufr
Mensagens: 12  |  Subsets: 18543
 - ADPSFC: 18543 subsets

Resumo:
  Header(s) incompleto(s): 0
  Subsets sem eventos finais (TEVN/WEVN/etc.): 0
  Subsets fora da janela de tempo (±6.0h): 12
  Subsets com unidades/valores implausíveis: 7
  Subsets faltando erros de observação: 5
  Subsets ADPSFC sem PRSS/PWO/CAT (ou fora de faixa): 0

CSV escrito em: relatorio.csv

Status: ATENÇÃO
```

---

## Regras de Validação (checker)

### Estruturais

* **Header mínimo** por subset: `SID (opcional) XOB YOB DHR TYP ELV SAID T29`
* **Eventos**: presença de pelo menos 1 **evento final** (`TEVN`, `WEVN`, ou `PEVN/ZEVN/QEVN`), checando `TQM`/`WQM` quando aplicável
* **Janela de tempo**: `|DHR| ≤ twind_h` (padrão ±3 h)

### Unidades / Intervalos (plausibilidade)

* `XOB` (°): [-180, 180], `YOB` (°): [-90, 90]
* `TOB` (K): [180, 340]
* `POB` (hPa): [1, 1100]
* `ZOB` (m): [-500, 60000]
* `UOB/VOB` (m s⁻¹): [-200, 200]
* `QOB` (kg kg⁻¹): [0, 0.1]
* (ADPSFC) `PRSS` (hPa): [100, 1100], `PWO` (~mm): [0, 100]

> Esses **não** são limites WMO; são **faixas operacionais**. Ajuste conforme o seu domínio (alta montanha, eventos severos, etc.).

### Erros de Observação (opcional `--check-oberrs`)

Exige `TOE/WOE/POE/ZOE/QOE` **se** a variável correspondente existir no subset.

### Campos ADPSFC (opcional `--kind adpsfc`)

Exige `PRSS`, `PWO`, `CAT` (e valida faixa de `PRSS/PWO`) **quando** header e eventos estiverem OK.

### Aprovado x Reprovado

Um subset é **aprovado** quando **todas** as checagens aplicáveis passam.
O CSV traz **% aprovado** por `msg_type` para varrer lotes grandes.

---

## CSV de Saída

Colunas:

* `msg_type` — tipo de mensagem BUFR (ex.: `ADPSFC`, `ADPUPA`)
* `total`, `approved`, `approved_pct`
* `bad_hdr`, `no_events`, `time_oow`, `bad_units`
* `missing_oberrs`, `missing_adpsfc_fields`

Exemplo (trecho):

```
msg_type,total,approved,approved_pct,bad_hdr,no_events,time_oow,bad_units,missing_oberrs,missing_adpsfc_fields
ADPSFC,18543,18291,98.64,0,0,12,7,5,0
```

---

## Dicas de Uso / Boas Práticas

* **DHR relativo**: garanta que `DHR = (t_obs - idate) [h]`. Se `DHR` não existir no BUFR, calcule com `YEAR/MNTH/DAYS/HOUR/MINU/SECO`.
* **Eventos**: o GSI usa o **último evento** em `TEVN/WEVN/…`. Se você só quer “passar adiante”, escreva 1 evento final por variável com `QM/PC/RC` (0 se neutro).
* **Missing**: use o missing da BUFRLIB (`1.0e10`), nunca `-9999`.
* **TYP/T29**: mantenha coerentes com sua DX/tipo de mensagem; se necessário, use `--typ/--t29` no gerador.

---

## Solução de Problemas (Troubleshooting)

* **GSI rejeita tudo**
  → Rode `check_prepbufr.py`. Se aparecer:

  * `Header(s) incompleto(s)` > 0: você não gravou algum de `XOB YOB DHR TYP ELV SAID T29`.
  * `Subsets sem eventos finais` > 0: faltam as sequências `TEVN/WEVN/...`.
  * `Fora da janela`: `DHR` errado (absoluto, não relativo).
  * `Unidades implausíveis`: conversão/escala incorreta.
* **Sem `SID`**
  → O GSI costuma aceitar, mas é recomendável escrever `SID` (IA5). Depende da função de string disponível na sua build do `nceplibs-bufr`.
* **DX incompatível**
  → Use a `prepobs_prep.bufrtable` correta e garanta que o encoder a embute no arquivo (`openbf(..., "OUT", ...)`).

---

## Desenvolvimento

### Estilo e lint

* Python 3.8+ compatível
* Recomendações: `ruff`, `black`, `mypy` (opcional)

### Testes locais rápidos

* Gere um arquivo toy:

  ```bash
  python scripts/write_prepbufr_min.py
  python scripts/check_prepbufr.py toy.prepbufr --progress counter
  ```

---

## Contribuição

Contribuições são bem-vindas! Abra uma *issue* com:

* contexto (tipos de dados, DX usada, versão do GSI),
* amostras (quando possível),
* logs do `check_prepbufr.py` e CSV.

Para *PRs*:

1. Explique a alteração e impacto.
2. Mantenha compatibilidade Python 3.8+.
3. Evite dependências grandes.

---

## FAQ

**Isso substitui o QC do GSI?**
Não. O GSI faz o QC físico/estatístico. Aqui garantimos **estrutura** e **metadados** corretos.

**Preciso de `SID`?**
Opcional, mas recomendado. O GSI geralmente funciona sem `SID`.

**E se meu fluxo usar outros tipos (SATWND, AIRCFT, SFCSHP…)?**
O gerador já escreve `TEVN/WEVN/…` de forma genérica. Ajuste `--type` e campos extras conforme necessidade (ou abra uma issue com um exemplo).

**Como ajustar os limites de plausibilidade?**
Edite o dicionário `RANGE` no `check_prepbufr.py` para as faixas do seu domínio.

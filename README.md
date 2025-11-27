# PREPBUFR Tools (CPTEC)

Ferramentas em Python para **ler**, **gerar** e **validar** arquivos **PREPBUFR** compatíveis com o **GSI** (NCEP/NCEPLIBS). Feitas para pipelines de assimilação no CPTEC, mas úteis para qualquer fluxo BUFR→PREPBUFR.

## Conteúdo

* `prepbufr_cptec_gsi.py` — **Gerador** de PREPBUFR a partir de 1+ BUFRs de entrada, montando **event stacks** (`TEVN/WEVN/PEVN/ZEVN/QEVN`) e cabeçalhos no formato que o GSI espera.
* `scripts/check_prepbufr.py` — **Validador** (“smoke test”) focado em **compatibilidade com o GSI** + diagnósticos:
  - Checagens **estruturais mínimas** (impactam o status final):
    - header mínimo (`XOB YOB DHR TYP`);
    - (opcional) **erros de observação**, se `--check-oberrs`;
    - (opcional) campos **ADPSFC** (`PRSS/PWO/CAT`), se `--kind adpsfc`.
  - Checagens **diagnósticas** (não mudam o status, só informam):
    - janela de tempo (`|DHR| ≤ twind`);
    - presença de eventos finais (`TEVN/WEVN/PEVN/ZEVN/QEVN`);
    - **unidades / intervalos plausíveis** para variáveis-chave.
  - Gera **CSV** por `msg_type` com **% de subsets aprovados**.
* (Opcional) `read_prepbufr_min.py` — leitor mínimo para inspeções rápidas.
* (Opcional) `write_prepbufr_min.py` — escritor mínimo de PREPBUFR (toy).

> **Importante:** Estas ferramentas **não** substituem o QC físico/estatístico do GSI.  
> O foco é a **estrutura BUFR/PREPBUFR**, metadados e flags esperados pelo GSI, com diagnósticos básicos de consistência.

---

## Instalação rápida

```bash
# Clone e entre no diretório
git clone https://github.com/joaogerd/prepbufr-tools.git
cd prepbufr-tools

# Crie o ambiente (via mamba ou conda)
mamba env create -f environment.yml

# Ative o ambiente
conda activate prepbufr-tools

# Execute qualquer ferramenta
python scripts/check_prepbufr.py arquivo.prepbufr
````

---

## Instalação

### Requisitos

* Python **3.8+**
* [conda-forge] `nceplibs-bufr` (fornece o módulo `ncepbufr`)
* `numpy`
* (opcional) `tqdm` para barra de progresso

```bash
# Ambiente recomendado (conda ou mamba)
conda create -n prepbufr python=3.11 -y
conda activate prepbufr
conda install -c conda-forge nceplibs-bufr numpy -y
pip install tqdm  # opcional (para --progress bar)
```

---

## Tabelas DX

Para **gerar** PREPBUFR, use a DX do PREPBUFR (tipicamente `prepobs_prep.bufrtable`).

Coloque a tabela no repo (ou referencie seu caminho local) e passe com `--dx`:

```bash
python prepbufr_cptec_gsi.py --dx /caminho/prepobs_prep.bufrtable ...
```

---

## Uso rápido

### 1) Gerar PREPBUFR (GSI-ready)

```bash
python prepbufr_cptec_gsi.py \
  --type adpsfc \
  --dx /caminho/prepobs_prep.bufrtable \
  --idate 2025091100 \
  --out saida.prepbufr \
  in1.bufr in2.bufr
```

* `--type`: por enquanto `adpsfc` ou `adpupa` (ativa alguns campos/eventos por tipo).
* `--idate`: data-centro (YYYYMMDDHH) usada para calcular `DHR` relativo.
* (opcional) `--typ 120 --t29 0` para forçar `TYP/T29` se necessário.

### 2) Validar PREPBUFR (compatibilidade + diagnóstico)

```bash
python scripts/check_prepbufr.py saida.prepbufr \
  --progress bar --twind 6 \
  --check-oberrs --kind adpsfc \
  --csv por_tipo.csv \
  --report-csv falhas_detalhe.csv \
  --vars-csv falhas_unidades.csv
```

Saída típica:

```text
Arquivo: saida.prepbufr
Unidade de pressão adotada: hpa
Mensagens: 12  |  Subsets: 18543
 - ADPSFC: 18543 subsets

Resumo (diagnóstico + estrutural):
  Header(s) incompleto(s): 0
  Subsets sem eventos finais (TEVN/WEVN/etc.): 0
  Subsets fora da janela de tempo (±6.0h): 12
  Subsets com unidades/valores implausíveis: 7
  Subsets faltando erros de observação: 5
  Subsets ADPSFC sem PRSS/PWO/CAT (ou fora de faixa): 0

CSV detalhado de falhas: falhas_detalhe.csv
CSV agregado (falhas de unidade por variável): falhas_unidades.csv
CSV (por msg_type) escrito em: por_tipo.csv

Status: ATENÇÃO
```

> No exemplo acima, o status **ATENÇÃO** vem dos 5 subsets faltando erros de observação (`--check-oberrs`).
> As checagens de **tempo, eventos e unidades** são apenas diagnósticas.

---

## Regras de Validação (checker)

### 1. Checagens **estruturais** (afetam o status PASSOU/ATENÇÃO)

Essas são as que realmente interessam para saber se o arquivo é **compatível com o GSI**:

* **Header mínimo** por subset:

  * campos requeridos: `XOB`, `YOB`, `DHR`, `TYP`
  * se qualquer um deles estiver ausente → `bad_header`.

* **Erros de observação** (opcional com `--check-oberrs`):

  * se ligado, exige erro para cada variável presente:

    * `TOB → TOE`
    * `UOB/VOB → WOE`
    * `POB → POE`
    * `ZOB → ZOE`
    * `QOB → QOE`
  * ausência de qualquer erro observado para uma variável presente → `missing_oberrs`.

* **Campos ADPSFC** (opcional com `--kind adpsfc`):

  * quando `--kind adpsfc` é usado, subset considerado “significativo” deve ter:

    * `PRSS` (pressão à superfície) dentro da faixa,
    * `PWO` (água precipitável) dentro da faixa,
    * `CAT` (categoria).
  * falhas entram como `missing_adpsfc_fields`.

> Se **qualquer uma** dessas checagens falhar (de acordo com as opções usadas), o status final será `ATENÇÃO`.

---

### 2. Checagens **diagnósticas** (não alteram o status)

Estas servem para investigar **qualidade dos dados** e **coerência interna**, mas **não impedem** que o GSI leia o arquivo:

* **Eventos finais**:

  * procura por `TEVN` (checando `TQM`), `WEVN` (checando `WQM`) ou presença de `PEVN/ZEVN/QEVN`.
  * se nenhum evento for encontrado → `no_events`.
  * isso é **apenas informativo**, já que muitos PREPBUFR modernos não usam mais todos esses eventos.

* **Janela de tempo**:

  * checa se `|DHR| ≤ twind` (padrão `±3 h`, configurável com `--twind`).
  * subsets fora dessa janela são contabilizados como `time_out_of_window`.
  * não afeta o status (é apenas um aviso).

* **Unidades / Intervalos (plausibilidade)**:

  * checagem de ranges **operacionais** (não WMO) para valores lidos,
    depois de aplicar a unidade de pressão escolhida:

    * `XOB` (°)     : `[-180, 180]`
    * `YOB` (°)     : `[-90, 90]`
    * `TOB` (K)     : `[180, 340]`
    * `POB` (hPa)   : `[1, 1100]`
    * `ZOB` (m)     : `[-500, 40000]`
    * `UOB` (m/s)   : `[-120, 120]`
    * `VOB` (m/s)   : `[-120, 120]`
    * `QOB` (kg/kg) : `[0, 0.05]`
    * `PRSS` (hPa)  : `[500, 1085]`
    * `PWO` (mm)    : `[0, 100]`
    * `DHR` (h)     : `[-6, 6]`
    * `HOVI` (m)    : `[0, 100000]`
    * `TDO` (K)     : `[173, 333]`
    * `MSST` (K)    : `[271, 313]`
  * valores fora dessas faixas são contabilizados como `bad_units` e podem ser inspecionados em detalhe.
  * não alteram o status final, servem como “sanity check”.

> Se o seu domínio tiver condições extremas (ex.: alta montanha, convectivo severo, trópico muito úmido), esses limites podem e devem ser ajustados diretamente no dicionário `RANGE` dentro do `check_prepbufr.py`.

---

### 3. Erros de observação (`--check-oberrs`)

Quando `--check-oberrs` é usado:

* Sempre que uma variável observada existir no subset, é esperado que o **erro de observação** correspondente também exista;
* Ausência de erro para qualquer variável observada faz com que o subset **não seja aprovado**;
* Isso entra diretamente no status final (`Status: ATENÇÃO` se houver falhas).

---

### 4. Campos ADPSFC (`--kind adpsfc`)

Quando `--kind adpsfc` é usado:

* Para subsets ADPSFC “válidos” (header mínimo ok), o script exige:

  * `PRSS` presente e em faixa (`[500, 1085]` hPa após conversão);
  * `PWO` presente e em faixa (`[0, 100]`);
  * `CAT` presente.
* Subsets que não atendem esses requisitos contam em `missing_adpsfc_fields` e
  influenciam o status final.

---

### 5. Aprovado x Reprovado

Um subset é considerado **aprovado** se:

1. O header mínimo (`XOB YOB DHR TYP`) está completo; **e**
2. Se `--check-oberrs` foi usado:

   * não há falhas de erros de observação; **e**
3. Se `--kind adpsfc` foi usado:

   * não há falhas ADPSFC (`PRSS/PWO/CAT`).

As checagens de **tempo**, **eventos** e **unidades** são **diagnósticas** e **não entram** no critério de aprovação.

O CSV por `msg_type` mostra:

* `total` — número total de subsets do tipo.
* `approved` — quantos foram aprovados segundo as regras estruturais.
* `approved_pct` — `%` aprovado (útil para varrer lotes grandes).

---

## CSV de saída

### CSV agregado por tipo (`--csv`)

Colunas:

* `msg_type` — tipo de mensagem BUFR (ex.: `ADPSFC`, `ADPUPA`)
* `total`, `approved`, `approved_pct`
* `bad_hdr`, `no_events`, `time_oow`, `bad_units`
* `missing_oberrs`, `missing_adpsfc_fields`

Exemplo:

```text
msg_type,total,approved,approved_pct,bad_hdr,no_events,time_oow,bad_units,missing_oberrs,missing_adpsfc_fields
ADPSFC,18543,18291,98.64,0,0,12,7,5,0
```

### CSV detalhado (`--report-csv`)

Um registro por falha:

* `src_file`, `pressure_unit`, `msg_type`, `subset_idx`, `fail_type`, `var`

onde `fail_type` ∈ {`bad_header`, `no_final_event`, `time_out_of_window`, `bad_units`, `missing_oberrs`, `missing_adpsfc`}.

### CSV por variável (`--vars-csv`)

Contagem de falhas de unidade (`bad_units`) por `(msg_type, var)`, útil para ver, por exemplo, quantas vezes `QOB` ficou fora da faixa em `ADPUPA`.

---

## Dicas de uso / boas práticas

* **DHR relativo**: garanta que `DHR = (t_obs - idate) [h]`.
  Se `DHR` não existir no BUFR original, calcule a partir de `YEAR/MNTH/DAYS/HOUR/MINU/SECO`.

* **Eventos**:

  * O GSI usa tipicamente o **último evento** em `TEVN/WEVN/…`;
  * Muitos PREPBUFR operacionais modernos não usam mais todas as sequências de evento, então `no_events > 0` não significa, por si só, erro.

* **Missing**:

  * Use sempre o missing da BUFRLIB (`1.0e10`), nunca `-9999`.

* **TYP/T29**:

  * Mantenha coerentes com sua DX e com o `convinfo` do GSI;
  * Se necessário, use `--typ/--t29` no gerador para forçar valores.

* **Pressão**:

  * Use `--pressure-unit auto` se tiver dúvidas sobre as unidades de `POB/PRSS` no seu fluxo;
  * O script autodetecta a unidade pelo valor mediano das pressões.

---

## Solução de problemas (Troubleshooting)

* **GSI rejeita tudo / não usa quase nada**

  1. Rode:

     ```bash
     python scripts/check_prepbufr.py seu.prepbufr --check-oberrs --kind adpsfc --where 5
     ```

  2. Interprete:

     * `Header(s) incompleto(s) > 0`
       → verifique se está escrevendo `XOB`, `YOB`, `DHR`, `TYP` corretamente.

     * `Subsets faltando erros de observação` (com `--check-oberrs`)
       → o arquivo não tem `TOE/WOE/POE/ZOE/QOE` para variáveis presentes; o GSI tende a ignorar essas observações.

     * `Subsets ADPSFC sem PRSS/PWO/CAT` (com `--kind adpsfc`)
       → faltam campos essenciais para superfície.

     * `Subsets com unidades/valores implausíveis`
       → suspeita forte de unidade errada (ex.: pressão em Pa em vez de hPa), mas isso é diagnóstico: ajuste DX ou a forma como você escreve as variáveis.

     * `Subsets fora da janela de tempo` ou `sem eventos finais`
       → são diagnósticos; avalie se condizem com sua configuração do GSI (`twind`, uso de eventos) mas não são, por si só, “erros duros”.

---

## Desenvolvimento

### Estilo e lint

* Compatível com Python 3.8+.
* Recomenda-se: `ruff`, `black`, `mypy` (opcional).

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
* logs do `check_prepbufr.py` e CSVs gerados.

Para *PRs*:

1. Explique a alteração e o impacto esperado.
2. Mantenha compatibilidade com Python 3.8+.
3. Evite dependências grandes e irrelevantes para HPC/assimilação.

---

## FAQ

**Isso substitui o QC do GSI?**
Não. O GSI continua responsável pelo QC físico/estatístico.
Aqui garantimos a **estrutura** e fornecemos diagnósticos que ajudam a encontrar erros de DX, unidades e metadados antes de rodar o GSI.

**Preciso de `SID`?**
Não é obrigatório para o validador estrutural, mas é recomendado.
O GSI geralmente funciona sem `SID`, dependendo da configuração.

**E se meu fluxo usar outros tipos (SATWND, AIRCFT, SFCSHP…)?**
O gerador pode ser estendido para outros tipos.
O `check_prepbufr.py` já contabiliza todos os `msg_type` que encontrar e aplica as checagens estruturais mínimas (header / oberrs / ADPSFC).

**Como ajustar os limites de plausibilidade?**
Edite o dicionário `RANGE` dentro do `scripts/check_prepbufr.py` para refletir as faixas adequadas ao seu domínio (trópico úmido, alta montanha, etc.).



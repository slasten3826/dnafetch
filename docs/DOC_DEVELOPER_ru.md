# DNAfetch v3 — Документация для разработчиков

Архитектура, структуры данных, API, расширение.

---

## Обзор

DNAfetch — CLI-инструмент из 4 файлов на чистом Python 3.10+ без внешних зависимостей. Три слоя обработки: FETCH → TRANSLATE → DESCRIBE. Вход — FASTA, выход — .dpl (текстовый формат).

```
dnafetch.py      CLI, argparse, три команды (single/batch/compare)
fetch.py         Парсинг FASTA, извлечение CDS, трансляция в аминокислоты
translate.py     Маппинг нуклеотидов в операторы, вычисление метрик
describe.py      Генерация .dpl, классификация по порогам
```

---

## Структуры данных

### Gene (fetch.py)

```python
@dataclass
class Gene:
    name: str              # Имя гена (GRIN2A, INS, TP53...)
    organism: str          # Организм (Homo sapiens)
    transcript_id: str     # Accession (NM_000833.5)
    sequence: str          # Полная последовательность из FASTA
    cds: str               # Только CDS (от ATG до стоп-кодона включительно)
    codons: List[str]      # Список троек: ['ATG', 'GCC', 'TGA']
    amino_acids: List[str] # Список аминокислот: ['Met', 'Ala', '*']
    
    # Вычисляемые:
    cds_length: int        # len(cds)
    codon_count: int       # len(codons)
    protein_length: int    # codon_count - 1 (без стоп-кодона)
    stop_codon: str        # Последний кодон (TAA/TAG/TGA)
```

### PositionProfile (translate.py)

```python
@dataclass
class PositionProfile:
    frequencies: Dict[str, float]  # {'FLOW': 0.25, 'OBSERVE': 0.22, ...}
    dominant: str                  # 'FLOW'
    dominant_pct: float            # 0.25
    entropy: float                 # Нормированная энтропия [0, 1]
```

Энтропия: `H = -Σ(p * log2(p)) / log2(4)`. Нормирована на [0, 1]. Значение 1.0 — полностью равномерное распределение. 0.0 — один оператор занимает 100%.

### ComposePair (translate.py)

```python
@dataclass
class ComposePair:
    label: str        # 'FLOW→FLOW'
    op1: str          # 'FLOW'
    op2: str          # 'FLOW'
    frequency: float  # 0.117
```

COMPOSE — переход нуклеотида позиции 1 к нуклеотиду позиции 2 внутри кодона. 16 возможных пар (4×4). Частоты нормированы на 1.0.

### TranslationResult (translate.py)

```python
@dataclass
class TranslationResult:
    gene: Gene
    
    # Профили по позициям
    primary: PositionProfile     # Позиция 1
    modifier: PositionProfile    # Позиция 2
    wobble: PositionProfile      # Позиция 3
    overall: PositionProfile     # Все позиции вместе
    
    # COMPOSE
    compose_pairs: List[ComposePair]  # Отсортированы по frequency desc
    
    # Метрики
    gc_content: float            # Общий GC-контент CDS
    wobble_gc: float             # GC только в позиции 3
    complementary_balance: float # min(A,T)/max(A,T) * min(G,C)/max(G,C)
    top2_gap: float              # |freq_top1 - freq_top2| (overall)
    codon_diversity: float       # Уникальные кодоны / ожидаемые
    
    # Аминокислоты
    aa_distribution: Dict[str, float]  # Отсортирован по freq desc
    aa_std: float                      # std dev частот
    hydrophobic_ratio: float           # Доля гидрофобных AA
    charge_profile: Dict[str, float]   # {'positive': 0.13, 'negative': 0.11, 'neutral': 0.76}
```

---

## Маппинг

```python
# translate.py
NUC_TO_OP = {'A': 'FLOW', 'T': 'OBSERVE', 'G': 'LOGIC', 'C': 'CHOOSE'}
OP_TO_NUC = {'FLOW': 'A', 'OBSERVE': 'T', 'LOGIC': 'G', 'CHOOSE': 'C'}
```

Маппинг фиксированный, хардкожен. Нет возможности переопределения через конфигурацию (by design — иначе .dpl файлы от разных конфигураций несовместимы).

---

## API

### fetch.py

```python
def fetch_single(filepath: str) -> Gene
```
Парсит один FASTA файл. Если последовательность не начинается с ATG, ищет самый длинный ORF. Валидирует длину (кратность 3), стоп-кодон.

```python
def fetch_batch(filepath: str, limit: int = 0) -> Iterator[Gene]
```
Генератор. Читает multi-FASTA (поддерживает .gz). `limit=0` — без ограничений.

```python
def parse_gene(name: str, organism: str, transcript_id: str, sequence: str) -> Gene
```
Низкоуровневая функция. Принимает уже распарсенные данные.

### translate.py

```python
def translate_gene(gene: Gene) -> TranslationResult
```
Единственная публичная функция. Принимает Gene, возвращает TranslationResult со всеми метриками.

### describe.py

```python
def generate_dpl(result: TranslationResult) -> str
```
Генерирует текст .dpl файла. Возвращает строку.

```python
def write_dpl(result: TranslationResult, output_path: str) -> None
```
Генерирует и записывает .dpl в файл.

---

## Формат .dpl

Текстовый файл. Секции в квадратных скобках. Ключ-значение через ` = `. Комментарии через `#`. Списки в `[ ]`.

```
# Комментарий
[SECTION_NAME]
key = value
key = 0.1234
list_key = [
  item1: 0.1234  # inline comment
  item2: 0.5678
]
```

Секции: IDENTITY, FLOW, CONNECT, DISSOLVE, ENCODE, CHOOSE, OBSERVE, CYCLE, LOGIC, RUNTIME, MANIFEST, SIGNATURE. Порядок фиксированный.

### Парсинг .dpl

Для чтения .dpl достаточно простого парсера:

```python
def parse_dpl(content: str) -> dict:
    """Парсит .dpl в словарь секция → {ключ: значение}."""
    sections = {}
    current = None
    for line in content.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('[') and line.endswith(']'):
            current = line[1:-1]
            sections[current] = {}
            continue
        if '=' in line and current:
            key, _, val = line.partition('=')
            key = key.strip()
            val = val.strip()
            if '#' in val:
                val = val[:val.index('#')].strip()
            try:
                sections[current][key] = float(val)
            except ValueError:
                sections[current][key] = val
    return sections
```

---

## CLI

```bash
# Один ген → один .dpl
python dnafetch.py single gene.fasta
python dnafetch.py single gene.fasta -o output/gene.dpl
python dnafetch.py single gene.fasta --stdout

# Batch обработка (поддержка .gz)
python dnafetch.py batch gencode.fa.gz -o output/
python dnafetch.py batch gencode.fa.gz -o output/ -l 100   # лимит
python dnafetch.py batch gencode.fa.gz -o output/ -v        # verbose ошибки

# Сравнение двух .dpl
python dnafetch.py compare gene1.dpl gene2.dpl
```

### batch

Прогресс каждые 500 генов. Ошибочные гены пропускаются (невалидный CDS, нет ATG). С `-v` ошибки выводятся в stderr.

### compare

Выводит таблицу: числовые метрики с дельтами, строковые классификации с флагом совпадения (= / ≠).

---

## Константы и пороги

### Генетический код (fetch.py)

Стандартный генетический код (NCBI transl_table=1). 64 кодона → 20 аминокислот + стоп. Хардкожен в словаре `CODON_TABLE`.

### Свойства аминокислот (fetch.py)

```python
AA_PROPERTIES = {
    'Ala': {'class': 'nonpolar', 'hydrophobic': True, 'charge': 'neutral'},
    'Arg': {'class': 'positive', 'hydrophobic': False, 'charge': 'positive'},
    ...
}
```

20 аминокислот. Классификация: nonpolar, polar, positive, negative, aromatic, special.

### Пороги классификации (describe.py)

Все пороги эмпирические, калиброваны по контрольным точкам (Ebola, COVID-19, Pfizer BNT162b2, 6 генов человека).

**Encoding efficiency** (по wobble GC):
```
> 0.75 → engineered
> 0.60 → optimized
> 0.45 → moderate
> 0.30 → natural
≤ 0.30 → AT_rich
```

**Observation level** (по top2 gap):
```
< 0.005 → intense
< 0.02  → strong
< 0.05  → deep
< 0.08  → standard
< 0.12  → mild
≥ 0.12  → basic
```

**Cycle stability** (по top COMPOSE frequency):
```
> 0.12 → dominant
> 0.09 → strong
> 0.07 → moderate
≤ 0.07 → distributed
```

**Dissolution potential** (комбинация wobble entropy и complementary balance):
```python
score = wobble_entropy * 0.6 + (1.0 - comp_balance) * 0.4
> 0.5  → high
> 0.3  → moderate
> 0.15 → low
≤ 0.15 → crystallized
```

**Manifestation class** (по аминокислотному составу):
```python
if Gly > 15% and Pro > 10%:         structural
elif Cys > 4% and length < 200:     signaling_peptide
elif Leu > 12%:                      transmembrane
elif Ser > 10%:                      regulatory
else:                                enzymatic
```

---

## Расширение

### Добавление нового поля в .dpl

1. Вычислить метрику в `translate.py` (добавить поле в `TranslationResult`)
2. Добавить вывод в `describe.py` в нужную секцию
3. Если нужно для compare — добавить ключ в `compare_keys` или `str_keys` в `dnafetch.py`

### Добавление новой секции

Не рекомендуется. 10 секций = 10 команд ProcessLang. Добавление секции требует добавления команды в спецификацию ProcessLang.

### Добавление нового формата входных данных

Добавить парсер в `fetch.py`. Интерфейс: функция должна возвращать `Gene`. Пример для GenBank:

```python
def fetch_genbank(filepath: str) -> Gene:
    # Парсить GenBank, извлечь CDS features
    # Собрать sequence из CDS join()
    # Вернуть Gene(name=..., cds=..., ...)
    pass
```

### Изменение маппинга

Не делай этого. Все существующие .dpl файлы станут несовместимыми. Маппинг A→FLOW, T→OBSERVE, G→LOGIC, C→CHOOSE — часть спецификации DNAfetch. Если нужен другой маппинг — это другой инструмент.

---

## Обработка ошибок

**fetch.py:**
- Файл не найден → FileNotFoundError
- Нет ATG → ищет ORF; если ORF не найден → ValueError
- Длина CDS не кратна 3 → обрезает до ближайшей кратной
- Пустая последовательность → ValueError

**translate.py:**
- Неизвестный нуклеотид (N, R, Y) → пропускается при подсчёте частот
- Пустой список кодонов → ZeroDivisionError (не обработан, ожидается валидный вход от fetch)

**describe.py:**
- Пустые compose_pairs → `stability = none`, no top_compose line

**dnafetch.py batch:**
- Ошибка в отдельном гене → пропускается, счётчик errors++
- С `-v` ошибка выводится в stderr

---

## Производительность

На GENCODE v47 (~20,000 генов, 46 МБ .gz):
- Время: зависит от машины, порядка нескольких минут
- Память: O(1) для batch (генератор, один ген в памяти)
- Диск: ~20,000 .dpl файлов, каждый ~2-3 КБ = ~50-60 МБ

Узкое место — I/O (запись файлов), не вычисления.

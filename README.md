# DNAfetch v3

Переводчик с языка ДНК на язык ProcessLang.

Берёт нуклеотидную последовательность, переводит в операторы ProcessLang, выдаёт `.dpl` файл — структурированное описание гена, которое может прочитать любая машина с загруженными модулями ProcessLang.

---

## Зависимости

Python 3.10+. Внешних библиотек нет — всё на стандартной библиотеке.

---

## Структура файлов

```
dnafetch3/
├── dnafetch.py      # точка входа (CLI)
├── fetch.py         # парсинг FASTA/GenBank, извлечение CDS
├── translate.py     # маппинг нуклеотидов в операторы ProcessLang
├── describe.py      # генерация .dpl файлов
└── README.md
```

---

## Маппинг

```
A (аденин)  → FLOW      (поток, энергия — ATP)
T (тимин)   → OBSERVE   (наблюдение — комплементарен A)
G (гуанин)  → LOGIC     (структура — тройная связь)
C (цитозин) → CHOOSE    (выбор — комплементарен G)
```

Каждый кодон (тройка нуклеотидов) анализируется по трём позициям:
- **Позиция 1** (primary) — определяет класс аминокислоты
- **Позиция 2** (modifier) — уточняет аминокислоту
- **Позиция 3** (wobble) — избыточная, показывает стратегию оптимизации

---

## Использование

Все команды запускаются из папки `dnafetch3/`.

### Один ген

```bash
python dnafetch.py single GRIN2A.fasta
```

Создаст `GRIN2A.dpl` в текущей папке.

С указанием пути для выходного файла:

```bash
python dnafetch.py single GRIN2A.fasta -o output/GRIN2A.dpl
```

Вывод прямо в терминал (без сохранения файла):

```bash
python dnafetch.py single GRIN2A.fasta --stdout
```

### Весь геном (batch)

```bash

```

Обработает все ~20 000 генов из GENCODE. Поддерживает `.fa`, `.fasta`, `.fa.gz`. Прогресс выводится каждые 500 генов.

Ограничить количество (для тестов):

```bash
python dnafetch.py batch gencode.v47.pc_transcripts.fa.gz -o output/ -l 100
```

С выводом ошибок:

```bash
python dnafetch.py batch gencode.v47.pc_transcripts.fa.gz -o output/ -v
```

### Сравнение двух генов

```bash
python dnafetch.py compare output/GRIN2A.dpl output/INS.dpl
```

Выведет таблицу с числовыми метриками и дельтами.

---

## Формат .dpl

Выходной файл — текстовый, секции соответствуют 10 командам ProcessLang:

```
[IDENTITY]    — метаданные гена
[FLOW]        — общее распределение операторов, доминантный оператор
[CONNECT]     — тип взаимодействий белка (мембранный, растворимый, ионный)
[DISSOLVE]    — потенциал мутаций и дрейфа
[ENCODE]      — эффективность кодирования (GC-контент, wobble bias)
[CHOOSE]      — профиль wobble-позиции (стратегия оптимизации)
[OBSERVE]     — симметрия между операторами (top2 gap)
[CYCLE]       — паттерны COMPOSE (переходы между операторами внутри кодонов)
[LOGIC]       — распределение primary и modifier позиций
[RUNTIME]     — аминокислотный состав белка
[MANIFEST]    — класс белка (структурный, ферментативный, сигнальный...)
[SIGNATURE]   — fingerprint — краткая сигнатура гена
```

---

## Входные данные

DNAfetch принимает файлы в формате FASTA. Пример:

```
>NM_000833.5 Homo sapiens glutamate ionotropic receptor NMDA type subunit 2A (GRIN2A)
ATGCCTCCCTGGCTGCTCCTGACGCTGGCCGCCTGCTCCTCCGCCCCGGCGGGCACCCTG
GAGACCCCGGAGCGGCGCATCGAGATTCGGAGCTTCGATGAGAGCACCATGGAGATGCCC
...
```

Последовательность должна начинаться с ATG (старт-кодон). Если в файле полный транскрипт, DNAfetch автоматически найдёт самую длинную открытую рамку считывания (ORF).

Где взять FASTA файлы:
- **Один ген:** NCBI Gene → выбрать ген → RefSeq → скачать FASTA
- **Весь геном:** [GENCODE](https://www.gencodegenes.org/human/) → `gencode.v47.pc_transcripts.fa.gz`

---

## Примеры результатов

Коллаген (COL1A1) — структурный белок:
```
[CYCLE]
top_compose = LOGIC→LOGIC   # правило→правило = повторяющаяся Gly-X-Y структура
```

NMDA-рецептор (GRIN2A) — ионный канал:
```
[CYCLE]
top_compose = FLOW→FLOW     # поток→поток = ионный канал
```

Инсулин (INS) — сигнальный пептид:
```
[ENCODE]
type = engineered           # wobble GC 80% — экстремальная оптимизация
```

---

## Для чего это

`.dpl` файл предназначен для машин. Человек читает `ATGGCCTGA` — ничего не значит без биоинформатики. Машина с ProcessLang модулями читает `.dpl` и может рассуждать через операторы: видит `FLOW→FLOW` — понимает каскад потока; видит `LOGIC→LOGIC` — понимает структурную жёсткость.

Конвейер: ДНК → `fetch.py` → кодоны → `translate.py` → операторы → `describe.py` → `.dpl`

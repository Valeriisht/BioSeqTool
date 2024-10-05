# NASeqTool

- [**Русская версия**](#Биоинформатическая)
- [**English version**]

Биоинформатическая утилита NASeqTool для обработки последовательностей ДНК и РНК.

Репозиторий состоит из главного скрипта и дополнительных модулей, которые используются для обработки последовательностей ДНК/РНК и обработки данных секвенирования. 
Скрипты написаны на языке программирования Python. 

Структура репозитория имеет следующий вид: корневая папка с главным скриптом NASeqTool, где объявлены основные функции, и подпапка functions с импортируемыми модулями.
Утилита была разработана для работы на версии Python3.10 или новее.

### Установка

Клонируйте репозиторий на локальное устройство.

```sh
git clone git@github.com:Valeriisht/NASeqTool.git
```

### Содержание 

- #### NASeqTool
    - [**run_dna_rna_tools**](#run_dna_rna_tools)
    - [**filter_fastq_tool**](#filter_fastq_tool)


### run_dna_rna_tools

Утилита для манипуляции с молекулами ДНК/РНК. При исполнении скрипта подгружается модуль run_dna_rna_tool_module.py.

#### Использование
Вам понадобятся:

Последовательности ДНК/РНК, организованные в список (или в строку при выполнении операции над одной последовательностью).

#### Утилита может выполнять следующие операции: 

- #### transcribe 

Смысловая последовательность ДНК переводится в молекулу мРНК. Если на вход передана РНК, возвращается последовательность без изменений.

- #### reverse

Возвращает комплементарную цепь ДНК в направлении от 3' к 5'. Для молекулы РНК возвращается кДНК.

- #### reverse_complement 

Возвращает комплементарную цепь ДНК/кДНК в направлении от 5' к 3'.

- #### gc_content

Подсчитывается GC-состав последовательности и возвращается процент GC-нуклеотидов.

- #### is_na_sequence

Проверяется, является ли последовательность ДНК или РНК.

- #### protein_coding_sequence

Определяется длина кодирующей последовательности CDS.

### filter_fastq_tool

Утилита для обработки данных секвенирования в формате fastq. При исполнении скрипта подгружается модуль filter_fastq_module.py.

#### Использование
Вам понадобятся:

Последовательности ридов, организованные в словарь, где в качестве ключа выступают id последовательности, а в качестве значения кортеж с самой последовательностью и качеством последовательности в шкале phred33.

Пример структуры fastq файла:

![img.png](img.png)

Утилита выполняет следующие действия:

- Фильтрация последовательности относительно gc_состава. Подсчитывается gc- состав и определяется, находится ли значение в пределах пороговых.
По умолчанию пороговые значения определены от 0 до 100.

- Фильтрация последовательности относительно длины. Проверяется, находится ли длина рида в пределах заданных пороговых значений.
По умолчанию пороговые значения определены от 0 до 2**32.

- Фильтрация последовательности относительно качества рида. Проверяется, соответствует ли качество рида по шкале phred33 заданным пороговому значению. Качество рида переводится из шкалы Phred33 в ASCII.
По умолчанию пороговое значение равняется 0.

Результатом работы программы является словарь, содержащий название и нуклеотидную последовательность ридов, прошедших каждую из трех проверок.




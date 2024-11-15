- [**English version**](#NASeqTool)
- [**Русская версия**](#NASeqTool_)
  
# NASeqTool

<img align=right src="https://github.com/user-attachments/assets/70983c21-48f4-41d0-a7dc-9d5b63982eff" alt="#NASeqTool" width="100"/>

**NASeqTool** is a bioinformatics tool  for analyzing DNA and RNA sequences.

The repository consists of a main script and additional modules that are used for DNA/RNA sequence processing and sequencing data processing. The scripts are written in the Python programming language.

## Installation

Clone the repository to a local device.

```sh
git clone git@github.com:Valeriisht/NASeqTool.git
```

## Content

The repository structure is structured as the main folder with the master NASeqTool script, where the basic functions are declared, and the subfolder functions with imported modules.

<img width="400" alt="Rep Structure" src="https://github.com/user-attachments/assets/710383e6-8bbc-4b31-bc2e-e966fca2a58b">




The utility was designed to run on Python version 3.10 or newer.

- ### NASeqTool
    - [**run_dna_rna_tools**](#run_dna_rna_tools)
    - [**filter_fastq_tool**](#filter_fastq_tool)
    - [**bio_files_processor**](#bio_files_processor)

## run_dna_rna_tools

A utility for manipulating DNA/RNA molecules. The run_dna_rna_tool_module.py module is loaded when the script is executed. It is located in the NASeqTool.py script.



### Usage
You will need:

- DNA/RNA sequences are list (or string when performing an operation on a single sequence), which are required parameters for function calls.



### The utility can perform the following operations: 

- ### Function ```transcribe```.

The meaningful DNA sequence is translated into an mRNA molecule. If RNA is passed as input, the sequence is returned unchanged.

- ### Function ```reverse```

Returns the complementary DNA strand in the 3' to 5' direction. For an RNA molecule, cDNA is returned.

- ### Function ```reverse_complement```

Returns the complementary DNA/cDNA strand in the 5' to 3' direction.

- ### Function ```gc_content```

Counts the GC composition of the sequence and returns the percentage of GC nucleotides.

- ### Function ```is_na_sequence```.

Checks whether the sequence is DNA or RNA.

- ### Function ```length_cds```

Determines the length of the CDS coding sequence.

## filter_fastq_tool

The utility for analyzing sequencing data in fastq format. The filter_fastq_module.py module is loaded when the script is executed.It is located in the NASeqTool.py script.

### Usage

You will need:

- File in fastq format.

Example fastq file structure:

![image](https://github.com/user-attachments/assets/4c588d9b-6591-4adf-9ac7-85be25844346)

### The utility performs the steps below:

- Filtering the sequence with regard to gc_composition. The gc- composition is counted and determines if the value is within the thresholds.
By default, the threshold values are defined from 0 to 100.

- Filtering the sequence with regard to length. Checks whether the length of the sequence is within the specified threshold values.
By default, the threshold values are defined from 0 to 2**32.

- Sequence filtering with respect to sequence quality. Checks whether the quality on the phred33 scale is within the specified threshold values. The sequence quality is translated from phred33 scale to ASCII.
By default, the threshold value is 0.

The result of the program is an output_file containing only those reads that have passed each of the three filters.

## bio_files_processor

This utility is located in the root of the repository and is designed to work with large .gbk, .fasta, .txt files.

#### Usage

You will need:

- A .fasta file for the convert_multiline_fasta_to_oneline function.
- A .txt file for the parse_blast_output function.
- A .gbk file for the select_genes_from_gbk_to_fasta function.

### The utility performs the following actions:

- ### Function ```convert_multiline_fasta_to_oneline```

Converts fasta format from multiline to single line format

The function accepts 2 arguments (input_fasta and output_fasta - the default parameter). Reads the input fasta file, saves it to a new fasta file where each sequence is written to a single line. 

- ### Function ```parse_blast_output```

Parses the file resulting after alignment in the Blast program. Pulling the name of the best match with the database and saving the results in a single file.

The function takes 2 arguments (input_file, output_file) as input. The program reads the given txt file, for each QUERY best match against the base on the first line from the Description column. The set of obtained proteins is saved to a new file with one column sorted alphabetically.

- ### Function ```select_genes_from_gbk_to_fasta```

Detect a certain number of genes before and after each gene of interest and save their protein sequence into a new fasta file. The fasta-file can be sent to the input of BLAST explicitly. 

The function takes as input input_gbk (path to the input GBK file), genes (genes of interest near which neighbors are searched), n_before, n_after (number of genes before and after (>=1)), output_fasta (name of the output file).

# NASeqTool_

<img align=right src="https://github.com/user-attachments/assets/612014fa-fa0d-4102-a23c-0ef233ae3b48" alt="#NASeqTool" width="100"/>

**NASeqTool** биоинформатическая утилита для обработки последовательностей ДНК и РНК.

Репозиторий состоит из главного скрипта и дополнительных модулей, которые используются для обработки последовательностей ДНК/РНК и обработки данных секвенирования. 
Скрипты написаны на языке программирования Python. 


## Установка

Клонируйте репозиторий на локальное устройство.

```sh
git clone git@github.com:Valeriisht/NASeqTool.git
```



## Содержание 

Структура репозитория имеет следующий вид: корневая папка с главным скриптом NASeqTool, где объявлены основные функции, и подпапка functions с импортируемыми модулями.

<img width="400" alt="Rep Structure" src="https://github.com/user-attachments/assets/710383e6-8bbc-4b31-bc2e-e966fca2a58b">


Утилита была разработана для работы на версии Python3.10 или новее.



- ### NASeqTool
    - [**run_dna_rna_tools**](#run_dna_rna_tools_)
    - [**filter_fastq_tool**](#filter_fastq_tool_)
    - [**bio_files_processor**](#bio_files_processor_)

      


## run_dna_rna_tools_

Утилита для манипуляции с молекулами ДНК/РНК. При исполнении скрипта подгружается модуль run_dna_rna_tool_module.py. Находится в скрипте NASeqTool.py



### Использование
Вам понадобятся:

- Последовательности ДНК/РНК, организованные в список (или в строку при выполнении операции над одной последовательностью), которые являются обязательными параметрами для вызова функций.



### Утилита может выполнять следующие операции: 

- ### Функция ```transcribe ```

Смысловая последовательность ДНК переводится в молекулу мРНК. Если на вход передана РНК, возвращается последовательность без изменений.

- ### Функция ```reverse```

Возвращает комплементарную цепь ДНК в направлении от 3' к 5'. Для молекулы РНК возвращается кДНК.

- ### Функция ```reverse_complement```

Возвращает комплементарную цепь ДНК/кДНК в направлении от 5' к 3'.

- ### Функция ```gc_content```

Подсчитывается GC-состав последовательности и возвращается процент GC-нуклеотидов.

- ### Функция ```is_na_sequence```

Проверяется, является ли последовательность ДНК или РНК.

- ### Функция ```length_cds```

Определяется длина кодирующей последовательности CDS.

## filter_fastq_tool_

Утилита для обработки данных секвенирования в формате fastq. При исполнении скрипта подгружается модуль filter_fastq_module.py. Находится в скрипте NASeqTool.py

### Использование

Вам понадобятся:

- Файл в формате fastq.

Пример структуры fastq файла:

![image](https://github.com/user-attachments/assets/3cbf162c-ddb6-4b33-b8c5-4747d5873ebc)

### Утилита выполняет следующие действия:

- Фильтрация последовательности относительно gc_состава. Подсчитывается gc-состав и определяется, находится ли значение в пределах пороговых.
По умолчанию пороговые значения определены от 0 до 100.

- Фильтрация последовательности относительно длины. Проверяется, находится ли длина рида в пределах заданных пороговых значений.
По умолчанию пороговые значения определены от 0 до 2**32.

- Фильтрация последовательности относительно качества рида. Проверяется, соответствует ли качество рида по шкале phred33 заданным пороговому значению. Качество рида переводится из шкалы Phred33 в ASCII.
По умолчанию пороговое значение равняется 0.

Результатом работы программы является output_file, содержащий только те риды, прошедших каждую из трех фильтраций.


## bio_files_processor_

Утилита находится в корне репозитория и расчитана на работу с большими файлами формата .gbk, .fasta, .txt.

### Использование

Вам понадобятся:

- Файл в формата .fasta для функции convert_multiline_fasta_to_oneline.
-  Файл в формата .txt для функции parse_blast_output.
filter_fastq_update
- Файл в формата .gbk для функции select_genes_from_gbk_to_fasta.


### Утилита выполняет следующие действия:

- ### Функция ```convert_multiline_fasta_to_oneline```

Преобразование формата fasta формата из много строчного в однострочный

Функция принимает на вход 2 аргумента (input_fasta и output_fasta - параметр хадан по умолчанию). Читает поданный на вход fasta-файл, сохраняет в новый fasta-файл, в котором каждая последовательность записана в одну строку. 

- ### Функция ```parse_blast_output```

Парсинг файла, полученного после выравнивания в программе Blast. Вытаскивание названия наилучшего совпадения с базой и сохранение  результатов в одном файле.

Функция принимает на вход 2 аргумента (input_file, output_file). Программа читает заданный txt файл, для каждого запроса QUERY наилучшее совпадение по базе по первой строке из столбца Description. Набор полученных белков сохраняется в новый файл одним столбцом отсортированным по алфавиту.

- ### Функция ```select_genes_from_gbk_to_fasta```

Обнаружение определненного количества генов до и после каждого из гена интереса и сохранение их белковой последовательности  в новый fasta-файл. Fasta-файл может быть отправлен прямо в явном виде на вход BLAST. 

Функция ринимает на вход  input_gbk (путь до входного GBK файла), genes (гены интереса, рядом с которыми ищутся соседи), n_before, n_after (количество генов до и послe (>=1)), output_fasta (название выходного файла).










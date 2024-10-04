import sys
sys.path.append("functions")

from functions.filter_fastq_module import (
    is_good_quality_threshold,
    is_good_gc_content,
    is_good_length,
)
from functions.run_dna_rna_tools_module import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    gc_content,
    protein_coding_sequence,
    is_na_sequence,
)


def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Запускает выбранную из словаря функцию.
    Функции импортируются из модуля run_dna_rna_tools.

    :param args:
    :type args: str

    :raises ValueError: если function_name не определена в dict_functions
    :raises ValueError: если seq не является ДНК/РНК

    :rtype: list[any] | str
    :return: результат  работы применяемой функции в виде списка или строки
    """
    *sequences, name_function = args
    dict_functions = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "protein_coding_sequence": protein_coding_sequence,
    }
    result = []
    if name_function not in dict_functions.keys():
        raise ValueError(f"Function {name_function} is not defined")
    for seq in sequences:
        if not is_na_sequence(seq):
            raise ValueError(f"{seq=} is not DNA or RNA")
        result.append(str(dict_functions[name_function](seq)))
    return result if len(sequences) > 1 else result[0]


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | float = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """Функция принимает на вход словарь,
    где в качестве ключа выступает наименование рида,
    в качестве значения последовательность рида и его качество.
    Далее к каждому риду применяет функцию из модуля filter_fastq_module.
    В итоге, возвращается словарь с ридами,
    удовлетворяющими заданным пороговым значениям для фильтрации.

    :param seqs:
    :type seqs: dict[str, tuple[str, str]]
    :param gc_bounds:
    :type gc_bounds: tuple | float
    :param length_bounds:
    :type length_bounds:  tuple | float
    :param quality_threshold:
    :type quality_threshold: float

    :rtype: dict[str, tuple[str, str]]
    :return: dict с отфильтрованными ридами
    """
    good_seqs = {}
    for name_seq, seq_data in seqs.items():
        seq, quality = seq_data
        if all(
            [
                is_good_gc_content(seq, gc_bounds),
                is_good_length(seq, length_bounds),
                is_good_quality_threshold(quality, quality_threshold),
            ]
        ):
            good_seqs[name_seq] = seq_data
    return good_seqs

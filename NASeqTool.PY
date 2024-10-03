import sys
sys.path.insert(0, 'functions')

from functions import filter_fastq_module as fastq
from functions import run_dna_rna_tools_module as tools


def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Запускает выбранную из словаря функцию. Функции импортируются из модуля run_dna_rna_tools.

    :param args:
    :type args: str

    :raises ValueError: если function_name не определена в dict_functions
    :raises ValueError: если seq не является ДНК/РНК

    :rtype: list[any] | str
    :return: результат  работы применяемой функции в виде списка или строки
    """
    *sequences, name_function = args
    dict_functions = {
        "transcribe": tools.transcribe,
        "reverse": tools.reverse,
        "complement": tools.complement,
        "reverse_complement": tools.reverse_complement,
        "gc_content": tools.gc_content,
        "protein_coding_sequence": tools.protein_coding_sequence
    }
    result = []
    if name_function not in dict_functions.keys():
        raise ValueError(f"Функция {name_function} не определена")
    for seq in sequences:
        if not tools.is_na_sequence(seq):
            raise ValueError(f"Ошибка ввода данных: {seq} не ДНК или РНК")
        result.append(str(dict_functions[name_function](seq)))
    if len(sequences) == 1:
        return result[0]
    return result

def filter_fastq(seqs:dict[str, tuple], gc_bounds:tuple=(0,100), length_bounds:float=0.2**32, quality_threshold:float=0) -> dict[str, tuple]:
    """ Функция принимает на вход словарь со значением последовательности рида и его качеством.
    Далее к каждому риду применяет функцию из модуля filter_fastq_module.
    В итоге, возвращается словарь с ридами, удовлетворяющие стандарту качества.

    :param seqs:
    :type seqs: Dict[str, tuple]
    :param gc_bounds:
    :type gc_bounds: tuple
    :param length_bounds:
    :type length_bounds: float
    :param quality_threshold:
    :type quality_threshold: float

    :raises ValueError:

    :rtype: dict[str, tuple]
    :return Словарь с качественными ридами
    """
    good_seqs = {}
    for name_seq, seq in seqs.items():
        if all([fastq.is_good_gc_bounds(seq[0], gc_bounds),fastq.is_good_length_bounds(seq[0], length_bounds), fastq.is_good_quality_threshold(seq[1], quality_threshold)]):
            good_seqs[name_seq] = (seq[0], seq[1])
    return good_seqs


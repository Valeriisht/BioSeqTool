import sys
import os
sys.path.append("functions")

from functions.filter_fastq_module import (
    is_good_quality,
    is_good_gc_content,
    is_good_length,
)
from functions.run_dna_rna_tools_module import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    gc_content,
    length_cds,
    is_na_sequence,
)


def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Runs the function selected from the dictionary.
    Functions are imported from the run_dna_rna_tools module.

    :param args: DNA/RNA sequences
    :type args: str

    :raises ValueError: if function_name is not defined in dict_functions
    :raises ValueError: if the seq is not DNA/RNA

    :rtype: list[any] | str
    :return: the result of the applied function in the form of a list or a string
    """
    *sequences, name_function = args
    dict_functions = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "length_cds": length_cds
    }
    result = []
    if name_function not in dict_functions.keys():
        raise ValueError(f"Function {name_function} is not defined")
    for seq in sequences:
        if not is_na_sequence(seq):
            raise ValueError(f"{seq=} is not DNA or RNA")
        result.append(str(dict_functions[name_function](seq)))
    return result if len(sequences) > 1 else result[0]


def  main(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | float = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> str:
    """The function takes a dictionary as input,
    where the key is the name of the reid,
    the value is the sequence of the reid and its quality.
    Then it applies the function from filter_fastq_module to each reid.
    As a result, it returns a dictionary of reids,
    satisfying the specified threshold values for filtering.

    :param input_fastq:
    :type input_fastq: str
    :param output_fastq:
    :type output_fastq: str
    :param gc_bounds: threshold for filtering by gc-content
    :type gc_bounds: tuple | float
    :param length_bounds: length filtering threshold
    :type length_bounds:  tuple | float
    :param quality_threshold: quality filtering threshold
    :type quality_threshold: float

    :rtype: dict[str, tuple[str, str]]
    :return: dict with filtered reids
    """
    if not os.path.exists(input_fastq):
        raise SystemError("Error: File does not exist")
    with open(input_fastq, "r") as read_file, open(output_fastq, "w") as write_file:
        while True: #читаем, пока не закончиться файл
            seq_data = [read_file.readline().strip() for _ in range(4)]
            if not seq_data:
                break
            name_seq, seq, comment, quality = seq_data #проверка на отсутствие 4 строк ????
            if all(
            [
                is_good_gc_content(seq, gc_bounds),
                is_good_length(seq, length_bounds),
                is_good_quality(quality, quality_threshold),
             ]):
                for line in seq_data:
                    write_file.write(f"{line}\n")
    return output_fastq

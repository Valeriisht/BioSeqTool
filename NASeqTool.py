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

def read_fastq_file(input_fastq:str) -> dict[str, tuple[str, str]]:
    if not os.path.exists(input_fastq):
        raise SystemError("Error: File does not exist")
    with open(input_fastq, "r") as fastq_file:
        seqs = {}
        one_seq = []
        for fastq_seq in fastq_file:
            if len(one_seq) < 5:
                one_seq.append(fastq_seq.strip())
            else:
                seqs[one_seq[0]] = (one_seq[-2], one_seq[-1])
                one_seq = []
    return seqs

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | float = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """The function takes a dictionary as input,
    where the key is the name of the reid,
    the value is the sequence of the reid and its quality.
    Then it applies the function from filter_fastq_module to each reid.
    As a result, it returns a dictionary of reids,
    satisfying the specified threshold values for filtering.

    :param input_fastq:
    :type input_fastq: str
    :param seqs: a dictionary with the DNA sequence and its quality
    :type seqs: dict[str, tuple[str, str]]
    :param gc_bounds: threshold for filtering by gc-content
    :type gc_bounds: tuple | float
    :param length_bounds: length filtering threshold
    :type length_bounds:  tuple | float
    :param quality_threshold: quality filtering threshold
    :type quality_threshold: float

    :rtype: dict[str, tuple[str, str]]
    :return: dict with filtered reids
    """
    seqs = read_fastq_file(input_fastq)
    good_seqs = {}
    for name_seq, seq_data in seqs.items():
        seq, quality = seq_data
        if all(
            [
                is_good_gc_content(seq, gc_bounds),
                is_good_length(seq, length_bounds),
                is_good_quality(quality, quality_threshold),
            ]
        ):
            good_seqs[name_seq] = seq_data
    return good_seqs

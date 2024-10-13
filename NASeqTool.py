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
    :return: the result of the applied function
    in the form of a list or a string
    """
    *sequences, name_function = args
    dict_functions = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "length_cds": length_cds,
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
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | float = (0, 2 ** 32),
    quality_threshold: float = 0,
):
    """The function takes a file or path to file as input.
    The function from filter_fastq_module is
    applied to each reed for filtering by the input values.
    The result of the program is an output_file containing
    only those reads that have passed each of the three filters.

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

    """
    if not os.path.exists(input_fastq):
        raise SystemError("File does not exist")
    if not os.path.isdir("filtered"):
        os.makedirs("filtered")
    output_path = os.path.join("filtered", output_fastq)
    with (open(input_fastq, "r") as read_file, open(output_path, "w") as write_file):
        seq_data = []
        count = 0
        for line in read_file:
            if count != 4:
                seq_data.append(line.strip())
                count += 1
            if len(seq_data) == 4 and seq_data[0].startswith(">"):
                name_seq, seq, comment, quality = seq_data
                if all(
                    [
                        is_good_gc_content(seq, gc_bounds),
                        is_good_length(seq, length_bounds),
                        is_good_quality(quality, quality_threshold),
                    ]
                ):
                    for sequence in seq_data:
                        write_file.write(f"{sequence}\n")
                    seq_data = []
                    count = 0

filter_fastq("C:\\Users\\valer\\OneDrive\\Desktop\\example_fastq.fastq", "filter_fastq_output.fastq")

import functions.filter_fastq_module as fastqc
import functions.run_dna_rna_tools_module as tools


def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Запускает выбранную из словаря функцию"""
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

def filter_fastq(seqs:dict[str, tuple], gc_bounds:tuple=(0,100), length_bounds:float=0.2**32, quality_threshold:float=0):
    good_seqs = {}
    for name_seq, seq, quality_seq in seqs.items():
        if all([fastqc.is_good_gc_bounds(seq, gc_bounds),fastqc.is_length_bounds(seq, length_bounds)]):
            good_seqs[name_seq] = (seq, quality_seq)
    return good_seqs
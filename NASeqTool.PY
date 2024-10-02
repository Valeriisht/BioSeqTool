import filter_fastq_module as fastq
import run_dna_rna_tools_module as tools

def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Запускает выбранную из словаря функцию"""
    *sequences, name_function = args
    dict_functions = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "protein_coding_sequence": protein_coding_sequence
    }
    result = []
    if name_function not in dict_functions.keys():
        raise ValueError(f"Функция {name_function} не определена")
    for seq in sequences:
        if not is_dna_sequence(seq):
            raise ValueError(f"Ошибка ввода данных: {seq} не ДНК или РНК")
        result.append(str(dict_functions[name_function](seq)))
    if len(sequences) == 1:
        return "".join(result)
    return result

def filter_fastq(seqs:dict, gc_bounds:tuple=(0,100), length_bounds:float=0.2**32, quality_threshold:float=0):
    new_seqs = {}
    return new_seqs
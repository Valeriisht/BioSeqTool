

def run_dna_rna_tools(*args: str) -> list[any] | str:
    """Запускает выбранную из словаря функцию"""
    name_function, *sequences = args[::-1]
    dict_functions = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "protein_coding_sequence": protein_coding_sequence
    }  # eval
    result = []
    if name_function not in dict_functions.keys():
        raise ValueError(f"Функция {name_function} не определена")
    for seq in sequences[::-1]:
        if not is_dna_sequence(seq):
            raise ValueError(f"Ошибка ввода данных: {seq} не ДНК или РНК")
        result.append(str(dict_functions[name_function](seq)))
    if len(sequences) == 1:
        return "".join(result)
    return result


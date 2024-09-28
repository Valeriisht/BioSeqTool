
def transcribe(seq: str) -> str:
    """Функция переводит смысловую молекулу ДНК в молекулу мРНК.\
    Если на вход передана РНК, возвращается  последовательность без изменений
    """
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq: str) -> str:
    """Возвращает последовательность ДНК или РНК от 3' к 5'"""
    return seq[::-1]


def complement(seq: str) -> str:
    """Возвращает комплементарную цепь ДНК/РНК в направлении от 3' к  5'"""
    dict_complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "U": "A",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "u": "a"
    }
    return "".join([dict_complement[nucleotide] for nucleotide in seq])


def reverse_complement(seq: str) -> str:
    """Возвращает комплементарную цепь ДНК/кДНК в направлении от 5' к 3'"""
    return complement(seq)[::-1]


def gc_content(seq: str) -> float:
    """Подсчитывает GC-состав последовательности \
     возвращает долю GC-нуклеотидов в последовательности"""
    count_g = seq.upper().count("G")
    count_c = seq.upper().count("C")
    return round((count_c + count_g) / len(seq), 2)


def is_coding_sequence(seq: str) -> bool:
    """Проверяет, может предположительно последовательность кодировать белок"""
    stop_codons = ("UAA", "UAG", "UGA")
    rna_seq = transcribe(seq)
    return (any([codon in rna_seq.upper() for codon in stop_codons])
            and "AUG" in rna_seq)


def protein_coding_sequence(seq: str) -> str:
    """Определяет длину кодирующей последовательности CDS"""
    stop_codons = ("UAA", "UAG", "UGA")
    rna_seq = transcribe(seq)
    if is_coding_sequence(seq):
        pos_start = rna_seq.upper().index("AUG")
        pos_stop = min([
            index_stop
            for index_stop in range(pos_start, len(rna_seq), 3)
            if rna_seq[index_stop:(index_stop + 3)].upper() in stop_codons
        ])
        len_seq = abs(pos_start - (pos_stop + 3))
        if len_seq > 0:
            return (f"Предположительно, {seq=} кодирующая последовательность\n"
                    f"Длина CDS составляет {len_seq} пар оснований")
        return f"Последовательность {seq=} не кодирует белок"
    return f"Последовательность {seq=} не кодирует белок"


def is_dna_sequence(seq: str) -> bool:
    """Проверяет, является ли последовательность ДНК или РНК"""
    dna = {"A", "T", "C", "G"}
    rna = {"A", "U", "C", "G"}
    return set(seq.upper()).issubset(dna) or set(seq.upper()).issubset(rna)

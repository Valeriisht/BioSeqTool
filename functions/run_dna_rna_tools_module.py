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
    "u": "a",
}


def transcribe(seq: str) -> str:
    """The function translates a meaningful DNA molecule into an mRNA molecule. \
    If RNA is given as input, the sequence is returned unchanged.
    """
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq: str) -> str:
    """Returns the DNA or RNA sequence from 3' to 5'"""
    return seq[::-1]


def complement(seq: str) -> str:
    """Returns the complementary DNA/cDNA strand in the 3' to 5' direction"""
    return "".join([dict_complement[nucleotide] for nucleotide in seq])


def reverse_complement(seq: str) -> str:
    """Returns the complementary DNA/cDNA strand in the 5' to 3' direction"""
    return complement(seq)[::-1]


def gc_content(seq: str) -> float | str:
    """Counts the GC content of the sequence \
     returns the percentage of GC nucleotides in the sequence"""
    count_g = seq.upper().count("G")
    count_c = seq.upper().count("C")
    return round((count_c + count_g) / len(seq), 2) * 100


def is_coding_sequence(seq: str) -> bool:
    """Checks to see if the putative sequence 
    can encode a protein"""
    stop_codons = ("UAA", "UAG", "UGA")
    rna_seq = transcribe(seq)
    return any([codon in rna_seq.upper()
                for codon in stop_codons]) and "AUG" in rna_seq


def length_cds(seq: str) -> str:
    """Defines the length of CDS"""
    stop_codons = ("UAA", "UAG", "UGA")
    rna_seq = transcribe(seq)
    if is_coding_sequence(seq):
        pos_start = rna_seq.upper().index("AUG")
        pos_stop = min(
            [
                index_stop
                for index_stop in range(pos_start, len(rna_seq), 3)
                if rna_seq[index_stop: (index_stop + 3)].upper() in stop_codons
            ]
        )
        len_seq = abs(pos_start - (pos_stop + 3))
        if len_seq > 0:
            return (
                f"Supposedly, {seq=} the coding sequence\n"
                f"The length of the CDS is {len_seq} base pairs"
            )
    return f"The {seq=} sequence does not encode a protein"


def is_na_sequence(seq: str) -> bool:
    """Checks whether the DNA or RNA sequence"""
    dna = {"A", "T", "C", "G"}
    rna = {"A", "U", "C", "G"}
    return set(seq.upper()).issubset(dna) or set(seq.upper()).issubset(rna)

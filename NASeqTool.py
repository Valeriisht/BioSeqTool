import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences

    Args:
        ABC (_type_): Abstract Base Class from abc module

    Methods:
        __len__(self): Returns the length of the biological sequence
        __getitem__(self, key): Allows accessing individual elements or slices of the sequence
        __repr__(self): Returns a string representation of the object
        correct_alphabet(sequence, alphabet): Static method to check  sequence
    """

    @abstractmethod
    def __len__(self):
        """Returns the length of the biological sequence

        Returns:
            int: Length of the biological sequence+
        """
        pass

    def __getitem__(self, key):
        """
        Accessing individual elements or slices of the sequence

        Args:
            key (int or slice): Index or slice of the sequence

        Returns:
            BiologicalSequence: Individual element or slice of the sequence
        """
        pass

    def __repr__(self):
        """
        Returns a string representation of the object

        Returns:
            str: String representation of the object
        """
        pass

    @staticmethod
    def correct_alphabet(sequence, alphabet):
        """
        Static method to check if a sequence contains only valid characters from a given alphabet.#+

        Args:
            sequence (str): Biological sequence
            alphabet (str): Valid contain for the biological sequence

        Raises:
            ValueError: If the sequence contains characters not present in the alphabet
        """
        pass


class PolymerSequence(BiologicalSequence):
    """_summary_

    Args:
        BiologicalSequence (_type_): _description_
    """

    def __init__(self, sequence):
        self.sequence = sequence
        self.alphabet = ""
        self.dict_complement = {}

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        if isinstance(key, int):
            if key < 0:
                key += len(self.sequence)
            elif key >= len(self.sequence):
                raise IndexError("The range out of thr list")
            return self.sequence[key]
        elif isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            start = 0 if start is None else start
            stop = len(self.sequence) - 1 if stop is None else stop
            step = 1 if step is None else step
            start = len(self.sequence) + start if start < 0 else start
            stop = len(self.sequence) + stop if stop < 0 else stop
            return self.sequence[start:stop:step]

    def __repr__(self):
        return f"{self.__class__.__name__}({self.sequence})"

    @staticmethod
    def correct_alphabet(sequence, alphabet):
        if not set(sequence).issubset(alphabet):
            raise ValueError("Incorrect entry")
        return True


class NucleicAcidSequence(PolymerSequence):
    """_summary_

    Args:
        PolymerSequence (_type_): _description_
    """

    def complement(self):
        """_summary_

        Args:
            dict_com (dict): _description_

        Returns:
            _type_: _description_
        """

        try:
            if PolymerSequence.correct_alphabet(self.sequence, self.alphabet):
                return "".join(
                    [self.dict_complement[nucleotide] for nucleotide in self.sequence]
                )
        except NotImplementedError:
            print("The method must be implemented in child classes")

    def reverse(self):
        """Returns the DNA or RNA sequence from 3' to 5'

        Returns:
            _type_: _description_
        """
        return self.sequence[::-1]

    def reverse_complement(self):
        """Returns the complementary DNA/cDNA strand in the 5' to 3' direction

        Returns:
            _type_: _description_
        """
        return self.complement()[::-1]

    def gc_content(self) -> float | str:
        """Counts the GC content of the sequence \
        returns the percentage of GC nucleotides in the sequence

        Returns:
            float | str: _description_
        """
        count_g = self.sequence.upper().count("G")
        count_c = self.sequence.upper().count("C")
        return round((count_c + count_g) / len(self.sequence) * 100, 2)


class AminoAcidSequence(PolymerSequence):
    """_summary_

    Args:
        PolymerSequence (_type_): _description_
    """

    def __init__(self, sequence):
        super().__init__(sequence)
        self.alphabet = "ACDEFGHIKLMNPQRSTVW"

    @staticmethod
    def correct_alphabet(sequence, alphabet):
        """Checks whether the AminoAcid sequence"""

        if not set(sequence).issubset(alphabet):
            raise ValueError("Incorrect entry")
        return True

    def find_motif(self, motif):
        """_summary_

        Args:
            motif (_type_): _description_

        Returns:
            _type_: _description_
        """

        position = self.sequence.find(motif)
        count = self.sequence.count(motif)

        return f"Motif {motif} is found on {position} position for {count} count"


class DNASequence(NucleicAcidSequence):
    """_summary_

    Args:
        NucleicAcidSequence (_type_): _description_
    """

    def __init__(self, sequence):
        super().__init__(sequence)
        self.alphabet = ("A", "G", "C", "T")
        self.dict_complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }

    def transcribe(self):
        """The function translates a meaningful DNA molecule into an mRNA molecule. \
        If RNA is given as input, the sequence is returned unchanged.
        

        Returns:
            _type_: _description_
        """
        return self.sequence.replace("T", "U").replace("t", "u")


class RNASequence(NucleicAcidSequence):
    """_summary_

    Args:
        NucleicAcidSequence (_type_): _description_
    """

    def __init__(self, sequence):
        super().__init__(sequence)
        self.alphabet = ("A", "G", "C", "U")
        self.dict_complement = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }

    def is_coding_sequence(self) -> bool:
        """Checks to see if the putative sequence
        can encode a protein

        Returns:
            bool: _description_
        """
        stop_codons = ("UAA", "UAG", "UGA")
        return (
            any([codon in self.sequence.upper() for codon in stop_codons])
            and "AUG" in self.sequence
        )

    def length_cds(self):
        """Defines the length of CDS

        Returns:
            str: _description_
        """
        stop_codons = ("UAA", "UAG", "UGA")
        if self.is_coding_sequence():
            pos_start = self.sequence.upper().index("AUG")
            pos_stop = min(
                [
                    index_stop
                    for index_stop in range(pos_start, len(self.sequence), 3)
                    if self.sequence[index_stop : (index_stop + 3)].upper()
                    in stop_codons
                ]
            )
            len_seq = abs(pos_start - (pos_stop + 3))
            if len_seq > 0:
                return (
                    f"Supposedly, {self.sequence} the coding sequence\n"
                    f"The length of the CDS is {len_seq} base pairs"
                )
        return f"The {self.sequence} sequence does not encode a protein"


# Функция Bio.SeqIO.parse() в Biopython используется для чтения последовательностей из файла или потока данных.
# Она возвращает итератор, который выдает объекты SeqRecord, позволяет обрабатывать последовательности по одной в порядке их следования в файле


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | float = (0, 2**32),
    quality_threshold: float = 0,
):
    """The function takes a file or path to file as input.
    The function from filter_fastq_module is
    applied to each read for filtering by the input values
    The result of the program is an output_file containing
    only those reads that have passed each of the three filters.

    Args:
        input_fastq (str): _description_
        output_fastq (str): _description_
        gc_bounds (tuple | float, optional): _description_. Defaults to (0, 100)
        length_bounds (tuple | float, optional): _description_. Defaults to (0, 2 ** 32)
        quality_threshold (float, optional): _description_. Defaults to 0
        input_fastq (str): The path to the input FASTQ file
        output_fastq (str): The name of the output FASTQ file
        gc_bounds (tuple | float, optional): The lower and upper bounds for GC content.#+
            Defaults to (0, 100)
        length_bounds (tuple | float, optional): The lower and upper bounds for sequence length.#+
            Defaults to (0, 2 ** 32)
        quality_threshold (float, optional): The minimum average quality score for a read.#+
            Defaults to 0.#+

    Raises:
        SystemError: _description_#-
        SystemError: If the input FASTQ file does not exist.#+
    """
    if not os.path.exists(input_fastq):
        raise SystemError("File does not exist")

    if not os.path.isdir("filtered"):
        os.makedirs("filtered")

    output_path = os.path.join("filtered", output_fastq)

    with open(input_fastq, "r") as read_file, open(output_path, "w") as write_file:

        for record in SeqIO.parse(read_file, "fastq"):
            seq = str(record.seq)
            quality = record.letter_annotations["phred_quality"]

            gc_content = GC(seq)

            gc_lower, gc_upper = (
                gc_bounds if isinstance(gc_bounds, tuple) else (0, gc_bounds)
            )
            is_good_gc = gc_lower <= gc_content <= gc_upper

            length_lower, length_upper = (
                length_bounds
                if isinstance(length_bounds, tuple)
                else (0, length_bounds)
            )
            is_good_length = length_lower <= len(seq) <= length_upper

            is_good_quality = sum(quality) / len(quality) >= quality_threshold

            filters = [is_good_gc, is_good_length, is_good_quality]

            if all(filters):
                SeqIO.write(record, write_file, "fastq")

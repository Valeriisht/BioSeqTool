import os

if not os.path.isdir("bio_files_output"):
    os.makedirs("bio_files_output")


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """The function takes as input a file or file path in .fasta format.
    The output file is returned,
    where multi-line protein sequences are converted into a single-line record.

    :param input_fasta: input_file .fasta format
    :type input_fasta: str
    :param output_fasta: output_file .fasta format
    :type output_fasta: str

    """
    if not os.path.exists(input_fasta):
        raise SystemError("File does not exist")
    if output_fasta is None:
        output_fasta = "multiline_fasta_output.fasta"  # Create a new file
    with open(input_fasta, "r") as read_file, open(
        os.path.join("bio_files_output", output_fasta), "w"
    ) as write_file:
        multi_line = ""
        for line in read_file:
            if line.strip().startswith(">"):
                if multi_line:
                    write_file.write(f"{multi_line}\n")
                write_file.write(f"{line}")
                multi_line = ""
            else:
                multi_line += line.strip()
        if multi_line:
            write_file.write(f"{multi_line}\n")


def parse_blast_output(input_file: str, output_file: str = None):
    """The function takes as input a file or path to a .txt file
    that is received after alignment in the blast program.
    As a result of the function operation,
    the description of the protein with
    the best database match is written to the output file.

    :param input_file: input_file .txt format
    :type input_file: str
    :param output_file: output_file .txt format
    :type output_file: str

    """
    if not os.path.exists(input_file):
        raise SystemError("File does not exist")
    if output_file is None:
        output_file = "parse_blast_output.txt"
    with (
        open(input_file, "r") as read_file,
        open(os.path.join("bio_files_output", output_file), "w") as write_file,
    ):
        protein_description = []
        flag_description_protein = False
        for line in read_file:
            # Pull only the protein descriptions from the file
            if line.startswith("Sequences producing significant alignments:"):
                flag_description_protein = True
            if line.startswith("Alignments"):
                flag_description_protein = False
            if (
                flag_description_protein
                and not line.strip().startswith("Description")
                and not line.strip().startswith("Scientific")
                and not line.strip().startswith("Sequences producing significant alignments")
            ):
                if "    " in line:
                    line = line.strip().split("    ")[0]
                elif "]" in line:
                    line = line.strip().split("]")[0] + "]"
                elif "   " in line:
                    line = line.strip().split("   ")[0]
                elif "..." in line:
                    line = line.strip().split("...")[0]
                elif "...." in line:
                    line = line.strip().split("....")[0]
                if line:
                    protein_description.append(line)
        if protein_description:
            protein_description = sorted(protein_description, key=str.lower)
            # Sort without case sensitivity
        for description in protein_description:
            if description and description != "\n":
                write_file.write(f"{description}\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: list,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = None
):
    """The function takes as input a file or path to a file in .gbk format.
    As a result of the function work, the genes
    (and their protein sequences) that are located
    next to the genes of interest are written
     to the output file in .fasta format.

    :param input_gbk: input_file .gbk format
    :type input_gbk: str
    :param output_fasta: output_file .fasta format
    :type output_fasta: str
    :param genes: genes of interest
    :type genes: list
    :param n_before: number
    :type n_before: int
    :param n_after: number
    :type n_after: int

    """
    if not os.path.exists(input_gbk):
        raise SystemError("File does not exist")
    if output_fasta is None:
        output_fasta = "select_genes_from_gbk_to_fasta_output.fasta"
    with open(input_gbk, "r", encoding="gbk") as read_file, open(
        os.path.join("bio_files_output", output_fasta), "w", encoding="utf-8"
    ) as write_file:
        flag_found_gene = False
        flag_protein_start = False
        list_gene_protein = []
        name_gene = ""
        for line in read_file:
            # Pull only the genes and their protein sequence from the file
            if line.strip().startswith("/gene"):
                flag_found_gene = True
                name_gene = line.strip('\\, \n,"').split('="')[1]
            if line.strip().startswith("/translation") and flag_found_gene:
                flag_protein_start = True
                protein_seq_len = line.strip().split('="')[1]
                line = read_file.readline()
            if line.strip().startswith("CDS") and name_gene:
                flag_found_gene = False
                flag_protein_start = False
                list_gene_protein.append(name_gene)
                list_gene_protein.append(protein_seq_len)
            if flag_found_gene and flag_protein_start:
                protein_seq_len += line.strip()

        def write_sequences(
            index: int, shift: int
        ):  # Function for writing sequences to file
            for _ in range(n_before, 0, -1):
                write_file.write(f">{list_gene_protein[index + shift]}\n")
                write_file.write(f"{list_gene_protein[index +shift + 1]}\n")
                index += shift

        for name_gene in list_gene_protein:
            for gene in genes:
                if gene in name_gene:
                    index_gene = list_gene_protein.index(name_gene)
                    write_sequences(index_gene, -2)
                    write_sequences(index_gene, 2)

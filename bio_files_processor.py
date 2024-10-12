import os


current_directory = os.getcwd()
if not os.path.isdir('bio_files_output'):
    os.makedirs("bio_files_output")

def convert_multiline_fasta_to_oneline(input_fasta:str, output_fasta:str=None)-> str:
    """ The function takes a file or path to file as input,
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

    :rtype: str
    :return: output_file
    """
    if not os.path.exists(input_fasta):
        raise SystemError("File does not exist")
    if output_fasta is None:
        output_fasta = input_fasta.replace('.fasta', '_out.fasta') #создаем на основе старого файла
    with open(input_fasta, "r") as read_file, open(os.path.join("bio_files_output", output_fasta), "w") as write_file:
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
    return output_fasta


convert_multiline_fasta_to_oneline("example_multiline_fasta.fasta")

def parse_blast_output(input_file:str, output_file:str=None)-> str:
    """The function takes a file or path to file as input,
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

    :rtype: str
    :return: output_file
    """
    if not os.path.exists(input_file):
        raise SystemError("File does not exist")
    if output_file is None:
        output_file = input_file.replace('.txt', '_out.txt')
    with open(os.path.join(current_directory, input_file), "r") as read_file, open(os.path.join("bio_files_output", output_file), "w") as write_file:
        protein = []
        while True:
            line = read_file.readline()
            if line.startswith("Sequences producing significant alignments:"):
                line = read_file.readline(3)
                while not line.startswith("Alignments"):
                    line = line.strip().split("    ")
                    protein.append(line[0])
                    line = read_file.readline()
            if not line:  # Проверка на конец файла
                break
        protein = sorted(protein)
        for description in protein:
            if description and description != "\n":
                write_file.write(f"{description}\n")
    return output_file



def select_genes_from_gbk_to_fasta(input_gbk:str, genes:list, n_before:int=1, n_after:int=1, output_fasta:str|None = None) -> str:
    """The function takes a file or path to file as input,
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

    :rtype: str
    :return: output_file
    """
    if not os.path.exists(input_gbk):
        raise SystemError("File does not exist")
    if output_fasta is None:
        output_fasta = input_gbk.replace('.gbk', '_gbk.fasta')
    with open(input_gbk, "r", encoding='gbk') as read_file, open(os.path.join("bio_files_output", output_fasta), "w") as write_file:
        count = 0
        current_genes_position = []
        while True:
            line = read_file.readline()
            if not line:  #Проверка на конец файла
                break
            if line.strip().startswith("/gene"):
                current_genes_position.append(read_file.tell())
                name_gene = line.strip().split("=")[1] #извлечение имени гена
                right_gene= [gene in name_gene for gene in genes] #Проверка, есть ли такой ген в списке генов интереса
                if any(right_gene):
                    current_position = read_file.tell()
                    read_file.seek(current_genes_position[-(n_before-)]) #перемещаемся назад
                    name_for_seq = line.strip().split("=")[1]1
                    for seq in read_file:
                        if count >= n_before:
                            count = 0
                            break
                        if seq.strip().startswith("/gene"):
                            name_for_seq = seq.strip().split("=")[1]
                        if seq.strip().startswith("/translation"):
                            protein_seq_len = seq.strip().split("=")[1]
                            while not seq.startswith("CDS"):
                                seq = read_file.readline().strip()
                                protein_seq_len += seq
                            write_file.write(f">{name_for_seq}\n")
                            write_file.write(f"{protein_seq_len}\n")
                            count += 1
                    read_file.seek(current_position)
                    name_for_seq = line.strip().split("=")[1]
                    for seq in read_file:
                        if count >= n_after:
                            count = 0
                            read_file.seek(current_position)
                            break
                        if seq.strip().startswith("/gene"):
                            name_for_seq = seq.strip().split('=')[1]
                        if seq.strip().startswith("/translation"):
                            protein_seq_len = seq.strip().split("=")[1]
                            while not seq.startswith("CDS") :
                                seq = read_file.readline().strip()
                                if not seq:
                                    break
                                protein_seq_len += seq
                            write_file.write(f">{name_for_seq}\n")
                            write_file.write(f"{protein_seq_len}\n")
                            count += 1
        return output_fasta
select_genes_from_gbk_to_fasta("ex_gbk.txt", ["dtpD"], 2, 1)

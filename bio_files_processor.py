import os
import sys

current_directory = os.getcwd()
if not os.path.isdir('bio_files_output'):
    os.makedirs("bio_files_output")

def convert_multiline_fasta_to_oneline(input_fasta:str, output_fasta:str=None)-> str:
    if not os.path.exists(input_fasta):
        raise SystemError("Error: File does not exist")
    if output_fasta is None:
        output_fasta = input_fasta.replace('.fasta', 'out.fasta') #создаем на основе старого
    with open(os.path.join(current_directory, input_fasta), "r") as read_file, open(os.path.join("bio_files_output", output_fasta), "w") as write_file:
        multi_line = ""
        for line in read_file:
            if line.strip().startswith(">"):
                if multi_line:
                    write_file.write(f"{multi_line}\n")
                write_file.write(f"{line}\n")
                multi_line = ""
            else:
                multi_line += line.strip()
        if multi_line:
            write_file.write(f"{multi_line}\n")
    return output_fasta


def parse_blast_output(input_file:str, output_file:str=None)-> str:
    if not os.path.exists(input_file):
        raise SystemError("Error: File does not exist")
    if output_file is None:
        output_file = "default_output.txt"
    with open(input_file, "r") as read_file, open(output_file, "w") as write_file:
        line = ""
        proteins =[]
        while not line.strip().startswith("Sequences producing significant alignments:"):
            line = read_file.readline()
        line = read_file.readline().strip()
        line = read_file.readline().strip() #Пропускаем строчку с названиями
        position = read_file.tell()
        read_file.seek(position)
        for line in read_file:
            write_file.write(line[0])
    return output_file


convert_multiline_fasta_to_oneline("example_multiline_fasta.fasta")


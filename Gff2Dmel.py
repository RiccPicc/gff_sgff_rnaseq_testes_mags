from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO


# Path to the input file containing Glossina fuscipes proteins in FASTA format
input_file = "/Users/riccardo/Library/CloudStorage/GoogleDrive-riccardo.piccinno@unipv.it/My Drive/Spiroplasma/VectorBase-63_GfuscipesIAEA2018_AnnotatedProteins.fasta"
# Path to the BLAST database containing Drosophila melanogaster proteins
blast_db = "/Users/riccardo/blast/db/Dmeldb"

output_file = "blast_results.txt"

blastp_cline = NcbiblastpCommandline(query=input_file, db=blast_db, out=output_file,
                                     outfmt="6 qseqid stitle evalue ",
                                     max_target_seqs=1)

stdout, stderr = blastp_cline()

new_output = "Gff2Dmel.txt"

new_lines = ""

with open(output_file, "r") as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.split("\t")
        line[0] = line[0].split(".")[0]
        line[1] = line[1].split("|")[2].split("=")[1].strip(" ")
        line = "\t".join(line)
        new_lines += line

with open(new_output, "w+") as handle:
    handle.write(new_lines)

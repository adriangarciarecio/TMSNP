import re

file1 = open("uniprot_tm.txt")
file2 = open("patho_clinvar.txt")
codi = "hola"
for line in file2.readlines():
    parts = line.split()
    gene = parts[3]
    if gene in open("uniprot_tm.txt").read():
        print(gene, parts[0], parts[1], parts[2], parts[4])

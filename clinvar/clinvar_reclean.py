import re

file1 = open("uniprot_tm.txt")
file2 = open("clinvar_clean.txt")
codi = "clinvar_clean.txt"
for line in file2.readlines():
    parts = line.split()  # split line into parts
    gene = parts[0]
    file1.seek(0, 0)
    for line in file1.readlines():
        parts2 = line.split()
        uniprot = parts2[0]
        gene2 = parts2[1]
        if gene2 == gene:
            codi = uniprot
            print(codi, parts[1], parts[2], parts[3], parts[4])
        if gene2 != gene:
            codi2 = uniprot


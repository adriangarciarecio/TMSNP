# Cheks that proteins are indeed membrane proteins
# Integrates all the information

import re

list1 = open("gnomad_mutations.txt")
list3 = open("gnomad_all.txt", "w")
for line in list1.readlines():
    uniprot = line[0:6]
    linia = line[0:60]
    if uniprot in open("membrane_proteins.txt").read():
        list3.writelines(linia)

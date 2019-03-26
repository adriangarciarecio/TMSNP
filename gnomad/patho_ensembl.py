# Reads list.txt contains all membrane proteins ENS code
# Generates a list of ENS cods for Uniprot protein containing
# pathological mutations

# Read csv (from SQL) and extracts pathological snps

import sqlalchemy as sql
import os
import pandas as pd


list1 = open("ensembl_list.txt")  # input
list3 = open("list2.txt", "w")


# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")

snps = pd.read_sql("SELECT * FROM snps", engine)
proteins = snps["acc"].unique()

for line in list1.readlines():
    uniprot = line[0:6]
    if uniprot in proteins:
        list3.writelines(line)

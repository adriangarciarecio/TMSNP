#!/usr/bin/env python3

import sqlalchemy
import os
import pandas as pd


# Removes proteins that do not contain any Pathogenic SNP

engine = sqlalchemy.create_engine(
    f'mysql+mysqlconnector://lmcdb:{os.getenv("LMCDB_PASS")}@alf03.uab.cat/tmsnp'
)

snps = pd.read_sql("SELECT * FROM snps", engine)
proteins = snps["acc"].unique()

print("Filtering non-pathogenic proteins")
print("Before filtering")
print("Proteins:", len(proteins))
print("SNPs:", len(snps))

non_pathogenic_prots = []
for protein in proteins:
    sum_pathogenic = snps[snps["acc"] == protein]["pathogenic"].sum()
    if sum_pathogenic == 0:
        # print(protein)
        non_pathogenic_prots.append(protein)

for protein in non_pathogenic_prots:
    # print(protein)
    snps.drop(snps[snps["acc"] == protein].index, inplace=True)

print("After filtering")
print("Proteins:", len(snps["acc"].unique()))
print("SNPs:", len(snps))

snps.to_sql("snps", engine, if_exists="replace", index=False, chunksize=10000)

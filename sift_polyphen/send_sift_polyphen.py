#!/usr/bin/env python

# Takes the TMs from the MySQL tabe

import mysql.connector
import sqlalchemy
import os
import pandas as pd
import numpy as np
import sqlalchemy as sql

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")

df = pd.read_sql("select acc, snp_pos, aa_ref, aa_mut from snps", engine)
df.to_csv("to_polyphen_sift.txt", sep=" ", index=False, header=False)

print("Go to the following pages and use to_polyphen_sift.txt")
print("- Polyphen2: http://genetics.bwh.harvard.edu/pph2/bgi.shtml")
print(
    "- SIFT: http://sift.jcvi.org/protein_batch_submit.php?species=human PROVEAN Protein Batch"
)
print("From Polyphen2 download the SHORT results file!")

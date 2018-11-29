#!/usr/bin/env python

# Takes the TMs from the MySQL tabe

import mysql.connector
import sqlalchemy
import os
import pandas as pd
import numpy as np

# Set the connection
conn = mysql.connector.connect(
    host="alf03.uab.cat",
    user="lmcdb",
    password=os.getenv("LMCDB_PASS"),
    database="tmsnp",
)
mycursor = conn.cursor()
mycursor.execute("select acc, snp_pos, aa_ref, aa_mut from snps")

table_rows = mycursor.fetchall()
df = pd.DataFrame(table_rows, columns=["acc", "snp_pos", "aa_ref", "aa_mut"])
df.to_csv("to_polyphen_sift.txt", sep=" ", index=False, header=False)

print("Go to the following pages and use to_polyphen_sift.txt")
print("- Polyphen2: http://genetics.bwh.harvard.edu/pph2/bgi.shtml")
print("- SIFT: http://sift.jcvi.org/protein_batch_submit.php?species=human")
print("From Polyphen2 download the SHORT results file!")

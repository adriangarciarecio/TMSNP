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
connection = engine.connect()

try:
    connection.execute("alter table snp_eval add column pph2 int(1)")
except:
    pass
try:
    connection.execute("alter table snp_eval add column sift int(1)")
except:
    pass


df = pd.read_csv("sift.tsv", sep="\t")
# df.columns = df.columns.str.strip() # trim spaces in column names
# df = df.applymap(lambda x: x.strip() if type(x) is str else x) # trim spaces in all fields
df = df.rename(columns={"PREDICTION (cutoff=-2.5)": "prediction"})
df["prediction"] = df["prediction"].str.replace("Deleterious", "1")
df["prediction"] = df["prediction"].str.replace("Neutral", "0")
df.prediction = df.prediction.fillna(9).astype(
    int
)  # 9 means null is a trick for the convertion to int
df["POSITION"] = df["POSITION"].fillna(0).astype(int)

for i in df.index:
    if i % 100 == 0:
        print(i)
    connection.execute(
        f"""update snp_eval set sift={df.loc[i]['prediction']} where acc="{df.loc[i]['PROTEIN_ID']}" and snp_pos={df.loc[i]['POSITION']} and aa_mut="{df.loc[i]['RESIDUE_ALT']}";"""
    )
connection.close()

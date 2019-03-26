#!/usr/bin/env python

# Processs GNOMAD data (after Mireia's scripts)
# Takes the TMs from the MySQL tabe

import mysql.connector
import sqlalchemy as sql
import os
import pandas as pd
import numpy as np

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

three2one = {
    "Cys": "C",
    "Asp": "D",
    "Ser": "S",
    "Gln": "Q",
    "Lys": "K",
    "Ile": "I",
    "Pro": "P",
    "Thr": "T",
    "Phe": "F",
    "Asn": "N",
    "Gly": "G",
    "His": "H",
    "Leu": "L",
    "Arg": "R",
    "Trp": "W",
    "Ala": "A",
    "Val": "V",
    "Glu": "E",
    "Tyr": "Y",
    "Met": "M",
}

df_tm = pd.read_sql("select * from tm_segments", engine)
df_clinvar = pd.read_csv(
    "clinvar_output.txt",
    sep=" ",
    names=["acc", "snp_pos", "aa_ref", "aa_mut", "pathogenic"],
)
df_clinvar["tm"] = 0
df_clinvar["aa_ref"].replace(three2one, inplace=True)
df_clinvar["aa_mut"].replace(three2one, inplace=True)
df_clinvar["id"] = np.nan
df_clinvar["snp_id"] = df_clinvar.index
df_clinvar["snp_id"] = "CLI_" + df_clinvar["snp_id"].astype(str)
df_clinvar["gene"] = np.nan
df_clinvar["snp_rs"] = np.nan
df_clinvar = df_clinvar[
    [
        "acc",
        "id",
        "gene",
        "snp_id",
        "snp_rs",
        "aa_ref",
        "aa_mut",
        "snp_pos",
        "pathogenic",
    ]
]

for i in df_clinvar.index:
    snp_data = df_clinvar.loc[i]
    df_prot_tms = df_tm[df_tm["acc"] == snp_data["acc"]]
    for j, row in df_prot_tms.iterrows():
        if row.tm_start <= snp_data["snp_pos"] and row.tm_final >= snp_data["snp_pos"]:
            # print(row,  snp_data['snp_pos'])
            df_clinvar.loc[i, "tm"] = 1
            break

# Filter TM only and update the database
df_clinvar_tm = df_clinvar[df_clinvar["tm"] == 1]
del df_clinvar_tm["tm"]
df_clinvar_tm.to_csv("clinvar_tm.csv")
df_clinvar_tm.to_sql("snps", engine, if_exists="append", index=False)
connection.close()

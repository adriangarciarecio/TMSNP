#!/usr/bin/env python

# Processs GNOMAD data (after Mireia's scripts)
# Takes the TMs from the MySQL tabe

import mysql.connector
import sqlalchemy as sql
import os
import pandas as pd
import numpy as np
import sys

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
select = int(sys.argv[1])  # Selection mode
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

if select == 1:
    print("> Creating snps_all!")
    try:
        connection.execute("drop table if exists snps_all;")
        connection.execute(
            "CREATE TABLE snps_all AS SELECT * FROM snps;"
        )  # Use the table snps as base
    except:
        print("> Error!")

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

df_tm = pd.read_sql("select * from tm_segments;", engine)
df_gnomad = pd.read_csv(
    "gnomad.txt", sep="\t", names=["acc", "snp_pos", "aa_ref", "aa_mut", "gnomad_freq"]
)
df_gnomad["tm"] = 0
df_gnomad["aa_ref"].replace(three2one, inplace=True)
df_gnomad["aa_mut"].replace(three2one, inplace=True)
df_gnomad["id"] = np.nan
df_gnomad["snp_id"] = df_gnomad.index

df_gnomad["snp_id"] = "BIG_" + df_gnomad["snp_id"].astype(str)
df_gnomad["gene"] = np.nan
df_gnomad["snp_rs"] = np.nan
df_gnomad["pathogenic"] = 0

total_mut = len(df_gnomad)
for i in df_gnomad.index:  # range(3000):  for testing
    if i % 1000 == 0:
        print(str(round(i / total_mut * 100, 2)) + " Code: " + str(df_gnomad["acc"][i]))
    snp_data = df_gnomad.loc[i]
    df_prot_tms = df_tm[df_tm["acc"] == snp_data["acc"]]
    for j, row in df_prot_tms.iterrows():
        if row.tm_start <= snp_data["snp_pos"] and row.tm_final >= snp_data["snp_pos"]:
            # print(row, snp_data["snp_pos"])
            df_gnomad.loc[i, "tm"] = 1
            break

# Filter TM only and update the database
df_gnomad_tm = df_gnomad[df_gnomad["tm"] == 1]
del df_gnomad_tm["tm"]

# Add pathogenic column
try:
    if select == 0:
        connection.execute(
            "alter table snps add column gnomad_freq float(20,10) after pathogenic;"
        )
    if select == 1:
        connection.execute(
            "alter table snps_all add column gnomad_freq float(20,10) after pathogenic;"
        )
except:
    pass

# df_gnomad_tm.to_csv('gnomad_tm.csv')

if select == 0:
    print("> Modifying snps table!")
    df_gnomad_tm.to_sql("snps", engine, if_exists="append", index=False, chunksize=1000)

if select == 1:
    print("> Modifying snps_all table!")
    df_gnomad_tm.to_sql(
        "snps_all", engine, if_exists="append", index=False, chunksize=1000
    )

# in case of problems:
# delete * from snps where gnomad_freq is not NULL;
connection.close()

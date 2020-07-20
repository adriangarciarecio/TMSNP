#!/usr/bin/env python

# Annotate binding site, metal binding sites... based on the Uniprot

# sites.txt is not included since it was empty

# IT ALSO FIXES -1000 values in freq_ref freq_mut and entropy in snp_eval

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
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()
select = int(sys.argv[1])  # Selection mode
if select == 0:
    data_name = "snp_eval"
if select == 1:
    data_name = "snp_eval_all"

annotations = [
    "act_site",
    "binding_site",
    "metal",
    "mod_res",
    "carbohyd",
    "lipid",
    "mutagen",
]
annotations_2res = ["disulfide", "ca_bind", "dna_bind"]


total = len(annotations + annotations_2res)
for ii, annotation in enumerate(annotations):
    print(f"{annotation} ({ii+1}/{total})")
    try:
        connection.execute(f"alter table {data_name} add column {annotation} int(1);")
    except:
        pass

    df = pd.read_csv(f"{annotation}.txt", sep="\t", header=None)
    df.columns = ["acc", "snp_pos"]
    for i in df.index:
        out = f"""update {data_name} set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos']};"""
        connection.execute(out)

for j, annotation in enumerate(annotations_2res):
    print(f"{annotation} ({j+ii+2}/{total})")
    try:
        connection.execute(f"alter table {data_name} add column {annotation} int(1)")
    except:
        pass

    df = pd.read_csv(f"{annotation}.txt", sep="\t", header=None)
    df.columns = ["acc", "snp_pos1", "snp_pos2"]
    for i in df.index:
        out = f"""update {data_name} set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos1']};"""
        connection.execute(out)
        out = f"""update {data_name} set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos2']};"""
        connection.execute(out)

for annotation in annotations + annotations_2res:
    out = f"""update {data_name} set {annotation}=0 where {annotation} is NULL"""
    connection.execute(out)


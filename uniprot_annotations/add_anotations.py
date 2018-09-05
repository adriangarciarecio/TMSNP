#!/usr/bin/env python

# Annotate binding site, metal binding sites... based on the Uniprot

# sites.txt is not included since it was empty

# IT ALSO FIXES -1000 values in freq_ref freq_mut and entropy in snp_eval

import mysql.connector
import sqlalchemy
import os
import pandas as pd
import numpy as np

# Set the connection
conn = mysql.connector.connect(host='alf03.uab.cat', user='lmcdb',
                               password=os.getenv('LMCDB_PASS'), database='tmsnp')
mycursor = conn.cursor()

annotations = ['act_site', 'binding_site', 'metal', 'mod_res', 'carbohyd', 'lipid', 'mutagen']
annotations_2res = ['disulfide', 'ca_bind', 'dna_bind']


total = len(annotations + annotations_2res)
for ii, annotation in enumerate(annotations):
    print(f'{annotation} ({ii+1}/{total})') 
    try:
        mycursor.execute(f'alter table snp_eval add column {annotation} int(1)')
    except:
        pass

    df = pd.read_csv(f'{annotation}.txt', sep='\t', header=None)
    df.columns = ['acc', 'snp_pos']
    for i in df.index:
        out = f"""update snp_eval set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos']};"""
        mycursor.execute(out)
        conn.commit()

for j, annotation in enumerate(annotations_2res):
    print(f'{annotation} ({j+ii+2}/{total})')
    try:
        mycursor.execute(f'alter table snp_eval add column {annotation} int(1)')
    except:
        pass

    df = pd.read_csv(f'{annotation}.txt', sep='\t', header=None)
    df.columns = ['acc', 'snp_pos1', 'snp_pos2']
    for i in df.index:
        out = f"""update snp_eval set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos1']};"""
        mycursor.execute(out)
        conn.commit()
        out = f"""update snp_eval set {annotation}=1 where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['snp_pos2']};"""
        mycursor.execute(out)
        conn.commit()

for annotation in annotations + annotations_2res:
        out = f"""update snp_eval set {annotation}=0 where {annotation} is NULL"""
        mycursor.execute(out)
        conn.commit()

# replace -1000 by NULL
#columns = ['freq_ref', 'freq_mut', 'entropy']
#for column in columns:
#        out = f"""update snp_eval set {column}=NULL where {column}=-1000"""
#        mycursor.execute(out)
#        conn.commit()



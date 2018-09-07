#!/usr/bin/env python

# Takes the TMs from the MySQL tabe

import mysql.connector
import sqlalchemy
import os
import pandas as pd
import numpy as np

# Set the connection
conn = mysql.connector.connect(host='alf03.uab.cat', user='lmcdb',
                               password=os.getenv('LMCDB_PASS'), database='tmsnp')
mycursor = conn.cursor()

try:
    mycursor.execute('alter table snp_eval add column pph2 int(1)')
except:
    pass
try:
    mycursor.execute('alter table snp_eval add column sift int(1)')
except:
    pass


df = pd.read_csv('pph2-short.txt', sep='\t')
df.columns = df.columns.str.strip() # trim spaces in column names
df = df.applymap(lambda x: x.strip() if type(x) is str else x) # trim spaces in all fields
df['prediction'] = df['prediction'].str.replace('probably damaging', '1')
df['prediction'] = df['prediction'].str.replace('possibly damaging', '1')
df['prediction'] = df['prediction'].str.replace('benign', '0')
df.prediction = df.prediction.fillna(9).astype(int) # 9 means null is a trick for the convertion to int
df.pos = df.pos.fillna(0).astype(int)

for i in df.index:
    if i % 100 == 0: print(i)
    out = f"""update snp_eval set pph2={df.loc[i]['prediction']} where acc="{df.loc[i]['acc']}" and snp_pos={df.loc[i]['pos']} and aa_mut="{df.loc[i]['aa2']}";"""
    #print(out)
    mycursor.execute(out)
    conn.commit()



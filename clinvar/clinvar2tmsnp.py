#!/usr/bin/env python

# Processs GNOMAD data (after Mireia's scripts)
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
mycursor.execute("select * from tm_segments")

three2one = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
    'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 'Gly': 'G',
    'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 'Ala': 'A', 'Val':'V',
    'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

table_rows = mycursor.fetchall()
df_tm = pd.DataFrame(table_rows, columns=['acc', 'ini', 'end'])
df_clinvar = pd.read_csv('clinvar_output.txt', sep=' ', names=['acc', 'snp_pos', 'aa_ref', 'aa_mut', 'pathogenic'])
df_clinvar['tm'] = 0
df_clinvar['aa_ref'].replace(three2one, inplace=True)
df_clinvar['aa_mut'].replace(three2one, inplace=True)
df_clinvar['id'] = np.nan
df_clinvar['snp_id'] = df_clinvar.index
df_clinvar['snp_id'] = 'CLI_' + df_clinvar['snp_id'].astype(str)
df_clinvar['gene'] = np.nan
df_clinvar['snp_rs'] = np.nan
df_clinvar = df_clinvar[['acc', 'id', 'gene', 'snp_id', 'snp_rs', 'aa_ref', 'aa_mut', 'snp_pos', 'pathogenic']]

for i in df_clinvar.index:
    snp_data = df_clinvar.loc[i]
    df_prot_tms = df_tm[df_tm['acc'] == snp_data['acc']]
    for j, row in df_prot_tms.iterrows():
        if row.ini <= snp_data['snp_pos'] and row.end >= snp_data['snp_pos']:
           #print(row,  snp_data['snp_pos'])
           df_clinvar.loc[i, 'tm'] = 1
           break

# Filter TM only and update the database
df_clinvar_tm = df_clinvar[df_clinvar['tm']==1]
del df_clinvar_tm['tm']
df_clinvar_tm.to_csv('clinvar_tm.csv')
engine = sqlalchemy.create_engine(f'mysql+mysqlconnector://lmcdb:{os.getenv("LMCDB_PASS")}@alf03.uab.cat/tmsnp')
df_clinvar_tm.to_sql('snps', engine, if_exists='append', index=False)


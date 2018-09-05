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

mycursor.execute("select acc, id, pfam, snp_id, snp_rs, aa_ref, aa_mut, snp_pos, subs_mat, freq_ref, freq_mut, entropy, pathogenic from snp_eval")


table_rows = mycursor.fetchall()
df = pd.DataFrame(table_rows, columns=['acc', 'id', 'pfam', 'snp_id', 'snp_rs', 'aa_ref', 'aa_mut', 'snp_pos', 'subs_mat', 'freq_ref', 'freq_mut', 'entropy', 'pathogenic'])
# Fix spaces in amino acids? This should go in previous scripts

# Sort; Uniprot should be the first as we want to keep it in case of duplicates
# IDs start with  Uniprot: VAR_, GNOMAD: GNO_, ClinVar: CLI:
df = df.sort_values(['acc', 'snp_pos', 'snp_id'], ascending=False)
df_unique = df.drop_duplicates(subset=['acc', 'snp_pos', 'aa_mut'], keep='first')
df_unique = df_unique.sort_values(['acc', 'snp_pos', 'snp_id'])

engine = sqlalchemy.create_engine(f'mysql+mysqlconnector://lmcdb:{os.getenv("LMCDB_PASS")}@alf03.uab.cat/tmsnp')
df_unique.to_sql('snp_eval', engine, if_exists='replace', index=False, chunksize=10000)

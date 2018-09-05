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
mycursor.execute("select acc, id, gene, snp_id, snp_rs, aa_ref, aa_mut, snp_pos, pathogenic from snps")

table_rows = mycursor.fetchall()
df = pd.DataFrame(table_rows, columns=['acc', 'id', 'gene', 'snp_id', 'snp_rs', 'aa_ref', 'aa_mut', 'snp_pos', 'pathogenic'])
# Fix spaces in amino acids? This should go in previous scripts
df['aa_ref'] = df['aa_ref'].str.strip()
df['aa_mut'] = df['aa_mut'].str.strip()
df = df[df['aa_mut'] !='dup']  # remove duplications
df = df[df['aa_ref'] !='MISS']
df['aa_mut'] = df['aa_mut'].str.replace('MISS', '-') # missing residues
df['aa_mut'] = df['aa_mut'].str.replace('del', '-') # indel
df['aa_mut'] = df['aa_mut'].str.replace('Ter', '/') # indel

# Sort; Uniprot should be the first as we want to keep it in case of duplicates
# IDs start with  Uniprot: VAR_, GNOMAD: BIG_ (BroadInstitute GNOMAD), ClinVar: CLI
# The order we want is Uniprot > Clinvar > GNOMAD (reverse order)
df = df.sort_values(['acc', 'snp_pos', 'snp_id'], ascending=False)
df_unique = df.drop_duplicates(subset=['acc', 'snp_pos', 'aa_mut'], keep='first')
df_unique = df_unique.sort_values(['acc', 'snp_pos', 'snp_id'])

engine = sqlalchemy.create_engine(f'mysql+mysqlconnector://lmcdb:{os.getenv("LMCDB_PASS")}@alf03.uab.cat/tmsnp')
df_unique.to_sql('snps', engine, if_exists='replace', index=False, chunksize=10000)

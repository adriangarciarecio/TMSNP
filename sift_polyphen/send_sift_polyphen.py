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
mycursor.execute("select acc, snp_pos, aa_ini, aa_fin from snps")

table_rows = mycursor.fetchall()
df = pd.DataFrame(table_rows, columns=['acc', 'snp_pos', 'aa_ini', 'aa_fin'])
df.to_csv('to_polyphen.txt', sep=' ', index=False, header=False)



#engine = sqlalchemy.create_engine(f'mysql+mysqlconnector://lmcdb:{os.getenv("LMCDB_PASS")}@alf03.uab.cat/tmsnp')
#df_unique.to_sql('snps', engine, if_exists='replace', index=False)
print('Polyphen2: http://genetics.bwh.harvard.edu/pph2/bgi.shtml')
print('SIFT: http://sift.jcvi.org/protein_batch_submit.php?species=human')

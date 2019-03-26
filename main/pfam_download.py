#!/usr/bin/env python3

import os
from prody import *
import sqlalchemy as sql
import pandas as pd
from iker_snp import *

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

path_pfam = "./pfam/pfam_download/"
l_ok = os.listdir(path_pfam)
l_fail = []
data = pd.read_sql("select distinct pfam from receptor_pfam;", engine)
l_pfam = data.pfam
i = 1
for pfam in l_pfam:
    if pfam not in l_ok:
        print(f"> Downloading pfam {pfam} ({((len(l_ok)+i)/len(l_pfam))*100} %)!")
        try:
            fetchPfamMSA(pfam, alignment="full", format="fasta", folder=path_pfam)
            os.rename(f"{path_pfam}{pfam}_full.fasta", f"{path_pfam}{pfam}")
        except:
            print(f"Download fail")
            l_fail.append(pfam)
        i += 1
print(l_fail)

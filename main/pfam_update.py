import glob
import numpy as np
import math
import pandas as pd
import mysql.connector
import sqlalchemy as sql
import os
from bs4 import BeautifulSoup
import requests
import re

# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )

engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")

data = pd.read_sql("select * from receptor_pfam;", engine)
acc = data.acc
pfam = data.pfam
dic_ap = {}
new_acc = []
new_pfam = []
acc = data.acc.unique()
print("> Getting PFAM codes for accession!")
for i, a in enumerate(acc):
    print("> Progress: " + str((i * 100) / len(acc)) + "%")
    url = "https://pfam.xfam.org/protein/" + a + "#tabview=tab0"
    r = requests.get(url)
    soup = BeautifulSoup(r.text, "html.parser")
    b = soup.find_all("td", text="Pfam")
    for element in b:
        ele = str(element).split('"')
        ele = ele[1].split("_")
        if a not in dic_ap.keys():
            dic_ap[a] = [ele[1]]
        if ele[1] not in dic_ap[a]:
            dic_ap[a].append(ele[1])
# We have dictionary with each accession with all pfam contained in, now we need to extract two list ACC/PFAM
for key in dic_ap.keys():
    pfams = dic_ap[key]
    for p in pfams:
        new_acc.append(key)
        new_pfam.append(p)

pfam_codes = pd.DataFrame({"acc": new_acc, "pfam": new_pfam})

connection = engine.connect()

print("> Updating PFAM codes on MySQL")
pfam_codes.to_sql("receptor_pfam", engine, if_exists="replace", index=False)
connection.close()

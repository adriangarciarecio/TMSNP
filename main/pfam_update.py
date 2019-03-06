import glob
import numpy as np
import math
import pandas as pd
import mysql.connector
import sqlalchemy
import os
from bs4 import BeautifulSoup
import requests
import re

engine = sqlalchemy.create_engine(
    f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
)

data = pd.read_sql("select * from receptor_pfam;", engine)
acc = data.acc
pfam = data.pfam
new_pfam = list(pfam)
new_acc = list(acc)
acc = data.acc.unique()
print("> Getting PFAM codes for accession!")
for i, a in enumerate(acc):
    print("> Progress: " + str((i * 100) / len(acc)) + "%")
    url = "https://pfam.xfam.org/protein/" + a + "#tabview=tab0"
    r = requests.get(url)
    soup = BeautifulSoup(r.text, "html.parser")
    b = soup.find_all("td", text="Pfam")
    if len(b) != 1:
        for element in b:
            ele = str(element).split('"')
            ele = ele[1].split("_")
            if ele[1] not in new_pfam:
                new_pfam.append(str(ele[1]))
                new_acc.append(a)

    else:
        b = str(b).split('"')
        b = b[1].split("_")
        if b[1] not in new_pfam:
            new_pfam.append(str(b[1]))
            new_acc.append(a)

pfam_codes = pd.DataFrame({"acc": new_acc, "pfam": new_pfam})

connection = engine.connect()

print("> Updating PFAM codes on MySQL")
pfam_codes.to_sql("receptor_pfam", engine, if_exists="replace", index=False)
connection.close()

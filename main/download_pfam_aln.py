#!/usr/bin/env python3

import requests
import mysql.connector
import os
import subprocess
import sqlalchemy as sql
import pandas as pd
from iker_snp import *

###seleccionar els diferents pfam de les taules de mysql#####################3

path_pfam = "./pfam/pfam_download/"

if not os.path.exists(path_pfam):
    os.makedirs(path_pfam)

path_pfam_sys = path_pfam[:-1]

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

print("Downloading data from PFAM...")

data = pd.read_sql("select distinct pfam from receptor_pfam;", engine)
l_pfam = data.pfam

# curl
total_pfams = len(l_pfam)
for i, pfam_key in enumerate(l_pfam):
    print(pfam_key, "remaining families:", total_pfams - i)
    link = (
        "http://pfam.xfam.org/family/"
        + pfam_key
        + "/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0"
    )
    run_curl(link, path_pfam + pfam_key)

print("... Done!!!")
################################################################################
# Manage errors

print("Trying to get alignments that were not downloaded...")
abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
l_pfam_err = ["control"]
c_down = 0

while len(l_pfam_err) != 0:
    if c_down < 20:
        pfam_list = os.listdir(path_pfam_sys)
        l_pfam_err = list()
        for pfam in pfam_list:
            l_data = file_to_lines(path_pfam + pfam)
            for align in l_data:
                if align[0:9] == "<!DOCTYPE":
                    print("Error: " + pfam)
                    l_pfam_err.append(pfam)
        print(l_pfam_err)

        for pfam_key in l_pfam_err:
            os.remove(path_pfam + pfam_key)
            link = (
                "http://pfam.xfam.org/family/"
                + pfam_key
                + "/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0"
            )
            run_curl(link, path_pfam + pfam_key)

        c_down += 1
        print(c_down)
    else:
        for bad_pfam in l_pfam_err:
            print("filename:", path_pfam + bad_pfam)
            # Write a file with an error message instead of removing the file
            # outf1 = open(path_pfam + bad_pfam, 'w')
            # print('ERROR_ERROR/69-117             ----------------------------------', file=outf1)
            # outf1.close()
            os.remove(path_pfam + bad_pfam)
            print("Removed by error while the pfam file is downloading: " + bad_pfam)
        break

l_pfam_err = ["control"]
c_down = 0

while len(l_pfam_err) != 0:
    if c_down < 5:
        pfam_list = os.listdir(path_pfam_sys)
        l_pfam_err = list()
        for pfam in pfam_list:
            l_data = file_to_lines(path_pfam + pfam)
            for align in l_data:
                if abc.find(align[0]) == -1:
                    if pfam not in l_pfam_err:
                        l_pfam_err.append(pfam)
                        print("Error: " + pfam)
        print(l_pfam_err)

        for pfam_key in l_pfam_err:
            os.remove(path_pfam + pfam_key)
            link = (
                "http://pfam.xfam.org/family/"
                + pfam_key
                + "/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0"
            )
            run_curl(link, path_pfam + pfam_key)

        c_down += 1
        print(c_down)
    else:
        for bad_pfam in l_pfam_err:
            os.remove(path_pfam + bad_pfam)
            outf1 = open(path_pfam + bad_pfam, "w")
            print(
                "ERROR_ERROR/69-117             ----------------------------------",
                file=outf1,
            )
            outf1.close()
            print(
                "Removend and created in blank. Error found while the file was downloading from server: "
                + bad_pfam
            )
        break
connection.close()

#!/usr/bin/env python3

# Generates the MySQL tables with the SNP and related information for the Membrane proteins in the Uniprot

# INFO: some receptors appear twice (or more) when they have multiple pfam domains assigned.

import requests
import sqlalchemy as sql
import pandas as pd
import os
from bs4 import BeautifulSoup
import glob
import numpy as np
import math
import re
from iker_snp import *

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

# Dirs where PFAM downloaded files will go
path_pfam = "./pfam/pfam_download/"
path_raw = "./pfam/path_raw/"

dirs = ["./pfam/pfam_download/", "./pfam/path_raw/"]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

# paths without the ending /
path_pfam_sys = path_pfam[:-1]
path_raw_sys = path_raw[:-1]

################### MySQL TABLES CREATOR ##############
print("Preparing MySQL tables based on UniProt data ...")
connection.execute("DROP TABLE IF EXISTS receptor_pfam")
connection.execute("DROP TABLE IF EXISTS tm_segments")
connection.execute("DROP TABLE IF EXISTS snps")

connection.execute("CREATE TABLE receptor_pfam (acc varchar(20), pfam varchar(20));")
connection.execute(
    """CREATE TABLE snps (acc varchar(20), id varchar(30), gene varchar(20), snp_id varchar(30), 
                    snp_rs varchar(30), aa_ref varchar(30), aa_mut varchar(30), snp_pos int, tm_pos int,
                     pathogenic int, PRIMARY KEY (snp_id));"""
)
connection.execute(
    "CREATE TABLE tm_segments (acc varchar(20),  tm_start int, tm_final int);"
)

# All Human TM proteins tagged as reviewed
url_uniprot = "http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atransmem%29+AND+organism%3A%22Homo+sapiens+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score&format=txt"
req = requests.get(url_uniprot)
l_data = req.text.splitlines()

data_prot = extract_prot(l_data)
data_snp_loc = extract_snp_regions(l_data)
l_snps_transm = descart_no_transm(data_snp_loc)

for n_pos in range(len(data_prot)):
    for n_pfam in range(len(data_prot[n_pos][3])):
        acc_prot = data_prot[n_pos][1]
        #        id_protein = data_prot[n_pos][0]
        #        gen_protein = data_prot[n_pos][2]
        pfam_prot = data_prot[n_pos][3][n_pfam]
        connection.execute(
            ("INSERT INTO receptor_pfam VALUES" + str((acc_prot, pfam_prot)))
        )

for n_pos in range(len(l_snps_transm)):
    for n_snps in range(len(l_snps_transm[n_pos][2])):
        acc_prot = l_snps_transm[n_pos][0]
        var_snip = l_snps_transm[n_pos][2][n_snps]
        chang_snp = l_snps_transm[n_pos][3][n_snps]
        snp_ini_sym = chang_snp[: chang_snp.find("->")]
        snp_end_sym = chang_snp[chang_snp.find("->") + 3 :]
        rs_snip = l_snps_transm[n_pos][4][n_snps]
        pos_snp = int(l_snps_transm[n_pos][5][n_snps])
        rm = l_snps_transm[n_pos][6][n_snps]
        dis = l_snps_transm[n_pos][7][n_snps]
        id = l_snps_transm[n_pos][8]
        gen = l_snps_transm[n_pos][9]
        connection.execute(
            (
                "INSERT INTO snps VALUES"
                + str(
                    (
                        acc_prot,
                        id,
                        gen,
                        var_snip,
                        rs_snip,
                        snp_ini_sym,
                        snp_end_sym,
                        pos_snp,
                        rm,
                        dis,
                    )
                )
            )
        )

for n_pos in range(len(l_snps_transm)):
    for n_reg in range(len(l_snps_transm[n_pos][1])):
        acc_prot = l_snps_transm[n_pos][0]
        reg_ini = int(l_snps_transm[n_pos][1][n_reg][:8])
        reg_end = int(l_snps_transm[n_pos][1][n_reg][9:])
        connection.execute(
            ("INSERT INTO tm_segments VALUES" + str((acc_prot, reg_ini, reg_end)))
        )

# Delete non-TM
connection.execute("delete from snps where tm_pos=0;")

# Delete tm column
connection.execute("alter table snps drop column tm_pos;")

connection.close()

print("... Done!!!")

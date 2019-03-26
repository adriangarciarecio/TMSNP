#!/usr/bin/env python3

################################################
import numpy as np
import numpy
import requests
import collections
from math import exp, expm1, log10, log
import subprocess
import mysql.connector
import os
import collections
import sqlalchemy as sql
import pandas as pd
from iker_snp import *

path_pfam = "./pfam/pfam_download/"

# Set the connection
# engine = sql.create_engine(
#     f"mysql://lmcdb:{os.getenv('LMCDB_PASS')}@alf03.uab.cat/tmsnp"
# )
engine = sql.create_engine(f"mysql://adrian:compmod5@localhost/tmsnp")
connection = engine.connect()

# Delete non-TM
connection.execute("drop table if exists snp_val;")

print("Evaluating SNPs based on the PFAM alignment")

data = pd.read_sql("select distinct acc from snps;", engine)
l_acc = data.acc

# print(l_acc)

###################################################################################
phat_matrix = np.array(
    [
        [5, -6, -2, -5, 1, -3, -5, 1, -3, 0, -1, -7, -1, -1, -3, 2, 0, -4, -3, 1],
        [-6, 9, -3, -7, -8, -2, -6, -5, -4, -6, -6, -1, -6, -7, -7, -6, -6, -7, -6, -7],
        [-2, -3, 11, 2, -2, 2, 0, -1, 4, -3, -3, -2, -2, -1, -4, 1, -1, -5, 2, -3],
        [-5, -7, 2, 12, -7, 0, 6, -2, -1, -5, -5, -5, -5, -5, -5, -4, -5, -7, -4, -5],
        [1, -8, -2, -7, 7, -5, -7, -2, -7, -3, -2, -10, -2, 0, -8, 1, -1, -4, -1, -2],
        [-3, -2, 2, 0, -5, 9, 1, -2, 2, -3, -3, -1, -1, -2, -3, -1, -3, 1, 0, -3],
        [-5, -6, 0, 6, -7, 1, 12, -3, -1, -5, -5, -4, -5, -5, -5, -3, -5, -7, -2, -5],
        [1, -5, -1, -2, -2, -2, -3, 9, -4, -2, -2, -5, -1, -2, -3, 1, -1, -5, -3, -2],
        [-3, -4, 4, -1, -7, 2, -1, -4, 11, -5, -4, -5, -4, -2, -6, -2, -4, -3, 3, -5],
        [0, -6, -3, -5, -3, -3, -5, -2, -5, 5, 2, -7, 3, 0, -4, -2, -1, -4, -3, 3],
        [-1, -6, -3, -5, -2, -3, -5, -2, -4, 2, 4, -7, 2, 1, -5, -2, -1, -3, -2, 1],
        [
            -7,
            -1,
            -2,
            -5,
            -10,
            -1,
            -4,
            -5,
            -5,
            -7,
            -7,
            5,
            -6,
            -7,
            -4,
            -5,
            -6,
            -8,
            -4,
            -8,
        ],
        [-1, -6, -2, -5, -2, -1, -5, -1, -4, 3, 2, -6, 6, 0, -5, -2, 0, -4, -2, 1],
        [-1, -7, -1, -5, 0, -2, -5, -2, -2, 0, 1, -7, 0, 6, -5, -2, -2, 0, 4, -1],
        [
            -3,
            -7,
            -4,
            -5,
            -8,
            -3,
            -5,
            -3,
            -6,
            -4,
            -5,
            -4,
            -5,
            -5,
            13,
            -3,
            -4,
            -6,
            -5,
            -4,
        ],
        [2, -6, 1, -4, 1, -1, -3, 1, -2, -2, -2, -5, -2, -2, -3, 6, 1, -5, -2, -2],
        [0, -6, -1, -5, -1, -3, -5, -1, -4, -1, -1, -6, 0, -2, -4, 1, 3, -7, -3, 0],
        [-4, -7, -5, -7, -4, 1, -7, -5, -3, -4, -3, -8, -4, 0, -6, -5, -7, 11, 1, -4],
        [-3, -6, 2, -4, -1, 0, -2, -3, 3, -3, -2, -4, -2, 4, -5, -2, -3, 1, 11, -3],
        [1, -7, -3, -5, -2, -3, -5, -2, -5, 3, 1, -8, 1, -1, -4, -2, 0, -4, -3, 4],
    ]
)

dict_matrix = {
    "A": 0,
    "R": 1,
    "N": 2,
    "D": 3,
    "C": 4,
    "Q": 5,
    "E": 6,
    "G": 7,
    "H": 8,
    "I": 9,
    "L": 10,
    "K": 11,
    "M": 12,
    "F": 13,
    "P": 14,
    "S": 15,
    "T": 16,
    "W": 17,
    "Y": 18,
    "V": 19,
}

################################################################################
# Not used until part II but run at start in case of connection problems
phat_max = True
# All Human TM proteins tagged as reviewed
print("... Getting list of TM proteins from the Uniprot")
url_uniprot = "http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atransmem%29+AND+organism%3A%22Homo+sapiens+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score&format=txt"
req = requests.get(url_uniprot)
print("... Got list of TM proteins!")
l_data = req.text.splitlines()  # before: l_data = file_to_lines(uniprot_filename)
total_acc = len(l_acc)
# l_acc = ["Q9Y653"]
################################################################################
if phat_max == True:
    ################################################################################
    print("\nPart I: Evaluating subsitutions based in the Phat matrix")

    connection.execute("DROP TABLE IF EXISTS snp_phat")
    connection.execute(
        """CREATE TABLE snp_phat (acc varchar(20), id varchar(30), gene varchar(20),
                        snp_id varchar(30), snp_rs varchar(30), aa_ref varchar(30),
                        aa_mut varchar(30), snp_pos int, pathogenic int,
                        gnomad_freq float(20,10), subs_mat int);"""
    )
    ################################################################################
    # for testing: truncate l_acc
    for i, acc in enumerate(l_acc):
        print(f"Processing {acc} ({i+1} / {total_acc})")
        data = pd.read_sql(
            "select tm_start, tm_final from tm_segments where acc = '" + acc + "';",
            engine,
        )
        reg_start = data.tm_start
        reg_start = list(reg_start)
        reg_end = data.tm_final
        reg_end = list(reg_end)

        data = pd.read_sql(
            "select distinct snp_id from snps where acc = '" + acc + "';", engine
        )
        l_var = data.snp_id
        l_var = list(l_var)

        for var in l_var:
            data = pd.read_sql(
                "select aa_ref from snps where snp_id = '" + var + "';", engine
            )
            aa_ref = data.aa_ref
            aa_ref = list(aa_ref)
            aa_ref = aa_ref[0]

            data = pd.read_sql(
                "select aa_mut from snps where snp_id = '" + var + "';", engine
            )
            aa_mut = data.aa_mut
            aa_mut = list(aa_mut)
            aa_mut = aa_mut[0]

            data = pd.read_sql(
                "select id from snps where snp_id = '" + var + "';", engine
            )
            id = data.id
            id = list(id)
            id = id[0]

            data = pd.read_sql(
                "select gene from snps where snp_id = '" + var + "';", engine
            )
            gene = data.gene
            gene = list(gene)
            gene = gene[0]

            data = pd.read_sql(
                "select snp_rs from snps where snp_id = '" + var + "';", engine
            )
            rs = data.snp_rs
            rs = list(rs)
            rs = rs[0]

            data = pd.read_sql(
                "select snp_pos from snps where snp_id = '" + var + "';", engine
            )
            snp_pos = data.snp_pos
            snp_pos = list(snp_pos)
            snp_pos = snp_pos[0]

            tm_okey = False
            for i, pos in enumerate(reg_start):
                if snp_pos >= pos and snp_pos <= reg_end[i]:
                    tm_okey = True
                else:
                    continue

            data = pd.read_sql(
                "select pathogenic from snps where snp_id = '" + var + "';", engine
            )
            patho = data.pathogenic
            patho = list(patho)
            patho = patho[0]

            data = pd.read_sql(
                "select gnomad_freq from snps where snp_id = '" + var + "';", engine
            )
            gno_freq = data.gnomad_freq
            gno_freq = list(gno_freq)
            gno_freq = gno_freq[0]

            if tm_okey == True:
                ## Changes here
                if len(aa_ref) >= 2 or len(aa_mut) >= 2 or aa_mut == "/":
                    score = None  # -1000
                else:
                    # print(aa_ref, aa_mut)
                    score = phat_matrix[dict_matrix[aa_ref]][dict_matrix[aa_mut]]

                if type(score) == numpy.int64:
                    score = int(score)
                connection.execute(
                    "INSERT INTO snp_phat VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                    (
                        acc,
                        id,
                        gene,
                        var,
                        rs,
                        aa_ref,
                        aa_mut,
                        snp_pos,
                        patho,
                        gno_freq,
                        score,
                    ),
                )

    ###############################################################################

print("\nPart II: Extracting frequencies from the PFAM alignment and computing entropy")

no_pfam = set()
# Entropy

connection.execute("DROP TABLE IF EXISTS snp_eval")
connection.execute(
    """CREATE TABLE snp_eval (acc varchar(20), id varchar(30),
    pfam varchar(20), snp_id varchar(30), snp_rs varchar(30),
    aa_ref varchar(30), aa_mut varchar(30), snp_pos int, gnomad_freq float(20,10), subs_mat int,
     freq_ref float, freq_mut float, entropy float, pathogenic int);"""
)

## entropia inicial a parti de les regions transmem del arxiu total d'uniprot

data_prot_entro = ext_prot_entro(l_data)
data_snp_reg_entro = create_l_prot_entro(data_prot_entro)

data_snp_reg = snp_reg_entro(l_data)
l_snp_tm_reg = rm_no_stm_entro(data_snp_reg)

adn_all = ""
for l_n_pos in range(len(l_snp_tm_reg)):
    for n_reg in range(len(l_snp_tm_reg[l_n_pos][1])):
        acc_prot = l_snp_tm_reg[l_n_pos][0]
        reg_ini = int(l_snp_tm_reg[l_n_pos][1][n_reg][:8])
        reg_end = int(l_snp_tm_reg[l_n_pos][1][n_reg][9:])
        for l_pos_adn in range(len(data_snp_reg_entro)):
            if data_snp_reg_entro[l_pos_adn][0] == acc_prot:
                adn_all += data_snp_reg_entro[l_pos_adn][1][reg_ini - 1 : reg_end]

d_freq_abs = dict(
    collections.Counter(adn_all)
)  # amb els 2 aminoacids sense res mes ni x ni Z ni B

d_freq_rel = dict()
all_aa = 0
for aa in d_freq_abs:
    all_aa += d_freq_abs[aa]

for aa in d_freq_abs:
    d_freq_rel[aa] = d_freq_abs[aa] / all_aa

suma = 0.0
for d_aa in d_freq_rel:
    suma += d_freq_rel[d_aa] * log(d_freq_rel[d_aa])

entro_ini = abs(suma)
## calcul de la frequencia i entropia de la posicio on cau l'snip###############
data = pd.read_sql("select distinct acc from snps;", engine)
total_acc_list = data.acc
total_acc_list = list(total_acc_list)

data = pd.read_sql("select distinct acc from snp_phat;", engine)
l_acc = data.acc
l_acc = list(l_acc)

for i, accs in enumerate(l_acc):
    print(f"Processing {accs} ({i+1} / {total_acc})")
    data = pd.read_sql(
        "select tm_start, tm_final from tm_segments where acc = '" + accs + "';", engine
    )
    reg_start = data.tm_start
    reg_start = list(reg_start)
    reg_end = data.tm_final
    reg_end = list(reg_end)

    data = pd.read_sql(
        "select distinct pfam from receptor_pfam where acc = '" + accs + "';", engine
    )
    l_pfam = data.pfam
    l_pfam = list(l_pfam)
    data = pd.read_sql(
        "select distinct id from snp_phat where acc = '" + accs + "';", engine
    )
    id = data.id
    id = list(id)
    id = id[0]

    data = pd.read_sql(
        "select distinct snp_id from snp_phat where acc = '" + accs + "';", engine
    )
    l_var = data.snp_id
    l_var = list(l_var)

    for snp_var in l_var:
        data = pd.read_sql(
            "select snp_rs from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        snp_rs = data.snp_rs
        snp_rs = list(snp_rs)
        snp_rs = snp_rs[0]

        data = pd.read_sql(
            "select snp_pos from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        snp_pos = data.snp_pos
        snp_pos = list(snp_pos)
        snp_pos = snp_pos[0]

        # Saber si es TM o no
        tm_okey = False
        for i, pos in enumerate(reg_start):
            if snp_pos >= pos and snp_pos <= reg_end[i]:
                tm_okey = True
            else:
                continue

        data = pd.read_sql(
            "select aa_mut from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        aa_mut = data.aa_mut
        aa_mut = list(aa_mut)
        aa_mut = aa_mut[0]

        data = pd.read_sql(
            "select aa_ref from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        aa_ref = data.aa_ref
        aa_ref = list(aa_ref)
        aa_ref = aa_ref[0]

        data = pd.read_sql(
            "select pathogenic from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        patho = data.pathogenic
        patho = list(patho)
        patho = patho[0]

        data = pd.read_sql(
            "select gnomad_freq from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        gno_freq = data.gnomad_freq
        gno_freq = list(gno_freq)
        gno_freq = gno_freq[0]

        data = pd.read_sql(
            "select subs_mat from snp_phat where snp_id = '" + snp_var + "';", engine
        )
        subs_mat = data.subs_mat
        subs_mat = list(subs_mat)
        subs_mat = subs_mat[0]

        if tm_okey == True:
            for pfam_ac in l_pfam:
                try:
                    l_data = file_to_lines(path_pfam + pfam_ac)
                except FileNotFoundError:
                    print(f"No {path_pfam}{pfam_ac}")
                    no_pfam.add(pfam_ac)
                l_prot_pos_aa = create_l_prot(l_data)
                l_prot_pos_total = create_l_variables(l_prot_pos_aa, id, l_var, snp_pos)
                l_prot_rs = create_l_rs(l_prot_pos_total, l_prot_pos_aa)
                freq_entro = cal_freq_entro(l_prot_rs, entro_ini, aa_mut)
                if aa_mut == "MISS":
                    freq_entro = None
                if freq_entro == None:
                    freq_ref = None
                    freq_mut = None
                    entropia = None
                else:
                    freq_ref = float(freq_entro[0])
                    freq_mut = float(freq_entro[1])
                    entropia = float(freq_entro[2])
                if freq_ref != None and freq_mut != None and entropia != None:
                    connection.execute(
                        "INSERT INTO snp_eval VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                        (
                            accs,
                            id,
                            pfam_ac,
                            snp_var,
                            snp_rs,
                            aa_ref,
                            aa_mut,
                            snp_pos,
                            gno_freq,
                            subs_mat,
                            freq_ref,
                            freq_mut,
                            entropia,
                            patho,
                        ),
                    )


# Delete non-TM
# connection.execute("drop table snp_phat;")

print("Proteins with no PFAM alignment:")
for el in no_pfam:
    print(el, end=" ")

print("... Done!")
connection.close()

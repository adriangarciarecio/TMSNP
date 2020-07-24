import os

os.environ["LD_LIBRAY_PATH"] = "/var/www/tmsnp/venv_tmsnp/lib"
import re
import math
import sys
import sqlalchemy
import mysql.connector

import numpy as np
import pandas as pd
import collections
from collections import namedtuple
from .predictor import *
from .iker_snp import *

from .serverside.serverside_table import ServerSideTable
from .serverside import table_schemas

import requests
import numpy

import json

# engine = sqlalchemy.create_engine(f"mysql+mysqlconnector://lmcdb:3usBXLxt5F@alf03.uab.cat/tmsnp")
engine = sqlalchemy.create_engine(
    f"mysql+mysqlconnector://adrian:D1m3rB0w!@localhost:3306/tmsnp"
)
connection = engine.connect()
from math import exp, log10, log
from collections import Counter
from flask import (
    Flask,
    render_template,
    request,
    send_file,
    send_from_directory,
    jsonify,
)

# from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)

# app.config[
#     "SQLALCHEMY_DATABASE_URI"
# ] = "mysql+mysqlconnector://adrian:D1m3rB0w!@localhost:3306/tmsnp"
# db = SQLAlchemy(app)


# class PROCESSING_DATA(db.Model):
#     __teblename__ = "tmsnp"
#     id = db.Column("TICKETNO", db.Integer, primary_key=True)
#     objectid = db.Column("ID", db.Integer)
#     city = db.Column("CITIES", db.String)
#     date = db.Column("DATE", db.DateTime)
#     severity = db.Column("CRITICAL", db.String(2))
#     status = db.Column("STATUS", db.Integer)
#     shortdesc = db.Column("DESC", db.String(20))


###########################################################
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

d_mut = {
    "A": "A_m",
    "R": "R_m",
    "N": "N_m",
    "D": "D_m",
    "C": "C_m",
    "Q": "Q_m",
    "E": "E_m",
    "G": "G_m",
    "H": "H_m",
    "I": "I_m",
    "L": "L_m",
    "K": "K_m",
    "M": "M_m",
    "F": "F_m",
    "P": "P_m",
    "S": "S_m",
    "T": "T_m",
    "W": "W_m",
    "Y": "Y_m",
    "V": "V_m",
}

###########################################################
@app.route("/")
def script_punct():
    return render_template("tmsnp.html")


@app.route("/processing", methods=["GET"])
def processing():
    nothing = pd.read_sql(f"select * from snp_eval_all limit 1;", engine)
    ori_col = nothing.columns
    draw = request.args.get("draw")
    start = request.args.get("start")
    length = request.args.get("length")
    search = request.args.get("search[value]")
    col = request.args.get(f"order[0][column]")
    order = request.args.get(f"order[0][dir]")
    if search == "":
        llista = pd.read_sql(
            f"select * from snp_eval_all ORDER BY {list(ori_col)[int(col)]} {order} LIMIT {length} OFFSET {start};",
            engine,
        )
        rec_filt = 199_025
    else:
        llista = pd.read_sql(
            f"select * from snp_eval_all WHERE {list(ori_col)[0]} REGEXP '^{search}' ORDER BY {list(ori_col)[int(col)]} {order};",
            engine,
        )
        rec_filt = len(llista)
        llista = pd.read_sql(
            f"select * from snp_eval_all WHERE {list(ori_col)[0]} REGEXP '^{search}' ORDER BY {list(ori_col)[int(col)]} {order} limit {length} OFFSET {start};",
            engine,
        )
    columns = [
        "Accession",
        "ID",
        "Code PFAM",
        "SNP ID",
        "Code rs",
        "Reference aminoacid",
        "Mutated aminoacid",
        "Position",
        "Gnomad frequency",
        "Score substitution matrix",
        "Reference frequency",
        "Mutated frequency",
        "Entropy",
        "Disease related",
        "Active site",
        "Binding site",
        "Metal",
        "mod_res",
        "Carbohyd",
        "Lipid",
        "Mutagen",
        "Disulfide",
        "Calcium bind",
        "DNA bind",
        "Polyphen prediction",
        "Sift prediction",
        "Muthtp prediction",
        "Muthtp score",
        "ID_pk",
    ]

    order_columns = [
        "Accession",
        "ID",
        "Code PFAM",
        "SNP ID",
        "Code rs",
        "Reference aminoacid",
        "Mutated aminoacid",
        "Position",
        "Disease related",
        "Gnomad frequency",
        "Score substitution matrix",
        "Reference frequency",
        "Mutated frequency",
        "Entropy",
        "Polyphen prediction",
        "Sift prediction",
        "Muthtp prediction",
        "Muthtp score",
        "ID_pk",
        "Active site",
        "Binding site",
        "Metal",
        "mod_res",
        "Carbohyd",
        "Lipid",
        "Mutagen",
        "Disulfide",
        "Calcium bind",
        "DNA bind",
    ]
    llista.columns = columns
    llista = llista[order_columns]
    # data = llista.to_json(orient="values")
    llista = llista.fillna("-")
    di = {0: "Neutral", 1: "Pathogenic", 2: "Likely pathogenic"}
    llista = llista.replace({"Disease related": di})
    data = llista.to_dict("split")
    data = {
        "draw": draw,
        "recordsTotal": 199_025,
        "recordsFiltered": rec_filt,
        "data": data["data"],
    }
    return jsonify(data)


@app.route("/tmsnpdb")
def tmsnpdb():
    return render_template(
        "tmsnpdb.html"
        # columns=columns,
        # rows=rows,
        # tabledata=total_llista.to_html(
        #     justify="center", classes="em_table", index=False, table_id="table"
        # ),
    )


@app.route("/results", methods=["GET", "POST"])
def results():
    input_seq = request.form["input_seq"].upper()
    if input_seq.find("_") != -1:
        total_llista = pd.read_sql(
            f"select acc, id, pfam, snp_rs, snp_pos, aa_ref, aa_mut, pathogenic from snp_eval_all where id = '{input_seq}';",
            engine,
        )
    else:
        total_llista = pd.read_sql(
            f"select acc, id, pfam, snp_rs, snp_pos, aa_ref, aa_mut, pathogenic from snp_eval_all where acc = '{input_seq}';",
            engine,
        )

    if len(total_llista) == 0:
        return render_template("error_no_acc.html", input_seq=input_seq)
    else:
        name_columns = [
            "Uniprot Accession",
            "Protein ID",
            "PFAM code",
            "SNP_rs",
            "Aminoacid SNP position",
            "Aminoacid reference",
            "Aminoacid mutated",
            "Pathogenic/Non-pathogenic",
        ]
        total_llista.columns = name_columns
        total_llista = total_llista.sort_values(
            by=["Pathogenic/Non-pathogenic"], ascending=False
        )
        total_llista["Pathogenic/Non-pathogenic"] = total_llista[
            "Pathogenic/Non-pathogenic"
        ].map({0: "Non-pathogenic", 1: "Pathogenic", 2: "Likely pathogenic"})
        total_llista.SNP_rs.fillna("-", inplace=True)
        return render_template(
            "results.html",
            input_seq=input_seq,
            tabledata=total_llista.to_html(
                justify="center", classes="em_table", index=False, table_id="table"
            ),
        )


@app.route("/results_pred", methods=["GET", "POST"])
def results_pred():
    dict_aa = {
        "phe": "F",
        "leu": "L",
        "ser": "S",
        "tyr": "Y",
        "cys": "C",
        "trp": "W",
        "pro": "P",
        "his": "H",
        "gln": "Q",
        "arg": "R",
        "ile": "I",
        "met": "M",
        "thr": "T",
        "asn": "N",
        "lys": "K",
        "val": "V",
        "ala": "A",
        "asp": "D",
        "glu": "E",
        "gly": "G",
    }
    path_pfam = "/var/www/tmsnp/tmsnp/pfam/"
    prot = request.form["input_seq_acc"].upper()
    pos = request.form["input_seq_pos"]
    aa_ref = request.form["input_seq_aai"].upper()
    aa_mut = request.form["input_seq_aaf"].upper()
    if aa_ref.lower() in dict_aa.keys():
        aa_ref = dict_aa[aa_ref.lower()]
    if aa_mut.lower() in dict_aa.keys():
        aa_mut = dict_aa[aa_mut.lower()]
    if prot == "" or pos == "" or aa_ref == "" or aa_mut == "":
        return render_template("error_empty.html")
    # Look if the user gives code or id uniprot
    if prot.find("_") != -1:
        data = pd.read_sql(
            "select acc from snp_eval_all where id = '" + prot + "';", engine
        )
    else:
        data = pd.read_sql(
            "select acc from snp_eval_all where acc = '" + prot + "';", engine
        )

    # First need to know if the position is on TM zone
    if data.empty:
        return render_template(
            "error_no_acc.html", input_seq=prot, aa_ref=aa_ref, pos=pos, aa_mut=aa_mut
        )
    else:
        # Patho or Nonpatho:
        if prot.find("_") != -1:
            data = pd.read_sql(
                "select acc from snp_eval where id = '" + prot + "';", engine
            )
        else:
            data = pd.read_sql(
                "select acc from snp_eval where acc = '" + prot + "';", engine
            )
        if data.empty:
            mode = "nonpathogenic"
            if prot.find("_") != -1:
                data = pd.read_sql(
                    "select acc from snp_eval_all where id = '" + prot + "';", engine
                )
            else:
                data = pd.read_sql(
                    "select acc from snp_eval_all where acc = '" + prot + "';", engine
                )
        else:
            mode = "pathogenic"
        if prot.find("_") != -1:
            data_tm = pd.read_sql(
                "select tm_start, tm_final from tm_segments where id = '" + prot + "';",
                engine,
            )
        else:
            data_tm = pd.read_sql(
                "select tm_start, tm_final from tm_segments where acc = '"
                + prot
                + "';",
                engine,
            )

        reg_start = data_tm.tm_start
        reg_start = list(reg_start)
        reg_end = data_tm.tm_final
        reg_end = list(reg_end)

        tm_okey = False
        for i, s_pos in enumerate(reg_start):
            if int(pos) >= s_pos and int(pos) <= reg_end[i]:
                tm_okey = True
            else:
                continue
        if tm_okey == False:
            return render_template(
                "error_no_aa_tm.html",
                input_seq=prot,
                aa_ref=aa_ref,
                pos=pos,
                reg_end=reg_end[-1],
                reg_start=reg_start[0],
            )

    if prot.find("_") != -1:
        data = pd.read_sql(
            "select aa_ref from snp_eval_all where id = '"
            + prot
            + "' and aa_ref = '"
            + aa_ref
            + "';",
            engine,
        )
    else:
        data = pd.read_sql(
            "select aa_ref from snp_eval_all where acc = '"
            + prot
            + "' and aa_ref = '"
            + aa_ref
            + "';",
            engine,
        )
    if data.empty:
        return render_template(
            "error_no_aa_ref.html",
            input_seq=prot,
            aa_ref=aa_ref,
            pos=pos,
            aa_mut=aa_mut,
        )

    if tm_okey == True:
        if prot.find("_") != -1:
            data = pd.read_sql(
                "select * from snp_eval_all where id = '"
                + prot
                + "' and aa_ref = '"
                + aa_ref
                + "' and snp_pos = '"
                + pos
                + "' and aa_mut = '"
                + aa_mut
                + "';",
                engine,
            )
            data_patho = pd.read_sql(
                "select * from snp_eval where id = '" + prot + "';", engine
            )
        else:
            data = pd.read_sql(
                "select * from snp_eval_all where acc = '"
                + prot
                + "' and aa_ref = '"
                + aa_ref
                + "' and snp_pos = '"
                + pos
                + "' and aa_mut = '"
                + aa_mut
                + "';",
                engine,
            )
            data_patho = pd.read_sql(
                "select * from snp_eval where acc = '" + prot + "';", engine
            )

        # PREDICTION OR NOT (IF SUBMIT EXIST ON DATABASE)
        if (
            data.empty
        ):  # First, if we have the info on the database we don't need perform the prediction. In case, we don't have the case introduce by the user perform the prediction
            # PATHOGENIC LIST2
            if data_patho.empty:
                patho = "Non-pathogenic"
            else:
                patho = "Pathogenic"
            # Calculate score
            subs_mat = phat_matrix[dict_matrix[aa_ref]][dict_matrix[aa_mut]]
            if type(subs_mat) == numpy.int64:
                subs_mat = int(subs_mat)

            # Calculate the entropy initial
            # url_uniprot = "http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atransmem%29+AND+organism%3A%22Homo+sapiens+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score&format=txt"
            # req = requests.get(url_uniprot)
            # l_data = req.text.splitlines()
            l_data = open("/var/www/tmsnp/tmsnp/uniprot_file.txt").readlines()
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
                            adn_all += data_snp_reg_entro[l_pos_adn][1][
                                reg_ini - 1 : reg_end
                            ]

            d_freq_abs = dict(collections.Counter(adn_all))

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

            # Calculate freq and entropy
            if prot.find("_") != -1:
                data = pd.read_sql(
                    "select distinct(acc), id, pfam from snp_eval_all where id = '"
                    + prot
                    + "' and aa_ref = '"
                    + aa_ref
                    + "';",
                    engine,
                )
            else:
                data = pd.read_sql(
                    "select distinct(acc), id, pfam from snp_eval_all where acc = '"
                    + prot
                    + "' and aa_ref = '"
                    + aa_ref
                    + "';",
                    engine,
                )

            id = data.id
            id = list(id)
            id = id[0]

            l_pfam = data.pfam
            l_pfam = list(l_pfam)

            # l_var = data.snp_id
            l_var = []

            l_freq_ref, l_freq_mut, l_entro, l_subs_mat, l_aa_ref, l_aa_mut, l_snp_pos = (
                [],
                [],
                [],
                [],
                [],
                [],
                [],
            )
            for pfam_ac in l_pfam:
                try:
                    l_data = file_to_lines(path_pfam + pfam_ac)
                    l_prot_pos_aa = create_l_prot(l_data)
                    l_prot_pos_total = create_l_variables(
                        l_prot_pos_aa, id, l_var, int(pos)
                    )
                    l_prot_rs = create_l_rs(l_prot_pos_total, l_prot_pos_aa)
                    freq_entro = cal_freq_entro(l_prot_rs, entro_ini, aa_mut)
                    freq_ref = float(freq_entro[0])
                    freq_mut = float(freq_entro[1])
                    entropia = float(freq_entro[2])
                    l_subs_mat.append(subs_mat)
                    l_freq_ref.append(freq_ref)
                    l_freq_mut.append(freq_mut)
                    l_entro.append(entropia)
                    l_aa_ref.append(aa_ref)
                    l_aa_mut.append(aa_mut)
                    l_snp_pos.append(pos)
                except FileNotFoundError:
                    continue
            # Complete table
            data["subs_mat"] = numpy.asarray(l_subs_mat)
            data["freq_ref"] = numpy.asarray(l_freq_ref)
            data["freq_mut"] = numpy.asarray(l_freq_mut)
            data["entropy"] = numpy.asarray(l_entro)
            data["aa_ref"] = numpy.asarray(l_aa_ref)
            data["aa_mut"] = numpy.asarray(l_aa_mut)
            data["snp_pos"] = numpy.asarray(l_snp_pos)
            if mode == "pathogenic":
                base_data = pd.read_csv(
                    "/var/www/tmsnp/tmsnp/models/test_ref_patho_v2.txt", sep="\t"
                )
                base_data[:] = 0
                base_data = base_data.iloc[:1]
                # Replace values:
                base_data.at[0, data["acc"][0]] = 1
                base_data.at[0, "subs_mat"] = data["subs_mat"]
                base_data.at[0, "freq_ref"] = data["freq_ref"]
                base_data.at[0, "freq_mut"] = data["freq_mut"]
                base_data.at[0, "entropy"] = data["entropy"]
                base_data.at[0, data["pfam"][0]] = 1
                base_data.at[0, data["aa_ref"][0]] = 1
                base_data.at[0, d_mut[data["aa_mut"][0]]] = 1
                entry = np.asarray([base_data.iloc[0][:-1]])
                # if mode == "nonpathogenic":
                #     base_data = pd.read_csv(
                #         "/var/www/tmsnp/tmsnp/models/test_nonpatho_v1.txt", sep="\t"
                #     )
                #     base_data[:] = 0
                #     base_data = base_data.iloc[:1]
                #     # Replace values:
                #     base_data.at[0, "subs_mat"] = data["subs_mat"]
                #     base_data.at[0, "freq_ref"] = data["freq_ref"]
                #     base_data.at[0, "freq_mut"] = data["freq_mut"]
                #     base_data.at[0, "entropy"] = data["entropy"]
                #     base_data.at[0, data["aa_ref"][0]] = 1
                #     base_data.at[0, d_mut[data["aa_mut"][0]]] = 1
                #     print(base_data.iloc[0][:-1])
                #     print(list(base_data.columns.values))
                #     entry = np.asarray([base_data.iloc[0][:-1]])
                # print(entry)
                pred_75 = get_predictions(entry, 0.25, mode)  # 1
                pred_80 = get_predictions(entry, 0.2, mode)  # 2
                pred_85 = get_predictions(entry, 0.15, mode)  # 3
                pred_90 = get_predictions(entry, 0.1, mode)  # 4
                pred_95 = get_predictions(entry, 0.05, mode)  # 5
                pred_99 = get_predictions(entry, 0.01, mode)  # 6

                # CONFIDENCE
                preds = [pred_75, pred_80, pred_85, pred_90, pred_95, pred_99]

                c = 0
                prot_disf = "Outside the domain of applicability"
                for p in preds:  #
                    if p != "Outside the domain of applicability":
                        if p == "Positive":
                            prot_disf = "Pathogenic"
                        if p == "Negative":
                            prot_disf = "Non-pathogenic"
                        c += 1
                if c == 0:
                    conf = ""
                if c == 1:
                    conf = " (75%)"
                if c == 2:
                    conf = " (80%)"
                if c == 3:
                    conf = " (85%)"
                if c == 4:
                    conf = " (90%)"
                if c == 5:
                    conf = " (95%)"
                if c == 6:
                    conf = " (99%)"

            # CLASSIFICATION:
            data = data.round({"freq_ref": 3, "freq_mut": 3, "entropy": 3})

            if patho == "Pathogenic":
                if prot_disf == "Pathogenic":
                    classi = "Pathogenic"
                elif prot_disf == "Non-pathogenic":
                    classi = "Non-pathogenic"
                else:
                    classi = "Outside the domain of applicability"
            if patho == "Non-pathogenic":
                classi = "Non-pathogenic"
            print(classi)

            # Extra info
            for name, aa in dict_aa.items():
                if aa == data["aa_ref"][0]:
                    aa_ref_3 = name.capitalize()
                if aa == data["aa_mut"][0]:
                    aa_mut_3 = name.capitalize()

            return render_template(
                "results_pred.html",
                acc=data["acc"][0],
                id=data["id"][0],
                pfam=data["pfam"][0],
                aa_ref=data["aa_ref"][0],
                aa_ref_3=aa_ref_3,
                aa_mut=data["aa_mut"][0],
                aa_mut_3=aa_mut_3,
                snp_pos=data["snp_pos"][0],
                subs_mat=data["subs_mat"][0],
                freq_ref=data["freq_ref"][0],
                freq_mut=data["freq_mut"][0],
                entropy=data["entropy"][0],
                conf=conf,
                classification=classi,
            )
        else:
            data = data.round({"freq_ref": 3, "freq_mut": 3, "entropy": 3})
            # CLASSIFICATION:
            if data.pathogenic.values[0] == 0:
                patho = "Non-pathogenic"
                classi = "Non-pathogenic"
            if data.pathogenic.values[0] == 1 or data.pathogenic.values[0] == 2:
                patho = "Pathogenic"
                classi = "Pathogenic"

            # Extra info
            for name, aa in dict_aa.items():
                if aa == data["aa_ref"][0]:
                    aa_ref_3 = name.capitalize()
                if aa == data["aa_mut"][0]:
                    aa_mut_3 = name.capitalize()
            if data["gnomad_freq"][0] == None:
                gnomad_freq = 0
            else:
                gnomad_freq = data["gnomad_freq"][0]
            return render_template(
                "results_pred_exist.html",
                acc=data["acc"][0],
                id=data["id"][0],
                pfam=data["pfam"][0],
                aa_ref=data["aa_ref"][0],
                aa_ref_3=aa_ref_3,
                aa_mut=data["aa_mut"][0],
                aa_mut_3=aa_mut_3,
                snp_pos=data["snp_pos"][0],
                subs_mat=data["subs_mat"][0],
                freq_ref=data["freq_ref"][0],
                freq_mut=data["freq_mut"][0],
                entropy=data["entropy"][0],
                pathogenic=patho,
                gnomad_freq=gnomad_freq,
            )


###########################################################
@app.route("/get-dataset/<filename>")
def show_original(filename):
    return send_from_directory(
        "/var/www/tmsnp/tmsnp/static", filename=filename, as_attachment=True
    )


@app.route("/569Vdataset", methods=["POST", "GET"])
def show_569V():
    # training_569 = pd.read_csv(
    #     "/var/www/tmsnp/tmsnp/static/training_patho_28052019.txt", sep="\t"
    # )
    # test_569 = pd.read_csv(
    #     "/var/www/tmsnp/tmsnp/static/training_patho_28052019.txt", sep="\t"
    # )
    # return render_template(
    #     "show_database.html",
    #     training=training_569.to_html(
    #         justify="center", classes="em_table", index=False, table_id="table_training"
    #     ),
    #     test=test_569.to_html(
    #         justify="center", classes="em_table", index=False, table_id="table_test"
    #     ),
    # )
    return send_file(
        "/var/www/tmsnp/tmsnp/static/569V_dataset.zip",
        attachment_filename="569V_dataset.zip",
    )


@app.route("/4Vdataset", methods=["POST", "GET"])
def show_4V():
    #     # training_569 = pd.read_csv(
    #     #     "/var/www/tmsnp/tmsnp/static/training_patho_28052019.txt", sep="\t"
    #     # )
    #     # test_569 = pd.read_csv(
    #     #     "/var/www/tmsnp/tmsnp/static/training_patho_28052019.txt", sep="\t"
    #     # )
    #     # return render_template(
    #     #     "show_database.html",
    #     #     training=training_569.to_html(
    #     #         justify="center", classes="em_table", index=False, table_id="table_training"
    #     #     ),
    #     #     test=test_569.to_html(
    #     #         justify="center", classes="em_table", index=False, table_id="table_test"
    #     #     ),
    #     # )
    return send_file(
        "/var/www/tmsnp/tmsnp/static/4V_dataset.zip",
        attachment_filename="4V_dataset.zip",
    )


@app.route("/about")
def script_about():
    return render_template("how_it_works.html")


# @app.route("/help")
# def script_help():
#     return render_template("help.html")


if __name__ == "__main__":
    app.run(debug=True)

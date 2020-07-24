# This module contains all the functions developped by Iker Reina for the TM-SNP project

import collections
from math import exp, expm1, log10, log
import subprocess
import sys


def file_to_lines(filename):
    """
    :param filename: arxiu a llegir
    :rtype: llista de cada linia darxiu
    """
    fil = open(filename)
    lines = fil.readlines()
    fil.close()
    return lines


def zerolistmaker(n):
    l_0 = [0] * n
    return l_0


def extract_prot(l_data):
    """
    :param l_data: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [ID uniprot, AC el code uniprot, GN nom gen, [PFAM]
    """
    pfam_loc = list()
    id = ""
    loc = ""
    code = ""
    gen = ""
    l1 = list()
    l_of_l = list()
    for data_str in l_data:
        if data_str[0:2] == "ID" and id == "":
            l_id = data_str.split("   ")
            id = l_id[1]
        elif data_str[0:2] == "AC" and code == "":
            l_ac = data_str.split(";")
            code = l_ac[0][5:]
        elif data_str[0:2] == "GN" and gen == "":
            l_gn = data_str.split("=")
            gen = l_gn[1][: l_gn[1].find(" ")].strip(";")
        elif data_str[0:2] == "DR":
            l_pfam = data_str.split(" ")
            if l_pfam[3] == "Pfam;":
                pfam_loc.append(l_pfam[4][:-1])
        elif data_str[0:2] == "FT":
            l_ft = data_str.split(" ")
            for n in l_ft:
                if n == "TRANSMEM":
                    loc = "transm_loc"
        elif data_str[0:2] == "//" and loc == "transm_loc":
            l1.append(id)
            l1.append(code)
            l1.append(gen)
            l1.append(pfam_loc)
            l_of_l.append(l1)
            pfam_loc = []
            l1 = []
            id = ""
            code = ""
            gen = ""
            loc = ""
        elif data_str[0:2] == "//" and loc == "":
            pfam_loc = []
            l1 = []
            id = ""
            code = ""
            gen = ""
        else:
            continue
    return l_of_l


def extract_snp_regions(l_data):
    """
    :param l_data: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR]] i descarta les proteines sense regio transmembrana
    """
    id = ""
    gen = ""
    id_var_loc = ""
    disease_match = ""
    loc = ""
    code = ""
    disease = list()
    sub_loc = list()
    var_loc = list()
    var_change = list()
    var_id = list()
    var_rs = list()
    l1 = list()
    l_of_l = list()
    for data_str in l_data:
        if data_str[0:2] == "ID" and id == "":
            l_id = data_str.split("   ")
            id = l_id[1]
        elif data_str[0:2] == "GN" and gen == "":
            l_gn = data_str.split("=")
            gen = l_gn[1][: l_gn[1].find(" ")].strip(";")
        elif data_str[0:2] == "AC" and code == "":
            l_ac = data_str.split(";")
            code = l_ac[0][5:]
        elif data_str[0:17] == "CC   -!- DISEASE:":
            disease_match = "trobat"
            if data_str.find("(") == -1:
                disease_match = "trobat_dolent"
            else:
                disease.append(data_str[data_str.find("(") + 1 : data_str.find(")")])
        elif data_str[0:9] == "CC       " and disease_match == "trobat_dolent":
            if data_str.find("(") != -1:
                disease_match = "trobat"
                disease.append(data_str[data_str.find("(") + 1 : data_str.find(")")])
            else:
                disease_match = "trobat_dolent"
        elif data_str[0:2] == "FT":
            if data_str[5:13] == "TRANSMEM":
                loc = "transm_loc"
                sub_loc.append(data_str[13:27])
            elif data_str[5:12] == "VARIANT":
                var_rs.append("-")
                id_var_loc = "trobat"
                var_loc.append(data_str[20:-2])
            if id_var_loc == "trobat" and data_str.find("->") != -1:
                var_change.append(
                    data_str[data_str.find("->") - 3 : data_str.find("->") + 5].strip(
                        " \n ."
                    )
                )
            if id_var_loc == "trobat" and data_str.find("Missing") != -1:
                var_change.append("MISS -> MISS")
            if id_var_loc == "trobat" and data_str.find("/FTId=") != -1:
                var_id.append(data_str[40:50])
                id_var_loc = ""
            if id_var_loc == "trobat" and data_str.find("dbSNP:rs") != -1:
                var_rs.pop(len(var_rs) - 1)
                var_rs.append(
                    data_str[
                        data_str.find("dbSNP:rs") + 6 : data_str.find("dbSNP:rs") + 17
                    ].rstrip(").\n")
                )
        elif data_str[0:2] == "//" and loc == "transm_loc":
            l1.append(code)
            l1.append(sub_loc)
            l1.append(var_loc)
            l1.append(var_change)
            l1.append(var_id)
            l1.append(var_rs)
            l1.append(disease)
            l1.append(id)
            l1.append(gen)
            l_of_l.append(l1)
            sub_loc = []
            var_loc = []
            var_id = []
            var_rs = []
            l1 = []
            var_change = []
            disease = []
            code = ""
            loc = ""
            id = ""
            gen = ""
        elif data_str[0:2] == "//" and loc == "":
            sub_loc = []
            var_loc = []
            var_id = []
            var_rs = []
            l1 = []
            var_change = []
            disease = []
            code = ""
            id = ""
            gen = ""
        else:
            continue
    return l_of_l


def descart_no_transm(data_loc):
    """
    :param data_locations:  una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR]]
    :return: una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR], [llista pos snp transmembraana]]
    """
    l_f, l_snps = [], []
    for n_prot in range(0, len(data_loc)):
        uniprot = data_loc[n_prot][0]
        pos_tm = data_loc[n_prot][1]
        snp_rs = data_loc[n_prot][2]
        mut = data_loc[n_prot][3]
        snp_var = data_loc[n_prot][4]
        snp_rs = data_loc[n_prot][5]
        disease = data_loc[n_prot][6]
        id = data_loc[n_prot][7]
        gen = data_loc[n_prot][8]
        if len(snp_rs) != len(snp_var):
            print("Problem with SNPs, some don't have code VAR_...")
            print("PROTEIN: " + id)
            sys.exit()
        if len(snp_var) != len(snp_rs):
            print("Problem with SNPs, the rs codes don't genereted well...")
            sys.exit()
        if len(snp_rs) != len(mut):
            print(snp_rs)
            print(mut)
            print(
                "hi ha un problema amb els canvis de la mutacio (aminoacid) no tots apareix el canvi....... "
            )
            sys.exit()
        l_snp_tm = zerolistmaker(len(snp_rs))
        l_snp_dis = zerolistmaker(len(snp_rs))
        for snp in snp_rs:
            l_snps.append(snp[:7].strip())
            for a_disease in disease:
                if snp.find(a_disease) != -1:
                    l_snp_dis.pop(snp_rs.index(snp))
                    l_snp_dis.insert(snp_rs.index(snp), 1)
        for coords in pos_tm:
            x_coord_tm = coords[:7].strip()
            y_coord_tm = coords[7:].strip()
            for snp in snp_rs:
                pos_snp = snp[:7].strip()
                if int(pos_snp) >= int(x_coord_tm) and int(pos_snp) <= int(y_coord_tm):
                    l_snp_tm.pop(len(l_snp_tm) - 1)
                    l_snp_tm.insert(snp_rs.index(snp), 1)
        if len(snp_var) != len(l_snp_dis):
            print("Different lenght")
        if len(l_snp_dis) != len(l_snp_tm):
            print("Different lenght")
        l_f.append(
            [
                uniprot,
                pos_tm,
                snp_var,
                mut,
                snp_rs,
                l_snps,
                l_snp_tm,
                l_snp_dis,
                id,
                gen,
            ]
        )
        l_snps = []
    return l_f


def run_curl(url, file_out):
    proc = subprocess.Popen(["curl", url], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    with open(file_out, "w") as fileout:
        fileout.write(out)


def create_l_prot(l_data):
    l1, l_prot_aa = [], []
    if ">" not in l_data[0]:
        for prot in l_data:
            l_prot = prot.split(" ")
            if l_prot[0].find("/") != -1:
                pos_barra = l_prot[0].find("/")
            else:
                print(
                    "Error in code of the protein. Don't found the space between name and position(/)"
                )
                print("Error in: " + prot)
                break
            l1.append(l_prot[0][:pos_barra])
            l1.append(l_prot[0][pos_barra + 1 :])
            for pos in l_prot:
                if pos.find(l_prot[0]) == -1 and pos != "":
                    l1.append(pos.rstrip("\n"))
            l_prot_aa.append(l1)
            l1 = []
        return l_prot_aa
    else:
        aminoacids = ""
        entry = 0
        for prot in l_data:
            if ">" in prot:
                if entry == 1:
                    l1.append(aminoacids)
                    aminoacids = ""
                    l_prot_aa.append(l1)
                    l1 = []
                tag = prot.split(">")
                name_pos = tag[1].split("/")
                name_prot = name_pos[0]
                coord_aa = name_pos[1][:-1]
                l1.append(name_prot)
                l1.append(coord_aa)
                entry = 1
            if ">" not in prot:
                aminoacids = aminoacids + prot[:-1]
                if prot == l_data[-1]:
                    l1.append(aminoacids)
                    l_prot_aa.append(l1)
        return l_prot_aa


def create_l_variables(l_prot_pos_aa, id_prot, snp_var, snp_pos):
    l_prot_pos_total = list()
    l1 = list()
    for n_prot in range(0, len(l_prot_pos_aa)):
        name_prot = l_prot_pos_aa[n_prot][0]
        coord_aa = l_prot_pos_aa[n_prot][1]
        aminoacids = l_prot_pos_aa[n_prot][2]
        if (
            name_prot == id_prot
            and snp_pos >= int(coord_aa[: coord_aa.find("-")])
            and snp_pos <= int(coord_aa[coord_aa.find("-") + 1 :])
        ):
            c_1 = int(coord_aa[: coord_aa.find("-")])
            c_2 = 0
            l1.append(name_prot)
            l1.append(snp_var)
            for caracter in aminoacids:
                if caracter != "-" and c_1 <= snp_pos:
                    c_1 = c_1 + 1
                    caracter_final = caracter
                if c_1 <= snp_pos:
                    c_2 = c_2 + 1
            l1.append(caracter_final)
            l1.append(c_2)
            l_prot_pos_total.append(l1)
            l1 = []
    return l_prot_pos_total


def create_l_rs(l_prot_pos_total, l_prot_pos_aa):
    l_prot, l1, l_rs = [], [], []
    for numero_prot in range(0, len(l_prot_pos_total)):
        name_prot = l_prot_pos_total[numero_prot][0]
        snip_rs = l_prot_pos_total[numero_prot][1]
        amino = l_prot_pos_total[numero_prot][2]
        pos = l_prot_pos_total[numero_prot][3]
        l1.append(name_prot)
        l1.append(snip_rs)
        l1.append(amino)
        for n_prot in range(0, len(l_prot_pos_aa)):
            aminoacids = l_prot_pos_aa[n_prot][2]
            l_prot.append(aminoacids[pos])
        l1.append(l_prot)
        l_rs.append(l1)
        l_prot = []
        l1 = []
    return l_rs


def cal_freq_entro(l_rs, entro_ini, aa):
    l_aa = [
        "Y",
        "G",
        "F",
        "M",
        "A",
        "S",
        "I",
        "L",
        "T",
        "V",
        "P",
        "K",
        "H",
        "Q",
        "E",
        "Z",
        "W",
        "R",
        "D",
        "N",
        "B",
        "C",
        "X",
    ]
    for n_reg_pfam in range(0, len(l_rs)):
        d_freq_abs = dict(collections.Counter(l_rs[n_reg_pfam][3]))
        if "-" not in d_freq_abs:
            gap = 0
        else:
            gap = d_freq_abs["-"]
            del d_freq_abs["-"]
        if "X" not in d_freq_abs:
            x = 0
        else:
            x = d_freq_abs["X"]
            del d_freq_abs["X"]
        if "Z" in d_freq_abs:
            print("hi ha una Z en laliniament")
        if "B" in d_freq_abs:
            print("hi ha una B en laliniament: ")
        len_no_gap = len(l_rs[n_reg_pfam][3]) - (
            gap + x
        )  # li trec els gaps i les x a la llargada total de les posns
        fx = d_freq_abs[l_rs[n_reg_pfam][2]] / len_no_gap
        d_freq_rel = dict()
        for aminoacid in d_freq_abs:
            d_freq_rel[aminoacid] = d_freq_abs[aminoacid] / len_no_gap
        if aa in d_freq_rel:
            fx_amino_nou = d_freq_rel[aa]
        else:
            fx_amino_nou = 0
        suma = 0.0
        for d_aa in d_freq_rel:
            suma += d_freq_rel[d_aa] * log(d_freq_rel[d_aa])
        entropia = entro_ini + suma
        return fx, fx_amino_nou, entropia


def ext_prot_entro(l_data):
    """
    :param l_data: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [ID uniprot, AC el code uniprot, GN nom gen, [PFAM]
    """
    dna_loc, l1, l_of_l = [], [], []
    loc = ""
    code = ""
    for data_str in l_data:
        if data_str[0:2] == "AC" and code == "":
            l_ac = data_str.split(";")
            code = l_ac[0][5:]
        elif data_str[0:2] == "  ":
            adn = data_str.strip("\n ")
            dna_loc.append(adn.replace(" ", ""))
        elif data_str[0:2] == "FT":
            l_ft = data_str.split(" ")
            for n in l_ft:
                if n == "TRANSMEM":
                    loc = "transm_loc"
        elif data_str[0:2] == "//" and loc == "transm_loc":
            l1.append(code)
            l1.append(dna_loc)
            l_of_l.append(l1)
            dna_loc, l1 = [], []
            code = ""
            loc = ""
        elif data_str[0:2] == "//" and loc == "":
            dna_loc, l1 = [], []
            code = ""
        else:
            continue
    return l_of_l


def create_l_prot_entro(data_loc):
    """
    :param data_locations:  una llista de l
    :return: una llista de llistes
    """
    l, l_f = [], []
    str_adn = ""
    for n_prot in range(0, len(data_loc)):
        uniprot = data_loc[n_prot][0]
        dna_full = data_loc[n_prot][1]
        for adn in dna_full:
            str_adn += adn
        l.append(uniprot)
        l.append(str_adn)
        l_f.append(l)
        str_adn = ""
        l = []
    return l_f


def snp_reg_entro(l_data):
    """
    :param l_data: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR]] i descarta les proteines sense regio transmembrana
    """
    id_var_loc = ""
    loc = ""
    code = ""
    sub_loc, var_loc, var_id, var_rs, l1, l_of_l = [], [], [], [], [], []
    for data_str in l_data:
        if data_str[0:2] == "AC" and code == "":
            l_ac = data_str.split(";")
            code = l_ac[0][5:]
        elif data_str[0:2] == "FT":
            if data_str[5:13] == "TRANSMEM":
                loc = "transm_loc"
                sub_loc.append(data_str[13:27])
            elif data_str[5:12] == "VARIANT":
                var_rs.append("-")
                id_var_loc = "trobat"
                var_loc.append(data_str[20:-2])
            if id_var_loc == "trobat" and data_str.find("/FTId=") != -1:
                var_id.append(data_str[40:50])
                id_var_loc = ""
            if id_var_loc == "trobat" and data_str.find("dbSNP:rs") != -1:
                var_rs.pop(len(var_rs) - 1)
                var_rs.append(
                    data_str[
                        data_str.find("dbSNP:rs") + 6 : data_str.find("dbSNP:rs") + 16
                    ].rstrip(").\n")
                )
        elif data_str[0:2] == "//" and loc == "transm_loc":
            l1.append(code)
            l1.append(sub_loc)
            l1.append(var_loc)
            l1.append(var_id)
            l1.append(var_rs)
            l_of_l.append(l1)
            sub_loc, var_loc, var_id, var_rs, l1 = [], [], [], [], []
            code = ""
            loc = ""
        elif data_str[0:2] == "//" and loc == "":
            sub_loc, var_loc, var_id, var_rs, l1 = [], [], [], [], []
            code = ""
        else:
            continue
    return l_of_l


def rm_no_stm_entro(data_loc):
    """
    :param data_locations:  una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR]]
    :return: una llista de llistes amb: [Accesion de el code uniprot, [ llista de posns transmembrana], [llista de SNPs RS]], [llista snp VAR], [llista pos snp transmembraana]]
    """
    l_f, coord_trans, l_snps = [], [], []
    for n_prot in range(0, len(data_loc)):
        uniprot = data_loc[n_prot][0]
        pos_tm = data_loc[n_prot][1]
        snp_rs = data_loc[n_prot][2]
        snp_var = data_loc[n_prot][3]
        snp_rs_s = data_loc[n_prot][4]
        if len(snp_rs) != len(snp_var):
            print("Problem with SNPs, some don't have code VAR_...")
            print("Protein: " + snp_var)
            sys.exit()
        if len(snp_var) != len(snp_rs_s):
            print("Problem with SNPs, the rs codes don't be created well...")
            sys.exit()
        l_snp_tm = zerolistmaker(len(snp_rs))
        for snp in snp_rs:
            l_snps.append(snp[:7].strip())
        for coords in pos_tm:
            coord_trans.append(coords)
            x_coord_tm = coords[:7].strip()
            y_coord_tm = coords[7:].strip()
            for snp in snp_rs:
                pos_snp = snp[:7].strip()
                if int(pos_snp) >= int(x_coord_tm) and int(pos_snp) <= int(y_coord_tm):
                    l_snp_tm.pop(len(l_snp_tm) - 1)
                    l_snp_tm.insert(snp_rs.index(snp), 1)
        l_f.append([uniprot, pos_tm, snp_var, snp_rs_s, l_snps, l_snp_tm])
        coord_trans, l_snps = [], []
    return l_f

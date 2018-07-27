#!/usr/bin/env python3

################################################
import numpy as np
import requests
import collections
from math import exp, expm1, log10, log
import subprocess
import mysql.connector
import os
import collections
from iker_snp import *


pfam_download_path = './pfam_download/'

# Set the connection
conn = mysql.connector.connect(host='alf03.uab.cat', user='lmcdb',
                               password=os.getenv('LMCDB_PASS'), database='tmsnp')
mycursor = conn.cursor()


print('Evaluating SNPs based on the PFAM alignment')

mycursor.execute("select distinct acc from snps;")
distinct_acc_tupple = mycursor.fetchall()
distinct_acc_list = list()
for tupple in distinct_acc_tupple:
    for acc in tupple:
        distinct_acc_list.append(acc)

print(distinct_acc_list)

###################################################################################


phat_matrix = np.array([[5, -6, -2, -5, 1, -3, -5, 1, -3, 0, -1, -7, -1, -1, -3, 2, 0, -4, -3, 1],
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
                        [-7, -1, -2, -5, -10, -1, -4, -5, -5, -7, -7, 5, -6, -7, -4, -5, -6, -8, -4, -8],
                        [-1, -6, -2, -5, -2, -1, -5, -1, -4, 3, 2, -6, 6, 0, -5, -2, 0, -4, -2, 1],
                        [-1, -7, -1, -5, 0, -2, -5, -2, -2, 0, 1, -7, 0, 6, -5, -2, -2, 0, 4, -1],
                        [-3, -7, -4, -5, -8, -3, -5, -3, -6, -4, -5, -4, -5, -5, 13, -3, -4, -6, -5, -4],
                        [2, -6, 1, -4, 1, -1, -3, 1, -2, -2, -2, -5, -2, -2, -3, 6, 1, -5, -2, -2],
                        [0, -6, -1, -5, -1, -3, -5, -1, -4, -1, -1, -6, 0, -2, -4, 1, 3, -7, -3, 0],
                        [-4, -7, -5, -7, -4, 1, -7, -5, -3, -4, -3, -8, -4, 0, -6, -5, -7, 11, 1, -4],
                        [-3, -6, 2, -4, -1, 0, -2, -3, 3, -3, -2, -4, -2, 4, -5, -2, -3, 1, 11, -3],
                        [1, -7, -3, -5, -2, -3, -5, -2, -5, 3, 1, -8, 1, -1, -4, -2, 0, -4, -3, 4]])

dict_matrix = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
               'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, }

###############################################################################################

mycursor.execute("DROP TABLE IF EXISTS snp_phat")
mycursor.execute("""CREATE TABLE snp_phat (acc varchar(20), id varchar(30), gene varchar(20),
                    snp_id varchar(30), snp_rs varchar(30), aa_ini varchar(30), 
                    aa_fin varchar(30), snp_pos int, pathogenic int, mat_punct int);""")

########################################################################################33


for distinct_acc in distinct_acc_list:
    print("S'ha analitzat i omplert la taula_matrix amb les dades de la matriu del: " + distinct_acc)
    mycursor.execute("select distinct snp_id from snps where acc = '" + distinct_acc + "';")
    distinct_var_tupple = mycursor.fetchall()
    distinct_var_list = list()
    for tupple in distinct_var_tupple:
        for var in tupple:
            distinct_var_list.append(var)

    for var in distinct_var_list:
        mycursor.execute("select aa_ini from snps where snp_id = '" + var + "';")
        distinct_aa_inici_tupple = mycursor.fetchall()
        for tupple in distinct_aa_inici_tupple:
            for aa_inici in tupple:
                distinct_aa_inici = aa_inici.strip()

        mycursor.execute("select aa_fin from snps where snp_id = '" + var + "';")
        distinct_aa_final_tupple = mycursor.fetchall()
        for tupple in distinct_aa_final_tupple:
            for aa_final in tupple:
                distinct_aa_final = aa_final.strip()

        mycursor.execute("select id from snps where snp_id = '" + var + "';")
        distinct_id_tupple = mycursor.fetchall()
        for tupple in distinct_id_tupple:
            for id in tupple:
                if type(id) == str:
                    distinct_id = id.strip()
                else:
                    distinct_id = id

        mycursor.execute("select gene from snps where snp_id = '" + var + "';")
        distinct_gene_tupple = mycursor.fetchall()
        for tupple in distinct_gene_tupple:
            for gene in tupple:
                if type(gene) == str:
                    distinct_gene = gene.strip()
                else:
                    distinct_gene = gene

        mycursor.execute("select snp_rs from snps where snp_id = '" + var + "';")
        distinct_rs_tupple = mycursor.fetchall()
        for tupple in distinct_rs_tupple:
            for rs in tupple:
                if type(rs) == str:
                    distinct_rs = rs.strip()
                else:
                    distinct_rs = rs

        mycursor.execute("select snp_pos from snps where snp_id = '" + var + "';")
        distinct_snp_pos_tupple = mycursor.fetchall()
        for tupple in distinct_snp_pos_tupple:
            for snp_pos in tupple:
                distinct_snp_pos = snp_pos  

        mycursor.execute("select pathogenic from snps where snp_id = '" + var + "';")
        distinct_pathogenic_tupple = mycursor.fetchall()
        for tupple in distinct_pathogenic_tupple:
            for pathogenic in tupple:
                distinct_pathogenic = pathogenic 

        if (len(distinct_aa_inici) >= 2 or len(distinct_aa_final) >= 2 or distinct_aa_final == '/'):
             puntuation = None # -1000
#            mycursor.execute(("INSERT INTO snp_phat VALUES" + str((distinct_acc, distinct_id, distinct_gene, var,
#                                                                      distinct_rs, distinct_aa_inici, distinct_aa_final,
#                                                                      distinct_snp_pos, 
#                                                                      distinct_pathogenic, puntuation))))
#            conn.commit()

        else:
            print(distinct_aa_inici, distinct_aa_final)
            puntuation = phat_matrix[dict_matrix[distinct_aa_inici]][dict_matrix[distinct_aa_final]]
#            mycursor.execute(("INSERT INTO snp_phat VALUES" + str((distinct_acc, distinct_id, distinct_gene, var,
#                                                                      distinct_rs, distinct_aa_inici, distinct_aa_final,
#                                                                      distinct_snp_pos, 
#                                                                      distinct_pathogenic, puntuation))))
            sql_str = f"""INSERT INTO snp_phat VALUES ('{distinct_acc}', '{distinct_id}', '{distinct_gene}', '{var}',
                                                                      '{distinct_rs}', '{distinct_aa_inici}', '{distinct_aa_final}',                                                                     {distinct_snp_pos}, {distinct_pathogenic}, {puntuation});"""
            #print(sql_str)
            mycursor.execute(sql_str)
            conn.commit()

#########################################################################################################################################################################



no_pfam = set()
# Entropy

mycursor.execute("DROP TABLE IF EXISTS snp_eval")
mycursor.execute(
    "CREATE TABLE snp_eval (acc varchar(20), id varchar(30), pfam varchar(20), snp_id varchar(30), snp_rs varchar(30), aa_ini varchar(30), aa_fin varchar(30), snp_pos int, fx_i float, fx_f float, entropy float, pathogenic int);")

######################################entropia inicial a parti de les regions transmem del arxiu total d'uniprot###########

# All Human TM proteins tagged as reviewed
url_uniprot = 'http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atransmem%29+AND+organism%3A%22Homo+sapiens+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score&format=txt'
req = requests.get(url_uniprot)
data_llista = req.text.splitlines()  #before: data_llista = file_to_lines(uniprot_filename)

data_taula_protines_entropia = extract_data_taula_proteines_entropia(data_llista)
data_locations_taula_snip_regions_entropia = crear_llista_prot_seqamino_entropia(data_taula_protines_entropia)

data_locations_taula_snip_regions = extract_data_taula_snip_regions_entropia(data_llista)
llista_snps_transmembrana_taula_snip_regions = descart_regio_no_transmembrana_entropia(
    data_locations_taula_snip_regions)

string_adn_total = ''
for num_posicio_llista in range(len(llista_snps_transmembrana_taula_snip_regions)):
    for numero_regions in range(len(llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1])):
        accesion_protein = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][0]
        region_protein_inici = int(
            llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1][numero_regions][:8])
        region_protein_final = int(
            llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1][numero_regions][9:])
        for num_posicio_llista_adn in range(len(data_locations_taula_snip_regions_entropia)):
            if (data_locations_taula_snip_regions_entropia[num_posicio_llista_adn][0] == accesion_protein):
                string_adn_total += data_locations_taula_snip_regions_entropia[num_posicio_llista_adn][1][
                                    region_protein_inici - 1:region_protein_final]

diccionario_freq_abs = dict(collections.Counter(string_adn_total))  # amb els 2 aminoacids sense res mes ni x ni Z ni B

diccionario_freq_relatives = dict()
total_aminoacids = 0
for aminoacid in diccionario_freq_abs:
    total_aminoacids += diccionario_freq_abs[aminoacid]

for aminoacid in diccionario_freq_abs:
    diccionario_freq_relatives[aminoacid] = diccionario_freq_abs[aminoacid] / total_aminoacids

suma = 0.0
for aminoacid_dict in diccionario_freq_relatives:
    suma += diccionario_freq_relatives[aminoacid_dict] * log(diccionario_freq_relatives[aminoacid_dict])

entropia_inicial = -suma
print(entropia_inicial)

##########################################################calcul de la frequencia i entropia de la posicio on cau l'snip#######################################

mycursor.execute("select distinct acc from snps;")
total_acc_tupple = mycursor.fetchall()
total_acc_list = list()
for tupple in total_acc_tupple:
    for tot in tupple:
        total_acc_list.append(tot)
print(total_acc_list)

mycursor.execute("select distinct acc from snp_phat;")
distinct_acc_tupple = mycursor.fetchall()
distinct_acc_list = list()
for tupple in distinct_acc_tupple:
    for acc in tupple:
        distinct_acc_list.append(acc)
print(distinct_acc_list)

for accses in distinct_acc_list:
    mycursor.execute("select distinct pfam from receptor_pfam where acc = '" + accses + "';")
    distinct_pfam_tupple = mycursor.fetchall()
    distinct_pfam_list = list()
    for tupple in distinct_pfam_tupple:
        for pfam in tupple:
            distinct_pfam_list.append(pfam)
    mycursor.execute("select distinct id from snp_phat where acc = '" + accses + "';")
    distinct_id_tupple = mycursor.fetchall()
    for tupple in distinct_id_tupple:
        for id in tupple:
            distinct_id = id
    #print(accses)
    sql_str = "select distinct snp_id from snp_phat where acc = '" + accses + "';"
    #print(sql_str)
    mycursor.execute(sql_str)
    distinct_var_tupple = mycursor.fetchall()
    distinct_var_list = list()
    for tupple in distinct_var_tupple:
        for var in tupple:
            distinct_var_list.append(var)
    for snip_var in distinct_var_list:
        mycursor.execute("select snp_rs from snp_phat where snp_id = '" + snip_var + "';")
        distinct_snip_rs_tupple = mycursor.fetchall()
        for tupple in distinct_snip_rs_tupple:
            for snip_rs in tupple:
                distinct_snip_rs = snip_rs
        mycursor.execute("select snp_pos from snp_phat where snp_id = '" + snip_var + "';")
        distinct_snip_tupple = mycursor.fetchall()
        for tupple in distinct_snip_tupple:
            for snip in tupple:
                distinct_snip = snip
        mycursor.execute("select aa_fin from snp_phat where snp_id = '" + snip_var + "';")
        distinct_aa_tupple = mycursor.fetchall()
        for tupple in distinct_aa_tupple:
            for aa in tupple:
                distinct_aa = aa
        mycursor.execute("select aa_ini from snp_phat where snp_id = '" + snip_var + "';")
        distinct_aa_i_tupple = mycursor.fetchall()
        for tupple in distinct_aa_i_tupple:
            for aa_i in tupple:
                distinct_aa_i = aa_i

        mycursor.execute("select pathogenic from snp_phat where snp_id = '" + snip_var + "';")
        distinct_path_tupple = mycursor.fetchall()
        for tupple in distinct_path_tupple:
            for path in tupple:
                distinct_path = path

        for pfam_ac in distinct_pfam_list:
            try:
                data_llista = file_to_lines(pfam_download_path + pfam_ac)
            except FileNotFoundError:
                print(f'No {pfam_download_path}{pfam_ac}')
                no_pfam.add(pfam_ac)
            llista_prot_posicions_amino = crear_llista_llistes_proteines(data_llista)
            llista_prot_posicio_total = crear_llista_idprotein_rs_amino_posiciototal(llista_prot_posicions_amino,
                                                                                     distinct_id, distinct_var_list,
                                                                                     distinct_snip)
            llista_prot_rs_amino_distinctprot = crear_llista_prot_rs_amino_distinctprot(llista_prot_posicio_total,
                                                                                        llista_prot_posicions_amino)
            freq_entropia = calcul_freq_entropia(llista_prot_rs_amino_distinctprot, entropia_inicial, distinct_aa)
            if distinct_aa == 'MISS':
                freq_entropia = None
            if (freq_entropia == None):
                fx_i = '-1000'
                fx_f = '-1000'
                entropia = '-1000'
                mycursor.execute(("INSERT INTO snp_eval VALUES" + str((accses, distinct_id, pfam_ac, snip_var,
                                                                        distinct_snip_rs, distinct_aa_i, distinct_aa,
                                                                        distinct_snip, fx_i, fx_f, entropia, 
                                                                        distinct_path))))
                conn.commit()
            else:
                fx_i = float(freq_entropia[0])
                fx_f = float(freq_entropia[1])
                entropia = float(freq_entropia[2])
                mycursor.execute(("INSERT INTO snp_eval VALUES" + str((accses, distinct_id, pfam_ac, snip_var,
                                                                        distinct_snip_rs, distinct_aa_i, distinct_aa,
                                                                        distinct_snip, fx_i, fx_f, entropia, 
                                                                        distinct_path))))
                conn.commit()

# Delete non-TM
mycursor.execute("drop table snp_phat;")
conn.commit()
print("Proteins with no PFAM alignment:")
for el in no_pfam:
    print(el, end=' ')

print("... Done!")
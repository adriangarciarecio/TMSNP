#!/usr/bin/env python3

# Generates the MySQL tables with the SNP and related information for the Membrane proteins in the Uniprot

# INFO: some receptors appear twice (or more) when they have multiple pfam domains assigned.

import requests
import mysql.connector
import os
from iker_snp import *

# Set the connection
conn = mysql.connector.connect(host='alf03.uab.cat', user='lmcdb',
                               password=os.getenv('LMCDB_PASS'), database='tmsnp')
mycursor = conn.cursor()


# Dirs where PFAM downloaded files will go
pfam_download_path = './pfam_download/'  # !!Cal posar el PATH correcte de la carpeta "TFM/pfam/pfam_download" que en aquest cas es: "/home/alumne/Escritorio/TFM/pfam/pfam_download"
pfam_download_raw = './pfam_download_raw/'  # !! CAl posar EL PATH DE LA CARPETA: "TFM/pfam/pfam_download_raw" QUE EN AQUEST CAS ES: "/home/alumne/Escritorio/TFM/pfam/pfam_download_raw"

create_dirs = ['pfam_download', 'pfam_download_raw']
for one_dir in create_dirs:
    if not os.path.exists(one_dir):
        os.makedirs(one_dir)

# paths without the ending /
pfam_download_path_sys = pfam_download_path[:-1]
pfam_download_raw_sys = pfam_download_raw[:-1]



###################taules mysql##############
print("Preparing MySQL tables based on UniProt data ...")
mycursor.execute("DROP TABLE IF EXISTS receptor_pfam")
mycursor.execute("DROP TABLE IF EXISTS tm_segments")
mycursor.execute("DROP TABLE IF EXISTS snps")

mycursor.execute("CREATE TABLE receptor_pfam (acc varchar(20), pfam varchar(20));")
mycursor.execute("""CREATE TABLE snps (acc varchar(20), id varchar(30), gene varchar(20), snp_id varchar(30), 
                    snp_rs varchar(30), aa_ini varchar(30), aa_fin varchar(30), snp_pos int, tm_pos int,
                     pathogenic int, PRIMARY KEY (snp_id));""")
mycursor.execute("CREATE TABLE tm_segments (acc varchar(20),  tm_start int, tm_final int);")

# All Human TM proteins tagged as reviewed

url_uniprot = 'http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atransmem%29+AND+organism%3A%22Homo+sapiens+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score&format=txt'
req = requests.get(url_uniprot)
data_llista = req.text.splitlines()  #before: data_llista = file_to_lines(uniprot_filename)


data_taula_protines = extract_data_taula_proteines(data_llista)
data_locations_taula_snip_regions = extract_data_taula_snip_regions(data_llista)
llista_snps_transmembrana_taula_snip_regions = descart_regio_no_transmembrana(data_locations_taula_snip_regions)

for num_posicio_llista in range(len(data_taula_protines)):
    for num_pfam in range(len(data_taula_protines[num_posicio_llista][3])):
        accesion_protein = data_taula_protines[num_posicio_llista][1]
        #        id_protein = data_taula_protines[num_posicio_llista][0]
        #        gen_protein = data_taula_protines[num_posicio_llista][2]
        pfam_protein = data_taula_protines[num_posicio_llista][3][num_pfam]
        mycursor.execute(("INSERT INTO receptor_pfam VALUES" + str((accesion_protein, pfam_protein))))
        conn.commit()

for num_posicio_llista in range(len(llista_snps_transmembrana_taula_snip_regions)):
    for numero_snips in range(len(llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][2])):
        accesion_protein = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][0]
        var_snip = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][2][numero_snips]
        change_snip = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][3][numero_snips]
        snip_inici_simbol = change_snip[:change_snip.find('->')]
        snip_final_simbol = change_snip[change_snip.find('->') + 3:]
        rs_snip = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][4][numero_snips]
        position_snip = int(llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][5][numero_snips])
        rm = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][6][numero_snips]
        disease = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][7][numero_snips]
        id = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][8]
        gen = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][9]
        #        print(str((accesion_protein, var_snip, rs_snip, position_snip, rm)))
        mycursor.execute(("INSERT INTO snps VALUES" + str((accesion_protein, id, gen, var_snip, rs_snip,
                                                               snip_inici_simbol, snip_final_simbol, position_snip, rm,
                                                               disease))))
        conn.commit()

for num_posicio_llista in range(len(llista_snps_transmembrana_taula_snip_regions)):
    for numero_regions in range(len(llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1])):
        accesion_protein = llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][0]
        region_protein_inici = int(
            llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1][numero_regions][:8])
        region_protein_final = int(
            llista_snps_transmembrana_taula_snip_regions[num_posicio_llista][1][numero_regions][9:])
        #        print(str((accesion_protein, region_protein_inici, region_protein_final)))
        mycursor.execute(
            ("INSERT INTO tm_segments VALUES" + str((accesion_protein, region_protein_inici, region_protein_final))))
        conn.commit()

# Delete non-TM
mycursor.execute("delete from snps where tm_pos=0;")

# Delete tm column
mycursor.execute("alter table snps drop column tm_pos;") 
conn.commit()

print("... Done!!!")

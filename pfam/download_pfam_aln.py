#!/usr/bin/env python3

import requests
import mysql.connector
import os
import subprocess
from iker_snp import *

########################################seleccionar els diferents pfam de les taules de mysql#####################3


pfam_download_path = './pfam_download/'
pfam_download_path_sys = pfam_download_path[:-1]

# Set the connection
conn = mysql.connector.connect(host='alf03.uab.cat', user='lmcdb',
                               password=os.getenv('LMCDB_PASS'), database='tmsnp')
mycursor = conn.cursor()


print("Downloading data from PFAM...")


mycursor.execute("select distinct acc from snps;")
distinct_acc_tupple = mycursor.fetchall()
distinct_acc_list = list()
for tupple in distinct_acc_tupple:
    for acc in tupple:
        distinct_acc_list.append(acc)


distinct_pfam_list = list()
for distinct_acc in distinct_acc_list:
    mycursor.execute("select distinct pfam from receptor_pfam where acc = '" + distinct_acc + "';")
    distinct_pfam_tupple = mycursor.fetchall()
    for tupple in distinct_pfam_tupple:
        for pfam in tupple:
            if pfam not in distinct_pfam_list:
                distinct_pfam_list.append(pfam)
            else:
                continue

print(distinct_pfam_list)

# curl
total_pfams = len(distinct_pfam_list)
for i, pfam_key in enumerate(distinct_pfam_list):
    print(pfam_key, 'remaining families:', total_pfams - i)
    link = 'http://pfam.xfam.org/family/' + pfam_key + '/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0'
    run_curl(link, pfam_download_path + pfam_key)

print("... Done!!!")
#########################################################################################################################################################################

# Manage errors

print("Trying to get alignments that were not downloaded...")
abecedari = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
llista_pfam_errors = ["control"]
download_counter = 0

while len(llista_pfam_errors) != 0:
    if download_counter < 20:
        pfam_list = os.listdir(pfam_download_path_sys)
        llista_pfam_errors = list()
        for pfam in pfam_list:
            data_llista = file_to_lines(pfam_download_path + pfam)
            for aliniament in data_llista:
                if aliniament[0:9] == '<!DOCTYPE':
                    print('hi ha error en el pfam: ' + pfam)
                    llista_pfam_errors.append(pfam)
        print(llista_pfam_errors)

        for pfam_key in llista_pfam_errors:
            os.remove(pfam_download_path + pfam_key)
            link = 'http://pfam.xfam.org/family/' + pfam_key + '/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0'
            run_curl(link, pfam_download_path + pfam_key)

        download_counter += 1
        print(download_counter)
    else:
        for pfam_dolent in llista_pfam_errors:
            os.remove(pfam_download_path + pfam_dolent)
            outf1 = open(pfam_download_path + pfam_dolent, 'w')
            print('ERROR_ERROR/69-117             ----------------------------------', file=outf1)
            outf1.close()
            print('eliminat per error en la descarrega del servidor (html) pfam el : ' + pfam_dolent)
        break

llista_pfam_errors = ["control"]
download_counter = 0

while len(llista_pfam_errors) != 0:
    if download_counter < 5:
        pfam_list = os.listdir(pfam_download_path_sys)
        llista_pfam_errors = list()
        for pfam in pfam_list:
            data_llista = file_to_lines(pfam_download_path + pfam)
            for aliniament in data_llista:
                if abecedari.find(aliniament[0]) == -1:
                    if pfam not in llista_pfam_errors:
                        llista_pfam_errors.append(pfam)
                        print('hi ha error en el pfam: ' + pfam)
        print(llista_pfam_errors)

        for pfam_key in llista_pfam_errors:
            os.remove(pfam_download_path + pfam_key)
            link = 'http://pfam.xfam.org/family/' + pfam_key +\
                   '/alignment/full/format?format=pfam&alnType=full&order=t&case=u&gaps=dashes&download=0'
            run_curl(link, pfam_download_path + pfam_key)

        download_counter += 1
        print(download_counter)
    else:
        for pfam_dolent in llista_pfam_errors:
            os.remove(pfam_download_path + pfam_dolent)
            outf1 = open(pfam_download_path + pfam_dolent, 'w')
            print('ERROR_ERROR/69-117             ----------------------------------', file=outf1)
            outf1.close()
            print('eliminat i creat en blanc per error en la descarrega del servidor (no completa la descarrega) pfam el : '
                  + pfam_dolent)
        break

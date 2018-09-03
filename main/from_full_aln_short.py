#!/usr/bin/env python3

import os
import gzip
from iker_snp import *

# Dirs where PFAM downloaded files will go
pfam_download_path = './pfam_download/'
pfam_download_path_sys = pfam_download_path[:-1]
pfam_download_raw = './pfam_download_raw/'
pfam_download_raw_sys = pfam_download_raw[:-1]
pfam_full = './Pfam-A.full.gz'

pfam_list_total_raw = os.listdir(pfam_download_raw_sys)

for dist_pfam in pfam_list_total_raw:
    llista_aliniament = file_to_lines(pfam_download_raw + dist_pfam)

    outf1 = open(pfam_download_path + dist_pfam, 'w')
    abecedari = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

    for linia in llista_aliniament:
        if abecedari.find(linia[0]) != -1:
            linia_majus = linia.upper()
            linia_guions = linia_majus.replace(".", "-")
            print(linia_guions.rstrip("\n"), file=outf1)
    outf1.close()

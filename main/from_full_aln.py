#!/usr/bin/env python3

import os
import gzip
from iker_snp import *

# Dirs where PFAM downloaded files will go
path_pfam = "./pfam/pfam_download/"
path_pfam_sys = path_pfam[:-1]

path_raw = "./pfam/path_raw/"
if not os.path.exists(path_raw):
    os.makedirs(path_raw)
path_raw_sys = path_raw[:-1]
pfam_full = "./pfam/Pfam-A.full.gz"

finder = ""
ok_files = os.listdir(path_raw_sys)

miss_pfam = [
    "PF12796",
    "PF00082",
    "PF00168",
    "PF04389",
    "PF00085",
    "PF00106",
    "PF07859",
    "PF03372",
    "PF13853",
    "PF00089",
    "PF00400",
    "PF00300",
    "PF00534",
    "PF00004",
    "PF00443",
    "PF00201",
    "PF16363",
    "PF00884",
    "PF01593",
    "PF04055",
    "PF15843",
    "PF14992",
    "PF15188",
    "PF00535",
    "PF00743",
    "PF00892",
    "PF00501",
    "PF13193",
    "PF01554",
    "PF00149",
    "PF00561",
    "PF12146",
    "PF00171",
    "PF07992",
    "PF13439",
    "PF12697",
    "PF03547",
    "PF00232",
    "PF00071",
    "PF01494",
    "PF00135",
    "PF01425",
    "PF00063",
    "PF00583",
]

for miss in miss_pfam:
    if not miss in ok_files:
        outf1 = open(path_raw + miss, "w")
        with gzip.open(pfam_full, "rb") as FileObj:
            for line in FileObj:
                line = line.decode("ISO-8859-1")
                if line.find("AC   " + miss) != -1:
                    finder = "trobat"
                    print("Found ", miss)
                if finder == "trobat":
                    print(line.rstrip("\n"), file=outf1)
                if line.find("//") != -1 and finder == "trobat":
                    finder = ""
                    print(line)
            outf1.close()

l_raw_total = os.listdir(path_raw_sys)

for dist_pfam in l_raw_total:
    l_align = file_to_lines(path_raw + dist_pfam)

    outf1 = open(path_pfam + dist_pfam, "w")
    abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

    for l in l_align:
        if abc.find(l[0]) != -1:
            major_line = l.upper()
            guide_line = major_line.replace(".", "-")
            print(guide_line.rstrip("\n"), file=outf1)
    outf1.close()

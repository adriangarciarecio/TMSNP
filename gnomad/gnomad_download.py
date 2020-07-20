# Download from GNOMAD non-pathogenic all missense mutations for membrane
#  proteins containing pathogenical mutations

import re
import requests
import urllib.request
from urllib.error import HTTPError, URLError
import sys

# list1 = open("ensembl_list_model.txt")
list1 = open("ensembl_list_all.txt")

gnomad_mut = open("gnomad_mutations.txt", "r")
codes_gnomad = []
for line in gnomad_mut.readlines():
    name = line.split()
    if name[0] not in codes_gnomad:
        codes_gnomad.append(name[0])
# sys.exit()
# RESTART
gnomad_mut = open("gnomad_mutations.txt", "a")
for line in list1.readlines():
    parts = line.split()
    if parts[0] not in codes_gnomad:
        try:
            file_name = parts[1].replace(";", "")
            url = "http://gnomad-old.broadinstitute.org/gene/" + file_name
            print(url)
            try:
                f = urllib.request.urlopen(url)
                data = f.read()
                data2 = str(data)
                data3 = data2.replace('"p', "ppp")
                words = re.findall(r"\b\w*ppp\.\w*\b", data3)
                for word in words:
                    pos_a = data3.find(word)
                    k = data3[pos_a : pos_a + 1000]
                    s = k.find("allele_num")
                    ss = k.find("allele_freq")
                    gnomad_mut.write(parts[0] + " " + word + k[ss + 13 : s - 3] + "\n")
                print(parts[0], parts[1])
            except IOError:
                print("Error")
            except ValueError:
                print("Error")
            except Exception as e:
                print("Error")
        except:
            print("Error code not found!")

# This module contains all the functions developped by Iker Reina for the TM-SNP project

import collections
from math import exp, expm1, log10, log
import subprocess
import sys

def file_to_lines(filename):
    '''
    :param filename: arxiu a llegir
    :rtype: llista de cada linia darxiu
    '''
    fil = open(filename)
    lines = fil.readlines()
    fil.close()
    return lines


def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros


def extract_data_taula_proteines(data_llista):
    '''
    :param data_llista: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [ID uniprot, AC el codi uniprot, GN nom gen, [PFAM]
    '''
    pfam_loc = list()
    general_id = ''
    loc = ''
    codi = ''
    gen = ''
    llista1 = list()
    llista_de_llistes = list()
    for data_string in data_llista:
        if (data_string[0:2] == 'ID' and general_id == ''):
            ID_list = data_string.split('   ')
            general_id = ID_list[1]
        elif (data_string[0:2] == 'AC' and codi == ''):
            AC_list = data_string.split(';')
            codi = AC_list[0][5:]
        elif data_string[0:2] == 'GN' and gen == '':
            GN_list = data_string.split('=')
            gen = GN_list[1][:GN_list[1].find(' ')].strip(';')
        elif data_string[0:2] == 'DR':
            PFAM_list = data_string.split(' ')
            if PFAM_list[3] == 'Pfam;':
                pfam_loc.append(PFAM_list[4][:-1])
        elif data_string[0:2] == 'FT':
            FT_list = data_string.split(' ')
            for n in FT_list:
                if n == 'TRANSMEM':
                    loc = 'transm_loc'
        elif (data_string[0:2] == '//' and loc == 'transm_loc'):
            llista1.append(general_id)
            llista1.append(codi)
            llista1.append(gen)
            llista1.append(pfam_loc)
            llista_de_llistes.append(llista1)
            pfam_loc = []
            llista1 = []
            general_id = ''
            codi = ''
            gen = ''
            loc = ''
        elif (data_string[0:2] == '//' and loc == ''):
            pfam_loc = []
            llista1 = []
            general_id = ''
            codi = ''
            gen = ''
        else:
            continue
    return llista_de_llistes


def extract_data_taula_snip_regions(data_llista):
    '''
    :param data_llista: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR]] i descarta les proteines sense regio transmembrana
    '''
    general_id = ''
    gen = ''
    id_var_loc = ''
    disease_match = ''
    loc = ''
    codi = ''
    disease = list()
    sub_loc = list()
    variant_loc = list()
    variant_change = list()
    variant_id = list()
    variant_rs = list()
    llista1 = list()
    llista_de_llistes = list()
    for data_string in data_llista:
        if (data_string[0:2] == 'ID' and general_id == ''):
            ID_list = data_string.split('   ')
            general_id = ID_list[1]
        elif data_string[0:2] == 'GN' and gen == '':
            GN_list = data_string.split('=')
            gen = GN_list[1][:GN_list[1].find(' ')].strip(';')
        elif (data_string[0:2] == 'AC' and codi == ''):
            AC_list = data_string.split(';')
            codi = AC_list[0][5:]
        elif (data_string[0:17] == 'CC   -!- DISEASE:'):
            disease_match = 'trobat'
            if data_string.find('(') == -1:
                disease_match = 'trobat_dolent'
            else:
                disease.append(data_string[data_string.find('(') + 1:data_string.find(')')])
        elif (data_string[0:9] == 'CC       ' and disease_match == 'trobat_dolent'):
            if data_string.find('(') != -1:
                disease_match = 'trobat'
                disease.append(data_string[data_string.find('(') + 1:data_string.find(')')])
            else:
                disease_match = 'trobat_dolent'
        elif data_string[0:2] == 'FT':
            if data_string[5:13] == 'TRANSMEM':
                loc = 'transm_loc'
                sub_loc.append(data_string[13:27])
            elif data_string[5:12] == 'VARIANT':
                variant_rs.append('-')
                id_var_loc = 'trobat'
                variant_loc.append(data_string[20:-2])
            if (id_var_loc == 'trobat' and data_string.find('->') != -1):
                variant_change.append(data_string[data_string.find('->') - 3:data_string.find('->') + 5].strip(' \n .'))
            if (id_var_loc == 'trobat' and data_string.find('Missing') != -1):
                variant_change.append('MISS -> MISS')
            if (id_var_loc == 'trobat' and data_string.find('/FTId=') != -1):
                variant_id.append(data_string[40:50])
                id_var_loc = ''
            if (id_var_loc == 'trobat' and data_string.find('dbSNP:rs') != -1):
                variant_rs.pop(len(variant_rs) - 1)
                variant_rs.append(
                    data_string[data_string.find('dbSNP:rs') + 6:data_string.find('dbSNP:rs') + 17].rstrip(').\n'))
        elif (data_string[0:2] == '//' and loc == 'transm_loc'):
            llista1.append(codi)
            llista1.append(sub_loc)
            llista1.append(variant_loc)
            llista1.append(variant_change)
            llista1.append(variant_id)
            llista1.append(variant_rs)
            llista1.append(disease)
            llista1.append(general_id)
            llista1.append(gen)
            llista_de_llistes.append(llista1)
            sub_loc = []
            variant_loc = []
            variant_id = []
            variant_rs = []
            llista1 = []
            variant_change = []
            disease = []
            codi = ''
            loc = ''
            general_id = ''
            gen = ''
        elif (data_string[0:2] == '//' and loc == ''):
            sub_loc = []
            variant_loc = []
            variant_id = []
            variant_rs = []
            llista1 = []
            variant_change = []
            disease = []
            codi = ''
            general_id = ''
            gen = ''
        else:
            continue
    return llista_de_llistes


def descart_regio_no_transmembrana(data_locations_net):
    '''
    :param data_locations:  una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR]]
    :return: una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR], [llista posicio snips transmembraana]]
    '''
    llist = list()
    final_list = list()
    snips_sols = list()
    for num_proteina in range(0, len(data_locations_net)):
        nom_uniprot = data_locations_net[num_proteina][0]
        posicions_transmem = data_locations_net[num_proteina][1]
        snp_rs = data_locations_net[num_proteina][2]
        canvi_mutacio = data_locations_net[num_proteina][3]
        snp_var = data_locations_net[num_proteina][4]
        snp_rs_sol = data_locations_net[num_proteina][5]
        disease = data_locations_net[num_proteina][6]
        id = data_locations_net[num_proteina][7]
        gen = data_locations_net[num_proteina][8]
        if (len(snp_rs) != len(snp_var)):
            print('hi ha un problema amb els snps, no tots tenen codi VAR_... ')
            print('la proteina es: ' + id)
            sys.exit()
        if (len(snp_var) != len(snp_rs_sol)):
            print('hi ha un problema amb els snps, no es crean be els codis rs dels snips... ')
            sys.exit()
        if (len(snp_rs_sol) != len(canvi_mutacio)):
            print(snp_rs_sol)
            print(canvi_mutacio)
            print('hi ha un problema amb els canvis de la mutacio (aminoacid) no tots apareix el canvi....... ')
            sys.exit()
        llista_snp_transmem = zerolistmaker(len(snp_rs))
        llista_snp_disease = zerolistmaker(len(snp_rs))
        for snips in snp_rs:
            snips_sols.append(snips[:7].strip())
            for a_disease in disease:
                if (snips.find(a_disease) != -1):
                    llista_snp_disease.pop(snp_rs.index(snips))
                    llista_snp_disease.insert(snp_rs.index(snips), 1)
        for coordenades in posicions_transmem:
            primera_coordenada_transmembrana = coordenades[:7].strip()
            segona_coordenada_transmembrana = coordenades[7:].strip()
            for snips in snp_rs:
                posicio_snp = snips[:7].strip()
                if (int(posicio_snp) >= int(primera_coordenada_transmembrana) and int(posicio_snp) <= int(
                        segona_coordenada_transmembrana)):
                    llista_snp_transmem.pop(len(llista_snp_transmem) - 1)
                    llista_snp_transmem.insert(snp_rs.index(snips), 1)
        if (len(snp_var) != len(llista_snp_disease)):
            print('no tenen mateixa llargada la llista 0')
        if (len(llista_snp_disease) != len(llista_snp_transmem)):
            print('no tenen mateixa llargada la llista 0')
        llist.append(nom_uniprot)
        llist.append(posicions_transmem)
        llist.append(snp_var)
        llist.append(canvi_mutacio)
        llist.append(snp_rs_sol)
        llist.append(snips_sols)
        llist.append(llista_snp_transmem)
        llist.append(llista_snp_disease)
        llist.append(id)
        llist.append(gen)
        final_list.append(llist)
        llist = []
        snips_sols = []
    return final_list


def run_curl(url, file_out):
    proc = subprocess.Popen(["curl", url], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    out = out.decode('utf-8')
    with open(file_out, 'w') as fileout:
        fileout.write(out)


def crear_llista_llistes_proteines(data_llista):
    llista1 = list()
    llista_prot_amino = list()
    for proteina in data_llista:
        llista_proteina = proteina.split(' ')
        if (llista_proteina[0].find('/') != -1):
            pos_barra = llista_proteina[0].find('/')
        else:
            print('hi ha un error en el codi de la proteina no es trova la separacio entre el nom i la posicio (/)')
            print('lerror esta en la proteina: ' + proteina)
            break
        llista1.append(llista_proteina[0][:pos_barra])
        llista1.append(llista_proteina[0][pos_barra + 1:])
        for posicio in llista_proteina:
            if (posicio.find(llista_proteina[0]) == -1 and posicio != ''):
                llista1.append(posicio.rstrip('\n'))
        llista_prot_amino.append(llista1)
        llista1 = []
    return llista_prot_amino


def crear_llista_idprotein_rs_amino_posiciototal(llista_prot_posicions_amino, id_protein, snip_var, snip_position):
    llista_prot_posicio_total = list()
    llista1 = list()
    for num_proteina in range(0, len(llista_prot_posicions_amino)):
        nom_proteina = llista_prot_posicions_amino[num_proteina][0]
        coord_aminoacids = llista_prot_posicions_amino[num_proteina][1]
        aminoacids = llista_prot_posicions_amino[num_proteina][2]
        if (nom_proteina == id_protein and snip_position >= int(
                coord_aminoacids[:coord_aminoacids.find('-')]) and snip_position <= int(
                coord_aminoacids[coord_aminoacids.find('-') + 1:])):
            counter_1 = int(coord_aminoacids[:coord_aminoacids.find('-')])
            counter_2 = 0
            llista1.append(nom_proteina)
            llista1.append(snip_var)
            for caracter in aminoacids:
                if (caracter != '-' and counter_1 <= snip_position):
                    counter_1 = counter_1 + 1
                    caracter_final = caracter
                if (counter_1 <= snip_position):
                    counter_2 = counter_2 + 1
            llista1.append(caracter_final)
            llista1.append(counter_2)
            llista_prot_posicio_total.append(llista1)
            llista1 = []
    return llista_prot_posicio_total


def crear_llista_prot_rs_amino_distinctprot(llista_prot_posicio_total, llista_prot_posicions_amino):
    llista_distinct_prot = list()
    llista1 = list()
    llista_prot_rs_amino_distinctprot = list()
    for numero_proteina in range(0, len(llista_prot_posicio_total)):
        nom_proteina = llista_prot_posicio_total[numero_proteina][0]
        snip_rs = llista_prot_posicio_total[numero_proteina][1]
        amino = llista_prot_posicio_total[numero_proteina][2]
        posicio = llista_prot_posicio_total[numero_proteina][3]
        llista1.append(nom_proteina)
        llista1.append(snip_rs)
        llista1.append(amino)
        for num_proteina in range(0, len(llista_prot_posicions_amino)):
            aminoacids = llista_prot_posicions_amino[num_proteina][2]
            llista_distinct_prot.append(aminoacids[posicio])
        llista1.append(llista_distinct_prot)
        llista_prot_rs_amino_distinctprot.append(llista1)
        llista_distinct_prot = []
        llista1 = []
    return llista_prot_rs_amino_distinctprot


def calcul_freq_entropia(llista_prot_rs_amino_distinctprot, entropia_inicial, distinct_aa):
    aminos_llista = ['Y', 'G', 'F', 'M', 'A', 'S', 'I', 'L', 'T', 'V', 'P', 'K', 'H', 'Q', 'E', 'Z', 'W', 'R', 'D', 'N',
                     'B', 'C', 'X']
    for numero_regio_pfam in range(0, len(llista_prot_rs_amino_distinctprot)):
        diccionario_freq_abs = dict(collections.Counter(llista_prot_rs_amino_distinctprot[numero_regio_pfam][3]))
        if '-' not in diccionario_freq_abs:
            gap = 0
        else:
            gap = diccionario_freq_abs['-']
            del diccionario_freq_abs['-']
        if 'X' not in diccionario_freq_abs:
            equis = 0
        else:
            equis = diccionario_freq_abs['X']
            del diccionario_freq_abs['X']
        if 'Z' in diccionario_freq_abs:
            print('hi ha una Z en laliniament')
        if 'B' in diccionario_freq_abs:
            print('hi ha una B en laliniament: ')
        llargada_sense_gap = len(llista_prot_rs_amino_distinctprot[numero_regio_pfam][3]) - (
                    gap + equis)  # li trec els gaps i les x a la llargada total de les posicions
        fx = diccionario_freq_abs[llista_prot_rs_amino_distinctprot[numero_regio_pfam][2]] / llargada_sense_gap
        #        print(fx)
        diccionario_freq_relatives = dict()
        for aminoacid in diccionario_freq_abs:
            diccionario_freq_relatives[aminoacid] = diccionario_freq_abs[aminoacid] / llargada_sense_gap
        #        print(diccionario_freq_relatives)
        if distinct_aa in diccionario_freq_relatives:
            fx_amino_nou = diccionario_freq_relatives[distinct_aa]
        else:
            fx_amino_nou = 0
        #        print(fx_amino_nou)
        suma = 0.0
        for aminoacid_dict in diccionario_freq_relatives:
            suma += diccionario_freq_relatives[aminoacid_dict] * log(diccionario_freq_relatives[aminoacid_dict])
        #        print(entropia_inicial)
        #        print(suma)
        entropia = entropia_inicial + suma
        #        print(entropia)
        return fx, fx_amino_nou, entropia


def extract_data_taula_proteines_entropia(data_llista):
    '''
    :param data_llista: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [ID uniprot, AC el codi uniprot, GN nom gen, [PFAM]
    '''
    dna_loc = list()
    loc = ''
    codi = ''
    llista1 = list()
    llista_de_llistes = list()
    for data_string in data_llista:
        if (data_string[0:2] == 'AC' and codi == ''):
            AC_list = data_string.split(';')
            codi = AC_list[0][5:]
        elif data_string[0:2] == '  ':
            adn = data_string.strip('\n ')
            dna_loc.append(adn.replace(" ", ""))
        elif data_string[0:2] == 'FT':
            FT_list = data_string.split(' ')
            for n in FT_list:
                if n == 'TRANSMEM':
                    loc = 'transm_loc'
        elif (data_string[0:2] == '//' and loc == 'transm_loc'):
            llista1.append(codi)
            llista1.append(dna_loc)
            llista_de_llistes.append(llista1)
            dna_loc = []
            llista1 = []
            codi = ''
            loc = ''
        elif (data_string[0:2] == '//' and loc == ''):
            dna_loc = []
            llista1 = []
            codi = ''
        else:
            continue
    return llista_de_llistes


def crear_llista_prot_seqamino_entropia(data_locations_net):
    '''
    :param data_locations:  una llista de llist
    :return: una llista de llistes
    '''
    llist = list()
    final_list = list()
    string_adn = ''
    for num_proteina in range(0, len(data_locations_net)):
        nom_uniprot = data_locations_net[num_proteina][0]
        dna_full = data_locations_net[num_proteina][1]
        for adn in dna_full:
            string_adn += adn
        llist.append(nom_uniprot)
        llist.append(string_adn)
        final_list.append(llist)
        string_adn = ''
        llist = []
    return final_list


def extract_data_taula_snip_regions_entropia(data_llista):
    '''
    :param data_llista: llista de cada linia del arxiu text de Uniprot
    :return: una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR]] i descarta les proteines sense regio transmembrana
    '''
    id_var_loc = ''
    loc = ''
    codi = ''
    sub_loc = list()
    variant_loc = list()
    variant_id = list()
    variant_rs = list()
    llista1 = list()
    llista_de_llistes = list()
    for data_string in data_llista:
        if (data_string[0:2] == 'AC' and codi == ''):
            AC_list = data_string.split(';')
            codi = AC_list[0][5:]
        elif data_string[0:2] == 'FT':
            if data_string[5:13] == 'TRANSMEM':
                loc = 'transm_loc'
                sub_loc.append(data_string[13:27])
            elif data_string[5:12] == 'VARIANT':
                variant_rs.append('-')
                id_var_loc = 'trobat'
                variant_loc.append(data_string[20:-2])
            if (id_var_loc == 'trobat' and data_string.find('/FTId=') != -1):
                variant_id.append(data_string[40:50])
                id_var_loc = ''
            if (id_var_loc == 'trobat' and data_string.find('dbSNP:rs') != -1):
                variant_rs.pop(len(variant_rs) - 1)
                variant_rs.append(
                    data_string[data_string.find('dbSNP:rs') + 6:data_string.find('dbSNP:rs') + 16].rstrip(').\n'))
        elif (data_string[0:2] == '//' and loc == 'transm_loc'):
            llista1.append(codi)
            llista1.append(sub_loc)
            llista1.append(variant_loc)
            llista1.append(variant_id)
            llista1.append(variant_rs)
            llista_de_llistes.append(llista1)
            sub_loc = []
            variant_loc = []
            variant_id = []
            variant_rs = []
            llista1 = []
            codi = ''
            loc = ''
        elif (data_string[0:2] == '//' and loc == ''):
            sub_loc = []
            variant_loc = []
            variant_id = []
            variant_rs = []
            llista1 = []
            codi = ''
        else:
            continue
    return llista_de_llistes


def descart_regio_no_transmembrana_entropia(data_locations_net):
    '''
    :param data_locations:  una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR]]
    :return: una llista de llistes amb: [Accesion de el codi uniprot, [ llista de posicions transmembrana], [llista de SNPs RS]], [llista snips VAR], [llista posicio snips transmembraana]]
    '''
    llist = list()
    final_list = list()
    coord_trans = list()
    snips_sols = list()
    for num_proteina in range(0, len(data_locations_net)):
        nom_uniprot = data_locations_net[num_proteina][0]
        posicions_transmem = data_locations_net[num_proteina][1]
        snp_rs = data_locations_net[num_proteina][2]
        snp_var = data_locations_net[num_proteina][3]
        snp_rs_sol = data_locations_net[num_proteina][4]
        if (len(snp_rs) != len(snp_var)):
            print('hi ha un problema amb els snps, no tots tenen codi VAR_... ')
            print('la proteina es: ' + snp_var)
            sys.exit()
        if (len(snp_var) != len(snp_rs_sol)):
            print('hi ha un problema amb els snps, no es crean be els codis rs dels snips... ')
            sys.exit()
        llista_snp_transmem = zerolistmaker(len(snp_rs))
        for snips in snp_rs:
            snips_sols.append(snips[:7].strip())
        for coordenades in posicions_transmem:
            coord_trans.append(coordenades)
            primera_coordenada_transmembrana = coordenades[:7].strip()
            segona_coordenada_transmembrana = coordenades[7:].strip()
            for snips in snp_rs:
                posicio_snp = snips[:7].strip()
                if (int(posicio_snp) >= int(primera_coordenada_transmembrana) and int(posicio_snp) <= int(
                        segona_coordenada_transmembrana)):
                    llista_snp_transmem.pop(len(llista_snp_transmem) - 1)
                    llista_snp_transmem.insert(snp_rs.index(snips), 1)
        llist.append(nom_uniprot)
        llist.append(posicions_transmem)
        llist.append(snp_var)
        llist.append(snp_rs_sol)
        llist.append(snips_sols)
        llist.append(llista_snp_transmem)
        final_list.append(llist)
        llist = []
        coord_trans = []
        snips_sols = []
    return final_list

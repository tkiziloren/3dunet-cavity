import os.path
from pathlib import Path
import logging
import sys

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

SRC_SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/src_data/scPDB"


def find_number_of_chains(protein_path):
    with open(protein_path) as fp:
        lines = fp.readlines()
        started = False
        chains = set()
        for item in range(len(lines)):
            line = lines[item]
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                started = True
                continue

            if started:
                if line.startswith("   "):
                    splitted = line.split()
                    chains.add(splitted[5])
                    if item + 1 == len(lines): # in case there is not any line after "@<TRIPOS>SUBSTRUCTURE" group. Eg check "2i35_1"
                        return chains
                else:
                    return chains


def get_list_of_proteins(txt_file_path):
    with open(txt_file_path) as fp:
        lines = fp.readlines()
        chains = set()
        for item in range(len(lines)):
            line = lines[item]
            chains.add(line.replace("\n", ""))

        return chains


def find_chains():
    dir_names = get_list_of_proteins("kept_suffixes_without_problematic.txt")
    count = 0

    dict_protein_groups = dict()

    for dir_name in dir_names:
        protein_path = os.path.join(SRC_SC_PDB_PATH, dir_name, "protein.mol2")
        protein_name = os.path.basename(dir_name)
        chains = find_number_of_chains(protein_path)

        number_of_chains = 0
        if chains is not None:
            number_of_chains = len(chains)

        chains_as_str = str(chains)
        count = count + 1
        if count % 100 == 0:
            print(count)

        if number_of_chains in dict_protein_groups:
            dict_protein_groups[number_of_chains].append(protein_name)
        else:
            dict_protein_groups[number_of_chains] = [protein_name]

        logging.info("{} {} {}".format(protein_name, number_of_chains, chains_as_str))

    for key in dict_protein_groups:

        with open("./groups/{}.txt".format(key), "a") as f:
            proteins = dict_protein_groups[key]
            for p in proteins:
                f.write(p + '\n')

        logging.info("There are {} proteins with {} chains:".format(len(dict_protein_groups[key]), key))


find_chains()

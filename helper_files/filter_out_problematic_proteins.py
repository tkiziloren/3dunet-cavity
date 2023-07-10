import os.path
from pathlib import Path
import logging
import sys

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

ALL_PROTEIN_LIST = "./kept_scPDB_suffixes.txt"
PROBLEMATIC_PROTEIN_LIST = "./problematic_proteins_for_gen_cache.txt"


def get_list_of_proteins(txt_file_path):
    with open(txt_file_path) as fp:
        lines = fp.readlines()
        chains = set()
        for item in range(len(lines)):
            line = lines[item]
            chains.add(line.replace("\n", ""))

        return chains


def run():
    all = get_list_of_proteins(ALL_PROTEIN_LIST)
    problematic = get_list_of_proteins(PROBLEMATIC_PROTEIN_LIST)
    filtered = all.difference(problematic)

    for item in filtered:
        print(item)

    logging.info("Size of all:{}, problematic:{}, filtered:{}".format(len(all), len(problematic), len(filtered)))




run()

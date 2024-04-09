import os.path
from pathlib import Path
import logging
import sys

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

ALL_PROTEIN_LIST = "../kept_suffixes_without_problematic.txt"
DRUGGABLE_PROTEIN_LIST = "../scPDB_druggable_without_duplicate_proteins.txt"


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
    druggable = get_list_of_proteins(DRUGGABLE_PROTEIN_LIST)

    count = 0
    x = set()
    duplicates = set()
    for text in druggable:
        for item in all:
            if text in item:
                if text in x:
                    duplicates.add(text)
                x.add(text)

    for item in duplicates:
        druggable.remove(item)

    x.clear()
    for text in druggable:
        for item in all:
            if text in item:
                count +=1
                x.add(text)
                #print("{}, {}".format(text, item))
                print("- {}".format(item))

    print(count)
    print(len(x))

    # print(len(x))
    # for item in x:
    #     print(item)

    # for item in filtered:
    #     print(item)
    #
    # logging.info("Size of all:{}, problematic:{}, filtered:{}".format(len(all), len(druggable), len(filtered)))
    #

run()

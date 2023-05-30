import os.path
import shutil
from tqdm import trange
import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("src", help="Source PDB Dataset Location")
parser.add_argument("dest", help="Destination PDB Dataset Location, where filtered items wil be stored.")
parser.add_argument("filter_file_path", help="Path of the file that includes the name of the proteins to be filtered.")
args = parser.parse_args()

SRC_SC_PDB_PATH = args["src"]
DST_SC_PDB_PATH = args["dest"]
FILTER_FILE_PATH = args["filter_file_path"]

with open(FILTER_FILE_PATH) as fp:
    lines = fp.readlines()
    for item in trange(len(lines)):
        line = lines[item]
        line = line.replace("\n", "")
        shutil.copytree(os.path.join(SRC_SC_PDB_PATH, line), os.path.join(DST_SC_PDB_PATH,line), dirs_exist_ok=True)
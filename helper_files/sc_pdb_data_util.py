import argparse
import os.path
import shutil
from pathlib import Path

from tqdm import trange

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-src", help="Source PDB Dataset Location")
parser.add_argument("-dest", help="Destination PDB Dataset Location, where filtered items wil be stored.")
parser.add_argument("-filter_file_path", help="Path of the file that includes the name of the proteins to be filtered.")
args = parser.parse_args()

SRC_SC_PDB_PATH = "/Users/tevfik/Sandbox/github/PHD/data/caches/scPDB_filtered"
DST_SC_PDB_PATH = "/Users/tevfik/Sandbox/github/PHD/data/caches/scPDB_filtered_done"

p = Path(SRC_SC_PDB_PATH)
dir_names = [f for f in p.iterdir() if f.is_dir()]

for item in trange(len(dir_names)):
    dir_path = dir_names[item]
    done_path = os.path.join(dir_path, "GenNumPyGrids.done")
    dir_name = os.path.basename(dir_path)
    if os.path.exists(done_path):
        src_path = dir_path
        dst_path = os.path.join(DST_SC_PDB_PATH, dir_name)
        # print(src_path)
        # print(dest_path)
        shutil.move(src_path, dst_path, copy_function = shutil.copytree)



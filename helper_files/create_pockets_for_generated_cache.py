# gen_cache doesn't create pocket information, this script calculates and generate pockets for given generated cache.
import os.path
from pathlib import Path
from tqdm import trange

import prody

SRC_SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/caches/scPDB_filtered"

p = Path(SRC_SC_PDB_PATH)
dir_names = [f for f in p.iterdir() if f.is_dir()]

for item in trange(len(dir_names)):
    dir_path = dir_names[item]
    protein_name = os.path.basename(dir_path)

    destPdbFile = os.path.join(dir_path, "protein_trans.pdb")
    destLigandPdbFile = os.path.join(dir_path, "ligand.pdb")
    destPocketFile = os.path.join(dir_path, "pocket.pdb")
    destSelectedPdbFile = os.path.join(dir_path, "selected.pdb")

    protein = prody.parsePDB(destPdbFile)
    ligand = prody.parsePDB(destLigandPdbFile)
    lresname = ligand.getResnames()[0]
    complx = ligand + protein
    complx = complx.select(f'same chain as exwithin 7 of resname {lresname}')

    structure = complx.select(f'protein and not resname {lresname}')
    prody.writePDB(destSelectedPdbFile, structure)
    # pdb2pqr fails to read pdbs with the one line header generated by ProDy...
    with open(destSelectedPdbFile, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(destSelectedPdbFile, 'w') as fout:
        fout.writelines(data[1:])
    selected = prody.parsePDB(destSelectedPdbFile)

    complx = ligand + selected
    lresname = ligand.getResnames()[0]
    pocket = complx.select(f'same residue as exwithin 4.5 of resname {lresname}')
    prody.writePDB(destPocketFile, pocket)
    os.remove(destSelectedPdbFile)


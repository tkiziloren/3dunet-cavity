import logging
import os.path
import sys
from openbabel import openbabel
from pathlib import Path

from prody import *

from pytorch3dunet.unet3d.metrics import *
from pytorch3dunet.unet3d.utils import atomgroup_to_molecule, molecule_to_atom_group

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

# RUN_NAME = "scPDB_nonRedundant_druggable_1600-200-200"
# RUNS_PATH = "../runs"
# SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/src_data/scPDB_converted"

RUN_NAME = "scPDB_druggable_kalasanty_1600-200-200"
RUNS_PATH = "../runs"
SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/src_data/scPDB_converted"


def find_dvo_and_pli():
    predictions_path = os.path.join(RUNS_PATH, RUN_NAME, "predictions")
    predictions_path = Path(predictions_path)
    dir_names = [f for f in predictions_path.iterdir() if f.is_dir()]
    result = list()
    total = 0
    count_dvo = 0
    count_pli = 0
    dvo_less_than4 = list()
    pli_less_than4 = list()
    for item in range(len(dir_names)):
        dir_path = dir_names[item]
        protein_name = os.path.basename(dir_path)

        prediction_path = os.path.join(dir_path, "pocket_pred.pdb")
        pocket_path = os.path.join(SC_PDB_PATH, protein_name, protein_name + "_pocket.pdb")
        ligand_mol2_path = os.path.join(SC_PDB_PATH, protein_name, protein_name + "_ligand.mol2")
        ligand = next(pybel.readfile("mol2", ligand_mol2_path))

        prediction = prody.parsePDB(prediction_path)
        pocket = prody.parsePDB(pocket_path)

        prediction_molecule = atomgroup_to_molecule(prediction)
        pocket_molecule = atomgroup_to_molecule(pocket)

        dvo = get_DVO(prediction_molecule, pocket_molecule)
        pli = get_PLI(ligand, prediction_molecule)
        centre_prediction = prody.calcCenter(prediction)
        centre_pocket = prody.calcCenter(pocket)

        distance = prody.calcDistance(centre_pocket, centre_prediction)
        if distance <= 4:
            dvo_less_than4.append((protein_name, dvo))
            pli_less_than4.append((protein_name, pli))

        result.append((protein_name, dvo))
        logging.info("Protein name:{} DVO:{}, PLI".format(protein_name, dvo))

    print("DVOs of the proteins with distance less than 4 angstroms")
    for item in dvo_less_than4:
        print(item)

    print("-----------------")

    print("PLIs of the proteins with distance less than 4 angstroms")
    for item in pli_less_than4:
        print(item)


    return result

res = find_dvo_and_pli()
# res.sort(key=lambda tup: tup[1])

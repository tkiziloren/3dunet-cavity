import logging
import os.path
import sys
from pathlib import Path

from prody import *

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

RUN_NAME = "scPDB_druggable_1600-200-200"
RUNS_PATH = "../runs"
SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/caches/scPDB_filtered"



# RUN_NAME = "scPDB_filtered800_100_100"
# RUNS_PATH = "../runs"
# SC_PDB_PATH = "/Volumes/Data/DiskYedek/DATA/caches/scPDB_filtered"


def find_dcc():
    predictions_path = os.path.join(RUNS_PATH, RUN_NAME, "predictions_old")
    predictions_path = Path(predictions_path)
    dir_names = [f for f in predictions_path.iterdir() if f.is_dir()]
    result = list()
    count_less_than_4 = 0
    for item in range(len(dir_names)):
        dir_path = dir_names[item]
        protein_name = os.path.basename(dir_path)
        prediction_path = os.path.join(dir_path, "pocket_pred.pdb")
        pocket_path = os.path.join(SC_PDB_PATH, protein_name, "pocket.pdb")
        prediction = prody.parsePDB(prediction_path)
        pocket = prody.parsePDB(pocket_path)
        centre_prediction = prody.calcCenter(prediction)
        centre_pocket = prody.calcCenter(pocket)

        distance = prody.calcDistance(centre_pocket, centre_prediction)

        result.append((protein_name, distance))
        logging.info("Protein name:{} Distance:{}".format(protein_name, distance))
        if distance <= 5:
            count_less_than_4 += 1

    dcc = 100 * (count_less_than_4 / len(result))
    logging.info("Count of the predictions_old with distance less than 4 angstroms:{}, Count Total:{}, DCC:{}".format(count_less_than_4, len(result), dcc))
    return result


res = find_dcc()

# res.sort(key=lambda tup: tup[1])

import prody
from Bio.PDB import *

SCPDB_PATH = '/Volumes/Data/Disk Yedek/DATA/src_data/scPDB'
SCPDB_CONVERTED_PATH = '/Volumes/Data/Disk Yedek/DATA/src_data/scPDB_converted'


proteins = ['1xdk_2', '2cjf_2', '2b3d_2', '1zxc_1', '1wc4_2']


for item in proteins:
    mol2_path = f"{SCPDB_PATH}/{item}/protein.mol2"
    pdb_path =  f"{SCPDB_CONVERTED_PATH}/{item}/{item}_protein.pdb"
    # structure = prody.parsePQR(pdb_path)
    polymer = prody.parsePDBHeader(pdb_path)

    pdbl = PDBList()
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(item, pdb_path)
    print("Done")



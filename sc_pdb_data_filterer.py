import os.path
import shutil
from pathlib import Path
from openbabel import openbabel
from tqdm import trange


SRC_SC_PDB_PATH = "/Users/tevfik/Sandbox/github/PHD/data/src_data/scPDB"
DST_SC_PDB_PATH = "/Users/tevfik/Sandbox/github/PHD/data/src_data/scPDB_converted"

p = Path(SRC_SC_PDB_PATH)

# All subdirectories in the current directory, not recursive.
dir_names = [f for f in p.iterdir() if f.is_dir()]
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("mol2", "pdb")

for item in trange(len(dir_names)):
    protein_path = dir_names[item]
    protein_name = Path(dir_names[item]).stem
    dst_path = os.path.join(DST_SC_PDB_PATH, protein_name)
    os.makedirs(dst_path, exist_ok=True)

    # Copy ligands mol2 and sdf if not exists
    dst_ligand_name = protein_name + "_ligand"
    src_ligand_sdf = os.path.join(protein_path, "ligand.sdf")
    dst_ligand_sdf = os.path.join(dst_path, dst_ligand_name + ".sdf")
    src_ligand_mol2 = os.path.join(protein_path, "ligand.mol2")
    dst_ligand_mol2 = os.path.join(dst_path, dst_ligand_name + ".mol2")

    if os.path.exists(src_ligand_mol2):
        if not os.path.exists(os.path.join(dst_path, dst_ligand_mol2)):
            shutil.copyfile(src_ligand_mol2, dst_ligand_mol2)
    else:
        print (f"Ligand mol2 file doesn't exist in the path:{src_ligand_mol2}. Skipping...")

    if os.path.exists(src_ligand_sdf):
        if not os.path.exists(os.path.join(dst_path, dst_ligand_sdf)):
            shutil.copyfile(src_ligand_sdf, dst_ligand_sdf)
    else:
        print(f"Ligand sdf file doesn't exist in the path:{src_ligand_sdf}. Skipping...")

    # Convert protein to pdb format from mol2 format and copy into destination
    dst_protein_name = protein_name + "_protein.pdb"
    src_protein_mol2 = os.path.join(protein_path, "protein.mol2")
    dst_protein_pdb = os.path.join(dst_path, dst_protein_name)
    if not os.path.exists(dst_protein_pdb):
        protein = openbabel.OBMol()
        obConversion.ReadFile(protein, src_protein_mol2)
        obConversion.WriteFile(protein, dst_protein_pdb)


    # Convert binding_site to pdb format from mol2 format and copy into destination
    dst_pocket_name = protein_name + "_pocket.pdb"
    src_pocket_mol2 = os.path.join(protein_path, "site.mol2")
    dst_pocket_pdb = os.path.join(dst_path, dst_pocket_name)
    if not os.path.exists(dst_pocket_pdb):
        pocket = openbabel.OBMol()
        obConversion.ReadFile(pocket, src_pocket_mol2)
        obConversion.WriteFile(pocket, dst_pocket_pdb)

    # Copy other helper files (in case they exist)
    other_file_names = ["cavity6.mol2", "cavityALL.mol2", "IFP.txt", "ints_M.mol2"]
    for file_name in other_file_names:
        src_file_path = os.path.join(protein_path, file_name)
        dst_file_path = os.path.join(dst_path, file_name)
        if os.path.exists(src_file_path) and not os.path.exists(dst_file_path):
            shutil.copyfile(src_file_path, dst_file_path)







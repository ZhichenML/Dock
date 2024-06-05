from rdkit import Chem
from rdkit.Chem import AllChem
import os

def smi2xyz(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
    AllChem.UFFOptimizeMolecule(mol)
    atoms = mol.GetAtoms()
    string = "\n"
    for i, atom in enumerate(atoms):
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
    string += "units angstrom\n"
    return string, mol


def read_ligand_smi(fpath= '/public/home/gongzhichen/code/Lingo3DMol/data/output/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/0_pred_3_GLP1_6xox_lily_pocket.pdb.mol'):
    with open(fpath, 'r') as f:
        smi = f.readlines()[0]
    return smi


def read_ligand_batch(fpath= '/public/home/gongzhichen/code/Lingo3DMol/data/output/Generated_Ligands/GLP1_6xox_lily_pocket.pdb'):
    ligand_files = os.listdir(fpath)
    ligands = []
    new_path = '/public/home/gongzhichen/code/Lingo3DMol/data/output/Generated_Ligands/GLP1_6xox_lily_pocket.pdb.revise'
    if os.path.exists(new_path):
        os.system('rm -rf {}'.format(new_path))
    os.makedirs(new_path, exist_ok=True)
    for ligand_file in ligand_files:
        if ligand_file.endswith('.mol'):

            smile = read_ligand_smi(os.path.join(fpath, ligand_file))
            mol = Chem.MolFromSmiles(smile, sanitize=True)
            mol = Chem.RWMol(mol)
            if mol is not None:
                for atom in mol.GetAtoms():
                    atom_type = atom.GetSymbol()
                    if atom_type == '*':
                        mol.ReplaceAtom(atom.GetIdx(),Chem.Atom('H'))
                
                # xyz, mol = smi2xyz(mol)
                # ligand_file = ligand_file.replace('.mol', '.xyz')
                new_ligand_file = os.path.join(new_path, ligand_file)
                # with open(new_ligand_file, 'w') as f:
                #     f.write(''.join(xyz))
                Chem.MolToMolFile(Chem.RemoveHs(mol), new_ligand_file)
            
    return 


# def pymol():
#     f"pymol -c -p <<_EOF                                                                                  
#         load 701_pred_5_GLP1_6xox_lily_pocket.pdb.mol 
#         save 701_pred_5_GLP1_6xox_lily_pocket.pymol.pdb, 701_pred_5_GLP1_6xox_lily_pocket.pdb
#         quit
#         _EOF"
    
def main():
    read_ligand_batch()

if __name__ == '__main__':
    main()
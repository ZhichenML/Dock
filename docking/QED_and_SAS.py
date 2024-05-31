
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski, QED, rdMolDescriptors, Descriptors
from rdkit.Chem import RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer



def calculate_mol_property(mol_file='/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6x1a_poc.pdb/9_pred_3_GLP1_6x1a_poc.pdb.mol', verbose=False):

    # Read the .mol file
    if mol_file.endswith('.pdb'):
        molecule = Chem.MolFromPDBFile(mol_file)
    else:
        molecule = Chem.MolFromMolFile(mol_file)
    logP = Descriptors.MolLogP(molecule)

    Chem.AssignStereochemistry(molecule, force=True, cleanIt=True)

    # Penalized LogP
    penalized_logp = Crippen.MolLogP(molecule, True)
    # Molecular weight
    molecular_weight_exact = rdMolDescriptors.CalcExactMolWt(molecule)
    molecular_weight = Descriptors.MolWt(molecule)

    # Quantitative Estimate of Drug-likeness (QED)
    qed = QED.qed(molecule)

    # Hydrogen bond donors (HBD) and acceptors (HBA)
    hbd = Lipinski.NumHDonors(molecule)
    hba = Lipinski.NumHAcceptors(molecule)
    # Number of rotatable bonds
    rotatable_bonds = Lipinski.NumRotatableBonds(molecule)

    num_atoms = molecule.GetNumAtoms()
    num_heavy_atoms = Descriptors.HeavyAtomCount(molecule)
    num_rings = Descriptors.RingCount(molecule)


    
    # Number of chiral centers
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(molecule)
    c_centers = len(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))

    # Total Polar Surface Area (TPSA)
    tpsa = rdMolDescriptors.CalcTPSA(molecule)
    tpsa_ = Descriptors.TPSA(molecule)

    # Bertz CT
    bertz_ct = Descriptors.BertzCT(molecule)



    # Synthetic Accessibility Score (SA score)
    sa_score = sascorer.calculateScore(molecule)

    property_dict = {}

    property_dict['logP'] = np.round(logP, 2)
    property_dict['qed'] = np.round(qed, 2)
    property_dict['sa_score'] = np.round(sa_score, 2)
    # property_dict['penalized_logp'] = penalized_logp
    property_dict['molecular_weight'] = np.round(molecular_weight, 2)
    # property_dict['molecular_weight_exact'] = molecular_weight_exact

    property_dict['rotatable_bonds'] = rotatable_bonds
    property_dict['num_atoms'] = num_atoms
    property_dict['num_heavy_atoms'] = num_heavy_atoms
    property_dict['num_rings'] = num_rings
    property_dict['chiral_centers'] = chiral_centers
    property_dict['hbd'] = hbd
    property_dict['hba'] = hba
    # property_dict['c_centers'] = c_centers
    property_dict['tpsa'] = np.round(tpsa, 2)
    # property_dict['tpsa_'] = tpsa_
    property_dict['bertz_ct'] = np.round(bertz_ct, 2)

    if verbose:
        print(property_dict)

        # Output the calculated properties
        print(f'LogP: {logP} (<5, lower is better)')
        print(f'QED: {qed} (higher is better)')
        print(f'SA Score: {sa_score} (lower is better)')
        print(f'Molecular Weight: {molecular_weight} {molecular_weight_exact} (<500)', )
        print(f'HBD: {hbd} (<5)')
        print(f'HBA: {hba} (<10)')
        print(f'Number of Rotatable Bonds: {rotatable_bonds} (<10)')

        print('Penalized LogP:', penalized_logp)

        print(f'Number of Chiral Centers: {chiral_centers} {c_centers} (less is better)')
        print(f'TPSA: {tpsa}')
        print(f'Bertz CT: {bertz_ct}')


    return property_dict




import os
import argparse
import pandas as pd


def main(args):
    # Read docking results CSV
    docking_results_df = pd.read_csv(args.docking_results_file)

    # Iterate over each ligand file and calculate properties
    property_dicts = []
    for index, row in docking_results_df.iterrows():
        
        ligand_file = row['Ligand File']
        docking_score = row['Minimized Affinity']  # Assuming the column name for docking score is 'Minimized Affinity'
        ligand_file_path = os.path.join(args.ligand_folder, ligand_file)
        property_dict = calculate_mol_property(ligand_file_path)
        property_dict['Ligand File'] = ligand_file
        property_dict['Docking Score'] = docking_score
        property_dicts.append(property_dict)

    
    # Convert property dictionaries to DataFrame
    properties_df = pd.DataFrame(property_dicts)
    
    # Reorder columns to have Ligand File and Docking Score as leftmost columns
    properties_df = properties_df[['Ligand File', 'Docking Score'] + [col for col in properties_df.columns if col not in ['Ligand File', 'Docking Score']]]

    # Save properties to CSV
    properties_df.to_csv(args.output_file, index=False)

    print(f"Properties saved to: {args.output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate molecular properties for ligands listed in docking results CSV.")
    parser.add_argument('--docking_results_file', default='/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/docking_results.csv', help='Path to the docking results CSV file')
    parser.add_argument('--ligand_folder', default='/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb', help='Path to the folder containing ligand files')
    parser.add_argument('--output_file', default='/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb//ligand_properties.csv', help='Path to save the ligand properties CSV file')
    args = parser.parse_args()

    main(args)


    




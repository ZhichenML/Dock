import pandas as pd
import numpy as np
import os
import sys
import subprocess
def run_docking(protein_file, ligand_file):
    try:
        result = subprocess.run(
            ['./docking_clean.sh', protein_file, ligand_file],
            check=True,
            capture_output=True,
            text=True
        )

        for line in result.stdout.splitlines():
            if "Minimized Affinity:" in line:
                minimized_affinity = line.split(":")[1].strip()
                return minimized_affinity
        return "Minimized Affinity not found"
    except subprocess.CalledProcessError as e:
        print(f"Error running docking: {e}")
        print(e.stdout)
        print(e.stderr)
        return None


if __name__ == '__main__':
    receptor_folder = '/public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1-6x1a'
    protein_file = os.path.join(receptor_folder, 'GLP1-6x1a-protein.pdb')
    pocket_file = os.path.join(receptor_folder, 'GLP1-6x1a-pocket.pdb')
    ligand_ref_file = os.path.join(receptor_folder, 'GLP1-6x1a-ligand.pdb')

    csv_file='/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6x1a_poc.pdb/docking_results.csv'
    out_csv_file = '/public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6x1a_poc.pdb/docking_results_new.csv'
    df = pd.read_csv(csv_file)
    minimized_affinity = run_docking(protein_file, ligand_ref_file)
    new_row = pd.DataFrame({'Ligand File': os.path.basename(ligand_ref_file), 'Minimized Affinity': [minimized_affinity]}, index=[0])
    df = pd.concat([new_row, df])
    df.to_csv(csv_file, index=False)




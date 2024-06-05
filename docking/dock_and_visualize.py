import subprocess
import os
import argparse
import tqdm
import time
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
os.chdir(os.path.dirname(os.path.abspath(__file__)))



def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def revise_ligand_batch(fpath= '/public/home/gongzhichen/code/Lingo3DMol/data/output/Generated_Ligands/GLP1_6xox_lily_pocket.pdb'):
    # fpath : ligand folder path
    ligand_files = os.listdir(fpath)
    
    new_path = fpath+'.revise'
    if os.path.exists(new_path):
        os.system('rm -rf {}'.format(new_path))
    os.makedirs(new_path, exist_ok=True)
    for ligand_file in ligand_files:
        if ligand_file.endswith('.mol'):
            with open(os.path.join(fpath, ligand_file), 'r') as f:
                smile = f.readlines()[0]
            
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

def run_docking(protein_file, ligand_file, ligand_ref, remove_tmp_files=True):
    
    try:
        result = subprocess.run(
            ['./docking_clean.sh', protein_file, ligand_file, ligand_ref,  str(remove_tmp_files).lower()],
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

def calculate_docking_scores(protein_file, ligand_folder, ligand_ref):
    docking_results = []

    all_ligand_files = os.listdir(ligand_folder)[:3]
    for ligand_file in tqdm.tqdm(all_ligand_files):
        ligand_path = os.path.join(ligand_folder, ligand_file)
        if os.path.isfile(ligand_path) and ligand_file.endswith('.mol'):
            minimized_affinity = run_docking(protein_file, ligand_path, ligand_ref)
            if minimized_affinity != "Minimized Affinity not found":
                docking_results.append((ligand_file, float(minimized_affinity)))
            else:
                print(f"Skipping {ligand_file} due to 'Minimized Affinity not found'")
    docking_results.sort(key=lambda x: x[1])
    return docking_results

def save_docking_results(docking_results, output_file):
    df = pd.DataFrame(docking_results, columns=['Ligand File', 'Minimized Affinity'])
    df.to_csv(output_file, index=False)

def create_pymol_script(protein_file, ligand_files, ligand_files_redock, pocket_file, ligand_ref_file, output_image_base):
    script_content = f"""
# PyMOL script to visualize docking results

# Load the protein and reference ligand files
load {protein_file}, protein
load {pocket_file}, pocket
load {ligand_ref_file}, ligand_ref

# Set visualization options
hide everything
show cartoon, protein
show cartoon, pocket
show sticks, ligand_ref
color blue, protein
color red, ligand_ref
color palegreen, pocket
"""
    for i, (ligand_file, ligand_file_redock) in enumerate(zip(ligand_files, ligand_files_redock), start=0):
        if i != 0:
            script_content += f"""
    # Load ligand {i}
    load {ligand_file}, ligand{i}
    show sticks, ligand{i}
    color yellow, ligand{i}

    # Load redocked pose of ligand {i}
    load {ligand_file_redock}, ligand_redock{i}
    show sticks, ligand_redock{i}
    color cyan, ligand_redock{i}
    """
        elif i == 0:
            script_content += f"""
    # Load redocked pose of ligand {i}
    load {ligand_file_redock}, ligand_ref_redock
    show sticks, ligand_ref_redock
    color cyan, ligand_ref_redock
    """
    script_content += f"""
# Create a selection that includes pocket, ligands, and ligand_ref
select pocket_ligand_combined, pocket or ligand* or ligand_ref

color atomic, (not elem C) 
# Zoom in on the combined selection
zoom pocket_ligand_combined

# Set the viewport size for rendering
viewport 800, 600

# Save the visualization as images in different formats
ray 800, 600

# Save as PNG
png {output_image_base}.png, ray=1

# Save the current session as a PSE file
save {output_image_base}.pse
"""
    script_file = 'visualize_docking.pml'
    with open(script_file, 'w') as file:
        file.write(script_content)
    return script_file


def run_pymol_script(script_file):
    try:
        result = subprocess.run(['pymol', '-cq', script_file], check=True, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running PyMOL script: {e}")
        print(e.stdout)
        print(e.stderr)


def convert_pdbqt_to_pdb(input_file, output_file):
    try:
        result = subprocess.run(['./pdbqt2pdb.sh', input_file, output_file], check=True, capture_output=True, text=True)
        # print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")






def main():



    compound_name = 'GLP1_6xox_lily'
    receptor_folder = '/public/home/gongzhichen/code/Lingo3DMol/data/input/receptors'
    receptor_folder = os.path.join(receptor_folder, compound_name)
    protein_file = os.path.join(receptor_folder, compound_name + '_protein.pdb')
    pocket_file = os.path.join(receptor_folder, compound_name + '_pocket.pdb')
    ligand_ref_file = os.path.join(receptor_folder, compound_name + '_ligand.pdb')

    ligand_folder = '/public/home/gongzhichen/code/Lingo3DMol/data/output/Generated_Ligands'
    ligand_folder = os.path.join(ligand_folder, compound_name + '_pocket.pdb')

    revise_ligand_batch(ligand_folder)
    ligand_folder = ligand_folder + '.revise'

    output_file = os.path.join(ligand_folder, 'docking_results.csv')
    output_image_base = os.path.join(receptor_folder, 'docking_result')

    parser = argparse.ArgumentParser(description="Calculate docking scores and visualize results with PyMOL.")
    parser.add_argument('--protein_file', default=protein_file, help='Path to the protein PDB file')
    parser.add_argument('--pocket_file', default=pocket_file, help='Path to the pocket PDB file')
    parser.add_argument('--ligand_ref_file', default=ligand_ref_file, help='Path to the reference ligand PDB file')
    parser.add_argument('--ligand_folder', default=ligand_folder, help='Folder containing ligand files')
    parser.add_argument('--output_file', default=output_file, help='Path to save the docking results CSV')
    parser.add_argument('--output_image_base', default=output_image_base, help='Base name for the output image files')

    args = parser.parse_args()

    if not os.path.isfile(args.protein_file):
        print(f"Protein file does not exist: {args.protein_file}")
        exit(1)
    if not os.path.isdir(args.ligand_folder):
        print(f"Ligand folder does not exist: {args.ligand_folder}")
        exit(1)
    if not os.path.isfile(args.pocket_file):
        print(f"Pocket file does not exist: {args.pocket_file}")
        exit(1)
    if not os.path.isfile(args.ligand_ref_file):
        print(f"Ligand reference file does not exist: {args.ligand_ref_file}")
        exit(1)

    start_time = time.time()

    

    print('Starting docking')
    
    # step 1: calculate docking scores and save results
    if not os.path.isfile(args.output_file):
        docking_results = calculate_docking_scores(args.protein_file, args.ligand_folder, args.pocket_file)
        
        minimized_affinity = run_docking(args.protein_file, args.ligand_ref_file)
        docking_results.insert(0, (os.path.basename(args.ligand_ref_file), float(minimized_affinity)))
        
        save_docking_results(docking_results, args.output_file)
    else:
        print("Using existing docking results file")
        df = pd.read_csv(args.output_file)
        docking_results = [(row['Ligand File'], row['Minimized Affinity']) for _, row in df.iterrows()]

    # step 2: use docking_results for visualization
    top_10_ligands = [os.path.join(args.ligand_folder, ligand) for ligand, _ in docking_results[:11]]
    
    # show redocking poses
    # mol->pdb->pdbqt->redock.pdbqt->pdb
    
    top_10_ligands_redock = [os.path.join('/public/home/gongzhichen/code/Lingo3DMol/docking/tmp', os.path.basename(ligand).replace('.mol', '-redock.pdbqt')) for ligand, _ in docking_results[:11]]
    top_10_ligands_redock[0] = top_10_ligands_redock[0].replace('.pdb', '-redock.pdbqt')
    print('Starting re-docking')
    ligand_pdb_files = []
    for ligand_file, ligand_file_redock in zip(top_10_ligands, top_10_ligands_redock):
        if not os.path.isfile(ligand_file_redock):
            print(f"Re-docking {os.path.basename(ligand_file)}")
            minimized_affinity = run_docking(args.protein_file, ligand_file, args.pocket_file, remove_tmp_files=False)
        
        output_file = os.path.join(os.path.dirname(ligand_file_redock), os.path.basename(ligand_file_redock).replace('.pdbqt', '.pdb'))
        ligand_pdb_files.append(output_file)
        convert_pdbqt_to_pdb(ligand_file_redock, output_file)
            
    script_file = create_pymol_script(args.protein_file, top_10_ligands, ligand_pdb_files, args.pocket_file, args.ligand_ref_file, args.output_image_base)
    
    print(f"Generated PyMOL script: {script_file}")
    print('Starting PyMOL visualization...')
    run_pymol_script(script_file)

    output_image_png = f"{args.output_image_base}.png"
    output_image_pse = f"{args.output_image_base}.pse"

    if os.path.isfile(output_image_png):
        print(f"Visualization saved as {output_image_png}")
    else:
        print(f"Failed to create the visualization: {output_image_png}")

    if os.path.isfile(output_image_pse):
        print(f"Session saved as {output_image_pse}")
    else:
        print(f"Failed to save the session: {output_image_pse}")

    print(f"Total time: {time.time() - start_time:.2f} seconds")

if __name__ == '__main__':
    main()




import subprocess
import os
import argparse
import tqdm
import time
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
from rdkit.Chem import rdMolTransforms
import os.path as osp
os.chdir(os.path.dirname(os.path.abspath(__file__)))
smina_bin = '/public/home/gongzhichen/code/Lingo3DMol/docking/smina.static'

import ligand2center

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def pdb_to_smiles(pdb_file):
    # Read the PDB file into an RDKit molecule
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=True)
    if mol is None:
        raise ValueError(f"Could not read molecule from PDB file: {pdb_file}")
    
    # Generate SMILES string from the molecule
    smiles = Chem.MolToSmiles(mol)
    return smiles

def smiles_to_3d_pdb(smiles, pdb_filename):
    # Convert SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    # Write to PDB file
    with open(pdb_filename, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))

# # Example usage
# smiles = "O=C(Nc1cccc(Cl)c1)c1ccc(C(=O)NCCN2CCOC(C(=O)Nc3ccc(F)cc3)n2)c(NCc2cccc3ccccc23)n1"  # Ethanol
# pdb_filename = "ethanol.pdb"
# smiles_to_3d_pdb(smiles, pdb_filename)


def smiles_to_3d_pdb_batch(fpath, new_path):
    # mol (smiles) to pdb
    # fpath : ligand folder path
    ligand_files = os.listdir(fpath)

    for ligand_file in ligand_files:
        if ligand_file.endswith('.mol'):
            with open(os.path.join(fpath, ligand_file), 'r') as f:
                smile = f.readline().strip()
            
            pdb_filename = os.path.join(new_path, ligand_file.strip('.mol')) 

            if not '*' in smile:
                smiles_to_3d_pdb(smile, pdb_filename)
            else:
                mol = Chem.MolFromSmiles(smile, sanitize=False)
                if mol is not None:
                    rw_mol = Chem.RWMol(mol)
                
                    for atom in rw_mol.GetAtoms():
                        atom_type = atom.GetSymbol()
                        if atom_type == '*':
                            rw_mol.ReplaceAtom(atom.GetIdx(),Chem.Atom('H'))
                    
                    # Add hydrogens
                    mol = rw_mol.GetMol()
                    Chem.SanitizeMol(mol)

                    mol = Chem.AddHs(mol)

                    # Generate 3D coordinates
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)

                    # Write to PDB file
                    with open(pdb_filename, 'w') as f:
                        f.write(Chem.MolToPDBBlock(mol))
        
    return

def prepare_ligand(ligand_file, out_file=None, verbose=0):
    
    ligand_name = osp.basename(ligand_file)
    if ligand_name.endswith('.pdb'):
        ligand_pdb_file = ligand_file
    else:
        raise ValueError('Unsupported ligand file format: {}'.format(ligand_file))
    
    if out_file is None:
        out_file = ligand_pdb_file + 'qt'
    if osp.exists(out_file):
        return out_file
    
    # prepare_target4.py cannot identify the absolute or relative path, so we need to cd to the ligand_mol2_dir and perform the prepare_liagnd
    ligand_dir, ligand_pdb_name = osp.dirname(ligand_pdb_file), osp.basename(ligand_pdb_file)
    out_file_name = osp.basename(out_file)

    command = f'cd {ligand_dir} && prepare_ligand4.py -l {ligand_pdb_name} -A hydrogens -o {out_file_name}'
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running docking: {e}")
        print(e.stdout)
        print(e.stderr)
 
    if verbose:
        if result.returncode == 0:
            print('prepare the ligand successfully:', out_file)
        else:
            print('prepare the ligand failed:', result.stderr)
    pdbqt_file = osp.join(ligand_dir, out_file_name)
    command = f'mv {pdbqt_file} {out_file}'
    subprocess.run(command, shell=True, capture_output=True, text=True)

    return out_file




def prepare_ligand_batch(ligand_folder, out_folder=None, verbose=0):
    if out_folder is None:
        out_folder = ligand_folder + '_pdbqt'
        if not osp.exists(out_folder):
            os.mkdir(out_folder)
    
    all_ligand_files = os.listdir(ligand_folder)
    for ligand_file in all_ligand_files:
        if ligand_file.endswith('.pdb'):
            
            ligand_pdb_file = os.path.join(ligand_folder, ligand_file)
            prepare_ligand(ligand_file=ligand_pdb_file, out_file=os.path.join(out_folder, ligand_file+ 'qt'), verbose=verbose)



    return out_folder

def prepare_target(protein_file, out_file=None, verbose=1):
    if out_file is None:
        out_file = protein_file + 'qt'
    if osp.exists(out_file):
        return out_file
    
    command = f'prepare_receptor4.py -r {protein_file} -o {out_file}'
    if osp.exists(protein_file+'qt'):
        return protein_file+'qt'
        
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if verbose:
        if result.returncode == 0:
            print('prepare the target successfully:', out_file)
        else:
            print('prepare the target failed:', result.stderr)

    return out_file

import statistics

def calculate_center_size(center_ligand):
    x_coords, y_coords, z_coords = [], [], []

    with open(center_ligand, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)

    if not x_coords or not y_coords or not z_coords:
        raise ValueError("No HETATM records found in the provided ligand file.")

    center_x = statistics.mean(x_coords)
    center_y = statistics.mean(y_coords)
    center_z = statistics.mean(z_coords)

    size_x = max(x_coords) - min(x_coords) + 8
    size_y = max(y_coords) - min(y_coords) + 8
    size_z = max(z_coords) - min(z_coords) + 8

    return center_x, center_y, center_z, size_x, size_y, size_z

# Example usage:
# center_x, center_y, center_z, size_x, size_y, size_z = calculate_center_size('center_ligand.pdb')
# print(center_x, center_y, center_z, size_x, size_y, size_z)


def run_docking(protein_file, ligand_file, pocket_file, center_x, center_y, center_z, size_x, size_y, size_z, remove_tmp_files=True):
    
    try:
        
        result = subprocess.run(
            ['./docking_clean.sh', 
             protein_file, 
             ligand_file, 
             pocket_file, 
             str(center_x),
             str(center_y),
             str(center_z),
             str(size_x),
             str(size_y),
             str(size_z), 
             str(remove_tmp_files).lower()],
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

def docking_with_smina(protein_pdbqt, lig_pdbqt, centroid, verbose=0, out_lig_sdf=None, save_pdbqt=True):
    
    cx, cy, cz = centroid

    protein_base = os.path.basename(protein_pdbqt)
    ligand_base = os.path.basename(lig_pdbqt)
    
    # Create a temporary directory for intermediate files
    tmp_dir = f"../data/tmp/tmp_{os.path.splitext(protein_base)[0]}"
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Define output filenames based on input filenames
    out_lig_pdbqt = os.path.join(tmp_dir, f"{os.path.splitext(ligand_base)[0]}_smina.pdbqt")
    
    if out_lig_sdf is None:
        out_lig_sdf = os.path.join(tmp_dir, ligand_base.replace('.pdbqt', '_smina.sdf')   )
    
    
    if not osp.exists(protein_pdbqt) or not osp.exists(lig_pdbqt):
        print(lig_pdbqt)
        raise NotImplementedError('no pdbqt file found')
    
    command = '''{smina_bin} \
        --receptor {receptor_pre} \
        --ligand {ligand_pre} \
        --center_x {centroid_x:.4f} \
        --center_y {centroid_y:.4f} \
        --center_z {centroid_z:.4f} \
        --size_x 30 --size_y 30 --size_z 30 \
        --cpu {cpu} \
        --out {out_lig_pdbqt} \
        --exhaustiveness {exhaust}
        obabel {out_lig_pdbqt} -O {out_lig_sdf} -h'''.format(smina_bin=smina_bin,
                                            receptor_pre = protein_pdbqt,
                                            ligand_pre = lig_pdbqt,
                                            centroid_x = cx,
                                            centroid_y = cy,
                                            centroid_z = cz,
                                            cpu = 64,
                                            out_lig_pdbqt = out_lig_pdbqt,
                                            exhaust = 8,
                                            out_lig_sdf = out_lig_sdf)
    try: 
        dock_result = subprocess.run(command, shell=True, capture_output=True, text=True)
        print(f"Docking with SMINA successful: {dock_result.stdout}")
        
                
    except subprocess.CalledProcessError as e:
        print(f"Error running docking: {e}")
        
        print(e.stderr.decode())
        return None

    # Extract minimizedAffinity from the output PDBQT file
    minimized_affinity = None
    try:
        with open(out_lig_pdbqt, 'r') as file:
            for line in file:
                if line.startswith('REMARK minimizedAffinity'):
                    minimized_affinity = line.split()[2]
                    break
    except FileNotFoundError:
        print(f"Output PDBQT file not found: {out_lig_pdbqt}")

    # Output the minimized affinity
    if minimized_affinity:
        print(f"Minimized Affinity: {minimized_affinity}")
    else:
        print("Minimized Affinity not found")    

    if verbose:
        if dock_result.returncode == 0:
            print('docking successfully:', out_lig_sdf)
        else:
            print('docking failed:', dock_result.stderr)
    
    if not save_pdbqt:
        os.remove(out_lig_pdbqt)
    process_sdf(out_lig_sdf)
    return minimized_affinity

def docking_smina_batch(protein_pdbqt, ligand_folder, center_x, center_y, center_z, size_x, size_y, size_z, verbose=0, save_pdbqt=False):
    all_ligand_files = os.listdir(ligand_folder)
    print(f"Calculating docking scores for {len(all_ligand_files)} ligands")
    centroid = (center_x, center_y, center_z)
    docking_results = []

    for ligand_file in tqdm.tqdm(all_ligand_files):
        if ligand_file.endswith('.pdbqt'):
            ligand_path = os.path.join(ligand_folder, ligand_file)
           
            minimized_affinity = docking_with_smina(protein_pdbqt, ligand_path, centroid, verbose=verbose, save_pdbqt=save_pdbqt)
            if minimized_affinity != None:
                docking_results.append((ligand_file, float(minimized_affinity)))
            else:
                print(f"Skipping {ligand_file} due to 'Minimized Affinity not found'")

        else:
            print(f"Skipping {ligand_file} due to invalid file format")
    
    docking_results.sort(key=lambda x: x[1])

    return docking_results

def calculate_docking_scores(protein_file, ligand_folder, pocket_file, center_x, center_y, center_z, size_x, size_y, size_z):
    docking_results = []

    all_ligand_files = os.listdir(ligand_folder)
    print(f"Calculating docking scores for {len(all_ligand_files)} ligands")
    for ligand_file in tqdm.tqdm(all_ligand_files):
        ligand_path = os.path.join(ligand_folder, ligand_file)
        if os.path.isfile(ligand_path) and ligand_file.endswith('.pdbqt'):
            
            minimized_affinity = run_docking(protein_file, ligand_path, pocket_file, center_x, center_y, center_z, size_x, size_y, size_z)
            if minimized_affinity != "Minimized Affinity not found":
                docking_results.append((ligand_file, float(minimized_affinity)))
            else:
                print(f"Skipping {ligand_file} due to 'Minimized Affinity not found'")
        else:
            print(f"Skipping {ligand_file} due to invalid file format")

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
        # if i != 0:
            script_content += f"""
    # Load ligand {i}
    load {ligand_file}, ligand{i}_{os.path.basename(ligand_file).split('.')[0]}
    show sticks, ligand{i}
    color yellow, ligand{i}

    # Load redocked pose of ligand {i}
    load {ligand_file_redock}, ligand_redock{i}
    show sticks, ligand_redock{i}
    color cyan, ligand_redock{i}
    """
    #     elif i == 0:
    #         script_content += f"""
    # # Load redocked pose of ligand {i}
    # load {ligand_file_redock}, ligand_ref_redock
    # show sticks, ligand_ref_redock
    # color blue, ligand_ref_redock
    # """
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

from chem import read_sdf, write_sdf, sdf2centroid, sdf2mol2, pdb2centroid
def process_sdf(sdf_file):
    mols = read_sdf(sdf_file)
    for mol in mols:
        affin = mol.GetProp('REMARK').splitlines()[0].split()[1]
        smina_remark = f'SMINA RESULT:      {affin}      0.000      0.000'
        mol.SetProp('REMARK', smina_remark)
    write_sdf(mols, sdf_file)







def main():

    compound_name = 'GLP1_6xox_lily'
    base_dir = '/public/home/gongzhichen/code/Lingo3DMol/data'

    receptor_folder = os.path.join(os.path.join(base_dir, 'input/receptors'), compound_name)
    protein_file = os.path.join(receptor_folder, compound_name + '_protein.pdb')
    pocket_file = os.path.join(receptor_folder, compound_name + '_pocket.pdb')
    ligand_ref_file = os.path.join(receptor_folder, compound_name + '_ligand.pdb')

    ligand_folder = os.path.join(os.path.join(base_dir, 'output/Generated_Ligands'), compound_name + '_pocket.pdb')

    new_path = ligand_folder + '_mol2pdb'
    if not os.path.isdir(new_path):
        
        if os.path.exists(new_path):
            os.system('rm -rf {}'.format(new_path))
        os.makedirs(new_path, exist_ok=True)
        os.system(f'cp {ligand_ref_file} {os.path.join(new_path, os.path.basename(ligand_ref_file))}')

        smiles_to_3d_pdb_batch(ligand_folder, new_path)
        print("Step 1: \nGenerating pdb ligand files from SMILES successful")
    
    ligand_folder = new_path

    pdbqt_folder = ligand_folder + '_pdbqt'
    if not os.path.isdir(pdbqt_folder):
        os.makedirs(pdbqt_folder, exist_ok=True)
        prepare_ligand_batch(ligand_folder, pdbqt_folder)
    
    prepare_target(protein_file)
    protein_file = protein_file + 'qt'
    print('Prepare the target and Ligand successfully')

    output_file = os.path.join(receptor_folder, 'docking_results.csv')
    output_image_base = os.path.join(receptor_folder, 'docking_result')

    parser = argparse.ArgumentParser(description="Calculate docking scores and visualize results with PyMOL.")
    parser.add_argument('--protein_file', default=protein_file, help='Path to the protein PDB file')
    parser.add_argument('--pocket_file', default=pocket_file, help='Path to the pocket PDB file')
    parser.add_argument('--ligand_ref_file', default=ligand_ref_file, help='Path to the reference ligand PDB file')
    parser.add_argument('--ligand_folder', default=pdbqt_folder, help='Folder containing ligand files')
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
    
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_center_size(args.ligand_ref_file)
    print(f"Center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f}), Size: ({size_x:.2f}, {size_y:.2f}, {size_z:.2f})")
    size_x = 30
    size_y = 30
    size_z = 30
    
    print('Step2: \nStarting docking')
    
    # step 1: calculate docking scores and save results
    # mol->pdb->pdbqt->redock.pdbqt->pdb
    
    if not os.path.isfile(args.output_file):
        
        docking_results = docking_smina_batch(args.protein_file, args.ligand_folder, center_x, center_y, center_z, size_x, size_y, size_z, verbose=1, save_pdbqt=True)
    
        save_docking_results(docking_results, args.output_file)
    else:
        print("Using existing docking results file")
        df = pd.read_csv(args.output_file)
        docking_results = [(row['Ligand File'].replace('.mol', '.pdbqt'), row['Minimized Affinity']) for _, row in df.iterrows()]

    # step 2: use docking_results for visualization
    top_10_ligands = [os.path.join(args.ligand_folder, ligand) for ligand, _ in docking_results[:21]]
    # top_10_ligands= [os.path.join(args.ligand_folder, "701_pred_5_GLP1_6xox_lily_pocket.pdb.mol")]
    # copy the reference ligand to the ligand directory if it does not exist
    
    # homemade_folder = os.path.join(os.path.join(receptor_folder, 'homemade'))
    # homemade_files = os.listdir(homemade_folder)
    
    # prepare_ligand_batch(homemade_folder, args.ligand_folder)
    # top_10_ligands += [os.path.join(args.ligand_folder, v+'qt') for v in homemade_files]

    # show redocking poses
    
    print('Starting re-docking for visualization')
    # check redock files exist
    tmp_file = '../data/tmp/tmp_' + os.path.splitext(os.path.basename(args.protein_file))[0]
    print(f"Saving temporary files to {tmp_file}")
    
    top_10_ligands_redock = [os.path.join(tmp_file, os.path.basename(ligand).replace('.pdbqt', '_smina.pdbqt')) 
                             for ligand in top_10_ligands]
    
    ligand_pdb_files = []
    for ligand_file, ligand_file_redock in tqdm.tqdm(zip(top_10_ligands, top_10_ligands_redock)):
        
        if not os.path.isfile(ligand_file_redock):
            
            print(f"Re-docking {os.path.basename(ligand_file)} to {os.path.basename(ligand_file_redock)}")
            docking_with_smina(args.protein_file, ligand_file, (center_x, center_y, center_z))
            # run_docking(args.protein_file, ligand_file, args.pocket_file, center_x, center_y, center_z, size_x, size_y, size_z, remove_tmp_files=False)
        
        output_file = os.path.join(os.path.dirname(ligand_file_redock), os.path.basename(ligand_file_redock).replace('.pdbqt', '.pdb'))
        ligand_pdb_files.append(output_file)
        command = f'obabel -ipdbqt {ligand_file_redock} -opdb -O {output_file} -h'
        subprocess.run(command, shell=True, capture_output=True, text=True)
        # convert_pdbqt_to_pdb(ligand_file_redock, output_file)
    
    # for ligand_file in top_10_ligands:
    #     ligand2center.move2center(args.pocket_file, ligand_file, ligand_file)

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




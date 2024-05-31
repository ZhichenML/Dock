import subprocess
import csv
import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def run_docking(protein_file, ligand_file):
    try:
        # Call the shell script with the provided arguments
        result = subprocess.run(
            ['./docking_clean.sh', protein_file, ligand_file],
            check=True,
            capture_output=True,
            text=True
        )

        # Parse the output to find the minimized affinity
        for line in result.stdout.splitlines():
            if "Minimized Affinity:" in line:
                minimized_affinity = line.split(":")[1].strip()
                return minimized_affinity

        # If minimized affinity is not found
        return "Minimized Affinity not found"
    except subprocess.CalledProcessError as e:
        # Handle errors in running the shell script
        print(f"Error running docking: {e}")
        print(e.stdout)
        print(e.stderr)
        return None


def calculate_docking_scores(protein_file, ligand_folder):
    docking_results = []
    
    all_ligand_files = os.listdir(ligand_folder)
    
    # Iterate over all files in the ligand folder
    for ligand_file in all_ligand_files:
        ligand_path = os.path.join(ligand_folder, ligand_file)
        
        # Ensure it's a file and has the correct extension
        if os.path.isfile(ligand_path) and ligand_file.endswith('.mol'):
            minimized_affinity = run_docking(protein_file, ligand_path)
            if minimized_affinity != "Minimized Affinity not found":
                docking_results.append((ligand_file, minimized_affinity))
            else:
                print(f"Skipping {ligand_file} due to 'Minimized Affinity not found'")
    

    
    # Sort the results by minimized affinity (converted to float)
    docking_results.sort(key=lambda x: float(x[1]))

    return docking_results

def save_docking_results(docking_results, output_file):
    # Save the docking results to a CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Ligand File', 'Minimized Affinity'])
        writer.writerows(docking_results)


if __name__ == '__main__':
    # one on one 
    # protein_file = '/public/home/gongzhichen/code/Lingo3DMol/docking/GLP1-6x1a-protein.pdb'
    # ligand_file = '/public/home/gongzhichen/code/Lingo3DMol/docking/0_pred_3_GLP1_6x1a_poc.pdb.mol'

    # minimized_affinity = run_docking(protein_file, ligand_file)
    # print(f"Minimized Affinity: {minimized_affinity}")

    protein_file = '/public/home/gongzhichen/code/Lingo3DMol/docking/GLP1-6x1a-protein.pdb'
    ligand_folder = '/public/home/gongzhichen/code/Lingo3DMol/try_out/GLP1_6x1a_poc.pdb'
    output_file = os.path.join(ligand_folder, 'docking_results.csv')

    docking_results = calculate_docking_scores(protein_file, ligand_folder)

    # Save the results to a CSV file
    save_docking_results(docking_results, output_file)
    
    # Print the results
    for ligand_file, minimized_affinity in docking_results:
        print(f"Ligand: {ligand_file} - Minimized Affinity: {minimized_affinity}")




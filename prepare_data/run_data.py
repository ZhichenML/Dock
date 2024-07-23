import subprocess
import os

def run_confgen(input_file, output_file, num_conformers=10):
    """
    Run Schrödinger's ConfGen to generate conformers.

    Parameters:
    input_file (str): Path to the input structure file.
    output_file (str): Path to the output structure file.
    num_conformers (int): Number of conformers to generate.
    """
    # Construct the ConfGen command
    confgen_cmd = [
        "confgen",
        "-imae", input_file,
        "-omae", output_file,
        "-WAIT",
        "-num_confs", str(num_conformers)
    ]

    # Run the ConfGen command
    result = subprocess.run(confgen_cmd, capture_output=True, text=True)

    # Check for errors
    if result.returncode != 0:
        print("Error running ConfGen:")
        print(result.stderr)
    else:
        print("ConfGen completed successfully.")
        print(result.stdout)


def main():
    input_file = "input.mae"  # Replace with your input file path
    output_file = "output_conformers.mae"  # Replace with your output file path
    num_conformers = 10  # Number of conformers to generate

    run_confgen(input_file, output_file, num_conformers)

def main_sch():
    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdmolfiles import SDMolSupplier
    from schrodinger import structure
    from schrodinger.job import jobcontrol

    # # Define the input SMILES
    input_smiles = 'COc1c(C(=O)[O-])cccc1C(=O)N1CCN(c2nn(CCC(C)=O)c(=O)cc2C)CC1'

    # # Create an RDKit molecule from the SMILES string
    rdkit_molecule = Chem.MolFromSmiles(input_smiles)
    rdkit_molecule = Chem.AddHs(rdkit_molecule)  # Add hydrogens
    AllChem.EmbedMolecule(rdkit_molecule)  # Generate 3D coordinates

    # # Save the RDKit molecule to an SDF file
    input_sdf = 'input.sdf'
    if not os.path.exists(input_sdf):    
        with Chem.SDWriter(input_sdf) as writer:
            writer.write(rdkit_molecule)
    data=SDMolSupplier(fileName=input_sdf, sanitize=True, removeHs=False, strictParsing=True)
    
    for i,mol in enumerate([data]):
        if i == 1: break
        print(f'Processing molecule {i+1}')
        mol_name = mol.GetPropsAsDict()['ID']
        # Convert the SDF file to a Maestro file using Schrödinger's structure module
        #input_mae = 'input_structure.mae'
        output_maegz = mol_name + '-out.maegz'
        output_mae = mol_name + '-out.mae'
        output_sdf = mol_name + '-out.sdf'

        # with structure.StructureWriter(input_mae) as writer:
        #     for struct in structure.StructureReader(input_sdf):
        #         writer.append(struct)
        
        # Set up and run the ConfGen job
        confgen_command = [os.path.join(os.environ['SCHRODINGER'], 'confgen'), '-n', '10', '-t', '600', '-j', mol_name, input_sdf]
        
        job = jobcontrol.launch_job(confgen_command)

        # Wait for the job to complete
        job.wait()


        # Load the generated conformers
        conformers = structure.StructureReader(output_maegz)

        # Save the conformers to a file
        with structure.StructureWriter('conformers.mae') as writer:
            for conformer in conformers:
                writer.append(conformer)

        print(f'{len(list(conformers))} conformers generated and saved to conformers.mae')



        # Define the structconvert command
        structconvert_command = [os.path.join(os.environ['SCHRODINGER'], 'utilities/structconvert'),   'conformers.mae',  output_sdf]

        # Run the structconvert command
        result = subprocess.run(structconvert_command, capture_output=True, text=True)

        # Check for errors
        if result.returncode != 0:
            print("Error during conversion:", result.stderr)
        else:
            print(f"Conversion successful. {output_sdf} generated.")

        # Optionally, print the stdout and stderr for debugging
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)




    





# Example usage
if __name__ == "__main__":
    main_sch()

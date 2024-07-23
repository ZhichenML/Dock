import subprocess
import os

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdmolfiles import SDMolSupplier


def get_confgen_sdf(input_sdf):
    from schrodinger import structure
    from schrodinger.job import jobcontrol

    data=SDMolSupplier(fileName=input_sdf, sanitize=True, removeHs=False, strictParsing=True)
    
    mol_name = data[0].GetPropsAsDict()['ID']
   
    
    output_maegz = mol_name + '-out.maegz'
    output_mae = mol_name + '-out.mae'
    output_sdf = mol_name + '-out.sdf'

    # zhichen: not use mae file, use sdf file directly
    # Convert the SDF file to a Maestro file using Schrödinger's structure module
    # input_mae = 'input_structure.mae'
    # with structure.StructureWriter(input_mae) as writer:
    #     for struct in structure.StructureReader(input_sdf):
    #         writer.append(struct)
  
    confgen_command = [os.path.join(os.environ['SCHRODINGER'], 'confgen'), '-n', '10', '-t', '600', '-j', mol_name, input_sdf]
    job = jobcontrol.launch_job(confgen_command)
    job.wait()

    conformers = structure.StructureReader(output_maegz)

    # Save the conformers to a file
    with structure.StructureWriter(output_mae) as writer:
        for conformer in conformers:
            writer.append(conformer)


    structconvert_command = [os.path.join(os.environ['SCHRODINGER'], 'utilities/structconvert'),   output_mae,  output_sdf]

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

def smiles_to_sdf(input_smiles='COc1c(C(=O)[O-])cccc1C(=O)N1CCN(c2nn(CCC(C)=O)c(=O)cc2C)CC1', input_sdf = 'input.sdf'):
    # Smiles with no index -> rdkit mol -> order -> sdf

    rdkit_molecule = Chem.MolFromSmiles(input_smiles)
    
    for atom in rdkit_molecule.GetAtoms():atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))
    rdkit_molecule = Chem.AddHs(rdkit_molecule)  # Add hydrogens
    rdkit_molecule.SetProp('ID', 'test')  
    Draw.MolToFile(rdkit_molecule, 'input_mol_with_indices.png', size=(500, 500), kekulize=True, wedgeBonds=True)
    AllChem.EmbedMolecule(rdkit_molecule)  # Generate 3D coordinates
    
    with Chem.SDWriter(input_sdf) as writer:
        writer.write(rdkit_molecule)

    

    




def sdf_to_smile_with_order(input_sdf='test-out.sdf'):
    # Create a supplier object
    supplier = Chem.SDMolSupplier(input_sdf)

    # Initialize variables to store the lowest energy and corresponding molecule
    lowest_energy = float('inf')
    lowest_energy_mol = None

    # Iterate over molecules (conformations)
    for mol in supplier:
        if mol is None:
            continue  # Skip molecules that couldn't be read

        # Get the energy property (assume it's stored in a property called 'Energy')
        energy = float(mol.GetProp('r_f3d_energy'))

        # Check if this conformation has the lowest energy
        if energy < lowest_energy:
            lowest_energy = energy
            lowest_energy_mol = mol

    # Check if a conformation was found
    if lowest_energy_mol is None:
        print("No valid conformations found.")
    else:
        print(f"Lowest energy conformation has energy: {lowest_energy}")
        
        # Process the lowest energy conformation
        # For example, convert it to SMILES
        for atom in lowest_energy_mol.GetAtoms():atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))
        mol_with_indices = Chem.MolFromSmiles(Chem.MolToSmiles(lowest_energy_mol, canonical=False))
        Draw.MolToFile(mol_with_indices, 'lowest2d.png', size=(500, 500), kekulize=True, wedgeBonds=True)
        
        smiles = Chem.MolToSmiles(lowest_energy_mol, canonical=False)
        print(f"SMILES of the lowest energy conformation: {smiles}")
        
        # You can also save this conformation to a new SDF file
        # output_sdf_file = 'path/to/your/lowest_energy_conformation.sdf'
        # writer = Chem.SDWriter(output_sdf_file)
        # writer.write(lowest_energy_mol)
        # writer.close()
        # print(f"Lowest energy conformation saved to {output_sdf_file}")





# Example usage
if __name__ == "__main__":
    input_smiles = 'COc1c(C(=O)[O-])cccc1C(=O)N1CCN(c2nn(CCC(C)=O)c(=O)cc2C)CC1'
    input_sdf = 'input.sdf'
    
    # smiles_to_sdf(input_smiles, input_sdf)
    get_confgen_sdf(input_sdf)
    
    # sdf_to_smile_with_order()

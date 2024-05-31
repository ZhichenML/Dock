import subprocess
import os
import argparse

def create_pymol_script(protein_file, ligand_file, pocket_file, ligand_ref_file, output_image_base):
    script_content = f"""
# PyMOL script to visualize docking results

# Load the protein and ligand files
load {protein_file}, protein
load {pocket_file}, pocket
load {ligand_ref_file}, ligand_ref
load {ligand_file}, ligand

# Set visualization options
hide everything
show cartoon, protein
show cartoon, pocket
show sticks, ligand_ref
show sticks, ligand

# Color the protein and ligand
color blue, protein
color yellow, ligand
color red, ligand_ref
color palegreen, pocket

# Create a selection that includes pocket, ligand, and ligand_ref
select pocket_ligand_combined, pocket or ligand or ligand_ref

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

def main():
    protein_file =    '/public/home/gongzhichen/code/Lingo3DMol/docking/GLP1-6x1a-protein.pdb'
    pocket_file =     '/public/home/gongzhichen/code/Lingo3DMol/docking/GLP1-6x1a-poc.pdb'
    ligand_ref_file = '/public/home/gongzhichen/code/Lingo3DMol/docking/GLP1-6x1a-ligand.pdb'
    ligand_file =     '/public/home/gongzhichen/code/Lingo3DMol/try_out/GLP1_6x1a_poc.pdb/96_pred_4_GLP1_6x1a_poc.pdb.mol'
    output_image_base='docking_result'

    parser = argparse.ArgumentParser(description="Visualize docking results with PyMOL.")
    parser.add_argument('--protein_file', default=protein_file, help='Path to the protein PDB file')
    parser.add_argument('--ligand_file', default=ligand_file, help='Path to the ligand file (PDB/MOL)')
    parser.add_argument('--pocket_file', default=pocket_file, help='Path to the pocket PDB file')
    parser.add_argument('--ligand_ref_file', default=ligand_ref_file, help='Path to the reference ligand PDB file')
    parser.add_argument('--output_image_base', default=output_image_base, help='Base name for the output image files')

    args = parser.parse_args()
    
    # Check if the input files exist
    for file in [args.protein_file, args.ligand_file, args.pocket_file, args.ligand_ref_file]:
        if not os.path.isfile(file):
            print(f"File does not exist: {file}")
            exit(1)

    # Generate the PyMOL script
    script_file = create_pymol_script(args.protein_file, args.ligand_file, args.pocket_file, args.ligand_ref_file, args.output_image_base)
    print(f"Generated PyMOL script: {script_file}")

    # Run the PyMOL script to generate the images
    run_pymol_script(script_file)

    # Check if the output images were created successfully
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

if __name__ == '__main__':
    main()




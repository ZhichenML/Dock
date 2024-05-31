
# PyMOL script to visualize docking results

# Load the protein and reference ligand files
load /public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1_6xox_lily/GLP1_6xox_lily_protein.pdb, protein
load /public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1_6xox_lily/GLP1_6xox_lily_pocket.pdb, pocket
load /public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1_6xox_lily/GLP1_6xox_lily_ligand.pdb, ligand_ref

# Set visualization options
hide everything
show cartoon, protein
show cartoon, pocket
show sticks, ligand_ref
color blue, protein
color red, ligand_ref
color palegreen, pocket

# Load ligand 1
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/GLP1_6xox_lily_ligand.pdb, ligand1
show sticks, ligand1
color yellow, ligand1

# Load ligand 2
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/538_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand2
show sticks, ligand2
color yellow, ligand2

# Load ligand 3
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/272_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand3
show sticks, ligand3
color yellow, ligand3

# Load ligand 4
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/271_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand4
show sticks, ligand4
color yellow, ligand4

# Load ligand 5
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/536_pred_3_GLP1_6xox_lily_pocket.pdb.mol, ligand5
show sticks, ligand5
color yellow, ligand5

# Load ligand 6
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/281_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand6
show sticks, ligand6
color yellow, ligand6

# Load ligand 7
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/223_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand7
show sticks, ligand7
color yellow, ligand7

# Load ligand 8
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/161_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand8
show sticks, ligand8
color yellow, ligand8

# Load ligand 9
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/139_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand9
show sticks, ligand9
color yellow, ligand9

# Load ligand 10
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/153_pred_4_GLP1_6xox_lily_pocket.pdb.mol, ligand10
show sticks, ligand10
color yellow, ligand10

# Load ligand 11
load /public/home/gongzhichen/code/Lingo3DMol/datasets/Generated_Ligands/GLP1_6xox_lily_pocket.pdb/329_pred_3_GLP1_6xox_lily_pocket.pdb.mol, ligand11
show sticks, ligand11
color yellow, ligand11

# Create a selection that includes pocket, ligands, and ligand_ref
select pocket_ligand_combined, pocket or ligand* or ligand_ref

# Zoom in on the combined selection
zoom pocket_ligand_combined

# Set the viewport size for rendering
viewport 800, 600

# Save the visualization as images in different formats
ray 800, 600

# Save as PNG
png /public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1_6xox_lily/docking_result.png, ray=1

# Save the current session as a PSE file
save /public/home/gongzhichen/code/Lingo3DMol/datasets/receptors/GLP1_6xox_lily/docking_result.pse

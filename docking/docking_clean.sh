#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 9 ] || [ "$#" -gt 10 ]; then
    echo "Usage: $0 protein_file ligand_file center_ligand center_x center_y center_z size_x size_y size_z [remove_dir]"
    exit 1
fi

# Input files
protein_file=$1
ligand_file=$2
center_ligand=$3
center_x=$4
center_y=$5
center_z=$6
size_x=$7
size_y=$8
size_z=$9
remove_dir=${10:-true}  # Default value is true if not provided

# Get the base names of the input files
protein_base=$(basename "$protein_file")
ligand_base=$(basename "$ligand_file")
echo ligand_base

# Create a temporary directory for intermediate files
tmp_dir="./tmp_${protein_base%.*}" 
mkdir -p "$tmp_dir"

# Step 1: Process the protein file
# Remove water molecules (residue name HOH) and add hydrogens to the protein structure
# obabel "$protein_file" -d -o pdb -O "$tmp_dir/${protein_base%.*}_no_waters.pdb" --filter "resn!=HOH"
# obabel "$tmp_dir/${protein_base%.*}_no_waters.pdb" -h -xr -O "$tmp_dir/${protein_base%.*}.pdbqt"
obabel -ipdb "$protein_file" -opdbqt -xr -O "$tmp_dir/${protein_base%.*}.pdbqt" -p 7.4 

# Step 2: Process the ligand file
# Convert the ligand file to PDBQT format using OpenBabel
convert_ligand_to_pdbqt() {
    input_file=$1
    output_file=$2

    input_extension="${input_file##*.}"

    case $input_extension in
        mol)
            obabel "$input_file" -O "$output_file" --gen3D best -p 7.4 --minimize --ff MMFF94 
            ;;
        pdb)
            obabel "$input_file" -O "$output_file" --gen3D best -p 7.4 --minimize --ff MMFF94 
            ;;
        *)
            echo "Unsupported file format: $input_extension"
            return 1
            ;;
    esac

    if [ $? -eq 0 ]; then
        echo "Conversion successful: $output_file"
        return 0
    else
        echo "Conversion failed"
        return 1
    fi
}

convert_ligand_to_pdbqt "$ligand_file" "$tmp_dir/${ligand_base%.*}.pdbqt"
if [ $? -ne 0 ]; then
    echo "Conversion of ligand file to .pdbqt failed"
    [[ "$remove_dir" == "true" ]] && rm -r "$tmp_dir"
    exit 1
fi

# Define output filenames based on input filenames
protein_pdbqt="$tmp_dir/${protein_base%.*}.pdbqt"
ligand_pdbqt="$tmp_dir/${ligand_base%.*}.pdbqt"
output_pdbqt="$tmp_dir/${ligand_base%.*}-redock.pdbqt"

# Run smina with the specified parameters
./smina.static -r "$protein_pdbqt" -l "$ligand_pdbqt" --center_x "$center_x" --center_y "$center_y" --center_z "$center_z" --size_x "$size_x" --size_y "$size_y" --size_z "$size_z" --cpu 64 --exhaustiveness 16 -o "$output_pdbqt"

# Confirm the docking run
if [ $? -eq 0 ]; then
    echo "Docking run successful: $output_pdbqt"
else
    echo "Docking run failed"
    [[ "$remove_dir" == "true" ]] && rm -r "$tmp_dir"
    exit 1
fi

# Extract minimizedAffinity from the output PDBQT file
minimized_affinity=$(grep -m 1 'REMARK minimizedAffinity' "$output_pdbqt" | awk '{print $3}')

# Output the minimized affinity
if [ -n "$minimized_affinity" ]; then
    echo "Minimized Affinity: $minimized_affinity"
else
    echo "Minimized Affinity not found"
fi

# Clean up the temporary directory based on the remove_dir parameter
if [[ "$remove_dir" == "true" ]]; then
    rm -r "$tmp_dir"
fi

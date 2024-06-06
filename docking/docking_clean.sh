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
protein_base=$(basename $protein_file)
ligand_base=$(basename $ligand_file)

# Create a temporary directory for intermediate files
tmp_dir="./tmp_${protein_base%.*}" 
mkdir -p $tmp_dir

# Check and transform ligand file from .mol to .pdb if necessary
ligand_extension="${ligand_file##*.}"
if [ "$ligand_extension" = "mol" ]; then
  ligand_pdb="${ligand_base%.*}.pdb"
  obabel $ligand_file -O $tmp_dir/$ligand_pdb
  if [ $? -ne 0 ]; then
    echo "Transformation of ligand file from .mol to .pdb failed"
    [[ "$remove_dir" == "true" ]] && rm -r $tmp_dir
    exit 1
  fi
  ligand_file=$tmp_dir/$ligand_pdb
else
  cp $ligand_file $tmp_dir/
  ligand_file=$tmp_dir/$ligand_base
fi

# Convert protein and ligand files to PDBQT format
obabel $protein_file -xr -O $tmp_dir/${protein_base%.*}.pdbqt
obabel $ligand_file -O $tmp_dir/${ligand_base%.*}.pdbqt -h

# Define output filenames based on input filenames
protein_pdbqt=$tmp_dir/${protein_base%.*}.pdbqt
ligand_pdbqt=$tmp_dir/${ligand_base%.*}.pdbqt
output_pdbqt=$tmp_dir/${ligand_base%.*}-redock.pdbqt

# Run smina with the specified parameters
./smina.static -r $protein_pdbqt -l $ligand_pdbqt --center_x $center_x --center_y $center_y --center_z $center_z --size_x $size_x --size_y $size_y --size_z $size_z --cpu 64 --exhaustiveness 16 -o $output_pdbqt

# Confirm the docking run
if [ $? -eq 0 ]; then
  echo "Docking run successful: $output_pdbqt"
else
  echo "Docking run failed"
  [[ "$remove_dir" == "true" ]] && rm -r $tmp_dir
  exit 1
fi

# Extract minimizedAffinity from the output PDBQT file
minimized_affinity=$(grep -m 1 'REMARK minimizedAffinity' $output_pdbqt | awk '{print $3}')

# Output the minimized affinity
if [ -n "$minimized_affinity" ]; then
  echo "Minimized Affinity: $minimized_affinity"
else
  echo "Minimized Affinity not found"
fi

# Clean up the temporary directory based on the remove_dir parameter
if [[ "$remove_dir" == "true" ]]; then
  rm -r $tmp_dir
fi

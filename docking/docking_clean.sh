#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 protein_file ligand_file"
  exit 1
fi

# Input files
protein_file=$1
ligand_file=$2

# Create a temporary directory for intermediate files
tmp_dir="./tmp"
mkdir -p $tmp_dir

# Get the base names of the input files
protein_base=$(basename $protein_file)
ligand_base=$(basename $ligand_file)

# Check and transform ligand file from .mol to .pdb if necessary
ligand_extension="${ligand_file##*.}"
if [ "$ligand_extension" = "mol" ]; then
  ligand_pdb="${ligand_base%.*}.pdb"
  obabel $ligand_file -O $tmp_dir/$ligand_pdb
  if [ $? -ne 0 ]; then
    echo "Transformation of ligand file from .mol to .pdb failed"
    rm -r $tmp_dir
    exit 1
  fi
  ligand_file=$tmp_dir/$ligand_pdb
else
  cp $ligand_file $tmp_dir/
  ligand_file=$tmp_dir/$ligand_base
fi

# Convert protein and ligand files to PDBQT format
obabel $protein_file -xr -O $tmp_dir/${protein_base%.*}.pdbqt
obabel $ligand_file -O $tmp_dir/${ligand_base%.*}.pdbqt

# Define output filenames based on input filenames
protein_pdbqt=$tmp_dir/${protein_base%.*}.pdbqt
ligand_pdbqt=$tmp_dir/${ligand_base%.*}.pdbqt
output_pdbqt=$tmp_dir/${ligand_base%.*}-redock.pdbqt

# Run smina with the specified parameters
./smina.static -r $protein_pdbqt -l $ligand_pdbqt --autobox_ligand $ligand_pdbqt --autobox_add 8 --cpu 64 --exhaustiveness 16 -o $output_pdbqt

# Confirm the docking run
if [ $? -eq 0 ]; then
  echo "Docking run successful: $output_pdbqt"
else
  echo "Docking run failed"
  rm -r $tmp_dir
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

# Clean up the temporary directory
rm -r $tmp_dir

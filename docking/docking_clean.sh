#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
  echo "Usage: $0 protein_file ligand_file center_ligand [remove_dir]"
  exit 1
fi

# Input files
protein_file=$1
ligand_file=$2
center_ligand=$3
remove_dir=${4:-true}  # Default value is true if not provided

# Get the base names of the input files
protein_base=$(basename $protein_file)
ligand_base=$(basename $ligand_file)

# Create a temporary directory for intermediate files
tmp_dir="./tmp" 
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
obabel $ligand_file -O $tmp_dir/${ligand_base%.*}.pdbqt

# Define output filenames based on input filenames
protein_pdbqt=$tmp_dir/${protein_base%.*}.pdbqt
ligand_pdbqt=$tmp_dir/${ligand_base%.*}.pdbqt
output_pdbqt=$tmp_dir/${ligand_base%.*}-redock.pdbqt

# Extract bounding box coordinates from the pocket file
min_x=$(grep "^ATOM" $center_ligand | awk 'NR==1{min=$7} {if($7<min) min=$7} END {print min}')
max_x=$(grep "^ATOM" $center_ligand | awk 'NR==1{max=$7} {if($7>max) max=$7} END {print max}')
min_y=$(grep "^ATOM" $center_ligand | awk 'NR==1{min=$8} {if($8<min) min=$8} END {print min}')
max_y=$(grep "^ATOM" $center_ligand | awk 'NR==1{max=$8} {if($8>max) max=$8} END {print max}')
min_z=$(grep "^ATOM" $center_ligand | awk 'NR==1{min=$9} {if($9<min) min=$9} END {print min}')
max_z=$(grep "^ATOM" $center_ligand | awk 'NR==1{max=$9} {if($9>max) max=$9} END {print max}')

center_x=$(echo "($min_x + $max_x) / 2.0" | bc -l)
center_y=$(echo "($min_y + $max_y) / 2.0" | bc -l)
center_z=$(echo "($min_z + $max_z) / 2.0" | bc -l)
size_x=$(echo "$max_x - $min_x + 8" | bc -l)
size_y=$(echo "$max_y - $min_y + 8" | bc -l)
size_z=$(echo "$max_z - $min_z + 8" | bc -l)

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

This docking score calculation module is developed based on SMINA, which is a fast fork of Autodock vina.
To use


conda install -c conda-forge openbabel
conda install -c conda-forge py3dmol

wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina.static
or conda install -c conda-forge smina

chmod u+x smina.static


for visualization

apt-get install pymol 


docking_clean.sh do the one on one docking and clean temporary files


dock_and_visualize.py do the docking and visualization of the results

QED_and_SAS.py calculate QED and SAS descriptors for the ligands (generated)
#!/bin/bash

## an example of the FLEXR workflow with a grid optimization of parameters
## author: Tim Stachowski, PhD

# define paths to coot1 and python3
cootdir=/opt/homebrew/Cellar/coot/1.0.05/bin/coot
python3dir=/opt/homebrew/bin/python3

# make list of PDBs
ls ????.pdb > pdb.list

# loop over models in input list
for i in `cat pdb.list`
do

i=${i%.pdb*}

# remove previously built alt confs
phenix.pdbtools ${i}".pdb" remove_alt_confs=True

# create maps
phenix.maps ${i}"_modified.pdb" ${i}".mtz"

# run ringer
mmtbx.ringer ${i}"_modified.pdb" ${i}"_modified_map_coeffs.mtz"

# test different sigma cutoffs
for j in 0.30 0.35 0.40
do

# test different geometry tolerances
for k in 30 35 40
do

# run flexr-find
${python3dir} flexr_find.py -f ${i}"_modified_ringer.csv" -t ${j} -g ${k}

# run flexr-build
${cootdir} --script flexr_build.py file.${i}"_modified"

# rename output files with grid parameters
# output of peak finding
mv "peak_finder_output_"${i}"_modified_ringer.csv" "peak_finder_output_"${i}"_"${j}"_"${k}"_modified_ringer.csv"
# identified alt confs
mv ${i}"_modified_ringer_alts.csv" ${i}"_"${j}"_"${k}"_modified_ringer_alts.csv"
# FLEXR built model
mv ${i}"_modified_newconfs.pdb" ${i}"_"${j}"_"${k}"_modified_newconfs.pdb"

done
done
done

# ready for refinement....
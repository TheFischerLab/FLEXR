#!/bin/bash

## an example of the FLEXR workflow with a grid optimization of parameters
## author: Tim Stachowski, PhD

# make list of PDBs
ls ????.pdb > pdb.list


for i in `cat pdb.list`
do

i=${i%.*}

# remove previously built alt confs
phenix.pdbtools ${i}".pdb" remove_alt_confs=True
mv ${i}"_modified.pdb" ${i}".pdb"

# create maps
phenix.maps ${i}".pdb" ${i}".mtz"

# run ringer
mmtbx.ringer ${i}".pdb" ${i}"_map_coeffs.mtz"

# test different sigma cutoffs
for j in 0.30 0.35 0.40
do

# test different geometry tolerances
for k in 30 35 40
do

# run flexr-find
/opt/homebrew/bin/python3 flexr_find.py -f ${i}"_ringer.csv" -t ${j} -g ${k}
mv "peak_finder_output_"${i}"_ringer.csv" "peak_finder_output_"${i}"_"${j}"_"${k}"_ringer.csv"
cp $i".pdb" ${i}"_"${j}"_"${k}".pdb"
cp $i".mtz" ${i}"_"${j}"_"${k}".mtz"
mv ${i}"_ringer_alts.csv" ${i}"_"${j}"_"${k}"_ringer_alts.csv"

# run flexr-build
/opt/homebrew/Cellar/coot/1.0.05/bin/coot --script flexr_build.py file-${i}"_"${j}"_"${k} #load_maps-True score-True
rm ${i}"_"${j}"_"${k}".pdb"


done
done
done

# ready for refinement....
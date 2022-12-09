#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#    RINGER-REFINE:
#    This script detects conformational changes from Ringer measurements and
#    automatically builds them into a multiconformer model for refinement.
#
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import time
from itertools import product
import pandas as pd
import numpy as np
from ringer_tools import create_parser, ringer_parser, peak_find, match_and_build


CLI = create_parser()
ARGS = CLI.parse_args()

A = ARGS.filename
B = ARGS.geotolerance
C = ARGS.sigmathreshold
D = ARGS.plot
E = ARGS.peakheight
F = ARGS.peakprominence
G = ARGS.peakwidth
H = ARGS.peakdistance
S = ARGS.step

## Intro messages
print(' ')
print('Welcome to FLEXR!')
print('A program for automated multi-conformer model building using unbiased electron density map sampling.')
time.sleep(1)
print(' ')
print('Brought to you by the Fischer Lab at St. Jude Children\'s Research Hospital')
print('Copyright reserved')
print('Please cite: ')
print(' ')
time.sleep(1)
#print('This program finds peaks in Ringer plots and integrates them.')
#print('Flipper detects gain/losses in side chain conformers between ')
#print('two datasets and cases where the predominate peak changes ')
#print('are FLIPS!')
time.sleep(1)
print(' ')
print(' ')
print('Type -h for info about options')
print(' ')
time.sleep(5)
print('Let\'s get started')
print('Checking dependencies...')
time.sleep(5)

#disable printing pandas warnings
pd.options.mode.chained_assignment = None

## check if rotamer library can be found
try:
    ## define location of rotamer library
    print('Loading ideal rotamer library...')
    library = pd.read_csv('rotamer_library_coot.csv',header=0)
    chi_labels = ['chi1_mean','chi2_mean','chi3_mean','chi4_mean']
    for i in chi_labels:
        library[i] = library[i].apply(lambda x : x+360 if x<0 else x)
except NameError:
    print('rotamer_library_coot.csv needs to be in the same directory')

## assign what residues will always have peaks at certain chi angles
chi1_allowed = ['SER','GLN','ASN','GLU','ASP','ARG','LYS','MET','CYS']
chi1_branched = ['THR','VAL','ILE']
chi2_branched = ['LEU']
chi2_allowed = ['ILE']
rings = ['PHE','TYR','TRP','HIS']

print(' ')
print('Loading the first dataset...')
print('Finding peaks....')
if D:
    print('Plotting...')

## load data, find peaks, and build single dataframe
chis = ['chi1','chi2','chi3','chi4']
df = pd.DataFrame([])
for i in chis:
    parsed,rot_angles = ringer_parser(A,S)
    tmp = peak_find(A,parsed,i,C,D,E,F,G,H,rot_angles)
    df = pd.concat([df,tmp])
df = df.sort_values(by=['chain','res_n','chi'])
df.to_csv('peak_finder_output_'+A,header=True,index=False)

print(' ')
print('Matching Ringer peaks to ideal rotamers...')
## loop through all residues an extract info
## compare chi angles from Ringer to ideal angles in library
alt_confs = []
for a,b,c in df[['res_n','res_type','chain']].drop_duplicates().values:
    angles = df[(df['res_n']==a)&(df['res_type']==b)&(df['chain']==c)]['peak_angles']
    res = b
    library_tmp = library[library['res_type']==res]
    angles = [x[:] for x in np.array(angles)]
    angles = [[x for x in sublist] for sublist in angles]
    if len(angles) > 1:
        alts = list(product(*angles))
    if len(angles) == 1:
        alts = list(angles)
    if len(alts) > 0:
        ## test which residues have certain number of peaks based on branching vs unbrached vs ring
        if res in chi1_allowed:
            # test if any of chis contain more than one peak
            more_than_one = [x>1 for x in [len(sublist) for sublist in angles]]
            if True in more_than_one:
                # for linear side chains with multiple chi angles
                if len(alts) > 1:
                    for k in alts:
                        if len(k) > 1:
                            res2build = match_and_build(k,res,library_tmp,B)
                            alt_confs.append((a,b,c,res2build))
                # for serine that only has chi1
                elif len(alts) == 1:
                    for k in alts[0]:
                        k = [k]
                        res2build = match_and_build(k,res,library_tmp,B)
                        alt_confs.append((a,b,c,res2build))
        if res in chi1_branched:
            if len(alts[0]) > 2:
                for k in alts[0]:
                    k = [k]
                    res2build = match_and_build(k,res,library_tmp,B)
                    alt_confs.append((a,b,c,res2build))
        if res in chi2_branched:
            if (len(angles[0]) > 1) or (len(angles[1]) > 2):
                for k in alts:
                    if len(k) > 1:
                        res2build = match_and_build(k,res,library_tmp,B)
                        alt_confs.append((a,b,c,res2build))
        if res in chi2_allowed:
            if (len(angles[0]) > 2) or (len(angles[1]) > 1):
                for k in alts:
                    if len(k) > 1:
                        res2build = match_and_build(k,res,library_tmp,B)
                        alt_confs.append((a,b,c,res2build))
        if res in rings:
            if len(angles[0]) > 1:
                for k in alts:
                    if len(k) > 1:
                        res2build = match_and_build(k,res,library_tmp,B)
                        alt_confs.append((a,b,c,res2build))

## organize output
try:
    print(' ')
    alt_confs = pd.DataFrame(alt_confs)
    alt_confs.columns = ['res_n','res_type','chain','rotamer']
    alt_confs = alt_confs.replace('',np.nan)
    alt_confs = alt_confs.dropna()
    ## continue with residues that have more than one conf
    more_than_one = alt_confs[['res_n','chain']].value_counts().loc[lambda x: x>1].index
    alt_confs = alt_confs[alt_confs[['res_n','chain']].apply(tuple,1).isin(list(more_than_one))]
    alt_confs = alt_confs.sort_values(by=['chain','res_n'])
    alt_confs.to_csv(A[:-4]+'_alts.csv',header=True,index=False)
    print('Matching rotamers saved to: ringer_alts.csv...')
    print('Ready to build with Coot!!')
    print(' ')
except IndexError:
    print('Sorry, no matching rotamers found')
    print(' ')

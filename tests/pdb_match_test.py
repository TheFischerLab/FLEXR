#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 14:15:36 2022

@author: Tim Stachowski, PhD
Fischer Laboratory
St. Jude Children's Research Hospital

Test if two PDBs are the same
"""

import numpy as np
import pandas as pd
import argparse
from glob import glob
import subprocess

CLI = argparse.ArgumentParser(description='A simple test to compare two PDBs')

CLI.add_argument(
    '-f1',
    '--file1',
    nargs="?",
    type=str,
    default=None,
    help='input file1'
)

CLI.add_argument(
    '-f2',
    '--file2',
    nargs="?",
    type=str,
    default=None,
    help='input file2'
)

CLI.add_argument(
    '-d',
    '--details',
    nargs="?",
    type=bool,
    default=False,
    help='show specifics on mismatches. Default = False'
)

CLI.add_argument(
    '-ia',
    '--ignore_atom_num',
    nargs="?",
    type=bool,
    default=True,
    help='Ignore atom numbering. Default = True'
)

ARGS = CLI.parse_args()

A = ARGS.file1
B = ARGS.file2
C = ARGS.details
D = ARGS.ignore_atom_num


def get_atom_info(file):
    if file.endswith(".pdb"):
        ## radius
        pdb_info = []
        with open(file) as data:
            read_data = data.read()
            lines = read_data.splitlines()
            for line in lines:
                if line.startswith('ATOM'):
                    atom_type = line[12:16].strip()
                    #remove hydrogens
                    if 'H' not in atom_type:
                        atom_het = line[0:7].strip()
                        atom_num = line[7:13].strip()
                        residue_name = line[17:20].strip()
                        chain = line[21].strip()
                        res_seq = int(line[22:26].strip())
                        alt_loc = line[16:17]
                        # occupancy
                        occ = float(line[56:60].strip())
                        b_fac = float(line[61:65].strip())
                        x,y,z = line[30:38].strip(),line[38:46].strip(),line[46:54].strip()
                        x,y,z = float(x),float(y),float(z)
                        pdb_info.append((file,atom_het,atom_num,atom_type,alt_loc,occ,residue_name,chain,res_seq,x,y,z,b_fac))
    pdb_info = pd.DataFrame(pdb_info,columns = ['model','atom_het','atom_num','atom_type','alt_loc','occupancy','res_type','chain','res_num','x','y','z','b_factor'])
    return pdb_info

print()
print('Comparing '+A+' with '+B)
print()

file1 = get_atom_info(A)
file2 = get_atom_info(B)

match = True
for c in file1['chain'].drop_duplicates():
    chain = file1[file1['chain']==c]
    for r in chain['res_num'].drop_duplicates():
        f1 = file1[(file1['chain']==c)&(file1['res_num']==r)]
        del f1['model']
        if D:
            del f1['atom_num']
        f2 = file2[(file2['chain']==c)&(file2['res_num']==r)]
        del f2['model']
        if D:
            del f2['atom_num']
        if np.array_equal(f1.values,f2.values):
            continue
        else:
            match=False
            print('Mismatch: ',c,r)
            if (C):
                print(f1)
                print(f2)
if match == True:
    print('Files match')
if (C == False) & (match == False):
    print()
    print('Use -d True to find specifics')








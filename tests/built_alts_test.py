#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 13:26:41 2022

@author: Tim Stachowski, PhD
Fischer Laboratory
St. Jude Children's Research Hospital

Test if alt locs are built into a multi conf pdb

"""

import numpy as np
import pandas as pd
import argparse
from glob import glob
import subprocess


CLI = argparse.ArgumentParser(description='Test if Ringer alt locs are in a PDB. Requires Phenix.')

CLI.add_argument(
    '-f',
    '--pdb',
    nargs="?",
    type=str,
    default=None,
    help='input file'
)

CLI.add_argument(
    '-r',
    '--alts',
    nargs="?",
    type=str,
    default=None,
    help='input file'
)

ARGS = CLI.parse_args()

A = ARGS.pdb
B = ARGS.alts

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

def rotalyze(file):
    call = '/Applications/phenix-1.20.1-4487/build/bin/phenix.rotalyze '+A+' >'+A[:-4]+'_rotalyze.txt'
    subprocess.call(call,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    rotamers = pd.read_csv(A[:-4]+'_rotalyze.txt',sep=':',header=0,skiprows=-1)
    ## extract residue number
    res_ns = []
    for j in rotamers['residue']:
        res_n = ''.join(i for i in j if i.isdigit())
        if res_n != '':
            res_ns.append(int(res_n))
        else:
            res_ns.append('NA')
    rotamers['res_n'] = res_ns
    rotamers = rotamers[rotamers['res_n']!='NA']
    rotamers['res_n'] = rotamers['res_n'].astype(int)
    rotamers['chain'] = [x.strip()[0] for x in rotamers['residue']]
    return rotamers

def test(file,rotamers):
    match = False
    for c in alts['chain'].drop_duplicates():
        chain = alts[alts['chain']==c]
        for r in chain['res'].drop_duplicates():
            alt = alts[(alts['chain']==c)&(alts['res']==r)]['rotamer'].values
            alt = list(alt)
            rot = rotamers[(rotamers['chain']==c)&(rotamers['res_n']==r)]['rotamer'].values
            rot = list(rot)
            if bool(set(alt)==set(rot)):
                match = True
            else:
                match = False
                print('Mismatch: ',c,r,rot,alt)

print('Testing: '+A+' against '+B)
rotamers = rotalyze(A)
pdb = get_atom_info(A)
alts = pd.read_csv(B)
test(A,rotamers)
print('Done')










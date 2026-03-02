#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 12:22:24 2025

@author: Tim Stachowski, PhD
Fischer Laboratory
St. Jude Children's Research Hospital
"""


import os
import sys
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

pdb = sys.argv[1]

# for debugging
#pdb = '6DMH.pdb'

def get_atom_info(file):
    if file.endswith(".pdb") or file.endswith(".ent"):
        ## radius
        pdb_info = []
        with open(file) as data:
            read_data = data.read()
            lines = read_data.splitlines()
            for line in lines:
                if line.startswith('ATOM'):
                    atom_type = line[12:16].strip()
                    #remove hydrogens
                    #if 'H' not in atom_type:
                    if (True):
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
    pdb_info = pdb_info[['model','chain','res_num','res_type','alt_loc','atom_type','occupancy','x','y','z','b_factor','atom_het']]
    return pdb_info

tests = ["1. Occ<0.1","2. MismatchedOcc","3. Lonely_Hs","4. SumOcc!=1","5. Wrong Highest Alt ID","6. Wrong # of Alt IDs"]

def minimum_occupancy(df):
    # lowest occupancy in PDB
    #minocc = df[df['atom_het']=='ATOM']['occupancy'].min()
    residues = df[df['atom_het']=='ATOM']
    residues = residues[residues['alt_loc']!=' ']
    residues = residues[~(residues['atom_type'].str.contains('H'))]

    residues = residues.groupby(by=['chain','res_num','alt_loc'])['occupancy'].min().reset_index()
    lowestocc = residues.nsmallest(1,'occupancy')
    print('1. Lowest rotamer occupancy: ')
    print(lowestocc.to_string(index=False))
    print('')

def low_occupancy(df):
    # residues < 0.1 occ
    residues = df[df['atom_het']=='ATOM']
    residues = residues[residues['alt_loc']!=' ']
    residues = residues[~(residues['atom_type'].str.contains('H'))]

    residues = residues.groupby(by=['chain','res_num','alt_loc'])['occupancy'].min().reset_index()
    lessthanten = residues[residues['occupancy']<0.1]
    print('2. Rotamers with occupancy < 0.1: ')
    if len(lessthanten)>0:
        print(lessthanten.to_string(index=False))
        print('')
        lessthanten[tests[0]] = 'X'
        lessthanten = lessthanten[['chain','res_num','alt_loc']+[tests[0]]]
        df = pd.merge(df,lessthanten,on=['chain','res_num','alt_loc'],how='left')
    else:
        print(None)
        print('')
        lessthanten[tests[0]] = ''
        lessthanten = lessthanten[['chain','res_num']+[tests[0]]]
        df = pd.merge(df,lessthanten,on=['chain','res_num'],how='left')
    return df

def mismatched_occupancies(df):
    # rotamers where the atomic occupancies are not the same
    mismatches = df[df['atom_het']=='ATOM']

    # consider alts or any side-chain?
    #mismatches = mismatches[mismatches['alt_loc']!=' ']

    # consider all atoms or non-h atoms?
    #mismatches = mismatches[~(mismatches['atom_type'].str.contains('H'))]

    mismatches = mismatches.groupby(by=['chain','res_num','alt_loc']).occupancy.nunique().eq(1).reset_index()
    mismatches = mismatches[mismatches['occupancy']==False]
    print('3. Rotamers with varying atomic occupancies: ')
    if len(mismatches)>0:
        print(mismatches.to_string(index=False))
        print('')
        mismatches[tests[1]] = 'X'
        mismatches = mismatches[['chain','res_num','alt_loc']+[tests[1]]]
        df = pd.merge(df,mismatches,on=['chain','res_num','alt_loc'],how='left')
    else:
        print(None)
        print('')
        mismatches[tests[1]] = ''
        mismatches = mismatches[['chain','res_num']+[tests[1]]]
        df = pd.merge(df,mismatches,on=['chain','res_num'],how='left')
    return df

def lonely_hydrogens(df):
    # if hydrogens do not have any accompanying side-chain
    residues = df[df['atom_het']=='ATOM']

    residues = pd.pivot_table(residues,index=['chain','res_num','alt_loc'],values='atom_type',aggfunc=list).reset_index()
    alts = residues[residues['alt_loc']!=' ']
    alts['atom_type'] = alts['atom_type'].astype(str)

    hydrogens=["['HA']","['HB1']","['HB2']","['HB3']","['H']","['HG']","['HG2']",\
    "['HG3']","['HE21']","['HE22']","['HD2']","['HD3']","['HB']","['HG11']",\
    "['HG12']","['HG13']","['HG21']","['HG22']","['HG23']","['HG1']","['HE2']",\
    "['HE3']","['HZ1']","['HZ2']","['HZ3']","['HD11']","['HD12']","['HD13']","['HD21']",\
    "['HD22']","['HD23']","['HE']","['HH11']","['HH12']","['HH21']","['HH22']","['HD1']",\
    "['HE1']","['HZ']","['HA2']","['HA3']","['HH']","['HH1']","['HH2']"]

    #lonelies = alts[alts['atom_type']=="['H']"]
    lonelies = alts[alts['atom_type'].isin(hydrogens)]
    print('4. Hydrogens without side-chains: ')
    if len(lonelies)>0:
        print(lonelies)
        print('')
        lonelies[tests[2]] = 'X'
        lonelies = lonelies[['chain','res_num','alt_loc']+[tests[2]]]
        lonelies = lonelies.drop_duplicates()
        df = pd.merge(df,lonelies,on=['chain','res_num','alt_loc'],how='left')
    else:
        print(None)
        print('')
        lonelies[tests[2]] = ''
        lonelies = lonelies[['chain','res_num']+[tests[2]]]
        lonelies = lonelies.drop_duplicates()
        df = pd.merge(df,lonelies,on=['chain','res_num'],how='left')
    return df


def occupancies_sum(df):
    # occupancies do not equal 1.0
    totals = df[df['alt_loc']!=' ']
    totals = totals[totals['atom_het']=='ATOM']
    totals = totals[~(totals['atom_type'].str.contains('H'))]

    totals = totals.groupby(by=['chain','res_num','alt_loc']).occupancy.mean().reset_index()
    totals = totals.groupby(by=['chain','res_num']).occupancy.sum().round(1).eq(1).reset_index()
    totals = totals[totals['occupancy']==False]
    print('5. Total occupancy across rotamers != 1.0: ')
    if len(totals)>0:
        print(totals.to_string(index=False))
        print('')
        totals[tests[3]] = 'X'
        totals = totals[['chain','res_num']+[tests[3]]]
        df = pd.merge(df,totals,on=['chain','res_num'],how='left')
    else:
        print(None)
        print('')
        totals[tests[3]] = ''
        totals = totals[['chain','res_num']+[tests[3]]]
        df = pd.merge(df,totals,on=['chain','res_num'],how='left')
    return df

def wrong_occupancy_alt(df):
    # if A isn't highest occupancy residue
    highest = df[df['atom_het']=='ATOM']
    highest = highest[~(highest['atom_type'].str.contains('H'))]

    highest = highest.groupby(by=['chain','res_num','alt_loc']).occupancy.mean().reset_index()
    highest = highest[highest['alt_loc']!= ' ']
    maxocc = highest.groupby(by=['chain','res_num']).occupancy.max().reset_index()
    highest = pd.merge(maxocc,highest,on=['chain','res_num','occupancy'],how='left')
    highest = highest[highest['alt_loc']!='A']
    highest = highest[highest['occupancy']>0.5]
    print('6. Highest occupancy rotamer not labeled A? (e.g. Occ_A < Occ_B): ')
    if len(highest)>0:
        print(highest.to_string(index=False))
        print('')
        highest[tests[4]] = 'X'
        highest = highest[['chain','res_num','alt_loc']+[tests[4]]]
        df = pd.merge(df,highest,on=['chain','res_num','alt_loc'],how='left')
    else:
        print(None)
        print('')
        highest[tests[4]] = ''
        highest = highest[['chain','res_num']+[tests[4]]]
        df = pd.merge(df,highest,on=['chain','res_num'],how='left')
    return df

def number_of_locids(df):
    # check correct number of labels
    residues = df[df['atom_het']=='ATOM']
    residues = residues[~(residues['atom_type'].str.contains('H'))]

    residues = residues[['chain','res_num','alt_loc']].drop_duplicates()
    residues = residues[residues['alt_loc']!=' ']
    residues['convert'] = [ord(x)-64 for x in residues['alt_loc']]
    convert = residues.groupby(by=['chain','res_num'])['convert'].max().reset_index()
    maxalt = residues.groupby(by=['chain','res_num'])['alt_loc'].count().reset_index()
    combine = pd.merge(maxalt,convert,on=['chain','res_num'],how='left')
    combine = combine[combine['convert']!=combine['alt_loc']]
    combine.columns = ['chain','res_num','#alt_locs','max alt_loc id']
    print('7. Improper series alt_loc id assigned (e.g. A->C instead of A->B): ')
    if len(combine)>0:
        print(combine.to_string(index=False))
        print('')
        combine[tests[5]] = 'X'
        combine = combine[['chain','res_num']+[tests[5]]]
        df = pd.merge(df,combine,on=['chain','res_num'],how='left')
    else:
        print(None)
        print('')
        combine[tests[5]] = ''
        combine = combine[['chain','res_num']+[tests[5]]]
        df = pd.merge(df,combine,on=['chain','res_num'],how='left')
    return df

#################

def main():
    print('')
    print('######        FLEX-CHECK          ######')
    print('###### MULTICONF MODEL VALIDATION ######')
    print('######   AUTHOR: Tim Stachowski   ######')
    print('')
    print('Input file: '+os.path.basename(pdb))
    print('')

    df = get_atom_info(pdb)

    # only interested in side-chains
    backboneatoms = ['N', 'CA', 'C', 'O']
    df = df[~(df['atom_type'].isin(backboneatoms))]

    #print(df)

    if (df['alt_loc']!=' ').any():

        #validation functions
        minimum_occupancy(df)
        df = low_occupancy(df)
        df = mismatched_occupancies(df)
        df = lonely_hydrogens(df)
        df = occupancies_sum(df)
        df = wrong_occupancy_alt(df)
        df = number_of_locids(df)

        ## export ##
        # only output residues with errors
        df = df[df.eq('X').any(axis=1)]
        # optional
        df = df.drop(['atom_het','atom_type','x','y','z','b_factor'],axis=1)
        df = df.drop_duplicates()
        # group by mean occupancy - to reduce the number of lines
        df = df.fillna(0)
        df = df.groupby(by=['model', 'chain', 'res_num', 'res_type', 'alt_loc']+tests)['occupancy'].mean().reset_index()
        dforder = ['model','chain','res_num','res_type','alt_loc','occupancy']+tests
        df = df[dforder]

        df[tests] = df[tests].replace(0,pd.NA)

        df['occupancy'] = [round(x,2) for x in df['occupancy']]

        df.to_csv('multiconf_refinement_check_output.csv',header=True,index=False)

        print('')
        print('Done.')
        print('')

    else:
        print('No residues with alt_loc IDs found...')
        print('Exiting...')


if __name__ == "__main__":
    main()

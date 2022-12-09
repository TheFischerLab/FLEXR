"""

Part of Ringer-refine
Script for building in alternative conformations with Coot

"""
import time
import subprocess
import os
import sys
import json
import pandas as pd
import coot
import coot_utils

## Settings
variables={
    'ca_or_all' : 1,       # add alt conf starting at ca
                           # or whole residue? 1 for whole, 0 for CA
    'sigma' : None,        # map contour level #1
    'sigma2' : None,       # map contour level #2
    'load_files' : True,   # useful for a dry run
    'build' : True,        # useful for a dry run
    'image' : False,       # take images of each residue to compare model and map. Broken.
    'file' : None,         # input model to use for building
    'help' : False,        # print options and exit
    'ref_model' : False,   # load input model. for example, single conf model
                           # a second time as a reference
    'exit' : True          # exit coot at the end
    }

## Set settings from user command line input
for i in variables:
    for j in sys.argv[3:]:
        arg = j.split('.')
        #convert string to bool
        if arg[1] == 'True':
            arg[1] = True
        if arg[1] == 'False':
            arg[1] = False
        if arg[0] == i:
            variables[i] = arg[1]
        else:
            continue

# Print options if user asks for help or missing input file
if variables['help'] is True:
    print('')
    print('Welcome to FLEXR... the Coot part')
    print('Need help?')
    print('See our GitHub page or paper for full options documentation')
    print('https://github.com/TheFischerLab/FLEXR')
    print('')
    print(json.dumps(variables,sort_keys=True,indent=4))
    print('Exiting Coot...')
    print('')
    coot.coot_no_state_real_exit(0)

# If no input file defined: exit
if variables['file'] is None:
    print('Input file not defined...')
    print('Exiting Coot...')
    coot.coot_no_state_real_exit(0)

# Possible amino acid residue atoms
atom_list = ['N','CA','C','O','CB','CG','SD','CE','CD1','CD2','OG1','CG2','OD1','ND2',\
'ND1','CE1','NE2','CG1','OD2','CD','OG','CE2','CZ','OH','NE','NH1','NH2','OE1','OE2',\
'NE1','CE3','CZ2','CZ3','CH2','SG','NZ']

# List of possible alt locs - more than you should need
alt_loc = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']

# Load list of alt locs to build from peak finding part
alt_confs_file = variables['file']+'_ringer_alts.csv'
alt_confs = pd.read_csv(alt_confs_file,header=0).sort_values(by=['chain','res_n'])

# Load files
if variables['load_files'] is True:
    print('Loading files...')
    try:
        model = variables['file']+'.pdb'
        imol = coot.read_pdb(model)
        # if taking images or inspecting these options are relevant
        if (variables['image']) or (variables['exit'] is False):
            if variables['ref_model'] is True:
                imol2 = coot.read_pdb(model)
            mtz = variables['file']+'.mtz'
            imol_map = coot.make_and_draw_map(mtz, "FWT", "PHWT", "", 0, 0)
            imol_map2 = coot.make_and_draw_map(mtz, "FWT", "PHWT", "", 0, 0)
            #viz settings
            coot.set_background_colour(1,2,2)
            coot.set_zoom(20)
            coot.set_reorienting_next_residue_mode(1)
            if variables['sigma'] is not None:
                coot.set_contour_level_in_sigma(imol_map,float(variables['sigma']))
            else:
                coot.set_contour_level_in_sigma(imol_map,0.3)
            if variables['sigma2'] is not None:
                coot.set_contour_level_in_sigma(imol_map2,float(variables['sigma2']))
            else:
                coot.set_contour_level_in_sigma(imol_map2,1)
            coot.set_map_colour(imol_map2,.3,.3,.3)
    except FileNotFoundError:
        print('Files missing...')
        sys.exit(1)

# Build conformers
if variables['build'] is True:
    mol_num = 0
    print('Building alternative conformers....')
    for p in alt_confs[['chain','res_n']].drop_duplicates().values:
        print(p)
        chain=p[0]
        resno=int(p[1])
        for k in range(0,len(alt_confs[(alt_confs['res_n']==resno)&\
                            (alt_confs['chain']==chain)])):
            print('Building alt number: ',k)
            rotamer = alt_confs[alt_confs['res_n']==resno]['rotamer'].values[k]
            print('For: ',resno,rotamer,chain)
            coot.set_go_to_atom_chain_residue_atom_name(chain,resno,'CA')
            if k == 0:
                coot.set_residue_to_rotamer_name(mol_num,chain,resno,'','',rotamer)
            if k == 1:
                coot.altconf()
                coot.set_add_alt_conf_split_type_number(variables['ca_or_all'])
                coot_utils.with_auto_accept(\
                [coot.add_alt_conf_py,mol_num,chain,resno,'','',1])
                coot.set_residue_to_rotamer_name(mol_num,chain,resno,'','B',rotamer)
            if k > 1:
                coot.altconf()
                coot.set_add_alt_conf_split_type_number(variables['ca_or_all'])
                coot_utils.with_auto_accept(\
                [coot.add_alt_conf_py,mol_num,chain,resno,'','A',1])
                for atom in atom_list:
                    coot.set_atom_string_attribute(\
                    mol_num,chain,resno,'',atom,'','alt-conf',alt_loc[k])
                coot.set_residue_to_rotamer_name(mol_num,chain,resno,'',alt_loc[k],rotamer)
            # take snapshot of each residue
            if variables['image'] is True:
                filename = os.getcwd()+'/'+variables['file']+'_'+str(resno)+'_ringer.'
                image_format = ' -png '
                r3d_call = 'render'+image_format+' -labels '+\
                           filename+'png'+' < '+filename+'r3d'
                try:
                    #import raster3d
                    coot.raster3d(filename+'r3d')
                    subprocess.call(r3d_call,shell=True)
                except subprocess.CalledProcessError as e:
                    print('Having trouble creating images....')
            time.sleep(0)
    coot.set_go_to_atom_chain_residue_atom_name(chain,resno+1,'CA')
    coot.write_pdb_file(imol,variables['file']+'_newconfs.pdb')
    print('Building finished....')
else:
    print('Building option set to false.')

if variables['exit'] is True:
    coot.coot_no_state_real_exit(0)
    print('Exiting Coot...')

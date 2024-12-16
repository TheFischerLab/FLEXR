#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#    FLEXR-MSA:
#    This script produces Ringer plots of proteins with dissimilar amino acid sequences
#    and enables comparisons of conformational changes.
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

import argparse
import os
import subprocess
import sys
import time
from glob import glob
import numpy as np
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
from Bio import AlignIO

#chi = ARGS.chi
#chi = chi[0]
#colors = ARGS.colors
#colors = colors.strip().split(',')
#safety = ARGS.safety
#reload = ARGS.reload
#pearson = ARGS.pearson
#render = ARGS.render

def intro():
    ## Intro messages
    print(' ')
    print('Welcome to FLEXR-MSA!')
    print('')
    time.sleep(1)
    print('This program produces Ringer plots from proteins with')
    print('non-identical amino acid sequences to enable comparing')
    print('conformational changes')
    time.sleep(1)
    print('Let\'s get started')
    print('')
    time.sleep(5)

def check_dependencies():
    ##check for dependencies
    ## try to find matplotlib
    print('Checking dependencies...')
    try:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
    except ImportError:
        print('Matplotlib not found...')
        sys.exit(1)
    ## need muscle installed for alignment
    ## try: brew install muscle
    MUSCLE = subprocess.call("test -e '{}'".format('muscle'),shell=True)
    if MUSCLE == 0:
        print("Can't find MUSCLE, make sure a variable for MUSCLE is set in your path")
        sys.exit(1)
    else:
        print('Dependencies look good...')
        print(' ')

        ## load and extract sequences from ringer files
        print('Grabbing all the _ringer.csv files in the directory...')
        print('Working on them in alphabetical order...')
        print(' ')

        return plt

def find_files():

    ringer_outputs = glob('*_ringer.csv')
    ringer_outputs = [x for x in ringer_outputs if 'peak_finder' not in x]
    ringer_outputs = sorted(ringer_outputs)

    return ringer_outputs

def parse_ringers(ringer_outputs,chiopt,step):

    ## make tmp dir to store new ringer files
    try:
        subprocess.call(['mkdir','tmp'],stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print('Overwriting old data....')
        print('')

    newringers = []
    new_outputs = []

    #make columns
    try:
        step = int(step[0])
    except:
        step = step
    angles = np.arange(0,360,step)
    columns = ['res','map','chi','?']+[str(x) for x in angles]

    for ringer in ringer_outputs:

        #read in file
        df = pd.read_csv(ringer,names=None)

        try:
            df.columns = columns
        except:
            print('Error interpretting file header.')
            print('Try using --step option.')
            print('Done.')
            sys.exit()

        ## split by map type
        df = df[df['map']=='2mFo-DFc']
        del df['map']

        # split by chain
        chains = [x.split()[1][0] for x in df['res']]
        df['chain'] = [x.split()[1][0] for x in df['res']]

        chains = pd.Series(chains).drop_duplicates()
        chains = list(chains)

        ## add res# column
        df['resseq'] = [''.join(x for x in j if x.isdigit()) for j in df['res']]

        ## res type colum
        df['res_type'] = [str(x)[:3] for x in df['res']]

        #save new files
        for chain in chains:
            tmp = df[df['chain']==chain]
            filename = ringer.split('_ringer.csv')[0]
            filename = filename+'_'+chain+'_ringer.csv'
            tmp.to_csv('./tmp/'+filename,header=True,index=False)
            new_outputs.append(filename)

            ## take only residues with proper chi values
            #tmp = tmp[tmp['chi']==chiopt]
            newringers.append(tmp)

    columns = tmp.columns

    #new_outputs = glob('./tmp/*ringer.csv')
    #new_outputs = [x.split('/')[2] for x in new_outputs]

    return newringers,columns,new_outputs,angles

def normalization(normopt,ringers):
    ##normalize ringer profiles
        if (normopt):
            try:
                print('Attempting normalization...')
                normed_ringers = []
                for ringer in ringers:
                    norm_ringer = []
                    for chi in ['chi1','chi2','chi3','chi4']:
                        tmp = ringer[ringer['chi']==chi]
                        for res in tmp.res.drop_duplicates():
                            tmp2 = tmp[tmp['res']==res]
                            vals = tmp2.values[:][0]
                            sigmas = vals[3:-3]

                            norm_sigmas = [float(i)/sum(sigmas)*10 for i in sigmas]
                            output = list(vals[:3])+list(norm_sigmas)+list(vals[-3:])
                            norm_ringer.append(output)

                    norm_ringer = pd.DataFrame(norm_ringer)
                    norm_ringer.columns = ringer.columns
                    normed_ringers.append(norm_ringer)
                print('Normalization done.')
                print('')
                return normed_ringers
            except:
                print('Problem normalizing Ringer profiles.')
                print('')
                return ringers
        else:
            print('Not normalizing Ringer profiles.')
            print('')
            return ringers

def make_color_legend(ringer_outputs,coloropt):
    ## Assign colors if user doesn't define enough
    coloropt = coloropt.split(',')
    if len(ringer_outputs) != len(coloropt):
        print('Not enough colors defined.')
        print('Assigning random colors...')
        print(' ')
        try:
            from matplotlib import colors as mcolors
            colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
            by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                    for name, color in colors.items())
            sorted_names = [name for hsv, name in by_hsv]
            colors = np.random.choice(sorted_names,len(ringer_outputs))
            #colors = sorted_names[:len(ringer_outputs)]
            return colors
        except ImportError:
            print('Cannot find matplotlib...')
            sys.exit(1)
    else:
        return coloropt

def create_legend(ringer_outputs,colors):
    ## print pairs of colors and files
    print('Legend: ')
    print('Data Color')
    with open('plot_legend.txt','w') as f:
        color_matchings = list(zip(ringer_outputs,colors))
        for x in color_matchings:
            print(x[0],x[1])
            f.write(x[0]+' '+x[1])
            f.write('\n')
    print(' ')
    time.sleep(1)
    ## swap amino acid codes from one to three letter to match ringer output for alignment
    amino_acids = [["Alanine","A","ALA"],["Arginine","R","ARG"],
                   ["Asparagine","N","ASN"],["Aspartic Acid","D","ASP"],
                   ["Cysteine","C","CYS"],["Glutamic Acid","E","GLU"],
                   ["Glutamine","Q","GLN"],["Glycine","G","GLY"],
                   ["Histidine","H","HIS"],["Isoleucine","I","ILE"],
                   ["Leucine","L","LEU"],["Lysine","K","LYS"],
                   ["Methionine","M","MET"],["Phenylalanine","F","PHE"],
                   ["Proline","P","PRO"],["Serine","S","SER"],
                   ["Threonine","T","THR"],["Tryptophan","W","TRP"],
                   ["Tyrosine","Y","TYR"],["Valine","V","VAL"]]
    amino_acids = np.array(amino_acids)
    return amino_acids

def gen_fasta(ringer_outputs,ringers,amino_acids):
    ## write fasta for alignment
    print('Extracting sequences...')
    print('Sequences saved as .fasta files...')
    print(' ')
    with open('ringer_alignment.fasta','w') as f:
        for m,n in list(zip(ringer_outputs,ringers)):
            SEQ = n.drop_duplicates(subset='res')
            SEQ = list(SEQ['res_type'])
            for j in range(len(SEQ)):
                for k in amino_acids:
                    if SEQ[j] == k[2]:
                        SEQ[j] = k[1]
            SEQ = ''.join(SEQ)
            f.write('>'+m)
            f.write('\n')
            f.write(SEQ)
            f.write('\n')

def muscle_alignment(reload):
    ## use MUSCLE to do a multi-sequence alignment
    print('Aligning sequences...')
    print(' ')
    time.sleep(1)
    subprocess.check_output(['muscle','-align','ringer_alignment.fasta','-output',
                             'ringer_alignment_muscle.fasta'],text=True)

    print('')
    print('Saving alignment to ringer_alignment.fasta ...')
    time.sleep(1)

    ## reorganize the alignment into pandas dataframe
    alignment = AlignIO.read('ringer_alignment_muscle.fasta','fasta')

    ids = []
    for record in alignment:
        ids.append(record.id)

    length = len(alignment[0].seq)
    matchings = []
    for i in range(length):
        res = []
        res.append(i)
        for record in alignment:
            SEQ = record.seq
            res.append(SEQ[i])
        matchings.append(res)

    print('Saving alignment with new ID to alignment_new_index.csv ...')
    print(' ')
    print(' ')
    time.sleep(1)
    matchings = pd.DataFrame(matchings,columns=['res']+ids)
    matchings = matchings.reindex(sorted(matchings.columns), axis=1)

    if reload == False:
        matchings.to_csv('alignment_new_index.csv',header=True,index=False)
    else:
        matchings = pd.read_csv('alignment_new_index.csv',header=0)

    return matchings

def reindexing(matchings,ringer_outputs,ringers):
    ## re-index ringer outputs with numbering from alignment
    for i in range(len(matchings.columns)-1):
        new_index = matchings.iloc[:,i].replace('-',np.nan)
        new_index = new_index.dropna()
        old_index = ringers[i]['res'].drop_duplicates()
        if len(new_index) == len(old_index):
            reindex = np.column_stack((old_index,new_index.index.values))
            reindex = pd.DataFrame(reindex,columns=['res','new_resseq'])
            reindex['new_resseq'] = reindex['new_resseq'].astype(int)
            ringers[i] = pd.merge(ringers[i],reindex,how='left',on='res')
        else:
            print(ringer_outputs[i],i,'not equal length lists')
            sys.exit()
    return ringers

def safety(safeopt,chiopt,amino_acids):
    ## which residues have certain unbranched torsion angles
    if safeopt == 'True':
        print('Safety on. see -h to read more')
        print(' ')
        if chiopt == 'chi1':
            chi_allow = "Ser Gln Asn Glu Asp Arg Lys Met Cys Leu".upper().split()
        if chiopt == 'chi2':
            chi_allow = "Gln Glu Arg Lys Met Ile".upper().split()
        if chiopt == 'chi3':
            chi_allow = "Lys Arg Met".upper().split()
        if chiopt == 'chi4':
            chi_allow = "Lys Arg".upper().split()
    else:
        print('Safety off - be careful')
        print(' ')
        chi_allow = amino_acids[:,[2]].flatten()
    return chi_allow

def chi_divide(ringers,chi_allow):
    ## take only residues with proper chi values
    for i in range(len(ringers)):
        ringers[i] = ringers[i][ringers[i]['res_type'].isin(chi_allow)]
    return ringers

def assemble_output(ringer_outputs,ringers):
    ## assemble separate ringers into a single dataframe for easy plotting
    final_df = pd.DataFrame([])
    for i in range(len(ringers)):
        ringers[i]['file'] = ringer_outputs[i]
        final_df = pd.concat([final_df,ringers[i]])
    final_df['new_resseq'] = final_df['new_resseq'].astype(int)
    final_df.to_csv('final_ringer_dataframe.csv',header=True,index=False)
    return final_df

def plotting(plt,chiopt,final_df,ringers,angles,colors):
    ## plotting
    # loop over chis
    for chiopt in ['chi1','chi2','chi3','chi4']:
        tmpringers = []
        for ringer in ringers:
            ringer1 = ringer[ringer['chi']==chiopt]
            tmpringers.append(ringer1)
        tmp = final_df[final_df['chi']==chiopt]
        print('Plotting results for: ', chiopt)
        time.sleep(1)
        print('Plot images saved in a new dir: ./'+chiopt)
        print('')
        ## try to make directory to save images:
        try:
            subprocess.call(['mkdir',chiopt],stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print('Overwriting old data....')
            #print('Please delete the directory ',chiopt,' and re-run Ringer-delta.')
            #sys.exit(1)

        ## plotting
        for i in set(tmp['new_resseq']):
            plt.figure(figsize=(10,5))
            plt.subplots_adjust(right=0.3)
            for j in range(len(tmpringers)):
                if i in tmpringers[j]['new_resseq'].values:
                    LABEL = [str(x).strip() for x in tmpringers[j][tmpringers[j]['new_resseq']==i]\
                    .iloc[0,[0,-3,-2,-1]].values]
                    LABEL = ' '.join(LABEL)
                    plt.plot(angles,tmpringers[j][tmpringers[j]['new_resseq']==i].iloc[0,3:-5],\
                    color=colors[j],label=LABEL)
                    plt.hlines(y=[0,0.3,1.0],xmin=0,xmax=360,colors=['black','grey','black'],\
                    linestyle=['solid','dashed','dashed'])
                    plt.xlim(0,355)
                    plt.xticks((0,180,355), size = 12)
                    plt.yticks(size=12)
                    plt.xlabel("X$^%s$ (˚)" % chiopt[-1],fontsize=14)
                    plt.ylabel("2mFo-DFc ($\\sigma$)",fontsize=14)
                    plt.title(i)
                else:
                    continue
            #plt.legend(bbox_to_anchor=(1.1, 1.05),fontsize=8)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
            plt.savefig(os.getcwd()+'/'+chiopt+'/'+str(i)+chiopt+'.png')
            plt.close()

def pearson(pearsonopt,final_df,ringers,chiopt):
    if (pearsonopt):
        #try:
            import re
            from matplotlib import pyplot as plt
            from matplotlib.ticker import AutoMinorLocator
            from src.tools import ringer_tools
            from src.tools.ringer_tools import pearson_correlation
            cc_outputs = []
            for chiopt in ['chi1','chi2','chi3','chi4']:
                cc_output = pearson_correlation(final_df,ringers,chiopt)
                cc_outputs.append(cc_output)

                cc_output = cc_output[cc_output['file1']!=cc_output['file2']]
                cc_output['residues'] = [x.split()[2] for x in cc_output['res1']]
                cc_output['residues'] = [int(x[:]) for x in cc_output['residues']]
                cc_output = cc_output.sort_values(by='residues',ascending=True)

                ref = cc_output['file1'].values[0]
                tmp = cc_output[cc_output['file1']==ref]
                tmp = tmp.groupby(by='residues')['cc'].median().reset_index()
                tmp.columns = ['residues','cc']

                #look at CC across sequence
                fig,ax = plt.subplots()
                ax.plot(tmp.residues,tmp.cc)
                ax.scatter(tmp.residues,tmp.cc,s=3)
                plt.ylabel('Pearson CC')
                plt.xlabel('Alignment residue (n)')
                #plt.minorticks_on()
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                plt.savefig(chiopt+'_cc.png')
                plt.close()

            return cc_outputs
        #except ImportError:
            print('Needs rigner_tools.py script')

def pymol_render(pymolopt,pearsonopt,ringer_outputs,cc_outputs):
    if (pymolopt) and (pearsonopt):
        try:
            import re
            from src.tools import ringer_tools
            from src.tools.ringer_tools import pearson_correlation,render_by_attribute
        except ImportError:
            print('Needs ringer_tools.py script')
            sys.exit(1)
        i = 1
        for cc_output in cc_outputs:
            file = ringer_outputs[0][:-13]
            cc_output = cc_output[cc_output['file1']==ringer_outputs[0]]
            cc_output = cc_output.groupby(by='res1')['cc'].median().reset_index()
            cc_output['res_t'] = [x[:3].strip() for x in cc_output['res1']]
            cc_output['res_n'] = [re.findall(r'\d+',x)[0] for x in cc_output['res1']]
            render_by_attribute(file+'.pdb', file+'_chi'+str(i)+'_out.pdb',cc_output['res_t'],\
            cc_output['res_n'],cc_output['cc'])
            i += 1

def MSAmain(colorin,reload,chi,safeopt,pearsonopt,render,step,normopt):

        intro()

        plt = check_dependencies()

        ringer_outputs = find_files()

        ringers,columns,ringer_outputs,angles = parse_ringers(ringer_outputs,chi,step)

        ringers = normalization(normopt,ringers)

        colors = make_color_legend(ringer_outputs,colorin)

        amino_acids = create_legend(ringer_outputs,colors)

        gen_fasta(ringer_outputs,ringers,amino_acids)

        matchings = muscle_alignment(reload)

        ringers = reindexing(matchings,ringer_outputs,ringers)

        chi_allow = safety(safeopt,chi,amino_acids)

        ringers = chi_divide(ringers,chi_allow)

        final_df = assemble_output(ringer_outputs,ringers)

        plotting(plt,chi,final_df,ringers,angles,colors)

        cc_output = pearson(pearsonopt,final_df,ringers,chi)

        pymol_render(render,pearsonopt,ringer_outputs,cc_output)

        print('Done.')

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Tools for analyzing Ringer measurements.
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN.
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
#
import argparse
import time
import os
import subprocess
import sys
import numpy as np
import pandas as pd

def create_parser():
    """Arguments for flipper"""

    CLI = argparse.ArgumentParser()

    CLI.add_argument(
        '-f',
        '--filename',
        nargs="?",
        type=str,
        default=None,
        help='input file'
    )

    CLI.add_argument(
        '-f2',
        '--filename2',
        nargs="?",
        type=str,
        default=None,
        help='input file'
    )

    CLI.add_argument(
        '-g',
        '--geotolerance',
        nargs="?",
        type=int,
        default=30,
        help='Tolerance for match between Ringer measured chi and ideal \
        chi in library. Default = 30'
    )

    CLI.add_argument(
        '-t',
        '--sigmathreshold',
        nargs=1,
        type=float,
        default=0.3,
        help='sigma threshold for peak finding. default = 0.3'
    )

    CLI.add_argument(
        '-ph',
        '--peakheight',
        nargs=1,
        type=float,
        default=0.03,
        help='Required height of peaks. default = 0.03'
    )

    CLI.add_argument(
        '-pp',
        '--peakprominence',
        nargs=1,
        type=float,
        default=0.05,
        help='Required prominence of peaks. default = 0.05'
    )

    CLI.add_argument(
        '-pw',
        '--peakwidth',
        nargs=1,
        type=int,
        default=1,
        help='Required width of peaks. default = 1'
    )

    CLI.add_argument(
        '-pd',
        '--peakdistance',
        nargs=1,
        type=int,
        default=5,
        help='Required minimal horizontal distance between neighboring peaks. \
        Default = 5'
    )

    CLI.add_argument(
        '-p',
        '--plot',
        nargs=1,
        type=str,
        default=False,
        help='Save individual plots showing peak finding results? This is slow. \
        Default = False'
    )

    CLI.add_argument(
        '-s',
        '--step',
        nargs=1,
        type=int,
        default=5,
        help='step sized used for sigma measurements'
    )

    return CLI


def ringer_parser(file,step):
    """ Organize Ringer outputs from mmtbx """

    step=step
    angles = np.arange(0,360,step=step)

    dataframe = pd.read_csv(file,header=None)
    dataframe.columns = ['res','map','chi','peak']+[*angles]

    ## extract residue number
    res_ns = []
    for j in dataframe['res']:
        res_n = ''.join(i for i in j if i.isdigit())
        res_ns.append(int(res_n))
    dataframe['res_n'] = res_ns

    ## extract chain
    chain = []
    for j in dataframe['res']:
        chain_n = ''.join(i for i in j if not i.isdigit())
        chain_n = chain_n[-1]
        chain.append(chain_n)
    dataframe['chain'] = chain

    ## extract res_type
    dataframe['res_type'] = dataframe['res'].str[:3]

    ## map type
    dataframe = dataframe[dataframe['map']=='2mFo-DFc']

    return dataframe,angles


def peak_find(file,dataframe,chi,sigma_threshold,plot,height,prominence,width,distance,angles):
    """ Find peaks in Ringer plots """

    try:
        from scipy.signal import find_peaks
    except ImportError:
        print('Please install SciPy')

    if plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            D = False
            print('Matplotlib not found... not plotting')
            sys.exit(1)
        try:
            subprocess.check_output(['mkdir',file[:-4]+'_'+chi+'_peaks'])
        except subprocess.CalledProcessError as e:
            print('Please delete the directory ',chi+'_peaks',' and re-run.')
            sys.exit(1)

    output = []

    ## loop through residues
    dataframe = dataframe[dataframe['chi']==chi]
    for index,i in dataframe.iterrows():
        res_all = i['res']
        res_n = i['res_n']
        chain = i['chain']
        res = i['res_type']
        chi = chi
        sigma_raw = i[angles].values
        sigma_raw = [float(x) for x in sigma_raw]
        sigma_raw = np.array(sigma_raw)
        dat = np.column_stack((angles,sigma_raw))
        dat = [(x,y) for (x,y) in dat if y >= sigma_threshold]
        sigma_above_threshold = [y for (x,y) in dat]
        angles_cut = [x for (x,y) in dat]
        peaks, properties = find_peaks(sigma_above_threshold,
                                       height=height,
                                       prominence=prominence,
                                       width=width,
                                       distance=distance)
        peak_n = len(peaks)
        area = []
        lefts = []
        rights = []
        for peak in range(peak_n):
            left = properties['left_bases'][peak]
            lefts.append(left)
            right = properties['right_bases'][peak]
            rights.append(right)
            data_to_int = np.column_stack([angles_cut[left:right],
                                           sigma_above_threshold[left:right]])
            angles_int = [j for (j,y) in data_to_int]
            sigma_int = [y for (j,y) in data_to_int]
            area.append(np.trapz(sigma_int,angles_int))
        areas_norm = []
        if len(area) > 1:
            if (lefts[0] == lefts[-1]) or (rights[0] == rights[-1]):
                if area[0] > area[-1]:
                    area[0] = area[0]-area[-1]
                    area[0] = np.absolute(area[0])
                else:
                    area[-1] = area[-1] - area[0]
                    area[-1] = np.absolute(area [-1])
                sum_area = sum(area)
                for areai in area:
                    areas_norm.append(areai/sum_area)
            else:
                sum_area = sum(area)
                for areai in area:
                    areas_norm.append(areai/sum_area)
        if len(area) == 1:
            areas_norm = 1.0
        if len(area) == 0:
            areas_norm = 0
        if plot:
            plt.figure(figsize=(5,5))
            plt.plot(angles,sigma_raw,color='black')
            plt.plot(angles_cut,sigma_above_threshold,'red',linestyle='dotted')
            plt.plot([angles_cut[x] for x in peaks],
                     [sigma_threshold]*len(peaks), "x")
            plt.axhline(y=sigma_threshold, linestyle="--", color="gray")
            plt.legend(['raw data',
                        'data used for peak \n finding and integration',
                        'peak','sigma cutoff'],
                        bbox_to_anchor=(1.05, 1),
                        loc='upper left')
            plt.xlabel(chi+'  angle (˚)')
            plt.ylabel('sigma')
            plt.title(res)
            if len(peaks)>0:
                plt.ylim(np.min(sigma_raw),np.max(sigma_raw)+sigma_threshold)
            plt.tight_layout()
            plotname = str(res_n)+'_'+chain+'_'+file+'_'+chi+'_peaks.png'
            plt.savefig(plotname,dpi=100)
            plt.close()
            try:
                subprocess.check_output(['mv',plotname,file[:-4]+'_'+chi+'_peaks'])
            except subprocess.CalledProcessError as e:
                continue
        peaks = [angles_cut[x] for x in peaks]
        output.append((res_all,res_n,res,chi,peaks,peak_n,areas_norm,chain))
    output = pd.DataFrame(output)
    output.columns = ['res','res_n','res_type','chi','peak_angles','peaks_n','areas_norm','chain']
    return output


def find_flips(peaks1,peaks2,chi):
    """Find changes (flips) in major/minor conformation of a residue between
       two datafiles based on peak integration """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]

    # only consider ringer plots with 2 peaks
    dataframe_f1 = dataframe_f1[dataframe_f1['peaks_n']==2]
    dataframe_f2 = dataframe_f2[dataframe_f2['peaks_n']==2]

    # only consider one chi angle
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    # format
    dataframe_f1['areas_norm'] = dataframe_f1['areas_norm'].astype(str).str[1:-1]
    dataframe_f1['areas_norm'] = [x.split(',') for x in dataframe_f1['areas_norm']]
    dataframe_f1['areas_norm'] = \
    [(float(x),float(y)) for x,y in dataframe_f1['areas_norm']]
    dataframe_f2['areas_norm'] = dataframe_f2['areas_norm'].astype(str).str[1:-1]
    dataframe_f2['areas_norm'] = [x.split(',') for x in dataframe_f2['areas_norm']]
    dataframe_f2['areas_norm'] = \
    [(float(x),float(y)) for x,y in dataframe_f2['areas_norm']]

    #format
    dataframe_f1['res_n'] = dataframe_f1['res_n'].astype(int)
    dataframe_f2['res_n'] = dataframe_f2['res_n'].astype(int)

    flips_pair = []

    for j in range(len(dataframe_f1)):
        res_n = int(dataframe_f1.iloc[j,0])
        if res_n in dataframe_f2['res_n'].values:
            f1_p1 = dataframe_f1.iloc[j,5][0]
            f2_p1 = dataframe_f2[dataframe_f2['res_n']==res_n]['areas_norm'].values[0][0]
            if f1_p1 > 0.5 > f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
            if f1_p1 < 0.5 < f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
    flips = pd.DataFrame(flips_pair,columns=['res_n',chi])
    return flips

def find_gainloss(peaks1,peaks2,chi):
    """ Calculate the changes in number of peaks
    (conformations) between two datafiles
    """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]

    # only consider one chi angle
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    dataframe_f1 = dataframe_f1[['res_n','peaks_n','res']]
    dataframe_f1.columns = ['res_n','peak_cryo','res']
    dataframe_f2 = dataframe_f2[['res_n','peaks_n','res']]
    dataframe_f2.columns = ['res_n','peak_rt','res']

    dataframe = pd.merge(dataframe_f1,dataframe_f2,on='res_n',how='left')
    dataframe = dataframe.fillna(0)

    differences = dataframe.peak_cryo.values - dataframe.peak_rt.values
    residues = dataframe['res_n'].values[:]

    dataframe = np.column_stack((residues,differences))

    dataframe = pd.DataFrame(dataframe,columns = ['res','peak_gain_loss'])
    dataframe['res'] = dataframe['res'].astype(int)
    dataframe['chi'] = chi

    return dataframe

def match_and_build(k,res,library_tmp,rotamer_geometry):
    """ This is for matching peaks found from find_peaks to rotamer geometries """
    chi_labels = ['chi1_mean','chi2_mean','chi3_mean','chi4_mean']
    labels = chi_labels[:len(k)]
    conf = np.column_stack((labels,k))
    conf = pd.DataFrame(conf.T)
    conf.columns = conf.iloc[0]
    conf = conf.drop(conf.index[0])
    res_lib = library_tmp[conf.columns]
    match = np.isclose(conf.values.astype(float).astype(int),\
    res_lib.values.astype(int),atol=rotamer_geometry)
    matches = [x.all() for x in match]
    matches = [i for i, x in enumerate(matches) if x]
    if len(matches) != 0:
        ## if multiple matches to that combination of angles, take the most common one
        rot_index = library_tmp.iloc[matches,:].sort_values(by='frequency%')
        res2build = rot_index['rotamer'].values[0]
        print(res,res2build)
        return res2build

def pearson_correlation(final_df,ringers,chi):
    """ Produces matrices and a CSV values of pairwise Pearson
    CC calculations per residue.
    Needs an alignment from Ringer-delta
    """

    try:
        import seaborn as sns
    except ImportError:
        print('Please install Seaborn')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('Please install matplotlib')

    print(' ')
    print('Calculating Pearson CC values....')
    time.sleep(1)

    try:
        subprocess.check_output(['mkdir',chi+'_cc'])
    except subprocess.CalledProcessError as e:
        print('Please delete the directory ',chi+'_cc',' and re-run.')
        sys.exit(1)

    files = pd.read_csv('alignment_new_index.csv',header=0)
    files = files.columns
    output = []
    for i in set(final_df['new_resseq']):
        res_df = []
        LABELS = []
        for j in range(len(ringers)):
            file = files[j]
            if i in ringers[j]['new_resseq'].values:
                LABEL = [x.strip() for x in ringers[j][ringers[j]['new_resseq']==i]\
                .iloc[0,[0,-2,-3]].values]
                LABEL = ' '.join(LABEL)
                LABELS.append(file+','+LABEL)
                tmp = ringers[j][ringers[j]['new_resseq']==i].iloc[0,3:-4]
                res_df.append(tmp)
        res_df = pd.DataFrame(res_df,columns=None)
        res_df = res_df.T
        res_df.columns = LABELS
        corrs = []
        for m in range(len(LABELS)):
            x = res_df.iloc[:,m]
            for n in range(len(LABELS)):
                y = res_df.iloc [:,n]
                corr = x.corr(y,method='pearson')
                corrs.append(corr)
                output.append((LABELS[m].split(',')[0],LABELS[m].split(',')[1],\
                LABELS[n].split(',')[0],LABELS[n].split(',')[1],corr))
        corrs = np.array(corrs)
        matrix = corrs.reshape((len(LABELS),len(LABELS)))
        ax = sns.heatmap(matrix,xticklabels=LABELS,yticklabels=LABELS,vmax=1,\
        annot=True,linewidths=.5)
        ax.set(title='Pearson Correlation Coefficient:   '+str(i))
        plt.tight_layout()
        plt.savefig(os.getcwd()+'/'+chi+'_cc/'+'cc_'+str(i)+'.png',dpi=300)
        plt.close()
    output = pd.DataFrame(output)
    output.columns = ['file1','res1','file2','res2','cc']
    output.to_csv('cc_'+chi+'.csv',header=True,index=False)
    print('Values saved to: cc_'+chi+'.csv')
    print('Matrix plots saved to: ./'+chi+'_cc')
    print(' ')

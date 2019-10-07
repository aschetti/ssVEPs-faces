#!/usr/bin/python

from scipy.io import loadmat
import numpy as np
import pandas as pd
from os.path import join as opj
from os.path import isfile
import argparse

'''
    It's meant to be used as an executable script
    from the command line, check the parser arguments 
    to see what you need to pass as inputs.
    Compatible with python 2 and 3
'''

def getRT(resp, times):
    # make array as long as resp
    resp = np.repeat(resp.reshape([resp.shape[0], 1]), 
            times.shape[0], axis=1)
    rt = times - resp
    return rt[((rt >= 200.0) & (rt <= 800.0))]

def cleandata(data):
    data = np.sort(data.reshape([1, np.prod(data.shape)]))
    return data[data > 0.0]


def cond(data, cond_type):
    alldata = data['mat'][:, 0]
    return np.where(alldata == cond_type)


def computeBlock(resp, times, cond_idx):
    times = times[:, cond_idx]
    times = cleandata(times)
    return getRT(resp, times)


parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-d', '--datadir', metavar='PATH', required=True)
parser.add_argument('-s', '--subjects', type=int, required=True)
parser.add_argument('-b', '--nblocks', type=int, required=True)
parser.add_argument('-f', '--filename', type=str, required=True)
args = parser.parse_args()

# conds and labels could be made as inputs,
# but it would be a bit annoying to write them 
# on the command line, given that they are quite a few
# It's not actually an issue given that they are equal
# across pilot 2 and exp 1

conds = np.array([1, 2, 3, 30, 60, 90,\
        5, 10, 15, 150, 300, 450])

labels = ['a_up_f', 'n_up_f', 'i_up_f', 'a_dw_f', 'n_dw_f', 'i_dw_f',\
        'a_up_m', 'n_up_m', 'i_up_m', 'a_dw_m', 'n_dw_m', 'i_dw_m']


blocks = range(1, args.nblocks + 1)
subjects = range(1, args.subjects + 1)

labels = labels * args.subjects
sub_list = np.sort(np.repeat(subjects, len(labels)/args.subjects))

m = []
for subj in subjects:
    t = []
    r = []
    for b in blocks:
        fname = opj(args.datadir, 
                'sub {} block {}.mat'.format(subj, b))
        if isfile(fname): 
            data = loadmat(fname)

        #we clean resp here to avoid unbalanced array shapes
        response = cleandata(data['responses'])
        r.append(response)
        t.append(data['times_mat'])
        
    times = np.hstack(t)
    resp = np.hstack(r)

    for j in conds:
        m.append(np.median(computeBlock(resp,
                                    times,
                                    cond(data, j)))
                                    ) 
        
dict_ds = {'subjects': sub_list, 'rt': m, 'cond': labels}
pd.DataFrame(dict_ds).to_csv('{}.csv'
        .format(args.filename), 
        index=False,
        na_rep='nan')

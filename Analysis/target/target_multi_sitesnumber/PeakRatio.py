import glob # to load general files

import Target_cohesin_loading.util_visual as utvis
from Target_cohesin_loading.funcs import *
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import ndimage

import pandas as pd
import h5py 
import time
import sys
import warnings
import ast

window_size = 1
rep = 1 
mon = 1000
site = 10
ba_time = 2500000
hist_dict = {}
c=1
face = 1
life = 50
birth = 1
base_loading = 0.0001
sep = 100
paramdict_CTCF={
            'CTCF_facestall':[face, face],
            'LEF_lifetime':[life, life],
            'LEF_birth':[base_loading, birth],
            'LEF_separation':sep
            }
total_sites = 10000
# modifying density based on new loading on target
density_multiplier = 1 + (((birth-base_loading)/base_loading)/total_sites)
paramdict_CTCF['LEF_separation'] /= density_multiplier
LEFNum = int(paramdict_CTCF['LEF_separation'])

path_dict = {}

directory=dire = '/scratch1/rahmanin/target_loading_cohesin/target/target_mult_sites/sims/'
#'/scratch1/rahmanin/target_loading_cohesin/target/target_mult_sites/sims_chn_md/'

for fname in glob.glob(directory+'folder*'):
    path_dict[fname.split('sims/')[1][:]]= fname
path_dict = dict(sorted(path_dict.items()))

file = open('data/relative_ratio_target_targetnum_res_rev_further_15kbnews.csv','w')
file.write('birth,target_s_num,deltactcf,clife,cof,sep,face,ratio\n')
dire = '/scratch1/rahmanin/target_loading_cohesin/target/target_mult_sites/sims/'
for name in list(path_dict.keys())[:]:
    try:
        params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
        face, back, clife, cof, life, slife, birth, target_s_num, deltactcf, pause, sep, site, mon, rep, step, vel = params
        if deltactcf!= 1600: continue
        c+=1
        mapN=mon*site
        lefs = h5py.File(dire+name+'/LEFPositions.h5','r')["positions"]
        lef_lefts = lefs[ba_time:,:,0].flatten()
        lef_rights = lefs[ba_time:,:,1].flatten()
        ctcfrightlist = np.array(h5py.File(dire+name+'/LEFPositions.h5','r')['CTCF_sites_right'])
        ctcfleftlist = np.array(h5py.File(dire+name+'/LEFPositions.h5','r')['CTCF_sites_left'])
        lst = np.array(list(ctcfrightlist))# + list(ctcfleftlist))
        
        lef_lefts = lefs[ba_time:,:,0].flatten()
        lef_rights = lefs[ba_time:,:,1].flatten()
        lef_positions = np.hstack((lef_lefts,lef_rights))
    
        a,b = np.histogram(  np.mod( np.hstack((lef_lefts,lef_rights)) , mapN ), np.arange(0,mapN,1))
    
        ratio = (np.sum(a[4998:5002])/((np.sum(a[5080:5084])+np.sum(a[4916:4920]))/2))
        
        file.write('%s,%s,%s,%s,%s,%s,%s,%s\n'%(birth,target_s_num, deltactcf,clife,cof,sep,face,ratio))

        hist_dict[name] = a/np.sum(a)
    
    except Exception as e:
        print('file %s should be process'%name)
        
file.close()
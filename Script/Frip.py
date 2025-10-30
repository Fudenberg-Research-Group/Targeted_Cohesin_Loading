import Target_cohesin_loading.utils as utils
import time
import glob
import time
import numpy as np
import ast
import pandas as pd
import warnings
import h5py 
import sys


path_dict = {}
directory='/scratch1/rahmanin/target_loading_cohesin/target/sims_chn/'
for fname  in glob.glob(directory+'folder*'):
    path_dict[fname.split('/sims_chn/')[1][:]]= fname
path_dict = dict(sorted(path_dict.items()))


# Defining parameters for Frip calculation
window_size = 1
numx,numy = 1,len(path_dict)
rep = 1 
mon = 1000
site = 10
base_time = 500
hist_dict = {}
c=1
frip_file = open('frip_target_respoints_w1_main_tau17_all_density_changing_chnexample.csv','w') #Aug 23 is for 100 sep, and Aug 24 is for 1000 sep
frip_file.write('birth,life,deltactcf,clife,cof,sep,face,frip\n')
for name in list(path_dict.keys())[:10]:
    print(name)
    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, deltactcf, pause, sep, site, mon, rep, step, vel = params

    c+=1

    mapN=mon*site
    lefs = h5py.File(path_dict[name]+'/LEFPositions.h5','r')["positions"]
    print(len(lefs))
    ctcfrightlist = np.array(h5py.File(path_dict[name]+'/LEFPositions.h5','r')['CTCF_sites_right'])
    ctcfleftlist = np.array(h5py.File(path_dict[name]+'/LEFPositions.h5','r')['CTCF_sites_left'])
    lst = np.array(list(ctcfrightlist) + list(ctcfleftlist))
    ### list of boundary elements on all replications
    lst_t = []
    for i in range(rep):
        lst_t += list(np.array(lst))
    print(lst_t)
    
    
    lef_lefts = lefs[base_time:,:,0].flatten()
    lef_rights = lefs[base_time:,:,1].flatten()
    lef_positions = np.hstack((lef_lefts,lef_rights))

    peak_monomers = utils.peak_positions(lst_t,window_sizes=np.arange(-window_size,(window_size)+1) )
    frip = utils.FRiP(mapN * rep, lef_positions, peak_monomers)
    frip_file.write('%s,%s,%s,%s,%s,%s,%s,%s\n'%(birth,life,deltactcf,clife,cof,sep,face,frip))

frip_file.close()



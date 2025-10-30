import glob # for making directory of files
import numpy as np
import ast
import os
import polychrom
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import polychrom.contactmaps
import cooltools
import cooltools.lib.plotting

from polykit.analysis import contact_maps as cms
import polykit 
import cooler
import warnings
import h5py 
import sys


path_dict = {}

directory = '/scratch1/rahmanin/target_loading_cohesin/target/target_mult_sites/sims/'
for name  in glob.glob(directory+'/folder_*'):
    path_dict[name.split('/sims/')[1][:]]= name
path_dict = dict(sorted(path_dict.items()))


monomer_per_replica = 1000

mapN = 1 * monomer_per_replica #number of monomer to 
total_monomers = 1000
mapstarts = (np.arange(0,1000 , monomer_per_replica))
min_time = 600
max_time = 300000000
skip_freq = 1
refresh = True
if refresh== True:
    map_dict_eq = {}
i=1
erfile = open('files_error.txt','w')

for name in path_dict.keys():
    newname = name + '.npz'
    if os.path.exists(f'./mapsnew/{newname}'): continue
    try: 
        URIs = polychrom.hdf5_format.list_URIs(path_dict[name])
        URIs_eq = np.array(URIs)[(np.array([int(i.split("::")[-1]) for i in URIs]) > min_time)&(np.array([int(i.split("::")[-1]) for i in URIs]) < max_time)][::skip_freq]
        mrc = polychrom.contactmaps.monomerResolutionContactMapSubchains(
            URIs_eq,
            mapstarts,
            mapN,
            cutoff=2.3, n=8)
        map_dict_eq[name] = mrc
        np.savez_compressed('./mapsnew/%s.npz' % (name), mrc) 
    except Exception as e:
        erfile.write('%s\n'%name)
        print(f"An error occurred with {name}: {e}")

erfile.close()





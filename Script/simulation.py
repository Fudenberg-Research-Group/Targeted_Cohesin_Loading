# This file represent the simulation for the general case of boosted loading at the center of lattice, with association rate, and dynamic barriers (CTCFs)
# ===============================
# Standard Library Imports
# ===============================
import glob
import ast
import os
import sys
import time
import warnings
import numpy as np

# ===============================
# Configuration File
# ===============================
from config import * # contains simulation parameters, file paths, and flags

# ===============================
# Scientific Libraries
# ===============================
import h5py
import numpy as np
import pandas as pd

# ===============================
# Project-Specific Imports
# ===============================
from Target_cohesin_loading.lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary
# import the functions here
import Target_cohesin_loading.funcs as funcs
from funcs import *

# ===============================
# Cooltools for Hi-C Analysis
# ===============================
import cooltools
import cooltools.lib.plotting

# ===============================
# Polychrom for MD Simulations
# ===============================
import polychrom
from polychrom import forces, forcekits, polymerutils
import polychrom.contactmaps
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.lib.extrusion import bondUpdater
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic


# Get filename from command-line argument (the format is provided in the following lines)
filename = sys.argv[-1] if len(sys.argv) > 1 else default_filename
default_filename_ = 'folder_face_1.0_back_0_Clife_191.9_Cof_17.0_life_66.0_slife_66.0_birth_0.1_base_0.0001_targetsnum_16_deltactcf_4600_pause_0.0_ipause_0.985_sep_74_site_10_monomer_1000_replica_1_steps_200.0_vel_1'

print(f'This is file name: {filename}')

# Extract parameters from filename (if needed for logic)
params = [ast.literal_eval(i) for i in filename.split('folder_')[1].split('_')[1::2]]
face, back, clife, cof, life, slife, birth, base_loading, targetsnum, deltactcf, pause, ipause, sep, site, monomer, replica, steps, vel = params

# Define only the main paramdict_CTCF
paramdict_CTCF = {
    'CTCF_facestall': [face, face],
    'CTCF_backstall': [back, back],
    'CTCF_lifetime': [clife, clife],
    'CTCF_offtime': [cof, cof],
    'LEF_lifetime': [life, life],
    'LEF_stalled_lifetime': [slife, slife],
    'LEF_birth': [base_loading, birth],
    'targetsnum': targetsnum,
    'deltactcf': deltactcf,
    'LEF_ipause': [ipause, ipause],
    'LEF_pause': [pause, pause],
    'LEF_separation': sep,
    'sites_per_monomer': site,
    'monomers_per_replica': monomer,
    'number_of_replica': replica,
    'steps': steps,
    'velocity_multiplier': vel
}

# Use updated function to generate file/folder name
file_name = funcs.paramdict_to_filename(paramdict_CTCF)
folder_name = '/simulations/folder_' + file_name.split('file_')[1]
folder = os.getcwd() + folder_name

# Create folder if it doesn’t exist
if os.path.exists(folder):
    print("Folder already exists.")
else:
    os.mkdir(folder)

# Compute total number of lattice sites and their types
monomers_per_replica = paramdict_CTCF['monomers_per_replica']
sites_per_monomer = paramdict_CTCF['sites_per_monomer']
sites_per_replica = monomers_per_replica*sites_per_monomer
monomer_types = np.zeros(monomers_per_replica, dtype=int)
site_types = np.repeat(monomer_types, sites_per_monomer)
number_of_replica = paramdict_CTCF['number_of_replica']
total_sites = number_of_replica * sites_per_replica

# Adjust LEF density and update paramdict
paramdict_CTCF['LEF_separation'] = funcs.adjust_LEF_density(paramdict_CTCF)
print(paramdict_CTCF['LEF_separation'])

# Determining the strong and base CTCF regions
typedict = {'strong_CTCF':1, 'base_CTCF':0}
target_site_num = targetsnum # number of loading sites on the lattice
loadsite = total_sites//2 # The location of loading sites on lattice
site_types[total_sites//2 - target_site_num//2: total_sites//2+(1+target_site_num)//2] = typedict['strong_CTCF']
site_types[:total_sites//2 - target_site_num//2] = site_types[total_sites//2+(1+target_site_num)//2:] = typedict['base_CTCF']

# CTCF boundary sites configuration 
CTCF_right_positions = np.array([loadsite+(deltactcf//2)+1])
CTCF_left_positions = np.array([loadsite-(deltactcf//2)])
# ==========================================================
########### 1d simulation parameters for lattice ###########
# ==========================================================
bins = np.linspace(0, TRAJECTORY_LENGTH, BLOCKSTEPS, dtype=int)
N = (paramdict_CTCF['monomers_per_replica']*paramdict_CTCF['number_of_replica'])
LEFNum = N // paramdict_CTCF['LEF_separation']
print('LEFNum is', LEFNum)                
translocator = make_translocator(LEFTranslocatorDynamicBoundary, 
                                 site_types,
                                 CTCF_left_positions,
                                 CTCF_right_positions, 
                                 **paramdict_CTCF)

hist = []
c =1    
with h5py.File(folder+"/LEFPositions.h5", mode='w') as myfile:
    dset = myfile.create_dataset("positions", 
                                 shape=(TRAJECTORY_LENGTH, LEFNum, 2), #edited
                                 dtype=np.int32, 
                                 compression="gzip")
     # creating data sets for boundary elements possible sites
    dset_ctcf_sites_right = myfile.create_dataset("CTCF_sites_right",
                                                 shape = (len(CTCF_right_positions)), 
                                                 compression = "gzip", 
                                                 data=CTCF_right_positions.copy())

    dset_ctcf_sites_left = myfile.create_dataset("CTCF_sites_left",
                                                shape = len(CTCF_left_positions), 
                                                compression="gzip",
                                                data=CTCF_left_positions.copy())

    # creating data sets for boundary elements positions
    dset_ctcf_positions_right = myfile.create_dataset("CTCF_positions_right",
                                      shape = (TRAJECTORY_LENGTH, len(CTCF_right_positions), 1), 
                                     compression = "gzip")
    dset_ctcf_positions_left = myfile.create_dataset("CTCF_positions_left",
                                     shape = (TRAJECTORY_LENGTH, len(CTCF_left_positions), 1), 
                                     compression = "gzip")
    
    translocator.steps(0)
    
    for st, end in zip(bins[:-1], bins[1:]):
        cur = []
        ctcf_right_cur= []
        ctcf_left_cur = []
        for i in range(st, end):
            translocator.step() 
            hist.append(translocator.LEFs.copy())
            cur.append(translocator.LEFs.copy())
            ctcf_positions_right = (translocator.stallProbRight)[CTCF_right_positions]*1
            ctcf_positions_left = (translocator.stallProbLeft)[CTCF_left_positions]*1
            
            ctcf_right_cur.append(ctcf_positions_right.reshape(len(ctcf_positions_right),1))
            ctcf_left_cur.append(ctcf_positions_left.reshape(len(ctcf_positions_left),1))
        cur = np.array(cur)
        ctcf_right_cur = np.array(ctcf_right_cur)
        ctcf_left_cur = np.array(ctcf_left_cur)
        dset[st:end] = cur
        dset_ctcf_positions_right[st:end] = ctcf_right_cur
        dset_ctcf_positions_left[st:end] = ctcf_left_cur
    myfile.attrs["N"] = N * paramdict_CTCF['sites_per_monomer']
    myfile.attrs["LEFNum"] = LEFNum

# ==========================================================
######## Molecular dynamics simulaiton ########
# ==========================================================
myfile = h5py.File(folder + "/LEFPositions.h5", mode='r')
sites_per_monomer = paramdict_CTCF['sites_per_monomer']
N = myfile.attrs["N"] // sites_per_monomer
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"][::sites_per_monomer]// sites_per_monomer
Nframes = LEFpositions.shape[0]

# Md simulation characteristics
box = (N / PARAMDICT_md['dens']) ** 0.33  # density = 0.1.

# initialize positions
data = grow_cubic(N, int(box) - 2)  # creates a compact conformation 

# assertions for easy managing code below 
assert (Nframes % PARAMDICT_md['restartSimulationEveryBlocks']) == 0 
assert (PARAMDICT_md['restartSimulationEveryBlocks'] % PARAMDICT_md['saveEveryBlocks']) == 0

savesPerSim = PARAMDICT_md['restartSimulationEveryBlocks'] // PARAMDICT_md['saveEveryBlocks']
simInitsTotal  = (Nframes) // PARAMDICT_md['restartSimulationEveryBlocks'] 

milker = polychrom.lib.extrusion.bondUpdater(LEFpositions)

reporter = HDF5Reporter(folder=folder, max_data_length=100, overwrite=True, blocks_only=False)

for iteration in range(simInitsTotal):
    # simulation parameters are defined below 
    a = Simulation(
            platform="cuda",
            integrator='langevin',  timestep=PARAMDICT_md['timesteps'], collision_rate=PARAMDICT_md['tmst'],
            error_tol=0.01,  
            GPU="0",
            N = len(data),
            reporters=[reporter],
            PBCbox=[box, box, box],
            precision="mixed")  # timestep not necessary for variableLangevin
    ############################## New code ##############################
    a.set_data(data)  # loads a polymer, puts a center of mass at zero

    a.add_force(
        forcekits.polymer_chains(
            a,
            chains=[(0, None, 0)],

                # By default the library assumes you have one polymer chain
                # If you want to make it a ring, or more than one chain, use self.setChains
                # self.setChains([(0,50,1),(50,None,0)]) will set a 50-monomer ring and a chain from monomer 50 to the end

            bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={
                'bondLength':1.0,
                'bondWiggleDistance':0.1, # Bond distance will fluctuate +- 0.05 on average
             },

            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                'k':PARAMDICT_md['stiff']
            },

            nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={
                'trunc':1.5, # this will let chains cross sometimes
                'radiusMult':1.05, # this is from old code
                #'trunc':10.0, # this will resolve chain crossings and will not let chain cross anymore
            },
            except_bonds=True,
    ))
    # ------------ initializing milker; adding bonds ---------
    kbond = a.kbondScalingFactor / (PARAMDICT_md['smcBondWiggleDist'] ** 2)
    bondDist = PARAMDICT_md['smcBondDist'] * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)

    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=PARAMDICT_md['restartSimulationEveryBlocks'])
    #print(milker.allBonds[0])
    for t,l in enumerate(milker.allBonds):
        for b in l:
            if (b[0] == 11296) or (b[1] == 11296):
                print(t,b)
    # If your simulation does not start, consider using energy minimization below
    if iteration==0:
        a.local_energy_minimization() 
    else:
        a._apply_forces()

    for i in range(PARAMDICT_md['restartSimulationEveryBlocks']):        
       # print("restart#",i)
        if i % PARAMDICT_md['saveEveryBlocks'] == (PARAMDICT_md['saveEveryBlocks'] - 1):  
            a.do_block(steps=paramdict_CTCF['steps'])
        else:
            a.integrator.step(paramdict_CTCF['steps'])  # do steps without getting the positions from the GPU (faster)
        if i < PARAMDICT_md['restartSimulationEveryBlocks'] - 1: 
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a

    reporter.blocks_only = True  # Write output hdf5-files only for blocks

    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()

myfile.close()


import numpy as np
import matplotlib.pylab as plt

from lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary
from funcs import *
import funcs
import cooltools
import cooltools.lib.plotting
#%matplotlib inline
import pandas as pd
import h5py 
import glob
import os
import h5py 
import time
import sys

import warnings

import ast
# for MD simulations

import polychrom
from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.lib.extrusion import  bondUpdater
import polychrom.contactmaps
import ast



filename = sys.argv[-1]

print('this is file name %s'%filename)


params = [ast.literal_eval(i) for i in filename.split('folder_')[1].split('_')[1::2]]
face, back, clife, cof, life, slife, birth, deltactcf, pause, sep, site, monomer, replica, steps, vel = params

paramdict_CTCF={
            'CTCF_facestall':[face, face],
            'CTCF_backstall':[0, 0],
            'CTCF_lifetime':[clife, clife],
            'CTCF_offtime':[cof, cof],
            'LEF_lifetime':[life, life],
            'LEF_stalled_lifetime':[slife, slife],
            'LEF_birth':[0.0001, birth],
            'deltactcf':deltactcf,
            'LEF_pause':[pause, pause],
            'LEF_separation':sep,
            'sites_per_monomer':10,
            'monomers_per_replica':1000,
            'number_of_replica':1,
            'steps':200,
            'velocity_multiplier':vel
            }
paramdict={
            'CTCF_facestall':[face],
            'CTCF_backstall':[0],
            'CTCF_lifetime':[clife],
            'CTCF_offtime':[cof],
            'LEF_lifetime':[life],
            'LEF_stalled_lifetime':[slife],
            'LEF_birth':[birth],
            'deltactcf':deltactcf,
            'LEF_pause':[pause],
            'LEF_separation':sep,
            'sites_per_monomer':10,
            'monomers_per_replica':1000,
            'number_of_replica':1,
            'steps':200,
            'velocity_multiplier':vel
            }

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'CTCF_number':'Cnum',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'deltactcf':'deltactcf',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }
file_name = funcs.paramdict_to_filename(paramdict, paramdict_keys)
folder_name = '/sims_rev_md/'+'folder_' + file_name.split('file_')[1]
folder = os.getcwd() + folder_name
if os.path.exists(folder):
    print("already exist")
else:
    os.mkdir(folder)

monomers_per_replica = paramdict_CTCF['monomers_per_replica']
sites_per_monomer = paramdict_CTCF['sites_per_monomer']
sites_per_replica = monomers_per_replica*sites_per_monomer
monomer_types = np.zeros(monomers_per_replica, dtype=int)
site_types = np.repeat(monomer_types, sites_per_monomer)
#print(len(site_types))

# Let's make some strong and weak CTCF regions
typedict = {'strong_CTCF':1, 'weak_CTCF':0}
site_types[5001] = typedict['strong_CTCF']
site_types[:5001] = site_types[5002:] = typedict['weak_CTCF']


# LEF/CTCF properties in type A monomers may be obtained from the paramdict as follows
LEF_lifetime = paramdict_CTCF['LEF_lifetime'][1]
LEF_velocity = paramdict_CTCF['velocity_multiplier']
CTCF_facestall = paramdict_CTCF['CTCF_facestall']
CTCF_offtime = paramdict_CTCF['CTCF_offtime']
#print(CTCF_offtime[typedict['strong_CTCF']], CTCF_offtime[typedict['weak_CTCF']])

# Create some CTCF boundary sites
CTCF_right_positions = np.array([5001+(deltactcf//2)+1])
CTCF_left_positions = np.array([5001-(deltactcf//2)])
CTCF_sites_right = np.array([5001+(deltactcf//2)+1])
CTCF_sites_left = np.array([5001-(deltactcf//2)])


#file = open('loopsize_between_convctcf_life_%s_vel_%s_birth_%s_deltactcf_%s.csv'%(life,vel,birth,deltactcf),'w')
#file.write('life,vel,clife,cof,loopmean,loopstd\n')

########### 1d simulation parameters for lattice ###########
Trajn = 3100000 # trajectory length in monomer 
trajectory_length = Trajn * paramdict_CTCF['sites_per_monomer'] #trajectory length in lattice land
print('trajectory length is %s'%trajectory_length)
pause_multiplier = 1/(1-pause)
trajectory_length = trajectory_length * pause_multiplier
num_dummy_steps = trajectory_length // 150 #dummy steps in lattice land
blocksteps = 5 
bins = np.linspace(0, trajectory_length, blocksteps, dtype=int)
N = (paramdict_CTCF['monomers_per_replica']*paramdict_CTCF['number_of_replica'])
LEFNum = N // paramdict_CTCF['LEF_separation']
                
translocator = make_translocator(LEFTranslocatorDynamicBoundary, 
                                 site_types,
                                 CTCF_left_positions,
                                 CTCF_right_positions, 
                                 **paramdict_CTCF)

hist = []
c =1    
with h5py.File(folder+"/LEFPositions.h5", mode='w') as myfile:
    dset = myfile.create_dataset("positions", 
                                 shape=(trajectory_length, LEFNum, 2), #edited
                                 dtype=np.int32, 
                                 compression="gzip")
     # creating data sets for boundary elements possible sites
    dset_ctcf_sites_right = myfile.create_dataset("CTCF_sites_right",
                                                 shape = (len(CTCF_sites_right)), 
                                                 compression = "gzip", 
                                                 data=CTCF_sites_right.copy())


    dset_ctcf_sites_left = myfile.create_dataset("CTCF_sites_left",
                                                shape = len(CTCF_sites_left), 
                                                compression="gzip",
                                                data=CTCF_sites_left.copy())

    # creating data sets for boundary elements positions
    dset_ctcf_positions_right = myfile.create_dataset("CTCF_positions_right",
                                      shape = (trajectory_length, len(CTCF_sites_right), 1), 
                                     #dtype = np.bool, 
                                     compression = "gzip")
    dset_ctcf_positions_left = myfile.create_dataset("CTCF_positions_left",
                                     shape = (trajectory_length, len(CTCF_sites_left), 1), 
                                     #dtype = np.bool, 
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
            ctcf_positions_right = (translocator.stallProbRight)[CTCF_sites_right]*1
            ctcf_positions_left = (translocator.stallProbLeft)[CTCF_sites_left]*1
            
            ctcf_right_cur.append(ctcf_positions_right.reshape(len(ctcf_positions_right),1))
            ctcf_left_cur.append(ctcf_positions_left.reshape(len(ctcf_positions_left),1))
        cur = np.array(cur)
        ctcf_right_cur = np.array(ctcf_right_cur)
        #print(st,end,ctcf_right_cur)
        #print(cur)
        #print(np.shape(dset[st:end]),np.shape(cur))
        ctcf_left_cur = np.array(ctcf_left_cur)
        dset[st:end] = cur
        #print(np.shape(ctcf_right_cur),np.shape(dset_ctcf_positions_right[st:end]))
        dset_ctcf_positions_right[st:end] = ctcf_right_cur
        dset_ctcf_positions_left[st:end] = ctcf_left_cur
    myfile.attrs["N"] = N * paramdict_CTCF['sites_per_monomer']
    myfile.attrs["LEFNum"] = LEFNum


lefs = np.array(hist)
#print(lefs)
min_time = 0
lef_lefts = lefs[min_time:,:,0].flatten()
lef_rights = lefs[min_time:,:,1].flatten()
lef_positions = np.hstack((lef_lefts,lef_rights))
loopmean = np.mean((lef_rights-lef_lefts))
loopstd = np.std((lef_rights-lef_lefts))

### Molecular dynamics simulaiton ###
myfile = h5py.File(folder + "/LEFPositions.h5", mode='r')
sites_per_monomer = paramdict['sites_per_monomer']
N = myfile.attrs["N"] // sites_per_monomer
print(N)
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"][::sites_per_monomer]// sites_per_monomer
Nframes = LEFpositions.shape[0]

# Md simulation characteristics
stiff = 1
dens = 0.2
box = (N / dens) ** 0.33  # density = 0.1.

smcStepsPerBlock = 1  # now doing 1 SMC step per block 
# initialize positions
data = grow_cubic(N, int(box) - 2)  # creates a compact conformation 
block = 0  # starting block 
steps= paramdict['steps'] #//velocity_multiplier



# new parameters because some things changed 
saveEveryBlocks = 10   # save every 10 blocks (saving every block is now too much almost)
restartSimulationEveryBlocks = 100

# parameters for smc bonds
smcBondWiggleDist = 0.2
smcBondDist = 0.5

# assertions for easy managing code below 
assert (Nframes % restartSimulationEveryBlocks) == 0 
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal  = (Nframes) // restartSimulationEveryBlocks 


tstp = 70 # timestep for integrator in fs
tmst = 0.01 # thermostat for integrator

milker = polychrom.lib.extrusion.bondUpdater(LEFpositions)

reporter = HDF5Reporter(folder=folder, max_data_length=100, overwrite=True, blocks_only=False)

for iteration in range(simInitsTotal):
    # simulation parameters are defined below 
    a = Simulation(
            platform="cuda",
            integrator='langevin',  timestep=tstp, collision_rate=tmst,
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
                'k':1.5
                # K is more or less arbitrary, k=4 corresponds to presistence length of 4,
                # k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff
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
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)

    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryBlocks)
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

    for i in range(restartSimulationEveryBlocks):        
       # print("restart#",i)
        if i % saveEveryBlocks == (saveEveryBlocks - 1):  
            a.do_block(steps=steps)
        else:
            a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)
        if i < restartSimulationEveryBlocks - 1: 
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a

    reporter.blocks_only = True  # Write output hdf5-files only for blocks

    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()

myfile.close()

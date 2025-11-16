import h5py 
import time
import numpy as np
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
from polykit.analysis import contact_maps as cms
import cooler


def perform_md_simulation(lef_file_path, paramdict, paramdict_md):

    # Loading 1d paramters
    sites_per_monomer = paramdict['sites_per_monomer']
    lef_file = h5py.File(lef_file_path + "/LEFPositions.h5", mode='r')
    LEFpositions = lef_file["positions"][::sites_per_monomer]// sites_per_monomer
    number_of_monomers = lef_file.attrs["N"]//sites_per_monomer
    LEFNum = lef_file.attrs["LEFNum"]
    Nframes = LEFpositions.shape[0]
    
    # Md simulation characteristics
    stiff = paramdict_md['stiff']
    dens = paramdict_md['dens']
    box = (number_of_monomers / dens) ** 0.33 
    smcStepsPerBlock = 1  # now doing 1 SMC step per block 
    saveEveryBlocks = paramdict_md['saveEveryBlocks']   # save every 10 blocks (saving every block is now too much almost)
    restartSimulationEveryBlocks = paramdict_md['restartSimulationEveryBlocks']
    # assertions for easy managing code below 
    assert (Nframes % restartSimulationEveryBlocks) == 0 
    assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0
    
    # initialize positions
    data = grow_cubic(number_of_monomers, int(box) - 2)  # creates a compact conformation 
    steps= paramdict_md['steps'] # number of md steps between 1d updates
      
    # parameters for smc bonds
    smcBondWiggleDist = 0.2
    smcBondDist = 0.5
    
    savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
    simInitsTotal  = (Nframes) // restartSimulationEveryBlocks 
    
    tstp = 70 # timestep for integrator in fs
    tmst = 0.01 # thermostat for integrator
    
    milker = polychrom.lib.extrusion.bondUpdater(LEFpositions)
    
    reporter = HDF5Reporter(folder=lef_file_path, max_data_length=100, overwrite=True, blocks_only=False)
    
    for iteration in range(simInitsTotal):
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
    
                bond_force_func=forces.harmonic_bonds,
                bond_force_kwargs={
                    'bondLength':1.0,
                    'bondWiggleDistance':0.1, # Bond distance will fluctuate +- 0.05 on average
                 },
    
                angle_force_func=forces.angle_force,
                angle_force_kwargs={
                    'k':stiff
                },
    
                nonbonded_force_func=forces.polynomial_repulsive,
                nonbonded_force_kwargs={
                    'trunc':1.5, # this will let chains cross sometimes
                    'radiusMult':1.05, # this is from old code
                },
                except_bonds=True,
        ))
        # ------------ initializing milker; adding bonds ---------
        kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
        bondDist = smcBondDist * a.length_scale
    
        activeParams = {"length":bondDist,"k":kbond}
        inactiveParams = {"length":bondDist, "k":0}
        milker.setParams(activeParams, inactiveParams)
    
        # this step actually puts all bonds in and sets first bonds to be what they should be
        milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                    blocks=restartSimulationEveryBlocks)
        if iteration==0:
            a.local_energy_minimization() 
        else:
            a._apply_forces()
    
        for i in range(restartSimulationEveryBlocks):        
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


def make_md_contact_maps(mapN, mapstarts, freq, min_time, input_path, cool_uri, coarsen=True):
    URIs = polychrom.hdf5_format.list_URIs(input_path)
    URIs_eq = np.array(URIs)[np.array([int(i.split("::")[-1]) for i in URIs]) > min_time][::freq]
    cont_matrix = polychrom.contactmaps.monomerResolutionContactMapSubchains(
        URIs_eq,
        mapstarts,
        mapN,
        cutoff=2.3,
        n=8)
    
    #saving data as .cool files
    chunksize=10000000
    mrc_new = cms.coolify(cont_matrix,
            cool_uri,
            chrom_dict={},
            binsize=2500,
            chunksize=chunksize)
    
    clr = cooler.Cooler(cool_uri+'.2500.cool')
    base_uri = cool_uri+'.2500.cool'
    output_uri = cool_uri+'.10000.cool'
    if coarsen:
        factor = 4
        clr_10 = cooler.coarsen_cooler(base_uri,output_uri,factor,chunksize)
# Genomic region
REGION = 'chr1:30000000-37500000'

# File paths
import os

# Read values from environment variables or use defaults
REGION = os.getenv("REGION", "chr1:30000000-37500000")
RUN_INDEX = os.getenv("RUN_INDEX", "1")

# Paths
PARAMDICT_PATH = 'data/paramdict.json'
OUTPUT_PATH = f'simulations/sims_random_mod_run{RUN_INDEX}/'
MAP_OUTPUT_DIRECTORY = './data/contact_maps_mod/'
OUTPUT_FILE_READS = f'Rad21_vs_ctcf_reads_random_mod_a_w_{REGION}_run{RUN_INDEX}.csv'

# Ensure output directory exists
os.makedirs(OUTPUT_PATH, exist_ok=True)
os.makedirs(MAP_OUTPUT_DIRECTORY, exist_ok=True)


# Load parameters from JSON
import json
with open(PARAMDICT_PATH, 'r') as json_file:
    PARAMDICT = json.load(json_file)

# Experimental paramters
import OccupancyInputCTCF.utils.snippet as snp
import OccupancyInputCTCF.utils.experimental_path as exp
input_occupancy = exp.ctcf_peaks  # Input occupancy file
ctcf_bed_df = snp.get_region_snippet(exp.ctcf_peaks, exp.ctcf_motifs, REGION)

# Simulation parameters
TRAJECTORY_LENGTH = 3150
N_SIMULATIONS = 5
LATTICE_SITE = 250
WINDOW_SIZE = 0#1

# Contact map parameters
START_FRAMES = [100, 200, 300, 400, 500]
END_FRAME_OFFSET = 800
EVERY_FRAME = 500
MAX_DIST = 300
RES_CONVERT = 40
REPLICATION_NUMBER = PARAMDICT['number_of_replica']

# Flags
USE_PREDICTED_OCCUPANCY = True
RUN_MD_SIM = False
calculate_frip = False

# Md simulation parameters
PARAMDICT_md = {}
PARAMDICT_md['stiff']=1.5
PARAMDICT_md['dens']=0.2
PARAMDICT_md['saveEveryBlocks']=5
PARAMDICT_md['restartSimulationEveryBlocks']=100
PARAMDICT_md['steps'] = 20
import sys
import ast
# Get filename from command-line argument (the format is provided in the following lines)
filename = sys.argv[-1] if len(sys.argv) > 1 else default_filename
default_filename_ = 'folder_face_1.0_back_0.0_Clife_300_Cof_180_life_660.0_slife_660.0_birth_0.1_pause_0.9_sep_25_site_10_replica_10_vel_1_steps_20'
#print('this is file name %s'%filename)

params=[ast.literal_eval(i) for i in filename.split('folder_')[1].split('_')[1::2]]
face, back,clife, cof, life, slife, birth, pause, sep, site, rep, vel, steps = params

PARAMDICT['LEF_lifetime']=[life]
PARAMDICT['LEF_stalled_lifetime']=[life]
PARAMDICT['LEF_separation']=sep
PARAMDICT['velocity_multiplier']=vel

#print(PARAMDICT)

# Parameters for making maps from md simulations
import OccupancyInputCTCF.utils.convert as convert
import numpy as np
monomer_per_replica = convert.get_lattice_size(REGION, LATTICE_SITE)//10
mapN = 3 * monomer_per_replica #number of monomer to 
total_monomers = monomer_per_replica * REPLICATION_NUMBER
freq = 1 #frequent frames
mapstarts = (np.arange(0,total_monomers-2*monomer_per_replica , monomer_per_replica))
min_time = 200 # the number of steps to disregard when calculating contacts (for equilibration purpose)




import sys
import os
import json
from config import *
import pandas as pd
import h5py
from scipy.stats import truncnorm
import matplotlib.pylab as plt
import numpy as np
from cooltools.lib.numutils import adaptive_coarsegrain
import multiprocessing as mp
from functools import partial
# Add custom library path
sys.path.append('/home1/rahmanin/start/polychrom/projects/Site_wise_occupancy/OccupancyInputCTCF/')

# Import utility modules Targeted_cohesin_loading
import Targeted_cohesin_loading.plots as mplot
import Targeted_cohesin_loading.convert as convert
import Targeted_cohesin_loading.makeparams as params
import Targeted_cohesin_loading.One_d_simulation as simulation
import Targeted_cohesin_loading.md_simulation as mdsimulation
import Targeted_cohesin_loading.cmap_utils as utils_s
import warnings
warnings.filterwarnings('ignore')
import time
start = time.time()

# =================== WORKFLOW ===========================
print("Starting workflow...")

# Step 1: Preprocessing bed files in selected region
print("Step 1: Preprocessing...")
ctcf_bed_file = f'processing/{REGION}.csv'
ctcf_bed_df.to_csv(ctcf_bed_file, index=False)
print("Step 1 complete. Output:", ctcf_bed_df)

# Step 2: Predicting Occupancy Rate
#print("Step 2: Predicting Occupancy Rate...")
#if USE_PREDICTED_OCCUPANCY:
#    ctcf_occup_bed_df = ml.predict_ctcf_occupancy(ctcf_bed_file)
#else:
#    ctcf_occup_bed_df = convert.convert_ctcf_occupancy(ctcf_bed_df, ctcf_peaks)
#print("Step 2 complete. Output:", ctcf_occup_bed_df)

# Step 3: Refining Occupancy for Overlapping Sites
#print("Step 3: Refining Occupancy...")
#refined_occupancy = convert.get_refined_occupancy(ctcf_occup_bed_df, REGION)
#print("Step 3 complete. Output:", refined_occupancy)

# Step 4: Generating Barrier List with Occupancy and Lifetimes
#print("Step 4: Generating Barrier List...")
#CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list = convert.get_ctcf_list(
#    refined_occupancy, PARAMDICT, insert_on='bound_time')
#print("Step 4 complete. Right barriers:", CTCF_right_positions, "Left barriers:", CTCF_left_positions)

##### adding the random implementation ######
genome_length = 3*2.7 # Gbp , the total genome length with averagely 3 copies
tot_motifs = 55073 # Total number of CTCF motifs from JASPAR
occup_mean = 0.33182353 # Average of CTCF site occupancy from Sonmezer data
ctcf_tot_number = 217000 # Paper cited in manuscript, method section
bound_fraction = 0.5 # Paper cited in manuscript, method section
occup_mean_p = (ctcf_tot_number/3)*bound_fraction/tot_motifs
multiplication = occup_mean_p/occup_mean
ctcf_per_10_meg = (10*1e6 / (genome_length *1e9) * bound_fraction * ctcf_tot_number) / occup_mean
ctcf_per_10_meg

bind_sites_filepath = 'tutorials/data/sonmezer_dataset_CTCT_binding.sites.filtered.mm10.tsv'
bind_freqs_filepath = 'tutorials/data/binding.frequencies.tsv'
ref_genome_filepath = '/project/fudenber_735/genomes/mm10/mm10.fa'


#ref_genome = pysam.FastaFile(ref_genome_filepath)
bind_sites = pd.read_table(bind_sites_filepath)
bind_freqs = pd.read_table(bind_freqs_filepath)
df = pd.read_csv(bind_freqs_filepath, sep = '\t')
ctcf_occup_exp = df[df['biological.state']=='Bound']


# filtering sites with provided all three frequencies
site_filter = []
for tfbs in bind_freqs['TFBS_cluster'].unique():
       if len(bind_freqs[bind_freqs['TFBS_cluster'] == tfbs]) == 3:
            site_filter.append(tfbs)              
bind_sites = bind_sites[bind_sites.rownames.isin(site_filter)]
bind_sites_revised= bind_sites.rename(columns={'rownames':'TFBS_cluster'})
bind_freqs_dtf = bind_freqs[bind_freqs.TFBS_cluster.isin(site_filter)]
bind_freqs_pivot = bind_freqs_dtf.pivot_table(index='TFBS_cluster', columns='biological.state', values='Freqs', aggfunc='first')
bind_freqs_pivot.reset_index(inplace=True)
bind_freqs_pivot = bind_freqs_pivot[['TFBS_cluster','Accessible','Bound', 'Nucleosome.occupied']]
states = bind_freqs['biological.state'].unique()
#bind_freqs_pivot
binds_sites_freqs = pd.merge(bind_sites_revised, bind_freqs_pivot, on='TFBS_cluster', how='left')
binds_sites_freqs.to_csv('sites_with_freqs.tsv',sep='\t')
binds_sites_freqs
### making labels
frequency_columns = ['Accessible', 'Bound', 'Nucleosome.occupied']

# Convert the selected columns to a NumPy array
experiment_occup = binds_sites_freqs[['chrom','start','end','Bound']]
experiment_occup['true_occupancy']=experiment_occup['Bound']
experiment_occup
experiment_occup_mod = experiment_occup.query('true_occupancy >=0.05')
ctcf_occup_rand = np.random.choice( (experiment_occup_mod['true_occupancy'].values), (3000,), replace=True)


ctcf_random_list_num = np.random.choice(np.arange(0, 30000), size=5000, replace=False)
CTCF_right_positions_list = []
CTCF_left_positions_list = []
CTCF_occup_list = []
occup_sum = 0
i = 0 
while occup_sum < ctcf_per_10_meg/(3/2):
    CTCF_right_positions_list.append(ctcf_random_list_num[i])
    CTCF_occup_list.append(ctcf_occup_rand[i])
    occup_sum += ctcf_occup_rand[i]
    i += 1
    CTCF_left_positions_list.append(ctcf_random_list_num[i])
    CTCF_occup_list.append(ctcf_occup_rand[i])
    occup_sum += ctcf_occup_rand[i]
    i += 1
    
ctcf_loc_list = CTCF_right_positions_list + CTCF_left_positions_list
CTCF_right_positions = np.array(CTCF_right_positions_list)
CTCF_left_positions = np.array(CTCF_left_positions_list) 


### adjusting new bound mean value
occup_mean_rand = np.mean(CTCF_occup_list)
occup_mean_rand
occup_exp = 0.65
unbound_mean = 330*(1-occup_exp)/occup_exp
multiplication = occup_mean_rand/occup_exp
bound_mean = occup_mean_rand*unbound_mean/(1-occup_mean_rand)
print(bound_mean,unbound_mean)


CTCF_occup_list = np.array(CTCF_occup_list)
ctcf_lifetime_list_ = CTCF_occup_list * unbound_mean / (1 - CTCF_occup_list)#  + pseudo) for scenario with ctcf occupancy up to 1
s = []

for i in range(len(ctcf_lifetime_list_)): 
    mean = ctcf_lifetime_list_[i]
    std_dev = ctcf_lifetime_list_[i]/3
    upper = mean
    a = (3 - mean) / std_dev
    b =  upper / std_dev
    random_value = truncnorm.rvs(a, b, loc=mean, scale=std_dev)
    s.append(random_value)
ctcf_lifetime_list = ctcf_lifetime_list_ + np.array(s)
ctcf_offtime_list = ctcf_lifetime_list*((1-CTCF_occup_list)/CTCF_occup_list)







# Step 5: Running 1D Simulation
print("Step 5: Running 1D Simulation...")
file_name = params.paramdict_to_filename(PARAMDICT)
output_directory = OUTPUT_PATH + 'folder_' + file_name.split('file_')[1]
os.makedirs(output_directory, exist_ok=True)

PARAMDICT['monomers_per_replica'] = convert.get_lattice_size(REGION, lattice_site=250)//10
#ctcf_params = convert.get_ctcf_list(refined_occupancy, PARAMDICT)
ctcf_params = CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list

output_dirs = []
for sim_id in range(1, N_SIMULATIONS+1):
    file_name = f"simulation_{sim_id}"
    output_directory_partial = os.path.join(output_directory, f"{file_name}")
    os.makedirs(output_directory_partial, exist_ok=True)
    output_dirs.append(output_directory_partial)
    
def Perform_1d_simulation(output_directory, PARAMDICT, ctcf_params, TRAJECTORY_LENGTH):
    return simulation.Perform_1d_simulation(PARAMDICT, ctcf_params, TRAJECTORY_LENGTH, output_directory)
    
# use partial to set up some parameters of the function and leave output dir later
partially_set_Perform_1d_simulation = partial(Perform_1d_simulation, PARAMDICT=PARAMDICT, ctcf_params=ctcf_params, TRAJECTORY_LENGTH=TRAJECTORY_LENGTH)

def worker_init_fn():
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
with mp.Pool(processes=mp.cpu_count(), initializer=worker_init_fn) as pool:
    pool.map(
        partially_set_Perform_1d_simulation, output_dirs
    )
    
print("Step 5 complete. Output at", OUTPUT_PATH)

# Optional Step 6: MD Simulation
if RUN_MD_SIM:
    print("Step 6: Running optional MD simulation...")
    for sim_id in range(1, N_SIMULATIONS+1):
        file_name = f"simulation_{sim_id}"
        LEF_FILE_PATH = os.path.join(output_directory, f"{file_name}")
        mdsimulation.perform_md_simulation(LEF_FILE_PATH, PARAMDICT, PARAMDICT_md)
        cool_uri = output_directory + '/simulation_%s/map.mcool'%(sim_id)
        mdsimulation.make_md_contact_maps(mapN, mapstarts, freq, min_time, LEF_FILE_PATH, cool_uri)
    print("Step 6 complete.")
    
### extracting maps from each simulaiton
# Step 7: Processing Outputs
print("Step 7: Processing Outputs...")
MAP_OUTPUT_DIRECTORY = output_directory + '/contact_maps'
os.makedirs(MAP_OUTPUT_DIRECTORY, exist_ok=True)
mapN = PARAMDICT['monomers_per_replica']*PARAMDICT['sites_per_monomer']
total_sites = mapN * PARAMDICT['number_of_replica']
lef_array, wmap, wchip, wchip_ctcf = [], [], [], []

for sim_id in range(1, N_SIMULATIONS+1):
    file_name = f"simulation_{sim_id}"
    output_directory_partial = os.path.join(output_directory, f"{file_name}")
    lef_file_path = os.path.join(output_directory_partial, 'LEFPositions.h5')
    if os.path.exists(lef_file_path):
        print(f"LEFPositions.h5 file found at {lef_file_path}.")
        lefs = h5py.File(lef_file_path, 'r')['positions']
        lef_array.append(lefs)
        print('length of lefs is %s'%len(np.array(lefs)))
        chip = utils_s.chip_seq_from_lef(lefs, mapN)
        cmap = utils_s.contact_map_from_lefs(lefs[:], mapN)
        ctcfchip = utils_s.chip_seq_from_ctcf(lef_file_path,mapN)
        wmap.append(cmap)
        wchip.append(chip)
        wchip_ctcf.append(np.array(ctcfchip))
    else:
        print(f"Error: LEFPositions.h5 file not found at {lef_file_path}")
lefs_array = np.vstack(lef_array)

# Multiprocessing for contact maps from graph based distance
map_output_dirs = utils_s.create_contact_map_folders(N_SIMULATIONS, MAP_OUTPUT_DIRECTORY+'/')
print(map_output_dirs)

end_frame = len(lefs_array)-END_FRAME_OFFSET
def worker(map_output_dir, START_FRAMES):
    utils_s.calculate_contact_maps(total_sites, lefs_array, START_FRAMES, end_frame, EVERY_FRAME, MAX_DIST, RES_CONVERT, REPLICATION_NUMBER, map_output_dir)

with mp.Pool(processes=mp.cpu_count()) as pool:
    pool.starmap(worker, zip(map_output_dirs, START_FRAMES))

w_map = []
for dirs in map_output_dirs:
    file_path = os.path.join(dirs, 'contact_map.npz')
    with np.load(file_path) as data:
        cmap = data['contact_map']  
        w_map.append(cmap)
whole_map = np.sum(w_map, axis=0)
whole_chip = np.sum(wchip, axis=0)
whole_chip_ctcf = np.sum(wchip_ctcf, axis=0)    
chip = utils_s.chip_seq_from_lef(lefs, mapN)
cmap = whole_map
print("Step 7 complete. Output at", OUTPUT_PATH)


# step 8: plots 
print("Step8: making experimental and simulated maps ...")
output_file = 'outputs/plots_%s_a.pdf'%(REGION)
    
mplot.plot_chip_hic(REGION, whole_chip, whole_chip_ctcf, whole_map, res= 200000, output_file=output_file)
print("Step 8 complete. Output at", output_file)


def Calculate_EOC_reads(paramdict, lefs_array, lst, window_size):
    rep = paramdict['number_of_replica'] 
    mon = paramdict['monomers_per_replica']
    site = paramdict['sites_per_monomer']
    mapN = paramdict['monomers_per_replica']*paramdict['sites_per_monomer']
    lst_t = []
    for i in range(rep):
        lst_t += list(np.array(lst)+i*mon*site)
    lef_lefts = lefs_array[:,:,0].flatten()
    lef_rights = lefs_array[:,:,1].flatten()
    lef_positions = np.hstack((lef_lefts,lef_rights))
    peak_monomers = utils_s.peak_positions(lst_t, window_sizes = np.arange(-window_size, (window_size)+1))
    num_sites_t = mapN*rep
    hist,edges = np.histogram(  lef_positions  , np.arange(num_sites_t+1) )
    reads = np.sum(hist[peak_monomers])
    return reads

ctcf_occup = ctcf_lifetime_list/ (ctcf_lifetime_list+ctcf_offtime_list)
traj = lefs_array.shape[0]
ctcf_reads = ctcf_occup*traj
ctcf_reads


file = open(OUTPUT_FILE_READS,'w')
file.write('Rad21,ctcf\n')
for i in range(len(ctcf_loc_list)):
    lst_e = [ctcf_loc_list[i]]
    reads = Calculate_EOC_reads(PARAMDICT, lefs_array, lst_e, WINDOW_SIZE)
    reads_per_ctcf = reads/ctcf_reads[i]
    file.write('%s,%s\n'%(reads,ctcf_reads[i]))
file.close()











import chromoscores.maputils as chrmap
import chromoscores.snipping as chrsnip
import chromoscores.scorefunctions as chrscores
from scipy.stats import pearsonr
import numpy as np

def map_from_lefs(dset, sites_per_replica):
    
    ll = np.mod(dset.reshape((-1, 2)), sites_per_replica)
    ll = ll[ll[:,1] > ll[:,0]]
    
    lmap = np.histogram2d(ll[:,0], ll[:,1], np.arange(sites_per_replica))[0]
    
    return (lmap + lmap.T)

def FRiP(num_sites_t, lef_positions, peak_positions ):
    
    hist,edges = np.histogram(  lef_positions  , np.arange(num_sites_t+1) )
    return np.sum(hist[peak_positions] )/len(lef_positions)
    
def peak_positions(boundary_list, window_sizes=[1]):
    """
    Calculate peak positions based on a boundary_list within window_sizes.

    Args:
        boundary_list (list): List of boundary values.
        window_sizes (list, optional): List of window sizes. Defaults to [1].

    Returns:
        np.ndarray: Array containing peak positions.
    """
    peak_monomers = np.array([])

    for i in window_sizes:
        inds_to_add = [boundary + i for boundary in boundary_list]
        peak_monomers = np.hstack((peak_monomers, inds_to_add))

    return peak_monomers.astype(int)
def per_k_multiplier(multiplier, num_site=1, tot_site_per_k=4):
    per_k_multiplier_conv = (1*multiplier*num_site+(tot_site_per_k-1))/tot_site_per_k
    return per_k_multiplier_conv

    
def create_matrix(n, higher_value, lower_value):
    # Create an n x n matrix filled with the lower value
    matrix = np.full((n, n), lower_value)
    
    for i in range(n):
        matrix[i, n - i - 1] = higher_value
    
    return matrix

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
















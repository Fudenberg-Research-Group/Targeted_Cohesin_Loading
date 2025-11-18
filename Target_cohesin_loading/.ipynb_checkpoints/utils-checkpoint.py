import chromoscores.maputils as chrmap
import chromoscores.snipping as chrsnip
import chromoscores.scorefunctions as chrscores
from scipy.stats import pearsonr
import numpy as np

def map_from_lefs(dset, sites_per_replica):
    """
    Construct a symmetric 2D contact map from LEF positions.

    This function takes an array of LEF head pairs and generates a contact
    matrix representing how frequently left-right LEF heads span each pair
    of lattice sites. It reshapes the input, wraps positions using modulo
    (to account for replica structure), filters valid pairs, and produces
    a histogram-based contact map.

    Parameters
    ----------
    dset : array-like
        Array containing LEF head coordinates. It may be flattened; the
        function reshapes it into shape (N, 2), where each row is a
        (left_head, right_head) pair.
    sites_per_replica : int
        Total number of lattice sites in a single replica. Used for modulo
        wrapping and histogram binning.

    Returns
    -------
    contact_map : ndarray, shape (sites_per_replica, sites_per_replica)
        Symmetric array where entry (i, j) counts how many LEFs span the
        sites i and j. The result is lmap + lmap.T to ensure symmetry.

    Notes
    -----
    - Pairs with right < left are removed.
    - Histogram binning automatically aggregates repeated interactions.
    """
    
    ll = np.mod(dset.reshape((-1, 2)), sites_per_replica)
    ll = ll[ll[:,1] > ll[:,0]]
    
    lmap = np.histogram2d(ll[:,0], ll[:,1], np.arange(sites_per_replica))[0]
    
    return (lmap + lmap.T)


def FRiP(num_sites_t, lef_positions, peak_positions):
    """
    Compute the Fraction of Reads in Peaks (FRiP) from LEF head positions.

    FRiP measures how many LEF positions fall inside a defined peak region.
    It is the ratio between the number of LEF heads that lie within peak
    indices and the total number of observed LEF head positions.

    Parameters
    ----------
    num_sites_t : int
        Total number of lattice sites across all replicas.
        Used to define histogram bin edges.
    lef_positions : array-like
        Flattened list of left and right LEF head positions (indices on the lattice).
    peak_positions : array-like
        Indices representing the window or peak region around a site
        (e.g., ±window_size around CTCF). FRiP measures enrichment here.

    Returns
    -------
    frip_value : float
        Fraction of all LEF positions that fall inside the peak region.
        Computed as sum(hist[peak_positions]) / len(lef_positions).

    Notes
    -----
    - Uses histogram binning for efficient counting.
    - Higher FRiP indicates more LEF accumulation at peak regions.
    """
    
    hist,edges = np.histogram(lef_positions, np.arange(num_sites_t+1))
    return np.sum(hist[peak_positions]) / len(lef_positions)

    
def peak_positions(boundary_list, window_sizes=[1]):
    """
    Generate lattice positions around specified boundaries using a window range.

    This function expands each boundary position into a set of nearby indices,
    defined by the values in `window_sizes`. For every boundary in the list, the
    function adds (boundary + w) for each w in window_sizes, producing a combined
    list of peak/enrichment positions.

    Parameters
    ----------
    boundary_list : list or array-like
        List of base boundary positions (e.g., CTCF site locations).
    window_sizes : list or array-like, optional
        Offsets added to each boundary position. These typically represent
        distance from the boundary (e.g., [-window, ..., +window]).
        Defaults to [1].

    Returns
    -------
    np.ndarray
        Array of integer positions representing all expanded peak positions.

    Notes
    -----
    - No bounds checking is performed; caller must ensure valid positions.
    - Useful for defining local windows around CTCF sites when counting LEF reads.
    """

    peak_monomers = np.array([])

    for i in window_sizes:
        inds_to_add = [boundary + i for boundary in boundary_list]
        peak_monomers = np.hstack((peak_monomers, inds_to_add))

    return peak_monomers.astype(int)

    
def per_k_multiplier(multiplier, num_site=1, tot_site_per_k=4):
    """
    Compute an averaged multiplier value across a block of size K.

    This function converts a per-site multiplier into a per-K-site multiplier
    by distributing the effect across `tot_site_per_k` sites. Only `num_site`
    sites carry the multiplier, and the remaining sites contribute the default value 1.

    Parameters
    ----------
    multiplier : float
        The multiplier applied to the active sites.
    num_site : int, optional
        Number of sites within the K-block to which the multiplier applies.
        Defaults to 1.
    tot_site_per_k : int, optional
        The total number of sites considered in the K-block.
        Defaults to 4.

    Returns
    -------
    float
        The averaged effective multiplier across the entire K-block.
    """
    per_k_multiplier_conv = (1 * multiplier * num_site + (tot_site_per_k - 1)) / tot_site_per_k
    return per_k_multiplier_conv



def create_matrix(n, higher_value, lower_value):
    """
    Create an n×n matrix with a constant background and a highlighted diagonal.

    Builds a square matrix filled with `lower_value`, then sets the
    anti-diagonal (top-right to bottom-left) entries to `higher_value`.

    Parameters
    ----------
    n : int
        Size of the matrix (n rows × n columns).
    higher_value : float or int
        Value placed on the anti-diagonal.
    lower_value : float or int
        Value used to fill the rest of the matrix.

    Returns
    -------
    numpy.ndarray
        An n×n matrix with the specified diagonal pattern.
    """
    matrix = np.full((n, n), lower_value)
    
    for i in range(n):
        matrix[i, n - i - 1] = higher_value
    
    return matrix


def Calculate_EOC_Reads(paramdict, lefs_array, lst, window_size):
        """
    Computes the total number of LEF encounters (reads) around specified CTCF sites (Extruders On CTCF).

    This function takes a list of CTCF positions for a single replica, expands them
    across all replicas, and counts how many LEF positions (left or right heads)
    fall within a given window around those sites.

    Parameters
    ----------
    paramdict : dict
        Dictionary containing simulation parameters, including:
            'number_of_replica'      – number of lattice replicas
            'monomers_per_replica'   – monomers per replica
            'sites_per_monomer'      – sites per monomer
    lefs_array : ndarray, shape (replica, number_of_LEFs, 2)
        Array containing LEF head positions.
        lefs_array[:,:,0] = left heads
        lefs_array[:,:,1] = right heads
    lst : list or array
        Positions of CTCF sites in a single replica.
    window_size : int
        Number of lattice sites to include on each side of the CTCF site
        when counting LEF positions.

    Returns
    -------
    reads : int
        Total number of LEF heads located within ±window_size
        of the CTCF sites across all replicas.

    Notes
    -----
    - CTCF positions listed in `lst` are expanded across all replicas.
    - The function counts *both* left and right LEF heads.
    - Uses a histogram over the full lattice to efficiently count encounters.
    """
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
    peak_monomers = peak_positions(lst_t, window_sizes = np.arange(-window_size, (window_size)+1))
    num_sites_t = mapN*rep
    hist,edges = np.histogram(  lef_positions  , np.arange(num_sites_t+1) )
    reads = np.sum(hist[peak_monomers])
    return reads
















import chromoscores.maputils as chrmap
import chromoscores.snipping as chrsnip
import chromoscores.scorefunctions as chrscores
from scipy.stats import pearsonr
import numpy as np

def assign_lattice_positions(dataframe, region, lattice_size=250):
    """
    Extract and annotate genomic entries within a specified region.

    This function selects rows from a genomic dataframe that fall inside a
    region string (e.g., 'chr1:100000-200000'), computes the midpoint of each
    interval, and assigns each entry to a lattice/bin index based on a chosen
    lattice size. It is useful for mapping genomic coordinates onto simulation
    lattice coordinates.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Genomic table containing at least 'chrom', 'start', and 'end'
        columns. Additional columns are preserved but not modified.
    region : str
        Genomic region in the format 'chrom:start-end', such as
        'chr1:150000-250000'.
    lattice_size : int, optional
        Bin size used to convert physical genomic coordinates into
        lattice indices. Defaults to 250.

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe corresponding to the selected region, with two
        additional columns:
            'mid'         – midpoint of each genomic interval
            'lattice_loc' – integer lattice index computed as
                            floor((mid - region_start) / lattice_size)

    Notes
    -----
    - Uses bioframe.parse_region_string and bioframe.select for region parsing.
    - Returned indices are zero-based lattice coordinates.
    """
    region_start = bioframe.parse_region_string(region)[1]
    region_dataframe = bioframe.select(dataframe, region, cols=['chrom', 'start', 'end'])
    region_dataframe['mid'] = (region_dataframe.end + region_dataframe.start) / 2
    region_dataframe['lattice_loc'] = ((region_dataframe['mid'] - region_start) // lattice_size).astype(int)
    region_dataframe = region_dataframe.reset_index(drop=True)
    return region_dataframe


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
    
def make_region_occupancy(file):
    """
    Compute per-lattice occupancy for a genomic region from prediction output.

    This function loads a CSV file containing predicted occupancy values and
    aggregates them by (lattice_loc, strand). For each lattice position, the
    combined occupancy is computed using the probability that at least one
    footprint occurs:
        combined = 1 − Π(1 − predicted_occupancy_i)

    The function returns a tidy table with one row per lattice position and
    strand.

    Parameters
    ----------
    file : str or path-like
        Path to a CSV file containing columns:
        ['chrom', 'start', 'end', 'mid', 'strand', 'lattice_loc',
         'predicted_occupancy'].

    Returns
    -------
    pandas.DataFrame
        A dataframe with aggregated occupancy and key genomic coordinates.
        Columns include:
            'chrom', 'start', 'end', 'mid',
            'strand', 'lattice_loc', 'predicted_occupancy'.

    Notes
    -----
    - Aggregation uses: 1 - product(1 - x), the probability of at least one event.
    - Uses drop_duplicates to recover genomic coordinates after groupby merge.
    """
    df = pandas.read_csv(file)
    result_c = df.groupby(['lattice_loc', 'strand'])['predicted_occupancy'] \
                 .apply(lambda x: 1 - ((1 - x).prod())).reset_index()
    result = result_c.merge(df.drop_duplicates(['lattice_loc', 'strand']),
                            on=['lattice_loc', 'strand'], how='left')
    result = result.rename(columns={'predicted_occupancy_x': 'predicted_occupancy'})
    result = result[['chrom', 'start', 'end', 'mid', 'strand', 'lattice_loc', 'predicted_occupancy']]
    return result


def chip_seq_from_lef(lef_positions, site_number_per_replica, min_time=0):
    """
    Produce a simulated ChIP-seq–style coverage profile from LEF positions.

    This function extracts LEF left and right head locations across simulation
    time, flattens them, wraps them onto a single replica using modulo, and
    returns a histogram representing coverage. The resulting array mimics a
    simple ChIP-seq signal derived from LEF occupancy.

    Parameters
    ----------
    lef_positions : ndarray
        Array of LEF head positions with shape
        (time, n_lefs, 2), where the last axis is (left_head, right_head).
    site_number_per_replica : int
        Total number of sites on one lattice replica; used for modulo wrapping.
    min_time : int, optional
        Earliest simulation timestep to include. Defaults to 0.

    Returns
    -------
    numpy.ndarray
        Histogram of LEF occupancy across the replica lattice, representing
        per-site ChIP-like signal.

    Notes
    -----
    - Counts both left and right heads equally.
    - Uses modulo to merge multiple replicas into a single occupancy track.
    - The output array length equals site_number_per_replica.
    """
    lef_lefts = lef_positions[min_time:, :, 0].flatten()
    lef_rights = lef_positions[min_time:, :, 1].flatten()
    lef_positions_aray = np.hstack((lef_lefts, lef_rights))
    hist, hist_ary = np.histogram(
        np.mod(lef_positions_aray, site_number_per_replica),
        np.arange(0, site_number_per_replica, 1)
    )
    return hist



def chip_seq_from_ctcf(lef_file_path, site_number_per_replica):
    """
    Generate a simulated CTCF ChIP-seq–like occupancy track from LEF simulation output.

    This function loads CTCF binding positions from an HDF5 file produced by the
    loop-extrusion simulation. It extracts left-facing and right-facing CTCF
    binding events, multiplies positions by the corresponding site indicators,
    concatenates them, and bins them into a histogram representing per-site
    occupancy. If a CTCF site is listed in both left and right arrays
    (rare but possible for convergent or overlapping annotations),
    the count is divided by two to avoid double-counting.

    Parameters
    ----------
    lef_file_path : str or path-like
        Path to the HDF5 simulation output file containing datasets:
        'CTCF_positions_right', 'CTCF_positions_left',
        'CTCF_sites_right', and 'CTCF_sites_left'.
    site_number_per_replica : int
        Total number of lattice sites in one simulation replica. Used as the
        histogram domain.

    Returns
    -------
    numpy.ndarray
        Array of length `site_number_per_replica` containing CTCF occupancy
        counts per site, approximating a ChIP-seq–style signal.

    Notes
    -----
    - Flattened CTCF position arrays are multiplied by their corresponding
      site masks before histogramming.
    - Positions equal to zero (non-bound or inactive sites) are removed.
    - Sites appearing in both left- and right-facing lists are halved to
      avoid duplication.
    - Output provides a simple, interpretable per-site occupancy measure.
    """
    ctcf_array_right = np.array(h5py.File(lef_file_path, 'r')['CTCF_positions_right'])
    ctcf_array_left = np.array(h5py.File(lef_file_path, 'r')['CTCF_positions_left'])
    ctcf_array_right_sites = np.array(h5py.File(lef_file_path, 'r')['CTCF_sites_right'])
    ctcf_array_left_sites = np.array(h5py.File(lef_file_path, 'r')['CTCF_sites_left'])

    ctcfrightary = np.concatenate([arr.flatten() * ctcf_array_right_sites
                                   for arr in ctcf_array_right if arr.size > 0])
    ctcfleftary = np.concatenate([arr.flatten() * ctcf_array_left_sites
                                  for arr in ctcf_array_left if arr.size > 0])

    ctcfs = np.concatenate([ctcfrightary[ctcfrightary > 0],
                            ctcfleftary[ctcfleftary > 0]])

    ctcfhist, hist_array = np.histogram(ctcfs,
                                        np.arange(0, site_number_per_replica, 1))

    # Correct positions that appear in both left and right site arrays.
    common_list = np.intersect1d(ctcf_array_right_sites,
                                 ctcf_array_left_sites)
    for elements in common_list:
        ctcfhist[elements] = ctcfhist[elements] / 2

    return ctcfhist


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

def peak_ratio_around(num_sites_t, lef_positions, peak_positions, neighbor_size):
    """
    Compute the ratio of LEF occupancy at peak positions relative to nearby non-peak positions.

    This function builds a histogram of LEF positions across all sites and compares
    the total occupancy at specified peak positions to the occupancy in a region
    outside the peaks but within a neighboring window.

    Parameters
    ----------
    num_sites_t : int
        Total number of sites (used to define the histogram bins).

    lef_positions : array-like
        Positions of LEFs along the sites. Each value represents the site index
        where a LEF is located.

    peak_positions : array-like
        Indices of sites considered as peaks.

    neighbor_size : int
        Size of the neighboring region used to define outside-peak positions.

    Returns
    -------
    float
        The ratio of total LEF counts at peak positions to total LEF counts
        at outside-peak neighboring positions.
    """
    outside_peak_positions = np.arange(
        np.max(peak_positions) + 3,
        np.max(peak_positions + neighbor_size)
    )
    hist, edges = np.histogram(lef_positions, np.arange(num_sites_t + 1))
    return np.sum(hist[peak_positions]) / np.sum(hist[outside_peak_positions])

def chip_seq_from_ctcf(lef_file_path, site_number_per_replica):
    """
    Construct a ChIP-seq–like occupancy profile for CTCF sites from LEF simulation output.

    This function reads CTCF binding information (left- and right-oriented sites)
    from an HDF5 file, aggregates all occupied CTCF positions across replicas,
    and builds a histogram representing CTCF occupancy along the genomic sites.
    Sites that appear in both left- and right-oriented arrays are counted once
    by averaging their contributions.

    Parameters
    ----------
    lef_file_path : str
        Path to the HDF5 file containing CTCF position and site datasets
        (CTCF_positions_right, CTCF_positions_left, CTCF_sites_right,
        CTCF_sites_left).

    site_number_per_replica : int
        Total number of genomic sites per replica, used to define histogram bins.

    Returns
    -------
    numpy.ndarray
        One-dimensional array of length `site_number_per_replica` representing
        the aggregated CTCF occupancy (ChIP-seq–like signal) across sites.
    """
    ctcf_array_right = np.array(h5py.File(lef_file_path, 'r')['CTCF_positions_right'])
    ctcf_array_left = np.array(h5py.File(lef_file_path, 'r')['CTCF_positions_left'])
    ctcf_array_right_sites = np.array(h5py.File(lef_file_path, 'r')['CTCF_sites_right'])
    ctcf_array_left_sites = np.array(h5py.File(lef_file_path, 'r')['CTCF_sites_left'])

    ctcfrightary = np.concatenate(
        [arr.flatten() * ctcf_array_right_sites for arr in ctcf_array_right if arr.size > 0]
    )
    ctcfleftary = np.concatenate(
        [arr.flatten() * ctcf_array_left_sites for arr in ctcf_array_left if arr.size > 0]
    )

    ctcfs = np.concatenate(
        [ctcfrightary[ctcfrightary > 0], ctcfleftary[ctcfleftary > 0]]
    )

    ctcfhist, hist_array = np.histogram(
        ctcfs, np.arange(0, site_number_per_replica, 1)
    )

    common_list = np.intersect1d(ctcf_array_right_sites, ctcf_array_left_sites)
    for elements in common_list:
        ctcfhist[elements] = ctcfhist[elements] / 2

    return ctcfhist
def chip_seq_from_lef(lef_positions, site_number_per_replica, min_time=0):
    """
    Construct a ChIP-seq–like occupancy profile for LEFs from simulation trajectories.

    This function aggregates left and right LEF positions across replicas and time,
    starting from a specified minimum time point, and builds a histogram representing
    LEF occupancy along genomic sites. Positions are wrapped using modulo to map
    LEFs onto a single replica coordinate system.

    Parameters
    ----------
    lef_positions : numpy.ndarray
        Array of LEF positions with shape (time, replicas, 2), where the last
        dimension corresponds to left and right LEF positions.

    site_number_per_replica : int
        Total number of genomic sites per replica, used to define histogram bins.

    min_time : int, optional
        Time index from which LEF positions are included in the aggregation.
        This can be used to discard early transient dynamics (default is 0).

    Returns
    -------
    numpy.ndarray
        One-dimensional array of length `site_number_per_replica` representing
        the aggregated LEF occupancy (ChIP-seq–like signal) across sites.
    """
    lef_lefts = lef_positions[min_time:, :, 0].flatten()
    lef_rights = lef_positions[min_time:, :, 1].flatten()

    lef_positions_aray = np.hstack((lef_lefts, lef_rights))

    hist, hist_ary = np.histogram(
        np.mod(lef_positions_aray, site_number_per_replica),
        np.arange(0, site_number_per_replica, 1)
    )

    return hist
















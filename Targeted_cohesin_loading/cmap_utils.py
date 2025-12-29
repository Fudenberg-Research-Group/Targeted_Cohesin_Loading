import bioframe

import numpy as np
import multiprocessing as mp
import os
from functools import partial
import h5py
import networkx as nx
import time 




##### Required functions for making 1D contact maps
def create_lattice_graph(n, Lefs):
    """
    Create a one-dimensional lattice graph with additional loop (LEF) connections.

    This function builds an undirected graph with n nodes arranged in a linear
    lattice, where each node is connected to its immediate right neighbor.
    Additional edges are added to represent loop connections defined by LEF pairs.

    Parameters
    ----------
    n : int
        Total number of nodes in the lattice.

    Lefs : array-like of tuple(int, int)
        List of node index pairs representing loop (LEF) connections to be added
        on top of the regular lattice edges.

    Returns
    -------
    networkx.Graph
        An undirected graph representing the lattice with loop connections.
    """
    G = nx.Graph()

    for i in range(n):
        if i + 1 < n:
            G.add_edge(i, i + 1)

    for i, j in Lefs:
        G.add_edge(i, j)

    return G




def closest_distance(G, start, end):
    """
    Compute the shortest path distance between two nodes in a graph.

    This function returns the length of the shortest path between the specified
    start and end nodes in the given graph. If no path exists, it returns infinity.

    Parameters
    ----------
    G : networkx.Graph
        Graph in which the distance is computed.

    start : int
        Index of the starting node.

    end : int
        Index of the target node.

    Returns
    -------
    int or float
        Length of the shortest path between start and end. Returns infinity if
        no path exists.
    """
    try:
        return nx.shortest_path_length(G, source=start, target=end)
    except nx.NetworkXNoPath:
        return float('inf')



        
def calculate_contact_maps(
    total_sites,
    lefs,
    str_frame,
    end_frame,
    every_frame,
    max_dist,
    res_convert,
    replication_number,
    output_dir
):
    """
    Compute and save a contact map from LEF-mediated lattice graphs over simulation frames.

    This function constructs contact maps by converting LEF positions into loop
    connections on a one-dimensional lattice, computing graph-based distances
    between site pairs, and accumulating contact strengths across time frames
    and replicated domains. Contacts decay with graph distance following a
    power-law relationship.

    Parameters
    ----------
    total_sites : int
        Total number of lattice sites before resolution conversion.

    lefs : numpy.ndarray
        Array of LEF positions with shape (time, replicas, 2), where each entry
        represents a pair of loop anchor positions.

    str_frame : int
        Starting frame index (inclusive) for contact map accumulation.

    end_frame : int
        Ending frame index (exclusive) for contact map accumulation.

    every_frame : int
        Step size for frame sampling (e.g., use every Nth frame).

    max_dist : int
        Maximum genomic separation (in sites) for which contacts are computed.

    res_convert : int
        Resolution conversion factor used to coarse-grain lattice sites.

    replication_number : int
        Number of replicated domains tiled along the lattice.

    output_dir : str
        Directory where the computed contact map will be saved.

    Returns
    -------
    None
        The function saves the contact map as a compressed NumPy file
        ('contact_map.npz') in the specified output directory.
    """
    N = total_sites // res_convert
    mod_i_values = mod_j_values = np.mod(
        np.arange(N // replication_number), N // replication_number
    )
    sites_p_r = N // replication_number
    contact_matrix = np.zeros((sites_p_r, sites_p_r))

    for frame in range(str_frame, end_frame, every_frame):
        slice_1 = lefs[frame, :, :] // res_convert

        for dupl in range(replication_number):
            start_idx = dupl * (N // replication_number)
            end_idx = (dupl + 1) * (N // replication_number)

            mask = (slice_1 > start_idx) & (slice_1 < end_idx)
            pair_mask = np.all(mask, axis=1)
            filtered_pairs = slice_1[pair_mask]
            final_pairs = np.mod(filtered_pairs, N // replication_number)

            G = create_lattice_graph(N // replication_number, final_pairs)

            for i in range(N // replication_number):
                for j in range(i + 1, N // replication_number):
                    if j < i + max_dist:
                        dist = closest_distance(G, i, j)
                        contact = 1 / (dist) ** 1.5
                        contact_matrix[i, j] += contact
                        contact_matrix[j, i] += contact

    np.savez_compressed(
        os.path.join(output_dir, 'contact_map.npz'),
        contact_map=contact_matrix
    )




def calculate_contact_map_save(
    lefs,
    str_frame,
    end_frame,
    every_frame,
    max_dist,
    res_convert,
    replication_number,
    output_dir
):
    """
    Compute, accumulate, and save a contact map from LEF-mediated graph distances.

    This function iterates over selected simulation frames, constructs a lattice
    graph with LEF-induced loop connections, and computes contact strengths between
    site pairs based on graph shortest-path distances. Contacts are accumulated
    across replicated domains and summed over all processed frames before saving.

    Parameters
    ----------
    lefs : numpy.ndarray
        Array of LEF positions with shape (time, replicas, 2), where each entry
        represents a pair of loop anchor positions.

    str_frame : int
        Starting frame index (inclusive) for contact map computation.

    end_frame : int
        Ending frame index (exclusive) for contact map computation.

    every_frame : int
        Step size for frame sampling (e.g., use every Nth frame).

    max_dist : int
        Maximum genomic separation (in sites) for which contact distances are
        explicitly computed.

    res_convert : int
        Resolution conversion factor used to coarse-grain LEF positions and
        lattice sites.

    replication_number : int
        Number of replicated domains tiled along the lattice.

    output_dir : str
        Directory where the accumulated contact map will be saved.

    Returns
    -------
    None
        The function saves the accumulated contact map as a compressed NumPy file
        ('contact_map.npz') in the specified output directory.
    """
    contact_map = []
    N = np.max(lefs) // res_convert + 1
    mod_i_values = mod_j_values = np.mod(
        np.arange(N // replication_number), N // replication_number
    )
    sites_p_r = N // replication_number

    for frame in range(str_frame, end_frame, every_frame):
        contact_matrix = np.zeros((sites_p_r, sites_p_r))
        lefs_t = lefs[frame, :, :] // res_convert
        start = time.time()

        G = create_lattice_graph(N, lefs_t)
        end = time.time()
        print(f"Elapsed time for graph: {end - start} seconds")

        start = time.time()
        for i in range(sites_p_r):
            for j in range(i + max_dist, sites_p_r):
                contact_matrix[i, j] = contact_matrix[j, i] = (
                    replication_number * (1 / (j - i) ** 1.5)
                )

        for dupl in range(replication_number):
            start_idx = dupl * (N // replication_number)
            end_idx = (dupl + 1) * (N // replication_number)

            for i in range(start_idx, end_idx):
                for j in range(i + 1, end_idx):
                    if j < i + max_dist:
                        dist = closest_distance(G, i, j)
                        contact = 1 / (dist + 1) ** 1.5
                        contact_matrix[
                            mod_i_values[i - start_idx],
                            mod_j_values[j - start_idx]
                        ] += contact
                        contact_matrix[
                            mod_j_values[j - start_idx],
                            mod_i_values[i - start_idx]
                        ] += contact

        contact_map.append(contact_matrix)
        end = time.time()
        print(f"Elapsed time for contact map calculation: {end - start} seconds")

    np.savez_compressed(
        os.path.join(output_dir, 'contact_map.npz'),
        contact_map=np.sum(contact_map, axis=0)
    )





def create_contact_map_folders(n, output_directory):
    """
    Create and return directories for storing multiple contact map outputs.

    This function creates `n` subdirectories inside the specified output directory,
    each named sequentially as 'contactmap_1', 'contactmap_2', ..., 'contactmap_n'.
    Existing directories are preserved.

    Parameters
    ----------
    n : int
        Number of contact map directories to create.

    output_directory : str
        Path to the parent directory where contact map folders will be created.

    Returns
    -------
    list of str
        List of full paths to the created contact map directories.
    """
    output_dirs = []
    for contact_id in range(1, n + 1):
        file_name = f"contactmap_{contact_id}"
        output_directory_partial = os.path.join(output_directory, file_name)
        os.makedirs(output_directory_partial, exist_ok=True)
        output_dirs.append(output_directory_partial)
    return output_dirs



def contact_map_from_lefs(dset, sites_per_replica):
    """
    Construct a symmetric contact map directly from LEF loop positions.

    This function reshapes LEF position data into left–right pairs, maps positions
    onto a single replica using modulo arithmetic, filters valid loops where the
    right position is greater than the left position, and builds a 2D histogram
    representing contact frequencies between site pairs. The final contact map
    is symmetrized.

    Parameters
    ----------
    dset : numpy.ndarray
        Array of LEF positions that can be reshaped into pairs of loop anchors.
        The data is expected to contain left and right positions in alternating
        columns.

    sites_per_replica : int
        Total number of genomic sites per replica, used to define histogram bins.

    Returns
    -------
    numpy.ndarray
        A symmetric 2D array of shape (sites_per_replica, sites_per_replica)
        representing the LEF-derived contact map.
    """
    lef_array = np.mod(dset.reshape((-1, 2)), sites_per_replica)
    lef_array = lef_array[lef_array[:, 1] > lef_array[:, 0]]

    lef_map = np.histogram2d(
        lef_array[:, 0],
        lef_array[:, 1],
        np.arange(sites_per_replica)
    )[0]

    return lef_map + lef_map.T



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
    






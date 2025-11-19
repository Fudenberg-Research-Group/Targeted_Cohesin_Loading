import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import bioframe
import cooler
import pyBigWig
import OccupancyInputCTCF.utils.experimental_path as exp
import cooltools
from cooltools.lib.plotting import *


### importing Hi-C maps 
clr = cooler.Cooler(exp.bonev_file+'::resolutions/10000')
### importing RAD21 & CTCF peaks/motifs
RAD21=bioframe.read_table(exp.RAD21_path, schema='bed')
ctcf = bioframe.read_table(exp.ctcf_peaks, schema='bed')
motif=bioframe.read_table(exp.motif_directory)

ctcf['mid']=(ctcf.end+ctcf.start)/2
RAD21['mid']=(RAD21.end+RAD21.start)/2

bw_rad21 = pyBigWig.open(exp.bw_path_rad21)
bw_ctcf = pyBigWig.open(exp.bw_path_ctcf)

def set_axes_ary(region, lattice_size=250, res = 100000):
    region_ = bioframe.parse_region_string(region)
    ### array to show in the map
    mb = 1e6
    start_reg = region_[1]
    end_reg = region_[2]
    ary = list(np.arange(0, (end_reg-start_reg)//(lattice_size)+1, res//lattice_size)) 
    ary_str = [str((elements/mb)) for elements in np.arange(start_reg, end_reg+1,res)]
    return ary, ary_str

def plot_chip_hic(region, chip, chip_ctcf, binned_matrix, clr=clr, lattice_size = 0.25, res= 100000, output_file='plots.png'):
    """
    Hi-C and ChIP-seq plots along with simulated data visualizations.
    
    Parameters:
    - region (str): Genomic region to visualize (e.g., 'chr6:50,000,000-51,000,000').
    - binned_matrix (ndarray): Simulated Hi-C matrix.
    - chip (ndarray): Simulated 1D extruder data.
    - output_file (str): Path to save the output plot.
    """

    region_ = bioframe.parse_region_string(region)
    chrom, start_reg, end_reg = region_

    values_rad21 = np.nan_to_num(bw_rad21.values(chrom, start_reg, end_reg))
    values_ctcf = np.nan_to_num(bw_ctcf.values(chrom, start_reg, end_reg))

    # Set up plot layout
    plt_width = 4
    f, axs = plt.subplots(
        figsize=(plt_width + plt_width + 2, plt_width + 5),
        ncols=4,
        nrows=3,
        gridspec_kw={'height_ratios': [4, 1, 1], "wspace": 0.01, 'width_ratios': [1, 0.04, 1, 0.04]},
        constrained_layout=True
    )

    norm = LogNorm(vmax=0.1)

    # Raw Hi-C data
    ax = axs[0, 0]
    im = ax.matshow(
        clr.matrix().fetch(region),
        norm=norm,
        cmap='fall'
    )
    ax.set_title(f'{chrom}:{start_reg:,}-{end_reg:,}')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    cax = axs[0, 1]
    plt.colorbar(im, cax=cax, label='raw counts')

    # RAD21 ChIP-seq
    ax1 = axs[1, 0]
    ax1.plot(np.arange(start_reg, end_reg), values_rad21)
    ax1.set_xlim([start_reg, end_reg])
    ax1.set_xlabel('position, bins')
    ax1.set_title('RAD21')
    axs[1, 1].set_visible(False)

    # CTCF ChIP-seq
    ax1 = axs[2, 0]
    ax1.plot(np.arange(start_reg, end_reg), values_ctcf, color='green', label='CTCF ChIP-seq')
    ax1.set_xlim([start_reg, end_reg])
    ax1.set_xlabel('position, bins')
    ax1.set_title('CTCF')
    axs[2, 1].set_visible(False)

    # Simulated Hi-C map
    ax = axs[0, 2]
    im = ax.matshow(np.log10(binned_matrix),  cmap='fall')
    ax.set_title('1D Simulation, %s'%region_[0])
    cax = axs[0, 3]
    plt.colorbar(im, cax=cax, label='raw counts')

    # Simulated extruder data
    ax1 = axs[1, 2]
    ax1.plot(np.arange(len(chip)) , chip / np.sum(chip)) #* lattice_size
    ax1.set_xlabel('extruders position, mbp')
    ax1.set_title('Simulation, Extruder')
    ax1.set_xlim(0, len(chip))
    xticks, xticklabels = set_axes_ary(region, lattice_size=250, res = res)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticklabels, ha="right")  
    axs[1, 3].set_visible(False)


    ax1 = axs[2, 2]
    ax1.plot(np.arange(len(chip_ctcf)), chip_ctcf / np.sum(chip_ctcf), color='red')
    ax1.set_xlabel('barrier position')
    ax1.set_title('Simulation, CTCF')
    ax1.set_xlim(0, len(chip_ctcf))
    xticksa, xticklabelsa = set_axes_ary(region, lattice_size=250, res = res)
    ax1.set_xticks(xticksa)
    ax1.set_xticklabels(xticklabelsa, ha="right") 
    axs[2, 3].set_visible(False)

    # Save and show plot
    plt.savefig(output_file)
    plt.show()

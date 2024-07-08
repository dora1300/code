"""
@Author:            Dominique A Ramirez
@Date:              2024 06 24
@Title:             plot_ete_comparisons_ns.py

@Description:       This script plots the density profiles of centered slabs, taking the average across
three different replicates to make the final product.

@Updates:

"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.size'] = 14


parser = argparse.ArgumentParser()
parser.add_argument('-exid', help="Experimental ID that is the prefix for the directories that I need to get into"
                    " e.g. DAR3-70", required=True, type=str)
parser.add_argument('-subdir', help="The letter codes delineating the different replicates. Provide as a series "
                    "of letter e.g. `-subdir c d e`", nargs='+', required=True)
parser.add_argument('-protein', help="The type of protein, e.g. `c3l2`", required=True, type=str)
parser.add_argument('-prefix', 
    help="The prefix of all ETE associated files (output from analysis)", 
    type=str)
parser.add_argument('-outname', help="The common name/ID you want to give to all the output files.", type=str,
                    default="output")
parser.add_argument('-start', help="The start time of analysis i.e. equilibration time [ps]", required=True, type=int)


# Parse the arguments and set up some global variables
args = parser.parse_args()
EXPID = args.exid
PROT = args.protein
OUTPUT = args.outname
START = args.start
PREFIX = args.prefix

subdirlist = args.subdir

templist = [273, 298, 310, 333]


"""
Parse the files for the COILS only!!
"""
# Parse the files, and add the data to the respective arrays for the downstream data analysis
ar273coils = []; ar298coils = []; ar310coils = []; ar333coils = []

for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_273_ETE_by_frame_COILS.csv",
        delimiter=',')
    ar273coils.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_298_ETE_by_frame_COILS.csv",
        delimiter=',')
    ar298coils.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_310_ETE_by_frame_COILS.csv",
        delimiter=',')
    ar310coils.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_333_ETE_by_frame_COILS.csv",
        delimiter=',')
    ar333coils.append(load_ete[1])
    
    
    
"""
Parse the files for the LINKERS only!!
"""
# Parse the files, and add the data to the respective arrays for the downstream data analysis
ar273linkers = []; ar298linkers = []; ar310linkers = []; ar333linkers = []

for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_273_ETE_by_frame_LINKERS.csv",
        delimiter=',')
    ar273linkers.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_298_ETE_by_frame_LINKERS.csv",
        delimiter=',')
    ar298linkers.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_310_ETE_by_frame_LINKERS.csv",
        delimiter=',')
    ar310linkers.append(load_ete[1])
    
for i, val in enumerate(subdirlist):
    print(f"replicate: {i}")
    load_ete = np.loadtxt(f"./{EXPID}{val}/combined_ETEanalysis/{PREFIX}_333_ETE_by_frame_LINKERS.csv",
        delimiter=',')
    ar333linkers.append(load_ete[1])
    
time_array = np.loadtxt(f"./{EXPID}{subdirlist[0]}/combined_ETEanalysis/{PREFIX}_273_ETE_by_frame_LINKERS.csv",
        delimiter=',')[0]



"""
Find the index where the time array corresponds to the start time 
provided in the arguments
"""
begin = np.argwhere(time_array == args.start)[0][0]




"""
Now it's time for some heinous averages. The order of the averages is important (I believe)
"""
# first up, we're doing time based averaged of each single replicate
time_avgs_273_coil = []; time_avgs_298_coil = []
time_avgs_310_coil = []; time_avgs_333_coil = []

for i in range(0, 3):
    time_avgs_273_coil.append(
        np.average(ar273coils[i][begin:])
)
    time_avgs_298_coil.append(
        np.average(ar298coils[i][begin:])
)
    time_avgs_310_coil.append(
        np.average(ar310coils[i][begin:])
)
    time_avgs_333_coil.append(
        np.average(ar333coils[i][begin:])
)

time_avgs_273_linker = []; time_avgs_298_linker = []
time_avgs_310_linker = []; time_avgs_333_linker = []

for i in range(0, 3):
    time_avgs_273_linker.append(
        np.average(ar273linkers[i][begin:])
)
    time_avgs_298_linker.append(
        np.average(ar298linkers[i][begin:])
)
    time_avgs_310_linker.append(
        np.average(ar310linkers[i][begin:])
)
    time_avgs_333_linker.append(
        np.average(ar333linkers[i][begin:])
)


# now we can do replicate averages
rep_avg_273_coil = np.average(time_avgs_273_coil)
rep_std_273_coil = np.std(time_avgs_273_coil)
rep_avg_298_coil = np.average(time_avgs_298_coil)
rep_std_298_coil = np.std(time_avgs_298_coil)
rep_avg_310_coil = np.average(time_avgs_310_coil)
rep_std_310_coil = np.std(time_avgs_310_coil)
rep_avg_333_coil = np.average(time_avgs_333_coil)
rep_std_333_coil = np.std(time_avgs_333_coil)

rep_avg_273_linker = np.average(time_avgs_273_linker)
rep_std_273_linker = np.std(time_avgs_273_linker)
rep_avg_298_linker = np.average(time_avgs_298_linker)
rep_std_298_linker = np.std(time_avgs_298_linker)
rep_avg_310_linker = np.average(time_avgs_310_linker)
rep_std_310_linker = np.std(time_avgs_310_linker)
rep_avg_333_linker = np.average(time_avgs_333_linker)
rep_std_333_linker = np.std(time_avgs_333_linker)


"""
NOW I CAN DO SOME PLOTTING
"""
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharey=True)

fig.supylabel("Distance (nm)")

axs[0, 0].bar(0, rep_avg_273_coil, yerr=rep_std_273_coil, width=0.6, 
        color="grey", ecolor="black", capsize=15, 
        label=f"Coils ({rep_avg_273_coil:.3f})")
axs[0, 0].bar(1, rep_avg_273_linker, yerr=rep_std_273_linker, width=0.6, 
        color="goldenrod", ecolor="black", capsize=15,
        label=f"Linkers ({rep_avg_273_linker:.3f})")
axs[0, 0].legend()
axs[0, 0].set_title("T = 273 K")
axs[0, 0].set_xticks([0, 1.0], ["coils", "linkers"])

axs[0, 1].bar(0, rep_avg_298_coil, yerr=rep_std_298_coil, width=0.6, 
        color="grey", ecolor="black", capsize=15, 
        label=f"Coils ({rep_avg_298_coil:.3f})")
axs[0, 1].bar(1, rep_avg_298_linker, yerr=rep_std_298_linker, width=0.6, 
        color="goldenrod", ecolor="black", capsize=15,
        label=f"Linkers ({rep_avg_298_linker:.3f})")
axs[0, 1].legend()
axs[0, 1].set_title("T = 298 K")
axs[0, 1].set_xticks([0, 1.0], ["coils", "linkers"])

axs[1, 0].bar(0, rep_avg_310_coil, yerr=rep_std_310_coil, width=0.6, 
        color="grey", ecolor="black", capsize=15, 
        label=f"Coils ({rep_avg_310_coil:.3f})")
axs[1, 0].bar(1, rep_avg_310_linker, yerr=rep_std_310_linker, width=0.6, 
        color="goldenrod", ecolor="black", capsize=15,
        label=f"Linkers ({rep_avg_310_linker:.3f})")
axs[1, 0].legend()
axs[1, 0].set_title("T = 310 K")
axs[1, 0].set_xticks([0, 1.0], ["coils", "linkers"])

axs[1, 1].bar(0, rep_avg_333_coil, yerr=rep_std_333_coil, width=0.6, 
        color="grey", ecolor="black", capsize=15, 
        label=f"Coils ({rep_avg_333_coil:.3f})")
axs[1, 1].bar(1, rep_avg_333_linker, yerr=rep_std_333_linker, width=0.6, 
        color="goldenrod", ecolor="black", capsize=15,
        label=f"Linkers ({rep_avg_333_linker:.3f})")
axs[1, 1].legend()
axs[1, 1].set_title("T = 333 K")
axs[1, 1].set_xticks([0, 1.0], ["coils", "linkers"])


plt.savefig(f"{args.outname}_ETE_averaged_plot.png", dpi=600)

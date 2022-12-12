"""
@Title:             coil_partner_exchange_analysis.py
@Name:              Mando A Ramirez
@Date:              2022 04 20

@Description:       This script calculates a lot of things!
FIRST ---- This script calculates how many partners an individual coil associates with other the course of the
simulation. This will be plotted as a distribution and represents the partner exchange frequency. This script will also
analyze the fraction of free versus oligomerized coils throughout the simulation.
    This script only counts A-A interactions and does not include any other possibilities. Thus, this script needs to
know where the A-beads are, which is provided through the '-sA' flag.
    I am making the operational choice that a coil has to make at least 3 A-A interactions to count as an interaction.

SECOND --- This script calculates partner swapping. Essentially, it keeps track of how often a partner exchange occurs,
where a partner exchange is defined by a partner being replaced by another partner, or the complete dissociation of an
interaction.

@Updates:

"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd


""" Set up the arg parser """
parser = argparse.ArgumentParser(description="Partner Exchange analysis tool - this calculates the number of partners"
                                             " that a coil has throughout a simulation. This will create a plot for "
                                             "each individual coil as well as output statistics for further analysis.\n"
                                             "This is currently ONLY valid for simulations of coils ALONE.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-l", help="no. of beads in a coil. They must all have the same size.", required=True, type=int)
parser.add_argument("-Ai", help="the index of the first A-bead in the very first coil. This is important for counting"
                                " only A-A interactions (1-indexed)", required=True, type=int)
parser.add_argument("-f", help="the stride of analysis i.e. analyze every fth frame, default = 1", default=1, type=int)
parser.add_argument("-o", help="text that you want in the output file", default="output_partnerFreq")
parser.add_argument("-plot", help="Pass flag to turn on plotting. Otherwise, only calculations are made",
                    action="store_true")


args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
nsize = args.l          # this is also 1-indexed, essentially
output = args.o

traj = md.load(trajectory, top=topology)

"""
First thing to do - make a collection that stores information about which bead belongs to which coil
I will create range operators correspond to bead indices for each coil, this will help me identify the parent coil for
any given bead in the system.
"""
coil_catalog = []
for i in range(0, nmodels):
    start_bead = 0 + (nsize * i)
    end_bead = nsize + (nsize * i)
    coil_catalog.append(np.arange(start_bead, end_bead, dtype=int))

# Now make a list of all the A beads, separated by coil
# and make sure to 0-index it for mdtraj!
coilA_indices = []
current_coilA = []
for i in range(0, nmodels):
    Aa_ind = (args.Ai - 1) + (nsize * i)
    Ad_ind = (args.Ai - 1) + (nsize * i) + 3
    while Ad_ind < (nsize * nmodels) and Aa_ind < (nsize * nmodels):
        # this loop will add indices of the A-beads and once the last A-bead within the bounds of the coil are added,
        # then the loop breaks
        if Aa_ind >= (nsize * (i+1)):         # will only trigger once Aa extends beyond size of coil
            if Ad_ind < (nsize * (i+1)):
                current_coilA.append(Ad_ind)
            break
        if Ad_ind >= (nsize * (i+1)):         # will only trigger once Ad extends beyond size of coil
            if Aa_ind < (nsize * (i+1)):
                current_coilA.append(Aa_ind)
            break
        current_coilA.append(Aa_ind)
        current_coilA.append(Ad_ind)
        Aa_ind += 7
        Ad_ind += 7

    coilA_indices.append(np.array(current_coilA))
    current_coilA = []

## coilA_indices is a list whose entries are all arrays corresponding to the indices of A beads in each coil
## i.e. each entry corresponds to A beads in coil index+1 and Abead indices (0-indexed) are stored in an array
## e.g. coilA_indices = [{Abead indices array for coil 1}, {Abead indicies array for coil 2}, ...]

"""     WARNING WARNING WARNING
        2022 04 20
    This currently does not count the right number of A beads in the final coil. I am moving on however and will come 
back to solving this problem
"""

"""
The storage container for this analysis is per coil. Each coil will be assigned a numpy array whose length is the 
number of coils in the system + 1 for free coil. It will represent something like this:
coil1_array = [frequency bound coil1, frequency bound coil2, frequency bound coil3,...,frequency unbound]

Each of these arrays will be stored in a larger list and will be callable by index, i.e.
master_list = [{coil1 array}, {coil2 array}, {coil3 array}, ...]
"""
coilF_list = []     # coilF = coil frequency list
i = 0
while i < nmodels:
    coilF_list.append(np.zeros(nmodels+1, int))
    i += 1

# add in a sanity check to make sure I have done things correctly
# this checks to make sure that the list of A indices and the storage container are the right size
if len(coilA_indices) != nmodels or len(coilF_list) != nmodels:
    print("There is an error. Somehow, there is a mismatch between the size of the data containers and the number"
          " of coils provided in the argparser. Please correct this. Exiting now")
    exit(1)



"""
Begin analysis of the trajectory
"""
# I will be doing the analysis frame by frame because I don't want to sort through a wickedly huge array containing
# information for every single frame in the trajectory, though I realize this might slow down the analysis quite a bit

for frame in range(0, traj.n_frames, args.f):
    for coil_n, coil_As in enumerate(coilA_indices):
        # coil_n is the index of the coil in question (0-indexed)
        # coil_As is the indices of A-beads in coil_n

        # this is a temporary array that will allow me to count the number of interactions of coil_n in this single frame
        coil_interaction_count = np.zeros(20, dtype=int)

        """
        this computes the neighbors for the frame for A-beads in coil_n, BUT only searches for neighbors that are other
        A type beads. The haystack is designed to only contain other A-bead indices
        """
        # coil_As are the query_indices
        # haystack is the haystack indices
        # this works by making a haystack of all the A-bead indices that do NOT belong to the current coil coil_n
        haystack = np.array([], dtype=int)
        for index_a, coil_range in enumerate(coilA_indices):
            if index_a == coil_n:
                continue
            else:
                haystack = np.concatenate((haystack, coil_range), dtype=int)

        coil_n_neighbors = md.compute_neighbors(traj[frame], 0.70, coil_As, haystack_indices=haystack)[0]    # only need the first item

        # Now go through, identify if there are sufficient contacts to warrant an interaction, and then identify which
        # other coil(s) that coil is interacting with
        # if the number of neighbors is 0, then obviously there is no interaction and coil_n is unbound!
        if len(coil_n_neighbors) < 3:
            coilF_list[coil_n][-1] += 1
            continue                    # trigger the move onto the next coil_n in the loop

        for neighbor_i in coil_n_neighbors:
            for index_b, coil_i_range in enumerate(coil_catalog):
                if neighbor_i in coil_i_range:
                    # neighbor_i is part of the coil specified by coil_i_range whose index is index_b!
                    coil_interaction_count[index_b] += 1

        # at the end of this loop, then I will have my array, coil_interaction_count, that has attributed every neighbor
        # to every other coil. e.g. if coil_1 is interacting with coil_3 and coil_5, the array will look like this:
        # coil_interaction_count = [0, 0, 5, 0, 4, 0, 0, ...] which means that coil_1 interacts 5 times with coil_3 and
        # 4 times with coil_5. Now I have to add to the master frequency counting list, coilF_list
        for other_coil_index, interaction_count in enumerate(coil_interaction_count):
            if interaction_count >= 3:
                coilF_list[coil_n][other_coil_index] += 1

        # now that this is done, return to the top of the loop and go through all the steps again for the next coil
        # then, once all coils are cycled through for this frame, move onto the next frame
        # obviously, this is a lot of iterations and steps. There is room for improving the analysis


"""
Now time to save the frequency data for later use!
"""
pd.DataFrame(coilF_list).to_csv(f"{output}.csv")

"""
Now move onto the plotting!
"""
# for convenience, I am going to divide the frequencies of every entry in coilF_list by the sum of interactions so that
# i have fraction frequencies
# due to the large number of plots are produced, I am going to do some directory handling
x_axis_labels = []
for entry_i in range(len(coilF_list)):
    coilF_list[entry_i] = coilF_list[entry_i] / np.sum(coilF_list[entry_i])
    x_axis_labels.append(f"Coil {entry_i+1}")
x_axis_labels.append("Unbound")
pd.DataFrame(coilF_list).to_csv(f"{output}_normalized.csv")

if args.plot:
    if os.path.isdir("./partner_frequency_plots"):
        os.chdir("./partner_frequency_plots")
    else:
        os.mkdir("partner_frequency_plots")
        os.chdir("./partner_frequency_plots")

    for coil_number, entry in enumerate(coilF_list):
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.bar(x_axis_labels, entry, width=0.5, color="blue", alpha=0.5)
        plt.xticks(rotation=45)
        plt.xlabel("Coils in simulation")
        plt.ylabel("Fractional frequency")
        plt.title(f"Frequency of interactions between Coil-*{coil_number+1}* and \nevery other coil in the simulation")
        plt.savefig(f"partnerFrequencyPlot_coil{coil_number+1}.png", dpi=300)
        plt.close()
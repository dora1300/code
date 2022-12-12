"""
@Title:         ete_analysis.py
@Author:        Mando Ramirez
@Date:          20211027

@Description:   This code is for analyzing trajectories of coarse grained coil simulations to calculate the end-to-end
distances of each individual helix rod present in the simulation. The purpose of this is to see if the helices bend
back on themselves or exist as elongated helices/coils.

@Updates:
20211105 - removed the extension necessity in the output part of the argparse, and I'm letting the code handle it now.

20211209 - fixed some small length calculation errors. Also fixed the plotting so the labels can handle however many
models there are.

20211210 - added a "print" functionality so that I can just print the averaged results to stdout instead of resaving
files or making plots.

20220307 - adding a "distribution" plotting functionality, to instead plot distributions of distances for a simulation
instead of a time-series. Probably more helpful in the long run, really. Also will add KDE functionality in the future.
"""

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import scipy.stats as st

#   Global constants section
threshold = 5.17        # this is the defined threshold length of a helix rod. Anything less suggests its collapsing on
                        #    itself
buffer = 1.0            # the buffer tolerance of the threshold. Essentially, anything that is greater than
                        #    (threshold - buffer) will get counted as still a rod shape.

#   Argument Parser section
parser = argparse.ArgumentParser(description="Rod-length analysis for coiled-coil simulations. Only measures end-to-end"
                                             " distances of coil helices.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_ETEanalysis")
parser.add_argument("-plot_timeseries", help="Provide flag to turn ON time series plotting, one for each model.", action="store_true")
parser.add_argument("-plot_distribution", help="Provide flag to turn ON distribution style plotting, one for each model.", action="store_true")
parser.add_argument("-print", help="switch to choose to ONLY print to stdout, no saving files!", default=False, type=bool)


args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
lhelix = args.l         # this is also 1-indexed, essentially
sim_code = args.id
output = args.o + ".csv"

# Create a directory for ETE plots if I provide the flag for plotting.
if args.plot_timeseries or args.plot_distribution:
    if os.path.exists("./ETE_plots"):
        pass
    else:
        os.mkdir("./ETE_plots")


# load the trajectory and create the atom pairs for distance calculation based on the length of the helices and
# the number of models
traj = md.load(trajectory, top=topology)

apairs = []
for i in range(nmodels):
    a_start = 0 + (lhelix * i)                  #fixed to handle the helix length implicitly
    a_end = (lhelix - 1) + (lhelix * i)         # fixed to handle the helix length implicitly. Also keep in mind atom
                                                # positions are 0-index thus (lhelix-1)
    apairs.append([a_start, a_end])
# check to make sure I haven't done something stupid. The last atom in the apairs list of lists should match the no.
# atoms in the trajectory/topology
if apairs[-1][1] == (traj.n_atoms-1):
    pass
else:
    print("The number of models and length of helices does not match the number of atoms in the trajectory. Exiting now.")
    raise Exception("AtomMismatchError")
    exit()


# compute the distances of the provided atom pairs
# This will create a 2d array that takes the form [[dist apair1, dist apair2, dist apair 3], ...] where every entry
# corresponds to a frame, and then within each frame entry, the distances are listed in order of pairs provided to it
dist = md.compute_distances(traj, np.array(apairs))


# Open the output file and save the data to it. This saves the length for each model calculated throughout the trajectory
# Do this for every single invocation of this script.
outfile = open(output, 'w')
firstline = "Frame No."
for i in range(nmodels):
    firstline += f",helix_{i+1}_dist(nm)"
firstline += "\n"
outfile.write(firstline)

for i in range(len(dist)):
    nextline = f"{i}"
    for j in range(len(dist[i])):
        nextline += f",{dist[i][j]:.4f}"
    nextline += "\n"
    outfile.write(nextline)

outfile.close()


# # Analyze the distances and output them into another output file. In this analysis, I will get the average distance
# # over time, along with the proportion of frames where the distance is above the defined threshold of 5.17nm
## COMMENTED OUT AS OF 2022 03 07 -- MIGHT BRING BACK EVENTUALLY
# averages = []
# for i in range(nmodels):
#     averages.append(np.average(dist[:,i]))
# stddev = []
# for i in range(nmodels):
#     stddev.append(np.std(dist[:,i]))
#
# for i in range(len(averages)):
#     print(f"Model Number: {i+1}")
#     print(f"Avg helix ETE distance (nm): {averages[i]}")
#     print(f"Stdev helix ETE distance (nm): {stddev[i]}")
#     print()
#
# prop_count_list = []            # this list holds the raw count of how many frames a given coil is above the
#                                 #     threshold length
# proportion_list = []            # this is the prop_count_list but an actual proportion by dividing the raw count by the
#                                 #    simulation length (i.e. no. frames)
# for i in range(nmodels):
#     h_count = 0
#     for j in range(len(dist[:,i])):
#         if dist[:,i][j] >= (threshold - buffer):
#             h_count += 1
#         else:
#             pass
#     prop_count_list.append(h_count)
#     proportion_list.append(h_count/traj.n_frames)
#
# # now open the file and save the data to it. This saves the proportion data.
# if not args.print:
#     outfile_split = output.split(".")
#     outfile2_name = outfile_split[0] +"_proportions.csv"
#     outfile2 = open(outfile2_name, 'w')
#     outfile2.write("Helix no.,average_dist(nm),stddev_dist(nm),raw_counts_as_rod,proportion_as_rod\n")
#     for i in range(nmodels):
#         outfile2.write(f"{i+1},{averages[i]:.4f},{stddev[i]:.4f},{prop_count_list[i]},{proportion_list[i]:.4f}\n")
#     outfile2.close()


# Do the plotting for TIME SERIES -- plot each model individually.
if args.plot_timeseries:
    for mod in range(nmodels):
        fig, ax = plt.subplots()
        ax.plot(np.arange(0, traj.n_frames), dist[:,mod], linewidth=0.5, label=f"Coil model: {mod+1}", color="blue")
        ax.axhline(y=np.average(dist[:,mod]), color="red", label=f"Avg dist: {np.average(dist[:,mod]):.3f}",
                   linewidth=0.5)
        plt.xlabel("Simulation frame no.")
        plt.ylabel("Distance (nm)")
        plt.ylim(0, np.max(dist[:,mod])+0.25)
        plt.legend()
        plt.title(f"ETE time series for {sim_code}\nCoil model: {mod +1}")
        plt.savefig(f"./ETE_plots/{sim_code}_model{mod+1}_ETE_timeseries.png", dpi=300)
        plt.close("all")

if args.plot_distribution:
    for mod in range(nmodels):
        #ete_kde = st.gaussian_kde(dist[:,mod])
        fig, ax = plt.subplots()
        ax.hist(dist[:, mod], density=True, bins="sqrt", color="blue", alpha=0.5)
        plt.xlabel("Distance (nm)")
        plt.xlim(0, np.max(dist[:,mod])+0.25)
        plt.ylabel("Frequency")
        plt.title(f"ETE distribution of distances for {sim_code}\nCoil model: {mod +1}")
        plt.savefig(f"./ETE_plots/{sim_code}_model{mod+1}_ETE_distribution.png", dpi=300)
        plt.close("all")

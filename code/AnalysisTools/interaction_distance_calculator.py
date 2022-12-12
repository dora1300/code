"""
@Title:         interaction_distance_calculator.py
@Author:        Mando Ramirez
@Date:          20220204

@Description:   This code is for calculating specific interaction distances in simulations, i.e. 1-4 and 1-5 distances.
This is designed to be used only with a single type of segment, either coil or linker, and not necessarily combinations
of segments. BUT it can handle multiple models of the same *type* of segment.

@ToDos:
20220204 - come back later and add a full output feature (frame averages for every frame, for every model)

@Updates:
20220209 - added more clarity to the figure file name.
"""

import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import numpy as np


""" Set up Argument Parser and parse the args! """
parser = argparse.ArgumentParser(description="Calculate the distances of 1-3, 1-4, and 1-5 interactions")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_rama")
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-f", help="Stride - analyze every f-th frame. Default is 1.", default=1, type=int)
# parser.add_argument("-full", help="Switch to save the full, frame by frame output for each model. Produces"
#                                   " a lot of data and large files. Provide flag to set to True.", action="store_true")

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially
output = args.o + ".csv"


""" Load the trajectory! """
traj = md.load(trajectory, top=topology)



""" Set up the list of indices for 1-3, 1-4, and 1-5 interactions """
i3_indices = []; i4_indices = []; i5_indices = []
for i in range(lhelix-2):
    i3_indices.append([i, i+2])
for i in range(lhelix-3):
    i4_indices.append([i, i+3])
for i in range(lhelix-4):
    i5_indices.append([i, i+4])



""" Main distance calculation section """
# # if I want to save the full, complete output then I will initialize special lists just for this purpose. Each entry is
# # itself a list of the frame averages, for all the frames
# if args.full:
#     full_model_i3 = []; full_model_i4 = []; full_model_i5 = []

# These lists are for holding the final averages and std devs for each model's ensemble totals. Each entry is a *SINGLE*
# number and the index corresponds to the model.
i3_tot_avg = []; i4_tot_avg = []; i5_tot_avg = []
i3_tot_std = []; i4_tot_std = []; i5_tot_std = []


for n in range(nmodels):
    # these lists are for storing the average and std for each individual *model*
    avg_i3 = []; std_i3 = [];  avg_i4 = []; std_i4 = []; avg_i5 = []; std_i5  = []

    # first, make sure you scale the atom pair indices for the correct model
    if n == 0:
        pass
    else:
        for i in range(len(i3_indices)):
            i3_indices[i][0] = i3_indices[i][0] + (n*lhelix)
            i3_indices[i][1] = i3_indices[i][1] + (n*lhelix)
        for i in range(len(i4_indices)):
            i4_indices[i][0] = i4_indices[i][0] + (n*lhelix)
            i4_indices[i][1] = i4_indices[i][1] + (n*lhelix)
        for i in range(len(i5_indices)):
            i5_indices[i][0] = i5_indices[i][0] + (n*lhelix)
            i5_indices[i][1] = i5_indices[i][1] + (n*lhelix)

    # now run the distance calculation for each frame desired. I will get a whole numpy array for each group of interaction
    # indices and will calculate a frame average and std. I will store the frame average and std into the provided lists
    for frame in range(0, traj.n_frames, args.f):
        i3_dist = md.compute_distances(traj[frame], np.array(i3_indices))
        i4_dist = md.compute_distances(traj[frame], np.array(i4_indices))
        i5_dist = md.compute_distances(traj[frame], np.array(i5_indices))
        avg_i3.append(np.average(i3_dist)); std_i3.append(np.std(i3_dist))
        avg_i4.append(np.average(i4_dist)); std_i4.append(np.std(i4_dist))
        avg_i5.append(np.average(i5_dist)); std_i5.append(np.std(i5_dist))
    # at this point, the avg_iN lists now have all the average distances by frame

    # Now make plots for each model, and have all three distances on the same plot
    fig, ax = plt.subplots()
    ax.plot(np.arange(0, traj.n_frames, args.f), avg_i3, color='blue', linestyle='-', linewidth=0.75, marker="",
            alpha=0.5, label=f"1-3: {np.average(avg_i3):.4f}")
    ax.plot(np.arange(0, traj.n_frames, args.f), avg_i4, color='tab:orange', linestyle='-', linewidth=0.75, marker="",
            alpha=0.5, label=f"1-4: {np.average(avg_i4):.4f}")
    ax.plot(np.arange(0, traj.n_frames, args.f), avg_i5, color='red', linestyle='-', linewidth=0.75, marker="",
            alpha=0.5, label=f"1-5: {np.average(avg_i5):.4f}")
    ax.axhline(y=0.545, color="black", linestyle=(0, (3, 5, 1, 5)), linewidth=0.75, label="1-3 Reference")
    ax.axhline(y=0.508, color="black", linestyle="--", linewidth=0.75, label="1-4 Reference")
    ax.axhline(y=0.630, color="black", linestyle="-", linewidth=0.75, label="1-5 Reference")
    plt.ylim(0, 1.25)
    plt.legend()
    plt.grid(axis="both", alpha=0.5, color="grey")
    plt.xlabel("Simulation frames")
    plt.ylabel("Distance (nm)")
    plt.title(f"{sim_code} - Model {n+1}\n Interaction distances analysis")
    plt.savefig(f"{sim_code}_model{n+1}_interactionDistancePlot.png", dpi=300)
    plt.close("all")

    # if args.full:
    #     full_model_i3.append(avg_i3)
    #     full_model_i4.append(avg_i4)
    #     full_model_i5.append(avg_i5)

    i3_tot_avg.append(np.average(avg_i3))
    i4_tot_avg.append(np.average(avg_i4))
    i5_tot_avg.append(np.average(avg_i5))
    i3_tot_std.append(np.std(avg_i3))
    i4_tot_std.append(np.std(avg_i4))
    i5_tot_std.append(np.std(avg_i5))



""" Now write out the finished file product! """
# if args.full:
#     fullout = open(output.split(".")[0] + "_full.csv", 'w')
#     fullout.write("Frame No.,")
#     for n in range(nmodels):
#         fullout.write(f"model{n+1}_1-3,model{n+1}_1-4,model{n+1}_1-5,")
#     # to be continued

# This is the shorter output, included with every analysis

outfile = open(output, 'w')
for n in range(nmodels):
    outfile.write(f"model{n+1}_1-3,avg,model{n+1}_1-3_stDev,model{n+1}_1-4,avg,model{n+1}_1-4_stDev,model{n+1}_1-5,avg,model{n+1}_1-5_stDev,")
outfile.write("\n")
results_str = ""
for i in range(len(i3_tot_avg)):
    results_str += f"{i3_tot_avg[i]},{i3_tot_std[i]},{i4_tot_avg[i]},{i4_tot_std[i]},{i5_tot_avg[i]},{i5_tot_std[i]},"
outfile.write(results_str)
outfile.close()

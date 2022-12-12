"""
@Title:         pitch_surrogate_analysis.py
@Author:        Mando Ramirez
@Date:          20211203

@Description:   This code calculates a metric that serves as a surrogate for helix pitch for the analysis of my coils.
This script calculates the i-i+3, and i-i+4 distances and averages them to calculate the "i-i+~3.5" distance. The
rationale for this is that in ideal coiled-coils, every turn i.e. every pitch contains 3.5 amino acids. So the distance
can give me a surrogate for the pitch. This metric will not be accurate if the helix loosens or tightens significantly.

Please note that although I use the word "pitch" throughout for shorthand, the calculated pitches are not actually the
pitch of the helix.

@Updates:
20211204 - added a feature to print the mean pitch surrogate for each model to stdout
"""

import mdtraj as md
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np
import os


"""
Set up of argument parser, and parsing of arguments
"""
parser = argparse.ArgumentParser(description="Coil pitch analysis using the i-i+3.5 surrogate method. "
                                             "This will also work for single PDB files if you provide the same PDB"
                                             " file for both the trajectory and topology flags.")
parser.add_argument("-t", help="trajectory file with extension (.xtc/.pdb)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-f", help="analyze every f'th frame in the trajectory", default=1, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_pitch_surrogate")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="switch to control plotting", default=False, type=bool)

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially

short_output = args.o + ".csv"
full_output = args.o + "_fullData.csv"

#open the output files now, and format them according to how I want
out1 = open(short_output, 'w')
out1.write(f"Model No.,Total_avg_coil_pitch,Total_std_coil_pitch\n")

out2 = open(full_output, 'w')
outline = "Frame,"
for i in range(nmodels):
    outline = outline + f"Model {nmodels+1} avg pitch,Model{nmodels+1} std of avg\n"
out2.write(outline)

"""
Load the trajectory
"""
traj = md.load(trajectory, top=topology)

# Remember, mdtraj 0-INDEXES the atoms.
"""
Main distance computing section of the code. 
"""
full_coil_data = []
full_coil_data_std = []
for n in range(nmodels):
    pitchpairs = []
    for atom in range(lhelix-4):        # this generates the different atom pairs, and stores them in the list pitchpairs
        i = atom + (n*lhelix)
        pitchpairs.append([i, i+3]) ; pitchpairs.append([i, i+4])
    pitch_array = np.array(pitchpairs)  # convert pitchpairs into an np array for use in mdtraj
    # Now loop through every frame in the trajectory, or however many frames I want to analyze. Calculate the distances
    # for the entire pitchpairs list for each single frame. Calculate the average coil pitch/distance for each single
    # frame, and save it in the list coil_dist_avgs_total
    coil_dist_avgs_total = [] ; coil_dist_avgs_stdev = []
    for frame in range(0, traj.n_frames, args.f):
        dist = md.compute_distances(traj[frame], pitch_array)
        # dist now contains the distances for every pair in the pitch_array. I have to average the distances for every
        # 2 distances in the list i.e. average dist[0]/dist[1] and dist[2]/dist[3] and so on, because these correspond
        # to the [i,i+3][i,i+4] I have in the pitch_array. Save each of the atom's averages into a new list, which should
        # correspond to half the size of pitchpairs
        atom_dist_avg = []
        for d in range(0, len(dist[0])-1, 2):
            atom_i_avg_dist = np.mean(np.array([dist[0][d], dist[0][d+1]]))
            atom_dist_avg.append(atom_i_avg_dist)
        # at this point, make sure I haven't done anything stupid and check that the size of coil_dist_avg is the same
        # as half of pitchpairs
        if len(atom_dist_avg) == (len(pitchpairs)/2):
            pass
        else:
            print("The number of atoms with contacts that were analyzed does not match the number of atoms that are "
                  "expected to be in the model. Exiting now.")
            exit()
        coil_dist_avgs_total.append(np.mean(np.array(atom_dist_avg)))
        coil_dist_avgs_stdev.append(np.std(np.array(atom_dist_avg)))

    # reminder, coil_dist_avgs_total contains, in each index, the averaged coil pitch for model n+1. The length of this
    # list is the same size as the number of frames in the range() interator
    # coil_dist_avgs_stdev contains, in each index, the standard deviation for the average coil pitch for model n+1.
    full_coil_data.append(coil_dist_avgs_total) ; full_coil_data_std.append(coil_dist_avgs_stdev)

    # This is the plotting portion. If I choose, I will plot the average coil pitch in each frame
    if args.plot:
        fig,ax = plt.subplots()
        ax.plot(np.arange(0, traj.n_frames, args.f), np.array(coil_dist_avgs_total), 'b-',
                linewidth=0.5, label=f"Model: {n+1}")
        ax.axhline(y=np.average(np.array(coil_dist_avgs_total)), color="red", linestyle="--",
                   linewidth=1.0, label=f"Avg {np.mean(np.array(coil_dist_avgs_total)):.4f}")
        ax.axhline(y=0.51, color="forestgreen", linestyle="--", label="Ref. coil", linewidth=1.0)
        ax.axhline(y=0.55, color="darkgreen", linestyle="--", label=r"Ref. $\alpha$-helix", linewidth=1.0)
        plt.ylim(0.45, 1.0)
        plt.xlabel("Simulation frames")
        plt.ylabel("Distance (nm)")
        plt.legend()
        plt.title(f"{sim_code} average coil pitch (surrogate) throughout time\nModel: {n+1}")
        plt.savefig(f"{sim_code}_avgCoilPitchPlot.png", dpi=300)
        plt.close("all")

    out1.write(f"{n+1},{np.mean(np.array(atom_dist_avg)):.4f},{np.std(np.array(atom_dist_avg)):.4f}\n")
    print(f"Model: {n+1}")
    print(f"Average pitch surr. (nm): {np.average(np.array(coil_dist_avgs_total)):.4f}")
    print(f"Stddev pitch surr. (nm): {np.average(np.array(coil_dist_avgs_stdev)):.4f}")
    print()

# Write out the data for the full output, which will be the total average and the total standard deviation of averages
# for each coil model for each frame that was analyzed
for j in range(len(full_coil_data[0])):
    fr = j * args.f
    writeline = f"{fr},"
    for k in range(len(full_coil_data)):
        writeline = writeline + f"{full_coil_data[k][j]},{full_coil_data_std[k][j]}"
    writeline += "\n"
    out2.write(writeline)


# Don't forget to close the output files!
out1.close()
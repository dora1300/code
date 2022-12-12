"""
@Title:             rmsd_analysis.py
@Name:              Mando A Ramirez
@Date:              2021 11 02

@Description:       This script calculates the RMSD for a given trajectory. It will save out the RMSD into a plot and
a separate graph. It assumes the reference frame is the frame provided in the topology file used to load the trajectory.
It can also compute the pairwise (all-to-all) RMSD and plot a heat map, and save a corresponding file.

@Updates:
20211105 - changed the argparse for the output file to remove the extension and have the code handle it instead

20220112 - fixed the plotting switch. Also adding a short output that contains only the average and stddev of the RMSD
data.

20220307 - updated to handle a reference!

20220318 - removed the nmodels flag since that information isn't used for the script
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description="RMSD analysis of Gromacs trajectories. Can do 1D and 2D RMSD analysis. "
                                             "This file is specifically for coil model simulations.\n"
                                             "Currently, this does not handle individual models within the system, only"
                                             " the system as a whole.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
#parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-reft", help="reference structure to do RMSD against. NEEDS to be a single structure file,"
                                  " not another trajectory.")
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file NO extension", default="output_rmsd")
parser.add_argument("-plot_timeseries", help="Invoke to turn on time series specific RMSD data. Only "
                                             "meaningful for 1D RMSDs", action="store_true")
parser.add_argument("-plot_distribution", help="Invoke to plot distributions of RMSD values for the trajectory. "
                                             "Only meaningful for 1D RMSDs", action="store_true")
parser.add_argument("-rmsd", help="(opt) switch to control what type of analysis to do. 1D is ALWAYS performed. 1 = "
                                  "only 1D RMSD. 2 = 1D and 2D RMSD.", default=1, type=int)


args = parser.parse_args()

trajectory= args.t
topology = args.p
reference = args.reft
#nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
output = args.o + ".csv"
rmsd_control = args.rmsd


# Load the trajectory and the reference frame
if args.reft is None:
    ref = md.load(topology)  # reference structure is the topology of the trajectory in this instance
    traj = md.load(trajectory, top=topology)  # sample trajectory to compute RMSD for
    # TODO: fix this because this isn't correct (should be first frame, not just the topology)
else:
    ref = md.load(reference)        # reference structure is the provided reference in this case
    traj = md.load(trajectory, top=topology)    # sample trajectory to compute RMSD for


# Step 1 -- align trajectory to the REFERENCE. This is very important because if I'm using something that is NOT from
# the trajectory, then I am absolutely going to introduce some weird artifacts into the RMSD. I need to align the trajectory
# to the reference, and then calculate RMSD.
aligned = md.Trajectory.superpose(traj, ref)


# Step 2 -- compute the 1D RMSD for the given trajectory. This ALWAYS occurs.
# the reference can be specified by user input.
#    Units are in nm
rmsd_1d = md.rmsd(aligned, ref)


# Compute the 2D RMSD but only if the switch is given to do so. Otherwise continue along the script like nothing happened.
if rmsd_control != 1:
    # this code for doing the pairwise RMSD calculation was found at: https://mdtraj.org/1.9.4/examples/clustering.html
    matrix2d = np.empty((aligned.n_frames, aligned.n_frames))
    for i in range(aligned.n_frames):
        matrix2d[i] = md.rmsd(aligned, aligned, i)
else:
    pass

# Plot the results for the 1D and 2D RMSD if it's provided
if args.plot_timeseries:
    fig, ax = plt.subplots()
    ax.plot(traj.time, rmsd_1d, 'b-', linewidth=0.35, label="1D RMSD")
    plt.xlabel("Time (ps)")
    plt.ylabel("Distance (nm)")
    plt.ylim(0, np.max(rmsd_1d)+0.1)
    plt.title(f"RMSD time series of {sim_code} trajectory in reference to \n{reference}")
    plt.savefig(f"{sim_code}_rmsd_timeseries.png", dpi=300)
    plt.close("all")

if args.plot_distribution:
    fig, ax = plt.subplots()
    ax.hist(rmsd_1d, density=True, bins="sqrt", color="blue", alpha=0.5)
    plt.xlabel("Distance (nm)")
    plt.ylabel("Frequency (counts)")
    plt.xlim(0, np.max(rmsd_1d)+0.5)
    plt.title(f"RMSD distribution of {sim_code} trajectory in reference to \n{reference}")
    plt.savefig(f"{sim_code}_rmsd_distribution.png", dpi=300)
    plt.close("all")

if rmsd_control != 1:
    fig, ax = plt.subplots()
    plt.imshow(matrix2d, cmap='summer', origin='lower')
    plt.xlabel("Frame")
    plt.ylabel("Frame")
    plt.title(f"Pairwise RMSD of {sim_code} trajectory")
    plt.colorbar(orientation='vertical', label=r'RMSD (nm)')
    plt.savefig(f"{sim_code}_2d_rmsdplot.png", dpi=600)
    plt.close('all')

# Now its time to save the RMSD data to a file for later use, if required. As of 20211103 I am not saving out the 2D
#   RMSD since that will be a large file and will be kinda crazy to output. And I'm not sure I'd need that data for
#   any additional analysis.
with open(output, 'w') as outfile:
    outfile.write("Frame no., Time (ps), RMSD (nm)\n")
    for i in range(len(rmsd_1d)):
        outfile.write(f"{i},{traj.time[i]},{rmsd_1d[i]:.4f}\n")

avg_rmsd = np.average(rmsd_1d)
std_rmsd = np.std(rmsd_1d)
with open(args.o + "_short.csv", 'w') as shortout:
    shortout.write("Avg RMSD,Stddev RMSD\n")
    shortout.write(f"{avg_rmsd:.3f},{std_rmsd:.3f}\n")


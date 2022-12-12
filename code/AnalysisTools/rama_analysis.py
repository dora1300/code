"""
@Title:         rama_analysis.py
@Author:        Mando Ramirez
@Date:          20220207

@Description:   This code is for analyzing trajectories of coarse grained coil or linker simulations to calculate all the
 pseudo-dihedrals and angles for each model in the trajectory. Importantly, the models must all be the same length.
 This code is an offshoot of the code "coil_rama_analysis.py".
 I am removing the assessment of "deviation" from reference, because that type of analysis is better handled by the KDE
 script I wrote.

@Updates:
20220214 - Script now outputs the statistics for each model (i.e. average and median) into an output file for handling
the data elsewhere

20220301 - fixed the reference to reflect the ACTUAL correct angle and torsion reference values see DAR3-15
"""

import mdtraj as md
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(description="Ramachandran analysis of coiled-coil coarse-grained models.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_rama")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="(1 = only global plots) or (2 = all plots)", type=int, default=0)

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially
output = args.o + ".csv"

## Added 2021 11 15 - changed the dihedral and angle bounds, and made a generalized variable to make updating it easier
dih_ref = 50
angle_ref = 93
dih_upp_bound = dih_ref + 30
dih_low_bound = dih_ref - 30
angle_upp_bound = angle_ref + 30
angle_low_bound = angle_ref - 30

# Adding a section to check if there is a directory to save my plots. If there is, it will continue the script, if not
# then it will make the directory for me. This will only happen if the plot argument is passed as True
# added 2021 11 11
if args.plot == 1 or args.plot == 2:
    if os.path.exists("./rama_plots"):
        pass
    else:
        os.mkdir("./rama_plots")

if args.plot > 2:
    print(f"You gave the plotting switch code: {args.plot}")
    print("The plotting swtich code is invalid. Please try again.")
    exit(0)

# load the trajectory and gather bond
traj = md.load(trajectory, top=topology)

# Analyze the pseudodihedrals and pseudobond angles for each dihedral/angle pair by looping through the atoms
#   mdtraj 0-indexes the atoms, and my current definitions of angles are:
#       pseudodihedral = n-1 -- n -- n+1 -- n+2
#       pseudoangle    = n-1 -- n -- n+1
#   so I will start my range at 1 and go to (lhelix-2) inclusive, as this will give all the dihedrals and angles for a single coil
#   I will also include a model number term to know how many coil models and thus how many angles I need to analyze.

for n in range(1, nmodels+1):
    """ These next two lists are important for storing EVERY dihedral and bond angle angle for each individual model.
    But importantly they only store the angles for a single model. They get refreshed for every model."""
    total_dihedrals = []  # added 2021 11 23 - useful for storing the raw dihedrals
    total_bondangles = []  # added 2021 11 23 - useful for storing the raw angles

    for j in range(1, (lhelix-2)):
        i = j + (lhelix*(n-1))      #seems stupid, but this is my scaling factor for the number of models
                                # i.e. if I'm on my first model n=1 and j = 1, then i = 1 + (35 * 1-1) == i = 1 and
                                # the resutling dih_coords are [0, 1, 2, 3] which is correct for the start of model 1.
                                # but if I'm on my 2nd model n=2 and j=1, then i=1+ (35*(2-1)) == i = 36 which means the
                                # resulting dih_coords are [35, 36, 37, 38] and this is correct because for model 2 the
                                # starting atom is 35 BECAUSE MDTRAJ 0-INDEXES THE ATOMS
        dih_coords = [i-1, i, i+1, i+2]
        ang_coords = [i-1, i, i+1]

        # Compute the dihedrals for the given atom coordinates, grab the 0 column (0th or 1st, doesn't matter) and convert
        #    to degrees. Do the same for angles. These both produce numpy arrays
        dihedral_deg = md.compute_dihedrals(traj, [dih_coords, dih_coords])[:, 0] * (180 / math.pi)
        angle_deg = md.compute_angles(traj, [ang_coords, ang_coords])[:, 0] * (180 / math.pi)

        if args.plot == 2:
            # if the plotting switch is turned to True (default is False) then it will make plots for ALL the dihedral/angle
            #    pairs
            # Generate a plot for this specific dihedral/angle pair
            fig, ax = plt.subplots()
            plt.scatter(dihedral_deg, angle_deg, marker='x', s=0.75, c=traj.time)
            cbar = plt.colorbar()
            cbar.set_label("Time [ps]")
            plt.axvline(dih_ref, color='b', linestyle='--')
            plt.axhline(angle_ref, color='b', linestyle='--')
            plt.xlabel(r"Pseudo-dihedral $\alpha$ (deg)")
            plt.ylabel(r"Pseudo-bond angle $\theta$ (deg)")
            plt.xlim(-180, 180)
            plt.ylim(-180, 180)
            plt.title(f"{sim_code} angle plot. Model No: {str(n)}.\nDihedral: {str(dih_coords)}. Angle: {str(ang_coords)}")
            plt.savefig(f"./rama_plots/{sim_code}_angle_plot_{i}.png", dpi=300)
            plt.close("all")
        else:
            pass

        # Added 2021 11 23
        # Save all the raw dihedrals and bond angles for making a graph and plotting
        for index in range(len(dihedral_deg)):
            total_dihedrals.append(dihedral_deg[index])
            total_bondangles.append(angle_deg[index])


    # at this point, all the angles for a single model are stored in the appropriate lists. Now I will calculate stats
    # on them
    avg_dih = np.average(np.array(total_dihedrals))
    med_dih = np.median(np.array(total_dihedrals))
    std_dih = np.std(np.array(total_dihedrals))

    avg_ba = np.average(np.array(total_bondangles))
    med_ba = np.median(np.array(total_bondangles))
    std_ba = np.std(np.array(total_bondangles))


    # Make an output file to store in the information about each model and the angle statistics
    with open(output, 'w') as f:
        f.write("Model,Median Torsion,Average Torsion,StDev Torsion,Median Bond Angles,Average Bond Angles,StDev Bond Angles\n")
        f.write(f"{n},{med_dih:.3f},{avg_dih:.3f},{std_dih:.3f},{med_ba:.3f},{avg_ba:.3f},{std_ba:.3f}\n")

    # Added 2021 11 23
    # This is a plot just of the distribution of angles, *not* of the deviation of angles.
    # This needs to be handled for each model
    if args.plot == 1 or args.plot == 2:
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(15,15))
        ax1.hist(np.array(total_dihedrals), bins='sqrt', density=True, color='cornflowerblue', alpha=0.75)
        ax1.axvline(x=dih_ref, color='darkviolet',linestyle='--', label=f"Reference value: {dih_ref}")
        ax1.axvline(x=avg_dih, color='blue', linestyle='--',
                    label=f"Average value: {avg_dih:.2f}")
        ax1.axvline(x=med_dih, color='darkgreen', linestyle='--',
                    label=f"Median value: {med_dih:.2f}")
        ax1.set_xlim(-180, 180)
        ax1.legend()
        ax1.set_xlabel(r"Absolute $\alpha$ (deg)")
        ax1.set_ylabel("Counts")
        ax1.set_title(rf"Pseudo-dihedral $\alpha$. {sim_code}, Model No: {str(n)}, ALL pseudo-dihedrals")

        ax2.hist(np.array(total_bondangles), bins='sqrt', density=True, color='lightcoral', alpha=0.75)
        ax2.axvline(x=angle_ref, color='darkviolet',linestyle='--', label=f"Reference value: {angle_ref}")
        ax2.axvline(x=avg_ba, color='blue', linestyle='--',
                    label=f"Average value: {avg_ba:.2f}")
        ax2.axvline(x=med_ba, color='darkgreen', linestyle='--',
                    label=f"Median value: {med_ba:.2f}")
        ax2.set_xlim(-180, 180)
        ax2.legend()
        ax2.set_xlabel(r"Absolute $\theta$ (deg)")
        ax2.set_ylabel("Counts")
        ax2.set_title(rf"Pseudo-bondangle $\theta$. {sim_code}, Model No: {str(n)}, ALL pseudo-bond angles")
        plt.savefig(f"{sim_code}_model{n}_global_angle_histogram_plot.png", dpi=300)
        plt.close()

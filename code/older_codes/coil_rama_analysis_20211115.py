"""
@Title:         coil_rama_analysis_20211115.py
@Author:        Mando Ramirez
@Date:          20211021

@Description:   This code is for analyzing trajectories of coarse grained coil simulations to calculate all the
 pseudo-dihedrals and angles for each coil model in the trajectory.

@Updates:
20211028 - added a plot switch to turn of plotting if I want (there's a lot of things to plot otherwise). Also added
a feature to automatically calculate the average proportion for each helix plus the std. deviation of each helix proportion

20211103 - clarified the helical checking portion based on the permitted Ramachandran region from Tozzini 2006 paper.
This should capture more accurate permitted flexibility, and provides additional justification.

20211104 - removed the helical determination zone from yesterday, and replaced it with a simplified rectangle version
that is at least based in part on something meaningful. Also changing the helix proportion cut off to 0.75

20211105 - added another metric which measures the deviation of a given (alpha, theta) pair from the idealized helix values
and plots the deviations in a histogram. Also removed the extension on the output from the argparse so
that the code can handle that itself.

20211110 - changed the name to "coil_rama_analysis.py" so that I can use this specifically for analyzing coil rama
angles

20211112 - added a feature to add the rama_plots directory if it doesn't already exist
"""

import mdtraj as md
import matplotlib.pyplot as plt
import math
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(description="Ramachandran analysis of coiled models to determine helicity.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_rama")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="switch to control plotting. Default is False.", type=bool, default=False)

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially
output = args.o + ".txt"

global_dih_devs = []            # this is useful for storing ALL the dih/angle deviations for a final plotting
global_angle_devs = []          # this is useful for storing ALL the dih/angle deviations for a final plotting
total_angle_count = 0           # this and the following three are used for tracking angle statistics
not_helix_count = 0
not_dihedral_count = 0
not_bondangle_count = 0

# Adding a section to check if there is a directory to save my plots. If there is, it will continue the script, if not
# then it will make the directory for me. This will only happen if the plot argument is passsed as True
# added 2021 11 11
if args.plot:
    if os.path.exists("./rama_plots"):
        pass
    else:
        os.mkdir("./rama_plots")


# load the trajectory and gather bond
traj = md.load(trajectory, top=topology)

# Analyze the pseudodihedrals and pseudobond angles for each dihedral/angle pair. by looping through the atoms
#   mdtraj 0-indexes the atoms, and my current definitions of angles are:
#       pseudodihedral = n-1 -- n -- n+1 -- n+2
#       pseudoangle    = n-1 -- n -- n+1
#   so I will start my range at 1 and go to (lhelix-2) inclusive, as this will give all the dihedrals and angles for a single coil
#   I will also include a model number term to know how many coil models and thus how many angles I need to analyze.

model_proportions = []      # I will save the proportion helix for each dih/ang pair in this list. In this list will be
                            # a list for each model, and the specific dih/ang pairs will be put into the corresponding
                            # model list i.e. if there are two models then
                            #       model_proportions = [ [model1_dihang_props], [model2_dihang_props] ]
for n in range(1, nmodels+1):
    model_proportions.append([])
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

        if args.plot:
            # if the plotting switch is turned to True (default is False) then it will make plots for ALL the dihedral/angle
            #    pairs
            # Generate a plot for this specific dihedral/angle pair
            fig, ax = plt.subplots()
            plt.scatter(dihedral_deg, angle_deg, marker='x', s=0.75, c=traj.time)
            cbar = plt.colorbar()
            cbar.set_label("Time [ps]")
            plt.hlines(y=81, xmin=25, xmax=104, color='r')
            plt.hlines(y=120, xmin=25, xmax=104, color='r')
            plt.vlines(x=25, ymin=81, ymax=120, color='r')
            plt.vlines(x=104, ymin=81, ymax=120, color='r')
            plt.axvline(52, color='b', linestyle='--')
            plt.axhline(92, color='b', linestyle='--')
            plt.xlabel(r"Pseudo-dihedral $\alpha$ (deg)")
            plt.ylabel(r"Pseudo-bond angle $\theta$ (deg)")
            plt.xlim(-180, 180)
            plt.ylim(-180, 180)
            plt.title(f"{sim_code} angle plot. Model No: {str(n)}.\nDihedral: {str(dih_coords)}. Angle: {str(ang_coords)}")
            plt.savefig(f"./rama_plots/{sim_code}_angle_plot_{i}.png", dpi=300)
            plt.close("all")
        else:
            pass

        phelix = []             # this is where I will store the proportion information if the dihedral/angle is a helix
                                # or not. The index for each entry and be used to get the correspond time from traj.time


        # Loop through the dihedral and angle data and figure out if and when this specific angle breaks out of the
        #    acceptable helix range. Record the proportion of frames that the dihedral is within the acceptable range
        for k in range(len(dihedral_deg)):
            total_angle_count += 1
            if (dihedral_deg[k] >= 25 and dihedral_deg[k] < 104) and (angle_deg[k] >= 81 and angle_deg[k] <= 120):
                # This is the acceptable tolerance range for the given angles for me to still call it a helix.
                phelix.append(1)
            elif (dihedral_deg[k] >= 25 and dihedral_deg[k] < 104) and (angle_deg[k] < 81 or angle_deg[k] > 120):
                # this is when the dihedral is in an acceptable range, but the angle is not
                not_bondangle_count += 1
            elif (dihedral_deg[k] < 25 or dihedral_deg[k] > 104) and (angle_deg[k] >= 81 and angle_deg[k] <= 120):
                # this is when the bond angle is in an acceptable range, but the dihedral is not
                not_dihedral_count += 1
            else:
                not_helix_count += 1
                phelix.append(0)


        # Added 2021 11 05
        # Calculate the deviation of the dihedral/angle bonds from the idealized helix pair and plot the resultant
        #    histogram
        if args.plot:
            dihedral_devs = []
            angle_devs = []
            for l in range(len(dihedral_deg)):
                dihedral_devs.append(round(dihedral_deg[l] - 52, 2))
                global_dih_devs.append(round(dihedral_deg[l] - 52, 2))
                angle_devs.append(round(angle_deg[l] - 92, 2))
                global_angle_devs.append(round(angle_deg[l] - 92, 2))
            fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(15,15))
            ax1.hist(np.array(dihedral_devs), bins='sqrt', density=True, color='cornflowerblue')
            ax1.axvline(x=(104-52), color='k', linestyle='--')
            ax1.axvline(x=(25-52), color='k', linestyle='--')
            ax1.set_xlim(-250, 130)
            ax1.set_xlabel(r"Absolute deviation of $\alpha$ (deg)")
            ax1.set_ylabel("Counts")
            ax1.set_title(rf"Deviation of dihedral $\alpha$. {sim_code}, Model No: {str(n)}, Dihedral: {str(dih_coords)}")

            ax2.hist(np.array(angle_devs), bins='sqrt', density=True, color='lightcoral')
            ax2.axvline(x=(120-92), color='k', linestyle='--')
            ax2.axvline(x=(81-92), color='k', linestyle='--')
            ax2.set_xlim(-280, 90)
            ax2.set_xlabel(r"Absolute deviation of $\theta$ (deg)")
            ax2.set_ylabel("Counts")
            ax2.set_title(rf"Deviation of angle $\theta$. {sim_code}, Model No: {str(n)}, Angle: {str(ang_coords)}")
            plt.savefig(f"./rama_plots/{sim_code}_angle_deviation_plot_{i}.png", dpi=300)
            plt.close("all")

        proportion_helix = round(sum(phelix)/len(phelix), 4)
        model_proportions[n-1].append(proportion_helix)


# If the Plot = True, then I will make a final global histogram deviation plot for deviations of dihedral and angles
# for EVERY single dih/angle pair in the simulation.
# added 2021 11 05
if args.plot:
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(15, 15))
    ax1.hist(np.array(global_dih_devs), bins='sqrt', density=True, color='cornflowerblue')
    ax1.axvline(x=(104 - 52), color='k', linestyle='--')
    ax1.axvline(x=(25 - 52), color='k', linestyle='--')
    ax1.set_xlim(-250, 130)
    ax1.set_xlabel(r"Absolute deviation of $\alpha$ (deg)")
    ax1.set_ylabel("Counts")
    ax1.set_title(rf"Deviation of dihedral $\alpha$. {sim_code}, Model No: {str(n)}, ALL pseudo-dihedrals")

    ax2.hist(np.array(global_angle_devs), bins='sqrt', density=True, color='lightcoral')
    ax2.axvline(x=(120 - 92), color='k', linestyle='--')
    ax2.axvline(x=(81 - 92), color='k', linestyle='--')
    ax2.set_xlim(-280, 90)
    ax2.set_xlabel(r"Absolute deviation of $\theta$ (deg)")
    ax2.set_ylabel("Counts")
    ax2.set_title(rf"Deviation of angle $\theta$. {sim_code}, Model No: {str(n)}, ALL pseudo-bond angles")
    plt.savefig(f"{sim_code}_global_angle_deviation_plot.png", dpi=300)
    plt.close("all")

# open the output file
outfile = open(output, 'w')
outfile.write("Ramachandran analysis of coiled-coil model. Calculating proportion of dihedral/bond angles in the helix range.\n")
#    add a second results file for abbreviated results, instead of the full rama results that I otherwise plot.
#    added 20211028
outfile_split = output.split(".")
outfile2_name = outfile_split[0] + "_abbrev.csv"
outfile2 = open(outfile2_name, 'w')
outfile2.write("Total angle pairs analyzed,No. angles not helix,No. dihedral not in range,No. bas not in range,"
               "Fraction angles not helix,Fraction dihedral not in range,Fraction bas not in range\n")
frac_not_helix = not_helix_count / total_angle_count
frac_not_dihedral = not_dihedral_count / total_angle_count
frac_not_ba = not_bondangle_count / total_angle_count
outfile2.write(f"{total_angle_count},{not_helix_count},{not_dihedral_count},{not_bondangle_count},{frac_not_helix},"
               f"{frac_not_dihedral},{frac_not_ba}\n")
outfile2.write("Helix no.,helix_prop_avg,helix_prop_std,helix_prop_min,helix_prop_max,no_helices_as_helix\n")

# Save out the data for each dih/ang pair for each model in an output file for later analysis/processing.
for i in range(len(model_proportions)):
    outfile.write(f"MODEL: {i+1}\n")
    outfile.write(f"Dihedral\t\t\tAngle\t\t\tProportion of simulation in helix\n")
    count = 1
    index = count + (35 * i)    # for the same reasons above, this helps me make sure my atoms line up with my models
    helix_count = 0
    for j in range(len(model_proportions[i])):
        outfile.write(f"{index-1}--{index}--{index+1}--{index+2}\t\t")
        outfile.write(f"{index-1}--{index}--{index+1}\t\t")
        outfile.write(f"{model_proportions[i][j]}\n")
        index += 1
        if model_proportions[i][j] >= 0.75:
            helix_count += 1
        else:
            pass

    model_proportion_avg = np.average(np.array(model_proportions[i]))           # average & std calc added 20211028
    model_proportion_std = np.std(np.array(model_proportions[i]))
    outfile2.write(f"{i+1},{model_proportion_avg},{model_proportion_std},{min(model_proportions[i])},"
                   f"{max(model_proportions[i])},{helix_count}\n")

    outfile.write(f"Minimum helix p: {min(model_proportions[i])}\n")
    outfile.write(f"Maximum helix p: {max(model_proportions[i])}\n")
    outfile.write(f"Angle pairs with proportion > 0.75: {helix_count}\n")            # fixed from 0.5 20211104
    outfile.write(f"Range of helix proportions: {round(max(model_proportions[i]) - min(model_proportions[i]),4)}\n")
    outfile.write(f"Total frames considered for each dihedral/angle pair: {traj.n_frames}\n")

# Close the output
outfile.close()
outfile2.close()

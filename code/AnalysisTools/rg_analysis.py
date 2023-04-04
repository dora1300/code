"""
@Title:         rg_analysis.py
@Author:        Mando Ramirez
@Date:          20211021

@Description:   This code is for analyzing trajectories of coarse grained coil simulations to calculate the radius of
gyration for each simulation. It will handle

@Updates:
20211103 - updated to normalize some data to the first frame to show change in Rg instead of raw numbers

20211105 - fixed the argparse to remove the extension in the output file, and changed it so that the code handles that itself

20211110 - added the length option to the arg parse. It's not needed, but I have it here so that I can use the same
options for all my analysis tools essentially. This is because I'm lazy

20220103 - revised the output writing so that the short version is saved as a .csv instead
of .txt, which will make parsing much, much easier

20220118 - changed the plotting switch to be more accurate pythonese

20220222 - updating the script to handle model level Rg versus system level Rg.

20220425 - added a arg parse feature to control the starting FRAME to analyze from. This is useful
for ignoring the first nth number of frames if you know the system is not correlated at that point.
"""

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="Radius of gyration analysis for coiled-coil simulations.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", type=int)
parser.add_argument("-l", help="the length (no. atoms) in the linker. Each linker must be the same length.",
                    type=int)
parser.add_argument("-sys", help="SWITCH to control for calculating the Rg of the entire system, INSTEAD of "
                                 "individual models. Provide flag to turn ON system RG. This "
                                 "overwrites the model level specified by -n and -l", action="store_true")
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_rgyr")
parser.add_argument("-plot", help="switch for plotting. Provide the flag to set plotting to True",
                    action="store_true")
parser.add_argument("-startF", help="the STARTING FRAME to begin the analysis from. This is useful to "
                                "ignore the beginning frames if the system is not in equilibration",
                                type=int, default=1)


args = parser.parse_args()

# Determine if there is system level Rg analysis or if I need to check for model level details
if args.sys is None:
    args.sys = False
if (args.n is None or args.l is None) and not args.sys:
    print("If the system Rg is not specified, then you need to provide -n and -l flags. Exiting now.")
    exit()

trajectory= args.t
topology = args.p
sim_code = args.id

if args.plot:
    if os.path.exists("./Rg_plots"):
        pass
    else:
        os.mkdir("./Rg_plots")


""" Do this part only is system level Rg is specified (arg.sys = True)"""
if args.sys:
    # load the trajectory and gather bond
    traj = md.load(trajectory, top=topology)

    # calculate the radius of gyration for the entire simulation. Very simple. I don't think I can do the Rg of just one
    # coil but then again I'm not sure I need that, so whole system will suffice. Also compute the AUC
    rgArray = md.compute_rg(traj[args.startF:traj.n_frames:1])
    rg_average = np.average(rgArray)
    rg_stddev = np.std(rgArray)

    auc = np.trapz(rgArray, np.arange(args.startF, traj.n_frames,1))
                        # the data is integrated based on the x-axis, which in this case will be the number of frames

    # plot the time course Rg and save the figure. This is based on the args.plot switch. This will also plot the normalized
    #    Rg plot.
    if args.plot:
        fig, ax = plt.subplots()
        ax.plot(traj[args.startF:traj.n_frames:1].time, rgArray, 'r-', label="Rg", linewidth=0.5)
        plt.xlabel("Time (ps)")
        plt.ylabel("Radius of gyration (nm)")
        plt.xlim(0, np.amax(traj[args.startF:traj.n_frames:1].time))
        plt.ylim(0, np.amax(rgArray)+1)
        plt.title(f"Radius of gyration of simulation: {sim_code}")
        plt.savefig(f"./Rg_plots/{sim_code}_system_Rg_plot.png", dpi=300)
    else:
        pass


    # open the output files
    output = args.o + "_system.csv"
    shortout = open(output, 'w')
    output_split = output.split(".")
    fullname = f"{output_split[0]}_fullData.csv"
    fullout = open(fullname, 'w')

    # Save the ensemble averaged data into a
    shortout.write("Time averaged Rg calculations,\n")
    shortout.write(f"Average Rg (nm),{round(rg_average,4)}\n")
    shortout.write(f"StDev of Rg,{round(rg_stddev, 4)}\n")
    shortout.write(f"AUC of Rg (nm*ps),{round(auc,4)}\n")
    shortout.write(f"Number of frames,{traj.n_frames-args.startF}\n")

    shortout.close()

    fullout.write("Radius of gyration analysis of coiled coil model. Full trajectory Rg calculations.\n")
    fullout.write("Time (ps),Rg (nm)\n")
    for i in range(len(rgArray)):
        fullout.write(f"{traj[args.startF:traj.n_frames:1].time[i]},{rgArray[i]:.4f}\n")
    fullout.close()
else:
    """ Execute this code if model level Rg is specified i.e. args.sys = False """
    nmodels = args.n  # this is 1-indexed, essentially
    lhelix = args.l   # This is also 1-indexed
    # load the trajectory
    traj = md.load(trajectory, top=topology)
    atom_indices = [0, lhelix-1]

    all_model_Rgs = []      # This holds every models Rgs by frame
    model_avg_rgs = []      # this stores the AVERAGE Rg for each model i.e. each entry is for the respective models
    model_std_rgs = []      # this stores the StDev for each model

    for n in range(nmodels):
        # This is the main loop to analyze Rg based on different models. This will work by subsetting trajectories from
        # the main trajectory and analyzing the subsets, and then saving those data in a file
        if n > 0:
            atom_indices[0] += lhelix; atom_indices[1] += lhelix

        subset_traj = traj.atom_slice(np.arange(atom_indices[0], atom_indices[1]+1))
        model_rg = md.compute_rg(subset_traj[args.startF:subset_traj.n_frames:1])
        all_model_Rgs.append(model_rg)
        model_avg_rgs.append(np.average(model_rg))
        model_std_rgs.append(np.std(model_rg))

        if args.plot:
            # plot the time course Rg and save the figure. This is based on the args.plot switch. This will also plot the normalized
            #    Rg plot.
            if args.plot:
                fig, ax = plt.subplots()
                ax.plot(subset_traj[args.startF:subset_traj.n_frames:1].time / 1000, model_rg, 'b-', label="Rg", linewidth=0.5)
                plt.xlabel("Time (us)")
                plt.ylabel("Radius of gyration (nm)")
                plt.xlim(0, np.amax(subset_traj[args.startF:subset_traj.n_frames:1]/1000))
                plt.ylim(0, np.amax(model_rg) + 1)
                plt.title(f"Radius of gyration of simulation: {sim_code}")
                plt.savefig(f"./Rg_plots/{sim_code}_model_{n+1}_Rg_plot.png", dpi=300)
            else:
                pass

    # Now save the output
    # open the output files
    output = args.o + "_byModel.csv"
    shortout = open(output, 'w')
    output_split = output.split(".")
    fullname = f"{output_split[0]}_fullData.csv"
    fullout = open(fullname, 'w')

    shortout.write(f"{sim_code}\n")
    shortout.write("Model,Average Rg,StDev Rg\n")
    for nmod in range(len(model_avg_rgs)):
        shortout.write(f"{nmod+1},{model_avg_rgs[nmod]},{model_std_rgs[nmod]}\n")

    fullout.write(f"{sim_code}\n")
    fullout.write(f"Time (us)")
    for ind in range(len(all_model_Rgs)):
        fullout.write(f",Model_{ind+1}")
    fullout.write("\n")

    for i, time in enumerate(traj[args.startF:subset_traj.n_frames:1].time):
        fullout.write(f"{(time / 1000):.4f}")
        for model in range(len(all_model_Rgs)):
            fullout.write(f",{all_model_Rgs[model][i]}")
        fullout.write("\n")

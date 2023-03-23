"""
@Title:             coil_oligomerization analysis.py
@Name:              Mando A Ramirez
@Date:              2022 06 19

@Description:       This script analyzes the oligomeric state of coil simulations to produce statistics about the
different types of oligomers and higher order assemblies that form with coils in a box. I am writing this script with
slab simulations in mind so I can track the types of oligomers that form through simulation time. I will be calcualte
the following:
    - the types of each oligomers through time (per-frame)
    - the distribution of each oligomers for the whole simulation
    - the number-per-frame and distribution of self-vs-other interactions
    - the number-per-frame and distribution of free coils for the simulation

This script will first be written to handle only homomeric simulations of models, but can easily be modified to handle
heteromeric simulations.

@Updates:
2022 09 28 - added a feature to plot the distribution of oligomer sizes over time, and added a feature to choose the
starting time for the simulation or at least as close to it as I can.

@Updates:
2023 02 21 - added at timing feature so I can keep track of how long this analysis takes.

@Updates:
2023 03 23 - added a feature to handle analyzing just a single frame i.e. pdb file.
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import time


""" Set up the arg parser """
parser = argparse.ArgumentParser(description="Coil oligomerization tool - this determines the types of oligomers and"
                                             " higher order assemblies present in your coil simulations!")
parser.add_argument("-df", help="Definition file: this contains the information on how to build the model, with"
                                " alternating coil and linker segments.")
parser.add_argument("-nmodels", help="The total number of models in the simulation", required=True, type=int)
parser.add_argument("-N", help="The total number of coil segments in the simulation. Necessary for the distribution"
                               " plotting.", required=True, type=int)
parser.add_argument("-l", help="No. of beads in the ENTIRE MODEL. This is not the size of an individual coil. This is "
                               "the size of an entire individual model. As of right now, ONLY models that have the"
                               " same no. of beads can be analyzed.",
                    required=True, type=int)
parser.add_argument("-cutoff", help="The cut-off distance used to determine if two beads are close enough to be in "
                                    "an oligomer. [nm]", default=1.25, type=float)
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-start", help="the starting frame (in ps) to begin the analysis. Default = 0 [ps]", default=0,
                    type=int)
parser.add_argument("-f", help="the stride of analysis i.e. analyze every fth frame, default = 1", default=1, type=int)
parser.add_argument("-single", help="a switch to turn on if the analysis is only for a single structure file",
    action="store_true", default=False)
parser.add_argument("-name", help="name that you'd like to add to the analysis output files", default="")
parser.add_argument("-timeseries", help="pass the flag to turn on plotting of oligomer data as a time series",
                    action="store_true", default=False)
parser.add_argument("-distro", help="pass this flag to turn on plotting of average oligomer population, as a "
                                    "distribution", action="store_true", default=False)
parser.add_argument("-verbose", help="pass this flag to output extra information to stdout.", action="store_true",
                    default=False)

args = parser.parse_args()

trajectory= args.t
topology = args.p

""" Print a warning statement """
print("**********")
print("PLEASE BE ADVISED: This code automatically converts the information provided in the -df file into 0-INDEXED "
      "positions! Your input is 1-indexed, but this code will automatically change it to 0-indexed."
      ""
      "Also note that this code only needs to run once, since it outputs the data and then you can modify plots of that"
      "data without doing the entire analysis all over again.")
print("**********\n")

""" Begin the timing procedure, added 20230221 """
time_start = time.perf_counter()

# load the trajectory!
if args.single:
    traj = md.load(trajectory, top=topology)
    args.f = 1
else:
    traj_load = md.load(trajectory, top=topology)
    frame_start = 0
    while frame_start <= traj_load.n_frames:
        if int(frame_start * traj_load.timestep) >= args.start:
            break
        else:
            frame_start += 1
    traj = traj_load[frame_start::]
    remaining_frames = traj_load.n_frames - frame_start
    analyzed_frames = remaining_frames / args.f

"""
Useful global constants
"""
A_CUTOFF = args.cutoff                 # in nanometers!
NMODELS = args.nmodels
MSIZE = args.l

INTERACTION_CUTOFF = 0.5       # unitless. This says that 50% of beads between any two coils need to be interacting
                                # for there to be an oligomer



"""
Step 1 - read the file about the model description to find where all the coils and A beads are.
code taken from other things like model_maker.py for file reading
"""
# Read the file and save only the information about the coils. Save the start/end indices for the coils and the start
# of the A-bead positions.
# POSITIONS ARE 0-INDEXED
Coil_positions_df = []
coil_counter = 0                # this is how many coils are present in a model
with open(args.df, "r") as file:
    for line in file:
        line_split = line.rstrip("\n").split(",")
        if str(line_split[0]) == "coil":
            # information in line_split = [coil_start, coil_end, a_bead_start]
            # information originally 1-indexed
            # this information needs to be 0-indexed, so I'm converting it!!
            coil_counter += 1
            Coil_positions_df.append([int(line_split[1])-1, int(line_split[2])-1, int(line_split[3])-1])
        else:
            pass

print(f"Coil model file: {args.df}")
print(f"Parsed coil positions (0-indexed): {Coil_positions_df}")

# Now, for each coil, generate the A- and D- position A-beads FOR THE FIRST MODEL ONLY. I need to enumerate all the
# a-bead positions first before I can generate all the a-beads for all the models
# ALL POSITIONS ARE 0-INDEXED
Coil1_A_beads = []
for i in range(len(Coil_positions_df)):
    segStart = Coil_positions_df[i][0]
    segStop = Coil_positions_df[i][1]
    a_posits = []
    d_posits = []
    a_index = Coil_positions_df[i][2]
    d_index = Coil_positions_df[i][2] + 3
    while a_index <= segStop:
        a_posits.append(a_index)
        a_index += 7
    while d_index <= segStop:
        d_posits.append(d_index)
        d_index += 7
    coil_i_beads = a_posits + d_posits
    coil_i_beads.sort()
    Coil1_A_beads.append(coil_i_beads)

print(f"Enumerated coil beads: {Coil1_A_beads}")
print()

# Now enumerate all the a-bead indices for each model present in the simulation. The variable MODELS_A_BEADS will have
# the 0-BASED INDICES (!!!!) for every coil model in the system. THIS IS FOR ALL COILS
MODELS_A_BEADS = []
MODEL_RANGES = []
for M in range(0, NMODELS):
    # This information stores the 0-INDEXED (!!!) ranges for each of the models in the system!
    model_additive = M * MSIZE
    model_start = 0 + (M * MSIZE)
    model_end = MSIZE + (M * MSIZE)
    MODEL_RANGES.append(range(model_start, model_end))
    # then I store the information about the coil a-bead indices (0-INDEXED!!)
    for c1i in Coil1_A_beads:       # select out each coil from Coil1 in Coil1_A_beads
        temp_abead_posits = []
        for c1i_num in c1i:         # loop through every a-bead position in coil1i selected above
            temp_abead_posits.append(c1i_num + model_additive)          # add the model addifier to the number and put
                    # all those numbers in to the temporary list to add to the final container of every coil
        MODELS_A_BEADS.append(temp_abead_posits)

if args.verbose:
    print(f"All model coils a beads: {MODELS_A_BEADS}")
    print(f"Model ranges {MODEL_RANGES}")
    print()
# now we have all the index information we need to start doing the analysis


"""
Step 2 - do the trajectory analysis now! I will do the analysis BY FRAME
"""
# These are the TIME SERIES/Frame data for oligomers and self-vs-other interactions. Each entry will correspond to the
# total of that category, in order of frames
FREE_COILS = []         #
DIMERS = []; TRIMERS = []; TETRAMERS = []; HIGHER = []
SELF = []; OTHER = []

"""
Here is a general overview of the analysis algorithm
set master variables
Iterate through every frame{
    set the frame variables = 0
    iterate through every coil in the system = coil1{
        set the individual coil variables=0
        iterate through every *other* coil in the system = coili{
            calculate the neighbors between coil1 and coili
            determine self vs other interactions
            count interacting partners
        }
        determine what the oligomer is
        add the coil variables to the frame variables
    }
    add the frame variables to the master variables
    *continue through analysis
}
"""
# Begin the frame loop, analyze trajectory frame by frame
for frame in range(0, traj.n_frames, args.f):
    free_coil_counter = 0
    frame_dimer = 0; frame_trimer = 0; frame_tetramer = 0; frame_higherorder = 0
    frame_self = 0; frame_other = 0

    # begin the primary loop through every coil and set variables
    for ci in range(len(MODELS_A_BEADS)):
        COIL_MODEL = 0
        COIL = MODELS_A_BEADS[ci]
        oligo_counter = 0
        self_interact = 0
        other_interact = 0
        # determine which model I'm working with
        for ri in range(len(MODEL_RANGES)):
            if COIL[0] in MODEL_RANGES[ri]:
                COIL_MODEL = ri
            else:
                pass

        # now loop through every other coil and see if COIL is neighbors with it
        for ci2 in range(len(MODELS_A_BEADS)):
            if ci2 == ci:   # don't analyze current coil, duh
                continue
            else:
                COILi = MODELS_A_BEADS[ci2]
                # now find the neighbors!
                neighbors = md.compute_neighbors(traj[frame], A_CUTOFF, COIL, haystack_indices=COILi)[0]

                if len(neighbors) < np.floor(len(COIL) * INTERACTION_CUTOFF):
                    # no neighbors OR not sufficient coil interactions between COIL and COILi, so no oligomer
                    continue
                else:
                    # we have a bona fide oligomer, add it to the list
                    oligo_counter += 1

                    for ri in range(len(MODEL_RANGES)):
                        if neighbors[0] in MODEL_RANGES[ri]:    # only need to look at first entry
                            if ri == COIL_MODEL:
                                self_interact += 1
                            else:
                                other_interact += 1

        # now we're done checking if COIL has neighbors with any other coil in the system
        if oligo_counter == 0:      # no interactions, free-floating coil
            free_coil_counter += 1
        elif oligo_counter == 1:
            frame_dimer += 1
            frame_self += self_interact /2
            frame_other += other_interact /2
        elif oligo_counter == 2:
            frame_trimer += 1
            frame_self += self_interact / 2
            frame_other += other_interact / 2
        elif oligo_counter == 3:
            frame_tetramer += 1
            frame_self += self_interact / 2
            frame_other += other_interact / 2
        else:
            frame_higherorder += 1

    # now we're done with every single coil analysis for the given frame. Add the running totals to the lists
    FREE_COILS.append(free_coil_counter)
    DIMERS.append(frame_dimer/2)
    TRIMERS.append(frame_trimer/3)
    TETRAMERS.append(frame_tetramer/4)
    # these four arrays above contain the number of each oligomers throughout time.
    HIGHER.append(frame_higherorder)
    SELF.append(frame_self)
    OTHER.append(frame_other)
    
""" End the timing work. This is because the majority of the time will happen above with the 
actual calculation of the haystack """
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60) / 60 
# this reports the elapsed time in terms of hours
print("-----------------------------")
print(f"Elapsed time for the analysis fo haystack searching: {elapsed_time:.4f} hours")
print("-----------------------------")


# Handle the plotting of important information
if args.distro:
    free_average = np.average(FREE_COILS); free_std = np.std(FREE_COILS)
    norm_free = (free_average * 1.) / float(args.N)
    norm_free_std = (free_std * 1.) / float(args.N)
    dimer_average = np.average(DIMERS); dimer_std = np.std(DIMERS)
    norm_dimer = (dimer_average * 2.) / float(args.N)
    norm_dimer_std = (dimer_std * 2.) / float(args.N)
    trimer_average = np.average(TRIMERS); trimer_std = np.std(TRIMERS)
    norm_trimer = (trimer_average * 3.) / float(args.N)
    norm_trimer_std = (trimer_std * 3.) / float(args.N)
    tet_average = np.average(TETRAMERS); tet_std = np.std(TETRAMERS)
    norm_tet = (tet_average * 4.) / float(args.N)
    norm_tet_std = (tet_std * 4.) / float(args.N)

    fig, ax = plt.subplots()
    ax.errorbar([1, 2, 3, 4], [norm_free, norm_dimer, norm_trimer, norm_tet],
             yerr=[norm_free_std, norm_dimer_std, norm_trimer_std, norm_tet_std], color="black", linestyle="",
             markersize=10, marker=".", ecolor="black", elinewidth=1, capsize=5)
    plt.xticks(np.arange(0, 6))
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    plt.xlabel("N-mer")
    plt.ylabel("Counts")
    plt.grid(color="black", linestyle=":", alpha=0.5)
    plt.savefig(f"{args.name}_oligomeric_analysis_distribution.png", dpi=600)


if args.timeseries:
    arFrames = np.arange(0+args.start, (traj.n_frames*traj.timestep)+args.start, (args.f*traj.timestep))
    fig, ax = plt.subplots()
    ax.plot(arFrames, np.array(FREE_COILS), color="black", label = "Free coils",
              linewidth=1)
    ax.plot(arFrames, np.array(DIMERS), color="blue", label = "Dimers",
             linewidth=1)
    ax.plot(arFrames, np.array(TRIMERS), color="red", label = "Trimers",
             linewidth=1)
    ax.plot(arFrames, np.array(TETRAMERS), color="tab:orange", label = "Tetramers",
              linewidth=1)
    ax.plot(arFrames, np.array(HIGHER), color="tab:green", label = "Higher order oligos*",
              linewidth=1)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.ylabel("Counts")
    plt.xlabel("Simulation time (ps)")
    plt.savefig(f"{args.name}_oligomeric_analysis_timeseries.png", dpi=600)
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(arFrames, np.array(SELF), color="black", label = "Intra-model interactions",
              linewidth=1)
    ax.plot(arFrames, np.array(OTHER), color="blue", label = "Inter-model interactions",
             linewidth=1)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.ylabel("Counts")
    plt.xlabel("Simulation time (ps)")
    plt.savefig(f"{args.name}_selfvsother_analysis.png", dpi=600)
    plt.close()


# Now save out all the data so I can use it later:
outputdata = [arFrames, np.array(FREE_COILS), np.array(DIMERS), np.array(TRIMERS), np.array(TETRAMERS),
              np.array(HIGHER), np.array(SELF), np.array(OTHER)]
pd.DataFrame(np.array(outputdata)).to_csv(f"{args.name}_oligoAnalysis_outData.csv")

print(f"Number of frames in trajectory after slicing from the starting point: {remaining_frames}")
print(f"The total number of analyzed frames, based on the '-f' argument: {analyzed_frames}")

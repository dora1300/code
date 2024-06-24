"""
@Title:             ete_slab_analysis.py
@Name:              Mando A Ramirez
@Date:              2024 06 19

@Description:       This script calculates the ETE distances for every single coil and linker 
in a SLAB simulation and reports out the average +/- StDev ETE distance for coils and linkers 
for every analyzed time frame.


"""

"""
Import section
"""
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd

plt.rcParams['font.size'] = 14


""" 
Set up the arg parser 
"""
parser = argparse.ArgumentParser(description="Coil multimerization tool - this determines the types of multimers and"
                                             " higher order assemblies present in your coil simulations!")
parser.add_argument("-df", help="Definition file: this contains the information on how to build the model, with"
                                " alternating coil and linker segments.")
parser.add_argument("-nmodels", help="The total number of models in the simulation", required=True, type=int)
parser.add_argument("-l", help="No. of beads in an ENTIRE SINGLE PROTEIN. This is not the size of an individual coil. This is "
                               "the size of an entire individual protein. As of right now, ONLY proteins that have the"
                               " same no. of beads can be analyzed.",
                    required=True, type=int)
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-start", help="[ps] the start TIME [ps] to begin the analysis. Default = 0 [ps]", default=0,
                    type=int)
parser.add_argument("-stop", help="[ps] the stop TIME [ps] to end the analysis. Default is the entire simulation", type=int)
parser.add_argument("-f", help="the stride of analysis i.e. analyze every fth frame, default = 1", default=1, type=int)
parser.add_argument("-fstep", help="[ps] The time value of each frame, e.g. 10000 ps. No default given.",
                    required=True, type=int)
parser.add_argument("-verbose", help="pass this flag to output extra information to stdout.", action="store_true",
                    default=False)
parser.add_argument('-output', help="common name for output files. Do NOT include any type of extension. "
                    "Example of file name is (args.output)_plot_of_ete.png. Default = 'output'",
                    type=str, default='output')

args = parser.parse_args()

trajectory= args.t
topology = args.p


"""
Useful global constants
"""
# CUTOFF = args.cutoff                 # in nanometers!
NMODELS = args.nmodels
PSIZE = args.l
# MASS = args.mass


"""
Load the trajectory!
"""
# if args.single:
#     exit(0)
#     # traj = md.load(trajectory, top=topology)
#     # args.f = 1
#     # args.timeseries = False     # can't do timeseries for just 1 frame!
#     #     # this also adds a user sanity check to make sure that the code doesn't 
        
#     #     # try to do something that doesn't make sense.
#     # remaining_frames = 1
#     # analyzed_frames = 1
# else:
traj_load = md.load(trajectory, top=topology)
frame_start = 0
while frame_start <= traj_load.n_frames:
    if int(frame_start * traj_load.timestep) >= args.start:
        break
    else:
        frame_start += 1
frame_end = 0
if args.stop is None:
    frame_end = traj_load.n_frames
else:
    while frame_end <= traj_load.n_frames:
        if int(frame_end * traj_load.timestep) >= args.stop:
            break
        else:
            frame_end += 1
if frame_end < frame_start:
    raise ValueError("The specified end point is before the start point! Exiting because this doesn't make sense")
    exit(1)
traj = traj_load[frame_start:frame_end:]
remaining_frames = traj_load.n_frames - frame_start
analyzed_frames = remaining_frames / args.f



""" 
Print a warning statement
"""
print("**********")
print("PLEASE BE ADVISED: This code automatically converts the information provided in the -df file into 0-INDEXED "
      "positions! Your input is 1-indexed, but this code will automatically change it to 0-indexed.")
print("**********\n")




"""
Step 1 - read the file about the model description to find where all the coils and linkers are.
Start and stop positions are recorded!
code taken from other things like multimerization_analysis.py for file reading
"""
# Read the file and save only the information about the coils. Save the start/end indices for the coils and the start
# of the A-bead positions.
# POSITIONS ARE 0-INDEXED
Coil_positions_df = []
Linker_positions_df = []
coil_counter = 0                # this is how many coils are present in a model
linker_counter = 0              # this is how many linkers are present in a model
with open(args.df, "r") as file:
    for line in file:
        line_split = line.rstrip("\n").split(",")
        if str(line_split[0]) == "coil":
            # information in line_split = [coil_start, coil_end, a_bead_start]
            # information originally 1-indexed
            # this information needs to be 0-indexed, so I'm converting it!!
            coil_counter += 1
            Coil_positions_df.append([int(line_split[1])-1, int(line_split[2])-1])
        else:
            linker_counter += 1
            Linker_positions_df.append([int(line_split[1])-1, int(line_split[2])-1])

# Part of the sanity check ==> check that the size of each protein provided in the args
# matches what is in the protein model file.
prot_size = Coil_positions_df[-1][1]
try:
    (prot_size + 1) == PSIZE
except:
    print("The protein size (no. beads) does not match the size in the file. Exiting now!")
    print(f"Provided protein size in args: {PSIZE}")
    print(f"Protein size from model file: {prot_size + 1}")
    exit(1)

if args.verbose:
    print(f"Coil model file: {args.df}")
    print(f"Number of coils parsed: {coil_counter}")
    print(f"Parsed coil positions for a single coil (0-indexed): {Coil_positions_df}")
    print(f"Number of linkers parsed: {linker_counter}")
    print(f"Parsed linkers positions for a single coil (0-indexed): {Linker_positions_df}")


# Now enumerate all the a-bead indices for each model present in the simulation. The variable MODELS_A_BEADS will have
# the 0-BASED INDICES (!!!!) for every coil model in the system. THIS IS FOR ALL COILS
MODELS_COIL_POSITS = []
MODELS_LINK_POSITS = []
MODEL_RANGES = []
for M in range(0, NMODELS):
    # This information stores the 0-INDEXED (!!!) ranges for each of the models in the system!
    model_additive = M * PSIZE
    model_start = 0 + (M * PSIZE)
    model_end = PSIZE + (M * PSIZE)
    MODEL_RANGES.append(range(model_start, model_end))
    # then I store the information about the coil a-bead indices (0-INDEXED!!)
    for cli in Coil_positions_df:           # select out each coil from the coil_df list
        temp_coil_posits = []
        for cli_num in cli:                 # loop through every a-bead position in coil1i selected above
            temp_coil_posits.append(cli_num + model_additive)          # add the model addifier to the number and put
                    # all those numbers in to the temporary list to add to the final container of every coil
        MODELS_COIL_POSITS.append(temp_coil_posits)

    for lni in Linker_positions_df:           # select out each linker from the linker_df list
        temp_link_posits = []
        for lni_num in lni:                 # loop through every a-bead position in coil1i selected above
            temp_link_posits.append(lni_num + model_additive)          # add the model addifier to the number and put
                    # all those numbers in to the temporary list to add to the final container of every coil
        MODELS_LINK_POSITS.append(temp_link_posits)

if args.verbose:
    print(f"All model coils positions (0-indexed): {MODELS_COIL_POSITS}")
    print(f"All model linkers positions (0-indexed): {MODELS_LINK_POSITS}")
    print(f"Model ranges (0-indexed): {MODEL_RANGES}")
    print()
# now we have all the index information we need to start doing the analysis



"""
Step 3 -- now it's time to do the actual distance analysis
"""
# compute the distances of the provided atom pairs for coils and linkers
# This will create a 2d array that takes the form [[dist apair1, dist apair2, dist apair 3], ...] where every entry
# corresponds to a frame, and then within each frame entry, the distances are listed in order of pairs provided to it
dist_coils = md.compute_distances(traj, np.array(MODELS_COIL_POSITS))
dist_links = md.compute_distances(traj, np.array(MODELS_LINK_POSITS))

time_array = np.arange(frame_start, frame_end) * args.fstep



"""
Step 4 -- now make a graph of average distances over time so I can see what's going on at the
time level!
"""
coils_avg_by_time = np.average(dist_coils, axis=1)
coils_std_by_time = np.std(dist_coils, axis=1)
links_avg_by_time = np.average(dist_links, axis=1)
links_std_by_time = np.std(dist_links, axis=1)

fig, ax = plt.subplots()
plt.plot(time_array, coils_avg_by_time, marker="", linewidth=1.0, color='black',
         label="Coils")
plt.fill_between(time_array, coils_avg_by_time+coils_std_by_time,
                 coils_avg_by_time-coils_std_by_time,
                 color='black', alpha=0.20)
plt.plot(time_array, links_avg_by_time, marker="", linewidth=1.0, color='goldenrod',
         label="Linkers")
plt.fill_between(time_array, links_avg_by_time+links_std_by_time,
                 links_avg_by_time-links_std_by_time,
                 color='goldenrod', alpha=0.20)
plt.ylim(0, 6)
plt.ylabel("Distance (nm)")
plt.xlabel("Time (ps)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{args.output}_ETE_vs_time_plot.png", dpi=600)



"""
Step 5 -- now I can save out my data before doing anything different!
"""
np.savetxt(f"{args.output}_ETE_by_frame_COILS.csv", 
           np.array([time_array, coils_avg_by_time, coils_std_by_time]),
           delimiter=",")
np.savetxt(f"{args.output}_ETE_by_frame_LINKERS.csv", 
           np.array([time_array, links_avg_by_time, links_std_by_time]),
           delimiter=",")
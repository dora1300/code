"""
@Title:         coil_ncf_analysis.py
@Author:        Mando Ramirez
@Date:          20220221

@Description:   This code is a generalized NCF analysis script that I can use for coil segments. It is currently being
written for simulations of only coils as of the initial date.

@Description of NCF method:     The definition of NCF is HARDCODED and is the average distance of ONLY 1-4 interactions

@ToDos:

@Updates:
20220310 - fixed an error (lines 90, 91) where the indices are scaled according to the model number. I had incorrectly
put a fixed helix length in there when it should have been a variable of lhelix. That is now fixed.
"""

import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import numpy as np


""" Set up Argument Parser and parse the args! """
parser = argparse.ArgumentParser(description="Calculates Native Contact Fraction for a given model. Designed specifically"
                                             " for coil models. NCF in this code is defined ONLY by the 1-4 interaction distances.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_rama")
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-tol", help="Set the multiplicative tolerance value for the NCF calculation.", required=True,
                    type=float)
parser.add_argument("-f", help="Stride - analyze every f-th frame. Default is 10.", default=10, type=int)
parser.add_argument("-v", help="Print a statement about the definition of NCF. Exits after printing.", action="store_true")
parser.add_argument("-full", help="Switch to save the full, frame by frame output for each model. Produces"
                                  " a lot of data and large files. Provide flag to set to True.", action="store_true")

args = parser.parse_args()

if args.v:
    output_text = """Coil NCF analysis script\n
    This script defines average 1-4 interactions as the native structure. Average values
    came from the optimized coiled-coil dimer from CC Builder. Only 1-4 interactions are being calculated for the NCF.\n
    Averaged reference distance is hard coded in the script.\n
    This requires that all coil copies are the same length but can handle as many models as you want."""
    print(output_text)
    exit()
else:
    pass

""" Hard coded values for native contact reference as well as NCF position generation! """
# Set the constant values of reference distances. INT stands for interaction
# THESE VALUES ARE IN NANOMETERS
INT14 = 0.5084026047434603

# Set the tolerance value TOL
# this acts as the distance tolerance. Essentially, if the simulated distance is larger than reference
# distance by (ref*TOL) then it is not a native contact
TOL = args.tol

# Parse the arguments and put them into variables
trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially
output = args.o + f"_TOL{TOL}.csv"


# Define the array of indices for the 1-4 interactions
i4_indices = []
for i in range(lhelix-3):
    i4_indices.append([i, i+3])
i4_np = np.array(i4_indices)




""" Perform the NCF analysis for all models in the simulation. """
ncfs_by_models = []     # this is a list of lists, every entry is a model's NCFs for each frame analyzed

for n in range(nmodels):
    # make a switch to edit the indices so that they are for the right model
    if n == 0:
        pass
    else:
        for pair_i in range(len(i4_indices)):
            i4_indices[pair_i][0] = i4_indices[pair_i][0] + (n * lhelix)
            i4_indices[pair_i][1] = i4_indices[pair_i][1] + (n * lhelix)
        i4_np = np.array(i4_indices)

    model_ncf = []      # this holds the NCFs in each frame for the given model n
    traj = md.load(trajectory, top=topology)

    for frame in range(0, traj.n_frames, args.f):
        # this is to change the stride of analysis and make it a little faster
        # This output is in NANOMETERS!
        dists_14_traj = md.compute_distances(traj[frame], i4_np)
        native_contact_count = 0
        ncf_frame = 0

        for d in range(len(dists_14_traj[0])):
            # this selects all the distances calculated in the given frame
            # with the number of entries defined by the length of array i4_np
            if (dists_14_traj[0][d]) <= (INT14*TOL):
                native_contact_count += 1
        ncf_frame = native_contact_count / len(i4_indices)
        model_ncf.append(ncf_frame)

    ncfs_by_models.append(model_ncf)


""" Now do the final averaging and output the data to some files """
model_avgs = []; model_stds = []

for model in ncfs_by_models:
    model_avgs.append(np.average(np.array(model)))
    model_stds.append(np.std(np.array(model)))
print(model_avgs)
# Write out the short output with averages and std deviations
with open(output, 'w') as short:
    short.write(f"NCF analysis for: {sim_code}\n")
    short.write(f"Model,Avg NCF,StDev NCF\n")

    for i, val in enumerate(model_avgs):
        short.write(f"{i+1},{val},{model_stds[i]}\n")


# Make a plot comparing the NCFs of each model to each other
xtick_labels = []
for entry in range(nmodels):
    xtick_labels.append(f"Model {entry+1}")
fig, ax = plt.subplots()
plt.errorbar(np.arange(nmodels), np.array(model_avgs), yerr=np.array(model_stds), ms=3, linestyle="", mec="black",
             mfc="black", marker=".", elinewidth=1, capsize=1.5, ecolor="red")
ax.set_xticks(np.arange(nmodels))
ax.set_xticklabels(xtick_labels)
ax.set_ylim(0, 1.1)
ax.set_xlabel("Models")
ax.set_ylabel("Native Contact Fraction")
plt.title(f"NCF analysis for {sim_code}\n Average NCF per model")
plt.savefig(f"{sim_code}_ncf_plot_TOL_{TOL}.png", dpi=300)


# Save the FULL output if desired.
if args.full:
    with open(f"{args.o}_full_ncf_data_TOL{TOL}.csv", 'w') as long:
        long.write(f"Frame")
        for num in range(nmodels):
            long.write(f",Model {num+1}")
        long.write("\n")

        for pseudoframe in range(len(ncfs_by_models[0])):
            long.write(f"{pseudoframe * args.f}")
            for mod in range(nmodels):
                long.write(f",{ncfs_by_models[mod][pseudoframe]}")
            long.write("\n")


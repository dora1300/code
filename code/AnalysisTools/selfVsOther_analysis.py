"""
@Title:         selfVsOther_analysis.py
@Date:          2023 05 11
@Author:        Dominique A Ramirez

@Description:   This script uses the "contact_list" output from multimerization_analysis.py,
which is a list of multimer contacts in every analyzed frame, and calculates
things about intra- and inter-chain interactions.

This needs a contact map input file to be useful, which comes from the multimerization_analysis.py
script.

THIS ONLY WORKS FOR HOMOPOLYMERS WITH THE SAME LENGTH COILS.

"""
import numpy as np
import pandas as pd
import argparse
import time
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="This script analyzes coil multimer contacts to produce results about:"
                                             " (1) The fraction of intra-chain (self) vs"
                                             "   inter-chain (other) interactions, on a per"
                                             "   coil basis.")
parser.add_argument("-data", help="Data file from multimerization analysis -- list of coil multimer"
                                  " contacts, each line marks a new frame (.csv)."
                                  " Please provide path if not in working directory.", required=True)
parser.add_argument("-ncoils", help="The total number of coils in the simulation", required=True,
                    type=int)
parser.add_argument("-nprots", help="The total number of individual proteins in the simulation", required=True,
                    type=int)
parser.add_argument("-cpp", help="Not C++ lol! The number of coils *per protein*. This analysis"
                                 " assumes you're working with homopolymers only. Support for heteropolymers is not yet"
                                 " supported", type=int, required=True)
parser.add_argument("-output", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)

args = parser.parse_args()

# Important variables for this analysis. Remember simulations are homopolymers only right now
RESULTS = args.data
OUTPUT = args.output
NCOILS = args.ncoils  # total number of coils in the simulation
NPROTS = args.nprots  # total number of protein molecules in the simulation
COILpp = args.cpp  # the number of COILs Per Protein

# "Allocate" the arrays I will use to store my analyses results
coil_int_array = np.empty((NCOILS,), dtype=object)
# this will become an array of lists, where each entry coil_int_array[i] = [S, O] for
# coil_i, where S = the total number of self-interactions (meaning same protein)
# and O = total number of other-interactions (different protein)
protein_ranges_arr = np.empty((NPROTS,), dtype=object)
# this will essentially hold the identities of each of the different proteins in the sim
# and by identities, I mean which coil indicies belong to which protein
# The identities will be ranges OF COIL INDICES


for i, obj in enumerate(coil_int_array):
    coil_int_array[i] = [0, 0]  # the total counts for self- and other- start at 0, duh

coil_start = 0
coil_end = 0 + COILpp
for i, obj in enumerate(protein_ranges_arr):
    protein_ranges_arr[i] = range(coil_start, coil_end)
    coil_start += COILpp;
    coil_end += COILpp

# Here is a sanity check to make sure I haven't f---ed something up
if coil_start != NCOILS:
    # note! coil_start is TECHNICALLY 0-indexed, whereas NCOILS is TECHNICALLY
    # 1-indexed! That's okay though, because this check will work regardless
    print("The combination of number of coils per protein, number of proteins in the "
          "simulation, and total number of coils does not match. Check out your inputs "
          "and try again.")
    print(f"N-coils = {NCOILS}; Coils per protein * N-proteins = {COILpp * NPROTS}")
    exit(1)

"""
Begin the self-/intra- vs other-/inter-chain analysis!
"""
# It's always good to have a timer going so I see how long things take:
time_start = time.perf_counter()

with open(RESULTS, 'r') as FILE:
    for line in FILE:
        list_of_contacts = line.rstrip("  \n").split("  ")
        # contains the list of contacts in the given line i.e. frame

        for contacts in list_of_contacts:
            # the mapping code came from some online help:
            # <https://stackoverflow.com/questions/6429638/
            # how-to-split-a-string-of-space-separated-numbers-into-integers>
            try:
                indices = list(map(int, contacts.split(" ")))
            except:
                # If there's no interaction in the given frame, then just move onto
                # the next frame!
                continue

            for I, coilI in enumerate(indices):
                # check that my method for matching coil with protein works
                if coilI in protein_ranges_arr[int(coilI / COILpp)]:
                    # my method works
                    pass
                else:
                    # my method doesn't work, I can't find the right protein!
                    # this is a problem obviously because then I can't determine self from other
                    print("There is a problem with matching the given coil index to"
                          " the right protein. Exiting now.")
                    print(f"Coil index in question: {coilI}")
                    print(f"Protein range that was selected: {protein_ranges_arr[int(coilI / COILpp)]}")
                    exit(1)

                for J, coilJ in enumerate(indices):
                    if J == I:
                        # don't count interactions between the same coil, duh!
                        continue
                    else:
                        if coilJ in protein_ranges_arr[int(coilI / COILpp)]:
                            # this is a SELF interaction! CoilJ is in the same protein as
                            # coilI
                            coil_int_array[coilI][0] += 1
                        else:
                            # otherwise, it's an OTHER interaction
                            coil_int_array[coilI][1] += 1

# Stop the timing analysis and print out the data!
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60)
print("-----------------------------")
print(f"Elapsed time for doing analysis: {elapsed_time:.4f} minutes")
print("-----------------------------")

"""
Now we can do the fractional analysis and save out data
"""
coil_int_frac_array = np.empty((NCOILS,), dtype=object)
# this will hold the normalized data. For each coil in coil_int_array I am going to
# divide the self and the other counts by the total interaction counts to get the percentage
# of self- vs other- interactions that each coil makes, in each protein.
total_interactions = 0
for i, coil in enumerate(coil_int_array):
    total_interactions = coil[0] + coil[1]
    if total_interactions == 0:
        coil_int_frac_array[i] = [0, 0]
    else:
        coil_int_frac_array[i] = [coil[0] / total_interactions,
                                  coil[1] / total_interactions]

"""
# Now time for the averaged analysis for each coil ID, across all proteins
"""
self_avgByCoil = np.empty((COILpp,))
self_stdByCoil = np.empty((COILpp,))
other_avgByCoil = np.empty((COILpp,))
other_stdByCoil = np.empty((COILpp,))

reshaped = coil_int_frac_array.reshape(NPROTS, COILpp)
self_average = []
other_average = []
for i in range(COILpp):
    sameCoils_diffProts = reshaped[:, i]
    self_average.clear();
    other_average.clear()
    for k in range(len(sameCoils_diffProts)):
        self_average.append(sameCoils_diffProts[k][0])
        other_average.append(sameCoils_diffProts[k][1])
    self_avgByCoil[i] = np.average(self_average)
    self_stdByCoil[i] = np.std(self_average)
    other_avgByCoil[i] = np.average(other_average)
    other_stdByCoil[i] = np.std(other_average)

error_config = {'elinewidth': 2,
                'capthick': 1.5,
                'capsize': 5}
WIDTH = 0.6
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))
AX1 = ax1.bar(np.arange(1, COILpp + 1, 1), self_avgByCoil, width=WIDTH,
              yerr=self_stdByCoil, color="goldenrod", alpha=0.75, edgecolor="black",
              ecolor="black", error_kw=error_config, label="Intra-chain")
ax1.legend()
ax1.set_ylim(0, 1.5)
ax1.grid(color="black", linestyle=":", alpha=0.35)
AX2 = ax2.bar(np.arange(1, COILpp + 1, 1), other_avgByCoil, width=WIDTH,
              yerr=other_stdByCoil, color="grey", alpha=0.75, edgecolor="black",
              ecolor="black", error_kw=error_config, label="Inter-chain")
ax2.set_ylim(0, 1.5)
ax2.grid(color="black", linestyle=":", alpha=0.35)
ax2.legend()
plt.xticks(np.arange(1, COILpp + 1, 1))
plt.xlabel("Coil position in protein")
fig.supylabel("Fraction of interactions, averaged")
plt.tight_layout()
plt.savefig(f"{OUTPUT}_selfVsOther_averaged.png", dpi=600)

output_data = pd.DataFrame({
    "Coil ID in protein": np.arange(1, COILpp + 1, 1),
    "Avg self interaction": self_avgByCoil,
    "StDev self interaction": self_stdByCoil,
    "Avg other interaction": other_avgByCoil,
    "StDev other interaction": other_stdByCoil
})
output_data.to_csv(f"{OUTPUT}_selfVsOther_averaged_interactions.txt", index=False)

"""
@Title:         contact_analysis.py
@Date:          2023 05 01
@Author:        Dominique A Ramirez

@Description:   This script uses the output of the multimerization_analysis.py, which is a list of multimer
contacts in every analyzed frame, and performs some analysis of those data. In particular, I am analyzing:
    - Total number of interactions that an average coil makes throughout the simulation
    - Total number of unique interactions that a coil sees
    - Lifetime of coil interactions

Description of what I'm defining as an interaction:
    I have made a simplification in how I define an interaction, which makes the lifetime analysis much easier than
the alternative.

@Updates:
2023 05 02 - I was incorrectly counting total interactions and unique partners for coils. That is now fixed.

2023 08 28 - I was still incorrectly counting total interactions for individual coils. I *was* summing the values
in the interaction arrays for each coil, but what's in those arrays are the indices of every coil that a given 
coil_i interacts with. Summing that produces a dasterdly high and not correct number. This is now fixed by
calculating the len() of the interaction array for each coil.
    It is true that the number of interactions per coil can be larger than the number of frames analyzed, but it 
    shouldn't be > (multimerization_state * numer_of_frames)
"""

import numpy as np
import scipy.stats as st
import pandas as pd
import argparse
import time
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="This script analyzes coil multimer contacts to produce results about:"
                                             " (1) Number of unique interactions that a coil will see,"
                                             " (2) Total number of interactions an average coil will make, and"
                                             " (3) Lifetime (in ps) of coil multimer interactions.")
parser.add_argument("-data", help="Contacts list data file multimerization analysis -- list of coil multimer"
                                  " contacts, each line marks an ANALYZED frame (.csv)."
                                  " Please provide path if not in working directory.", required=True)
parser.add_argument("-ncoils", help="The total number of coils in the simulation", required=True,
                    type=int)
#parser.add_argument("-nprots", help="The total number of individual proteins in the simulation", required=True,
#                    type=int)
parser.add_argument("-ft", help="[ps] Frame time, i.e. the amount of time that each frame in the Contacts Data File "
                                "is worth in simulation. Integers only.",
                    required=True, type=int)
parser.add_argument("-output", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)
parser.add_argument("-lifetimes", help="Pass this flag to only calculate the lifetimes of interactions.",
                    action="store_true", default=False)

args = parser.parse_args()


"""
Warning message about the frame time -- have to make sure the user understands!
"""
print("""
CAUTION -- be careful that you have calculated the frame time (flag -ft) correctly, and make sure that
the units you've calculated match what the code expects. This analysis code has absolutely no way of knowing
how much simulation time corresponds to each analyzed frame in the _contacts.txt data file.
""")


# Assign the arguments to their rightful variables
NCOILS = args.ncoils
#NMODELS = args.nprots
FRAME_TIME = args.ft          # units of ps
RESULTS = args.data


total_partners_array = np.empty((NCOILS,), dtype=object)
    # this will become an array of lists. Each entry array[i] corresponds to a unique coil (0..n..NCOILS). The entry
    # is a list of ALL the partners that coil_n sees throughout the simulation.
    # THIS MEANS THE AMOUNT OF ENTRIES IN EACH ARRAY IS NOT THE SAME AS THE NUMBER OF FRAMES ANALYZED
    # for larger multimers, each coil will have MULTIPLE partners per frame. This data may or may not be useful
time_in_multimer_array = np.zeros((NCOILS,), dtype=object)
    # this is an array of length NCOILS, where every entry arr[i] = the amount of time that the given coil_i is
    # in a multimer. This means that the len(arr[i]) <= Number of frames or len(interactions_by_frame) NECESSARILY!!
unique_partners_array = np.empty((NCOILS,), dtype=object)
    # this will become an array of lists. Each entry array[i] corresponds to a unique coil (0..n..NCOILS).
    # Each entry is a list of UNIQUE partners that coil_n sees throughout the simulation. Thus, no repeats
interactions_by_frame = []
    # This is a simple counting container that holds the total interactions i.e. the number of COILS engaged
    # in a multimer
multimer_lifetimes = []
    # this is also a simple counting container. Each entry is the number of frames that any given interaction
    # exists for in simulation before the interaction changes.
    # PLEASE NOTE!!! --> I am defining interaction as a UNIQUE MULTIMER i.e. if the multimer identity changes, then
    # the INTERACITON CHANGES

for i, obj in enumerate(total_partners_array):
    total_partners_array[i] = []
    unique_partners_array[i] = []


# Good to have a timer for this
time_start = time.perf_counter()


# These lists are useful only for the lifetime analysis
existing_multimers = []
    # this only holds the multimers that exist in the PREVIOUS frame
list_of_multimers = []
    # Each entry is a LIST which represents a MULTIMER. Each entry in the LIST is the 0-Indexed coil number
    # e.g. list_of_multimers = [[0,18], [5,9]]. Multimers[i] is the interaction [0,18], which is a dimer between
    # coils 0 and 18
list_of_multimer_counters = []
    # This keeps track of the lifetime of each of the multimers in the above lists. Thus, THE INDICES OF THIS LIST
    # CORRESPOND TO THE INDICES OF THE ABOVE LIST. e.g. list_of_multimer_counters = [10, 2].
    # list_of_multimer_counters[0] = 10, which is the number of frames that list_of_multimers[0]=[0,18] has existed
multimers_accounted_for = []
no_longer_existing = []
lifetime = 0

"""
Begin the actual contact analysis!
"""
with open(RESULTS, 'r') as f:
    for line in f:
        list_of_contacts = line.rstrip("  \n").split("  ")
        # track the total number of interactions in the frame (easy peasy to track)
        interactions_by_frame.append(len(list_of_contacts))

        multimers_accounted_for.clear()
        no_longer_existing.clear()
        for contacts in list_of_contacts:
            # the mapping code came from some online help:
            # <https://stackoverflow.com/questions/6429638/
            # how-to-split-a-string-of-space-separated-numbers-into-integers>
            try:
                # This turns all the coil indices into integers, using the map() function
                # map() applies a function to all items in a list (which here is generated using
                # .split()) and then returns an iterator, which I put into a list using list()
                indices = list(map(int, contacts.split(" ")))
            except:
                # this will fail only if there are no interactions!
                # Thus, I need to reset all the multimer_interaction counters, and capture the lifetime data
                # and add it to the multimer_lifetimes list
                for j, ints in enumerate(existing_multimers):
                    if ints not in multimers_accounted_for:
                        dead_multimer_i = list_of_multimers.index(ints)
                            # this multimer is dead because it was present in the previous frame, but not
                            # the current frame
                        lifetime = list_of_multimer_counters[dead_multimer_i]
                        multimer_lifetimes.append(lifetime)
                        list_of_multimer_counters[dead_multimer_i] = 0
                existing_multimers[:] = [x for x in existing_multimers if (x in multimers_accounted_for)]
                continue


            # This code handles finding unique partners for interacting coils
            # as well as keeping track of all partners a coil sees
            # Keep in mind that the multimer analysis outputs data that isn't duplicated
            # but I HAVE to duplicate it here so that I correctly count all interactions
            # for ALL coils
            for coili in range(len(indices)):
                time_in_multimer_array[indices[coili]] += 1
                coilj = 0
                while coilj < len(indices):
                    if coili == coilj:
                        pass
                    else:
                        total_partners_array[indices[coili]].append(indices[coilj])
                        if indices[coilj] not in unique_partners_array[indices[coili]]:
                            unique_partners_array[indices[coili]].append(indices[coilj])
                    coilj += 1


            # This code handles the lifetime of multimer interaction calculations
            if indices not in existing_multimers:
                # the current interaction is not in memory of the last frame i.e.
                # brand new to this frame
                # I need to first figure out if this is a multimer I've seen before
                # or if it's totally brand new to the simulation
                existing_multimers.append(indices)

                try:
                    # have I seen this interaction before? if so, then I'll get an index
                    # and add 1 to the counter
                    i = list_of_multimers.index(indices)
                    list_of_multimer_counters[i] += 1
                except:
                    # if I haven't seen this interaction before, a ValueError will be raised
                    # So, I'll then add the interaction to the list
                    list_of_multimers.append(indices)
                    # IMPORTANT! I am also adding a number to the counter list WHICH
                    # HAS THE SAME INDEX AS THE INDICES IN THE ABOVE list
                    # this counter will serve as a counter for how long the multimer
                    # exists
                    list_of_multimer_counters.append(1)
                multimers_accounted_for.append(indices)
            else:
                # the current interaction IS in the existing multimers
                i = list_of_multimers.index(indices)
                list_of_multimer_counters[i] += 1
                multimers_accounted_for.append(indices)

        # Now, all the contacts for the current frame have been counted for in the lifetime
        # assessment. I need to find all the contacts in EXISTING_MULTIMERS that are not
        # accounted for i.e. no longer exist, capture their lifetime data and save it to the
        # global list, then reset the value of their counter
        for j, ints in enumerate(existing_multimers):
            if ints not in multimers_accounted_for:
                dead_multimer_i = list_of_multimers.index(ints)
                    # this multimer is dead because it was present in the previous frame, but not
                    # the current frame
                lifetime = list_of_multimer_counters[dead_multimer_i]
                multimer_lifetimes.append(lifetime)
                list_of_multimer_counters[dead_multimer_i] = 0
        existing_multimers[:] = [x for x in existing_multimers if (x in multimers_accounted_for)]

    # and then finally, I need to grab all the lifetimes for interactions that exist through to
    # the final frame
    for j, ints in enumerate(existing_multimers):
        dead_multimer_i = list_of_multimers.index(ints)
        lifetime = list_of_multimer_counters[dead_multimer_i]
        multimer_lifetimes.append(lifetime)

# Stop the timing analysis and print out the data!
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60)
print("-----------------------------")
print(f"Elapsed time for reading: {elapsed_time:.4f} minutes")
print("-----------------------------")


plt.rcParams['font.size'] = 14
"""
Do some analysis of the LIFETIME data, and save a plot. Also, save out the list of lifetimes in case
I want to do something else with it.
    Don't forget to convert the lifetimes [frames] to simulation time [ps]!!
"""
# first, handle the plot (which is in ps)
fig, ax = plt.subplots()
ax.hist((np.array(multimer_lifetimes)*FRAME_TIME), density=True, color="grey",
        bins=int(np.max(multimer_lifetimes*FRAME_TIME)/2))
plt.xlabel("Multimer lifetimes (ps)")
plt.ylabel("Counts")
plt.title("Lifetime of multimer interactions")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_multimer_lifetimes(ps).png", dpi=600)
plt.close()

# second, handle the plot (which is in ns) NOTICE THE DIFFERENT UNITS FOR CONVENIENCE!
fig, ax = plt.subplots()
ax.hist((np.array(multimer_lifetimes)*FRAME_TIME)/1000, density=True, color="grey",
        bins=int(np.max(multimer_lifetimes*FRAME_TIME)/2))
plt.xlabel("Multimer lifetimes (ns)")
plt.ylabel("Counts")
plt.title("Lifetime of multimer interactions")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_multimer_lifetimes(ns).png", dpi=600)
plt.close()

# save stats in ps only
avg_lifetime = np.average(np.array(multimer_lifetimes) * FRAME_TIME)
std_lifetime = np.std(np.array(multimer_lifetimes) * FRAME_TIME)
max_lifetime = np.max(np.array(multimer_lifetimes) * FRAME_TIME)

with open(f"{args.output}_multimerLifetimes_stats.csv", 'w') as f:
    f.write(f"Average lifetime (ps),{avg_lifetime}\n")
    f.write(f"StDev lifetime (ps),{std_lifetime}\n")
    f.write(f"Max lifetime (ps),{max_lifetime}\n")

np.savetxt(f"{args.output}_multimerLifetimes(ps).csv",
           (np.array(multimer_lifetimes) * FRAME_TIME), delimiter=",")

pd.DataFrame(np.array(multimer_lifetimes) * FRAME_TIME).to_csv(f"{args.output}_multimerLifetimes(ps).csv")


# Here I can break if I'm only calculating the lifetime data!
if args.lifetimes:
    exit(0)
else:
    pass

"""
Total interactions throughout the simulation for each coil.
This is plotted as percentage of simulation time that a coil is engaged in an interaction
Keep in mind that the "interactions_by_frame" data is identical to the multimer data that I get from the other script

Updated 2023 05 11 -- the code is fixed to use time_in_multimer_array data, which is correct for percentage
simulation time that a coil is in an interaction
"""
percentage_time_coil_interactions = []
for i, coil in enumerate(time_in_multimer_array):
    # coil is the number of frames coil_i is in a multimer in the simulation. Thus, there can be no more than
    # len(interactions_by_frame) partners for coil_i. len(interactions_by_frame) is just how many frames that were
    # analyzed
    percentage_time_coil_interactions.append(coil / len(interactions_by_frame) * 100)

avg_interactions = np.sum(interactions_by_frame) / len(interactions_by_frame)
print(f"Avg NUMBER of interactions by frame: {avg_interactions}")
avg_percentage_interactions = np.average(percentage_time_coil_interactions)
print(f"Avg percentage of time a coil is in a multimer (for all frames): {avg_percentage_interactions}")

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(np.arange(0, NCOILS), percentage_time_coil_interactions, color="grey", linewidth=1.25)
plt.hlines(y=avg_percentage_interactions, xmin=0, xmax=NCOILS+5, color="goldenrod",
           label=f"Average percentage \ninteractions: {avg_percentage_interactions:.4f}")
plt.xlim(0, NCOILS + 1)
plt.ylim(0, 105)
plt.xlabel("Coil number")
plt.ylabel("Percentage of simulation time")
plt.title("Time that each coil is in a multimer")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.legend()
plt.tight_layout()
plt.savefig(f"{args.output}_percentageTimeInMultimer_byCoil.png", dpi=600)
plt.close()

pd.DataFrame(percentage_time_coil_interactions).to_csv(f"{args.output}_percentageCoilInteractions_byCoil.csv",
                                                       header=["Percentage of time in multimer"])

# Just for shits and giggles I'll write out the total number of partners that a coil sees. Maybe it's useful, maybe not
total_summed_partners_array = []
for i, coil in enumerate(total_partners_array):
    total_summed_partners_array.append(len(coil))
pd.DataFrame(total_summed_partners_array).to_csv(f"{args.output}_totalNumberPartners_byCoil.csv",
                                                       header=["Total number of partners throughout simulation"])

"""
This calculates the unique partners that each coil sees throughout the simulation.
Also there's a plot of it too:
    1) distribution of unique coil partners
    2) plot of unique coil partners for each coil
"""
n_partners_list = []
for coil in unique_partners_array:
    n_partners_list.append(len(coil))

pd.DataFrame(n_partners_list).to_csv(f"{args.output}_uniquePartners_byCoil.csv",
                                     header=["No. unique partners"])

fig, ax = plt.subplots(figsize=(7, 5))
plt.plot(np.arange(0, NCOILS), n_partners_list, color="grey", linewidth=1.24)
plt.xlim(0, NCOILS + 1)
plt.ylim(0, max(n_partners_list)+1)
plt.xlabel("Coil number")
plt.ylabel("Count")
plt.title("Number of unique partners for each coil")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_uniquePartners_byCoil.png", dpi=600)
plt.close()


fig, ax = plt.subplots()
plt.hist(n_partners_list, density=True, bins='sqrt', color="grey")
plt.xlabel("No. unique partners")
plt.ylabel("Counts")
plt.title(f"distribution of unique partners \nin entire simulation")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_uniquePartners_distribution.png", dpi=600)


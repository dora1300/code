"""
@Title:         multimer_lifetimes_analysis.py
@Date:          2023 08 22
@Author:        Dominique A Ramirez

@Description:   This script uses the output of the multimerization_analysis.py, which is a list of multimer
contacts in every analyzed frame, and performs some analysis of those data.
    This code is a spin-off from contact_analysis.py and is specific for determining multimer lifetimes only.
This code is specifically meant to handle "flickering" of multimer interactions by incorporating a lag-time feature.

This allows me to make edits to multimer calculation code without messing up the code I've already used for previous
analyses.

"""

import numpy as np
import scipy.stats as st
import pandas as pd
import argparse
import time
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="This script analyzes coil multimer contacts to produce results about:"
                                             " (1) Lifetime (in ps and ns) of coil multimer interactions.")
parser.add_argument("-data", help="Contacts list data file multimerization analysis -- list of coil multimer"
                                  " contacts, each line marks an ANALYZED frame (.csv)."
                                  " Please provide path if not in working directory.", required=True)
parser.add_argument("-ncoils", help="The total number of coils in the simulation", required=True,
                    type=int)
parser.add_argument('-lag', help="The 'lag time' for lifetime analysis. 'Lag-time' in this case means the number of "
                    "interrupted frames that are permitted for a given multimer interaction to still be "
                    "counted as an interaction. Large values of lag will likely result in unphysical interpretations. "
                    "Default = 0, which means that interactions broken by just 1 frame are separate interactions. "
                    "Default behavior is identical to just calculating the multimer lifetimes without lag.",
                    default=0, type=int)
#parser.add_argument("-nprots", help="The total number of individual proteins in the simulation", required=True,
#                    type=int)
parser.add_argument("-ft", help="[ps] Frame time, i.e. the amount of time that each frame in the Contacts Data File "
                                "is worth in simulation. Integers only.",
                    required=True, type=int)
parser.add_argument("-output", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)


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
NCOILS      = args.ncoils
#NMODELS    = args.nprots
FRAME_TIME  = args.ft          # units of ps
RESULTS     = args.data
LAG         = args.lag

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


# These lists are useful only for the lifetime analysis
existing_multimers = []
    # this only holds the multimers that exist in the PREVIOUS frame
list_of_multimers = []
    # a list of lists
    # Each entry is a LIST which represents a MULTIMER. Each entry in the LIST is the 0-Indexed coil number
    # e.g. list_of_multimers = [[0,18], [5,9]]. Multimers[i] is the interaction [0,18], which is a dimer between
    # coils 0 and 18
list_of_multimer_counters = []
    # a list of lists !!!!
    # This keeps track of the lifetime of each of the multimers in the above lists. Thus, THE INDICES OF THIS LIST
    # CORRESPOND TO THE INDICES OF THE ABOVE LIST. e.g. list_of_multimer_counters = [10, 2].
    # NOTICE -- this is different than what is used for the contact_analysis.py script
    # Each entry is a list which basically contains two values. 1 = multimer exists in a given frame. 
    # NaN = multimer does not exist in the given frame.
    # Each entry, whose location list_of_multimer_counters[i] is identical to the index used in list_of_multimers[i],
    # will have a list that corresponds to frames BUT! len(list_of_multimer_counters[i]) <= N_Frames
    # Because len(list_of_multimer_counters[i]) may be less than N_Frames if list_of_multimers[i] doesn't exist until partway
    # through the simulation, then you **CANNOT** use the indices list_of_multimers[i][0..j] as a measure of the simulation
    # time.
multimers_accounted_for = []
no_longer_existing = []
lifetime = 0


# Good to have a timer for this
time_start = time.perf_counter()
"""
Start by parsing the contacts list and turning the information within that file into
something that is usable for calculating lifetimes.
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
                # Thus, I need to add NaNs to all the multimers in my list of counter-lists because no 
                # interactions exist!
                for j, multi in enumerate(list_of_multimers):
                    if multi not in multimers_accounted_for:
                        absent_multimer_i = list_of_multimers.index(multi)
                        list_of_multimer_counters[absent_multimer_i].append(float('nan'))
                        continue


            # This code handles the lifetime of multimer interaction calculations
            # Have I already seen the current multimer? ('indices')
            # Yes! Okay, find it's counter list and append a 1 to it
            try:
                i = list_of_multimers.index(indices)
                list_of_multimer_counters[i].append(1)
                # have to keep track of what I've seen and what I haven't
                multimers_accounted_for.append(indices)
            except:
                # No, I haven't seen this multimer before! 
                # this will raise a ValueError which I will handle
                # So, I'll then add the interaction to the list
                list_of_multimers.append(indices)
                # IMPORTANT! I am also adding a *brand new list* to the counter list WHICH
                # HAS THE SAME INDEX AS THE INDICES IN THE ABOVE list
                # this list will keep track of when the multimer exists and when it doesn't
                list_of_multimer_counters.append([1])
                # have to keep track of what I've seen and what I haven't
                multimers_accounted_for.append(indices)

        # Okay, now I've gone through every contact (multimer) in the current frame ('line') and now I need to add NaNs 
        # to all the multimers that I DIDN't see
        for j, multi in enumerate(list_of_multimers):
            if multi not in multimers_accounted_for:
                # this means I didn't see the multimer 'multi' this frame!
                absent_multimer_i = list_of_multimers.index(multi)
                list_of_multimer_counters[absent_multimer_i].append(float('nan'))


# Stop the timing analysis for parsing the data
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60)
print("-----------------------------")
print(f"Elapsed time for parsing the contact list: {elapsed_time:.4f} minutes")
print("-----------------------------")


"""
Now, it is time for actually calculating the lifetimes. This will incorporate the lag-time
feature.
"""
# Good to have another timer for this
time_start = time.perf_counter()
# this is the chief list to hold the lifetimes
multimer_lifetimes = []


for contact in list_of_multimer_counters:
    lifetime_counter = 0
    nancounter = 0
    j = 0
    while j <= len(contact):
        # case 0, handle the end of the list and the break
        if j == len(contact):
            if contact[j-1] == 1:
                multimer_lifetimes.append(lifetime_counter)
                break
            else:
                if lifetime_counter != 0:
                    multimer_lifetimes.append(lifetime_counter)
                    break
                else:
                    break

        # case 1, frame = 1 and no nans have been seen, increase the lifetime
        if contact[j] == 1 and nancounter == 0:
            lifetime_counter += 1
            nancounter = 0
            j += 1
            continue

        # case 2, frame = 1, but we've seen NaNs. However! The number of NaNs is less than the lag
        # so we are counting all of the NaNs as an interaction, so add the NaNs to the counter
        if contact[j] == 1 and nancounter > 0 and nancounter <= LAG:
            lifetime_counter += (1 + nancounter)
            nancounter = 0
            j += 1
            continue

        # case 3, we've encountered a NaN in the MIDDLE of an interaction potentially.
        # Increase the nancounter and make sure its less than the LAG
        if contact[j] != 1 and nancounter < LAG and lifetime_counter > 0:
            nancounter += 1
            j += 1
            continue

        # case 4, We've encountered a NaN while monitoring an interaction with counter > 0
        # BUT, we're larger than the lag. Add the counter to the lifetimes as the interaction is now
        # dead and reset values
        # This is the only way to add a lifetime to the list (except the end) because the only way a lifetime
        # ends is when we encounter a long enough nan stretch.
        if contact[j] != 1 and nancounter >= LAG and lifetime_counter != 0:
            multimer_lifetimes.append(lifetime_counter)
            lifetime_counter = 0
            nancounter = 0
            j += 1
            continue
            
        # cases 5, We've encountered a NaN while not monitoring an interaction, so just move on.
        if contact[j] != 1 and nancounter <= LAG and lifetime_counter == 0:
            j += 1
            continue
        if contact[j] != 1 and nancounter >= LAG and lifetime_counter == 0:
            j += 1
            continue


# Stop the timing analysis for calculating lifetimes
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60)
print("-----------------------------")
print(f"Elapsed time for calculating lifetimes: {elapsed_time:.4f} minutes")
print("-----------------------------")


multimer_lifetimes_for_plotting = []

for i in multimer_lifetimes:
    if i == 1:
        continue
    else:
        multimer_lifetimes_for_plotting.append(i)

plt.rcParams['font.size'] = 14
"""
Do some analysis of the LIFETIME data, and save a plot. Also, save out the list of lifetimes in case
I want to do something else with it.
    Don't forget to convert the lifetimes [frames] to simulation time [ps]!!
"""
# first, handle the plot (which is in ps)
fig, ax = plt.subplots()
ax.hist((np.array(multimer_lifetimes_for_plotting)*FRAME_TIME), density=True, color="grey",
        bins=10000)
plt.xlabel("Multimer lifetimes (ps)")
plt.ylabel("Counts")
plt.xlim(-100, 100000)
plt.title("Lifetime of multimer interactions")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_multimer_lifetimes(ps).png", dpi=600)
plt.close()

# second, handle the plot (which is in ns) NOTICE THE DIFFERENT UNITS FOR CONVENIENCE!
fig, ax = plt.subplots()
ax.hist((np.array(multimer_lifetimes_for_plotting)*FRAME_TIME)/1000, density=True, color="grey",
        bins=10000)
plt.xlabel("Multimer lifetimes (ns)")
plt.ylabel("Counts")
#plt.xlim(-10, 100)
plt.title("Lifetime of multimer interactions")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.tight_layout()
plt.savefig(f"{args.output}_multimer_lifetimes(ns).png", dpi=600)
plt.close()

np.savetxt(f"{args.output}_multimerLifetimes(ps).csv",
           (np.array(multimer_lifetimes_for_plotting) * FRAME_TIME), delimiter=",")

pd.DataFrame((np.array(multimer_lifetimes_for_plotting) * FRAME_TIME)/1000).to_csv(f"{args.output}_multimerLifetimes(ns).csv")

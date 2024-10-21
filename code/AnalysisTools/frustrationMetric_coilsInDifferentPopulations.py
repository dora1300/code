"""
@Title:             frustrationMetric_coilsInDifferentPopulations.py
@Name:              Mando A Ramirez
@Date:              2024 10 16

@Description:       This script calculates a frustration metric, or something can be used as a frustration
metric, for slab simulations of CC proteins (will also work on single molecule data or single time frame
data).

This script will run on the entire contacts.txt. There is no option to change which frames are analyzed.
You can choose what to do with this information at the end.
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import time


"""
Set up the argparser!
"""
parser = argparse.ArgumentParser(description="This script calculates one type of a Frustration Metric: "
                                 "The number of coils that are (1) unbound; (2) bound in inter-chain "
                                 "interactions; and (3) bound in intra-chain interactions.")
parser.add_argument("-contacts_data", help="Contacts list data file multimerization analysis -- list of coil multimer"
                                  " contacts, each line marks an ANALYZED frame (.csv)."
                                  " Please provide path if not in working directory.", required=True)
parser.add_argument("-coils_per_protein", help="The number of coil segments per protein. CURRENTLY, this analysis only"
                    " supports simulations that use proteins of the same type, i.e. all proteins are C6L5. "
                    "Mixed coils-per-protein types will be added later.", required=True, type=int)
parser.add_argument("-total_coils", help="The total number of coils in the simulation", required=True,
                    type=int)
# parser.add_argument("-ft", help="[ps] Frame time, i.e. the amount of time that each frame in the Contacts Data File "
#                                 "is worth in simulation. Integers only. Only necessary for making the final graph",
#                             type=int, default=None)
parser.add_argument("-output_name", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)


args = parser.parse_args()
CONTACTS_FILE = args.contacts_data
CPP = args.coils_per_protein
TOT_COILS = args.total_coils
output_name = args.output_name



"""
Where the magic happens!
Here I will parse through the _contacts.txt file and analyze the multimers in each desired time frame
and count the frustration metric from these data.
"""

coils_not_in_multimers_list = []
coils_in_intrachain_list = []
coils_in_interchain_list = []


with open(CONTACTS_FILE, 'r') as DATA:
    # here is the set up for the analysis...
    skip_row = 0
    for FRAME in DATA:
        # set up the variables I need for the frame-by-frame analysis
        ncoils_IN_multimer = 0
        ncoils_in_intrachain = 0
        ncoils_in_interchain = 0
        multimers_in_frame = FRAME.rstrip("\n").split("  ") # all multimers are separated by a double space

        for multimers in multimers_in_frame:
            try:
                # I am stealing the next bit from the contact_analysis.py script, to turn all the indices into 
                # integers for easy handling!
                # This turns all the coil indices into integers, using the map() function
                # map() applies a function to all items in a list (which here is generated using
                # .split()) and then returns an iterator, which I put into a list using list()
                coil_indices = list(map(int, multimers.split(" ")))
            except:
                # this will happen if there is no interaction in the given frame. This doesn't necessarily
                # mean the end of the file, though! So just keep iterating through
                continue

            # First, I want to count how many coils are in the selected multimer. I will do this for all the 
            # multimers in the frame
            ncoils_IN_multimer += len(coil_indices)
            
            # Now, this is the part that I will actually check what kind of interaction is happening for the coils
            for i, coil in enumerate(coil_indices):
                coil_i_in_intrachain = 0
                coil_i_in_interchain = 0
                for j, coil2 in enumerate(coil_indices):
                    if i == j:
                        # this is the case where the selected coils are the same in the multimer. It doesn't make
                        # any sense to evaluate a coil interaction when looking at the same coil!
                        continue
                    elif (coil_indices[i] % CPP == 0) and (coil_indices[j] % CPP == 0):
                        # this is an interchain interaction. If both coils modulo coils-per-protein = 0
                        # then it has to be interchain because that means the coils are the start of two different proteins.
                        coil_i_in_interchain = 1
                        continue
                    elif int(coil_indices[i]/CPP) == int(coil_indices[j]/CPP):
                        # this is a unique check. If two coil indicies are within the same protein, then 
                        # the int() of the division coil/CPP will produce the same integer. Otherwise, if the 
                        # two coils are on different proteins, then it will return false
                        coil_i_in_intrachain = 1
                        continue
                    else:
                        coil_i_in_interchain = 1
                
                # Here's how to make sense of this
                # the first loop cycles through the coils in the multimer. Then, the second loop cycles through EVERY
                # other coil *in the multimer* and asks is this interaction intra- or inter-chain for coil_i
                # But, the variables coil_i_in_intrachain and coil_i_in_interchain act as switches, so that 
                # they can only be set to 1 once in the multimer, because coil_i can only contribute itself
                # ONCE to the pool of inter and intrachain coils
                # clear as mud?
                ncoils_in_intrachain += coil_i_in_intrachain
                ncoils_in_interchain += coil_i_in_interchain

        coils_in_interchain_list.append(ncoils_in_interchain)
        coils_in_intrachain_list.append(ncoils_in_intrachain)
        coils_not_in_multimers_list.append(TOT_COILS - ncoils_IN_multimer)


"""
At this point, the analysis is finished. Congrats!
I'm going to save out all the data into a nice .csv file, which I can then do something with.
I will also make a cursory graph for easy inspection of the data.
"""

fraction_coils_in_interchain = np.array(coils_in_interchain_list)/TOT_COILS
fraction_coils_in_intrachain = np.array(coils_in_intrachain_list)/TOT_COILS
fraction_coils_NOT_multimer = np.array(coils_not_in_multimers_list)/TOT_COILS

# this is print debugging. Commented out but preserved in case there's more problems lol.
# print("intra", coils_in_intrachain_list)
# print("inter", coils_in_interchain_list)
# print(coils_not_in_multimers_list)
# print()
# print("intra", fraction_coils_in_intrachain)
# print("inter", fraction_coils_in_interchain)
# print(fraction_coils_NOT_multimer)

# exit()

# Save out the data into a .csv file!
# Importantly, I am also going to save the raw counts per frame, as well as save out
# the fraction of each population relative to the total number of coils in the simulation
np.savetxt(f"{output_name}_frustrationMetric_NUMBER_coils_in_different_populations_.csv", 
           np.array([coils_in_intrachain_list, coils_in_interchain_list, coils_not_in_multimers_list]).T,
           delimiter=",")
np.savetxt(f"{output_name}_frustrationMetric_FRACTION_coils_in_different_populations.csv", 
           np.array([fraction_coils_in_intrachain,
                     fraction_coils_in_interchain,
                     fraction_coils_NOT_multimer]).T,
           delimiter=",")



# Now, make a time course graph of everything!
plt.rcParams['font.size'] = 14
fig, ax = plt.subplots(figsize=(8, 5))
plt.plot(np.arange(1, len(fraction_coils_in_intrachain)+1), fraction_coils_in_intrachain,
         color="sienna", marker=".", linewidth=1, ms=2.5,
         label="Intra-chain coils")
plt.plot(np.arange(1, len(fraction_coils_in_interchain)+1), fraction_coils_in_interchain,
         color="seagreen", marker="^", linewidth=1, ms=2.5,
         label="Inter-chain coils")
plt.plot(np.arange(1, len(fraction_coils_NOT_multimer)+1), fraction_coils_NOT_multimer,
         color="goldenrod", marker="+", linewidth=1, ms=2.5,
         label="Coils not in a multimer")
plt.xlabel("Simulation frame")
plt.ylabel("Fraction of total coils")
plt.grid(linestyle=":", color="black", alpha=0.35)
plt.ylim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig(f"{output_name}_frustrationMetric_FRACTION_coils_in_different_populations.png")
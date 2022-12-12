"""
@Title:         ncf_analysis_v1.py
@Author:        Mando Ramirez
@Date:          20220103

@Description:   This code implements a native contact fraction analysis for coil helices to calculate the fraction of a
simulation when a coil helix is making or not making native contacts. Native contacts are defined based off of reference
distances calculated between heptad pairs, using the energy minimized coil as the reference.
    This code uses a precalculated list of distances for each heptad pair (i.e. A1-A2, B1-B2) that was produced by
calculating the distances in the energy minimized coil. These distances are hard coded.

    THIS CODE IS SPECIFICALLY FOR DAR3-11. THIS IS NOT FOR ANY OTHER ANALYSIS. THIS WILL BE EXTENDED INTO ANOTHER GENERAL
    SCRIPT AFTER DAR3-11 IS COMPLETE.

@Updates:
"""
import os
import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import numpy as np
import prody as pd
import math as m


"""
Hard-coded reference list of distances

These distances were calculated using all 72 energy minimized structures from sim_round_1 and sim_round_2 from 
DAR3-11. Each item corresponds to the distance between a specific heptad pair and goes in order from pairs starting 
at atom 0 onward. To see the numbers, reference file "round1-2_emDistance_calculations_20220103.xlsx"
"""
list_of_distances =[
    10.5794,
    10.7123,
    10.7263,
    10.572,
    10.6313,
    10.7594,
    10.6379,
    10.5696,
    10.7206,
    10.7162,
    10.5824,
    10.6319,
    10.7493,
    10.6368,
    10.58,
    10.7129,
    10.7147,
    10.5723,
    10.6414,
    10.7583,
    10.6477,
    10.5698,
    10.7211,
    10.7265,
    10.5819,
    10.6327,
    10.7499
]

"""
Set up of argument parser, and parsing of arguments
"""
parser = argparse.ArgumentParser(description="NCF analysis of coil simulations. This defines NC as the distance of "
                                             "heptad pairs (i.e. A1-A2, B1-B2) that are within acceptable variation "
                                             "range compared to a reference distance.")
parser.add_argument("-t", help="trajectory file with extension (.xtc/.pdb)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-f", help="analyze every f'th frame in the trajectory", default=10, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_NCF")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="Boolean switch to control plotting. Default is True.", action="store_true")

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially

short_output = args.o + ".csv"
full_output = args.o + "_full.csv"
short = open(short_output, 'w')
full = open(full_output, 'w')

# Open the trajectory and load it
traj = md.load(trajectory, top=topology)


# Create a list of atom pairs upon which to calculate distances. This is because MDTRAJ can do the compute_distances
# function on an array of atom indices
atom_pairs = []
for i in range(0, (lhelix-1)-7):
    # the complicated range is necessary to account for lhelix being 1 indexed whereas mdtraj does 0-indexing
    atom_pairs.append([i, i+7])


# Set the TOLERANCE_VALUE
# this acts as the distance tolerance. Essentially, if the simulated distance is larger than reference
# distance by TOLERANCE_VALUE then there is not a native contact.
TOLERANCE_VALUE = 0.200

"""
Main section - calculation of NCF
"""
trajectory_ncf = []      # this is a list to store the NCF for each frame analyzed

for n in range(nmodels):
    # this is to account for multiple models, which I might have in my simulations!
    # this is also where I will calculate the actual NCFs, in the following loop

    for frame in range(0, traj.n_frames, args.f):
        pairs_distances = md.compute_distances(traj[frame], np.array(atom_pairs))
        native_contact_count = 0; ncf_frame = 0

        for d in range(len(pairs_distances[0])):
            # this loops through the distances calculated for the given frame and counts the number of native contacts
            difference = m.abs(pairs_distances[0][d] - list_of_distances[d])
            if difference <= TOLERANCE_VALUE:
                native_contact_count += 1

        ncf_frame = native_contact_count / len(list_of_distances)
        trajectory_ncf.append(ncf_frame)



ensemble_average = np.average(np.array(trajectory_ncf))
ensemble_stddev = np.std(np.array(trajectory_ncf))

"""
Now, write all the data into the output files
"""
short.write("")


# Close the files!
short.close()
full.close()
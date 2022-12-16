"""
@Author:            Dominique A Ramirez
@Date:              2022 12 16
@Title:             fraction_free_coil_analysis_alpine.py

@Description:       This script calculates the "fraction of free coils" in a given simulation of coils and reports that value throughout time, and as a time averaged value. This is useful to determine how many unpaired coil segments exist throughout a simulation and is valid to use for NPT, NVT, or production simulation.

The script calculates cluster sizes using an index file which specifies beads belonging to a coil. It is important that the user provides an .ndx file that contains atoms **ONLY OF THE COIL SEGMENTS** and that all the coil segments use *the same number* of atoms for calculation.

This script is specifically for use on RMACC's Alpine super computer.

@Updates:

"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

parser = argparse.ArgumentParser(description="Tool to analyze fraction free coil globally for "
    "a coil simulation. This script uses 'gmx clustsize' and ALWAYS uses pbc correction. This "
    "feature cannot be turned off."    
)
parser.add_argument('-t', help="(.xtc) Trajectory file to perform analysis on. This trajectory "
    "does not have to be pbc corrected.")

parser.add_argument('-s', help="(.tpr) Gromacs run-file, used as the topology")

parser.add_argument('-ndx', help="(.ndx) Gromacs index file containing the atoms of the COIL "
    "A-BEADS that are to be analyzed. Be sure this is correct. There should only be 1 group "
    "in the index file.")
    
parser.add_argument('-n', help="The number of atoms per coil segment in the index file. This is "
    "critical to get right because the method depends on this value being correct.",
    required=True, type=int)
    
parser.add_argument('-num_coil', help="The total number of coil segments in the system, or if "
    "you specified differently in the ndx file, the total number of coil segments in the .ndx "
    "file.", required=True, type=int)

parser.add_argument('-dt', help="[ps] Analyze every dt ps frame in the trajectory. Default is every 10000 frames to speed up calculations.", 
    default=10000, type=int)
    
parser.add_argument('-b', help="[ps] Starting simulation time to analyze density. Default = 0 ps.", 
    default=0, type=int)
    
parser.add_argument('-e', help="[ps] Ending simulation time. No default.", 
    required=True, type=int)
    
parser.add_argument('-cutoff', help="[nm] The cut-off used to calculate neighbors.", 
    required=True, type=float)

parser.add_argument('-output', help="Name of the output analysis file. Please provide .csv extension", default="output.csv")

parser.add_argument('-path', help="The full path to the xvg2csv.py script. This can be done "
    "better but for now this is what I have. Provide the xvg2csv.py as well.", required=True)


args = parser.parse_args()

xtc = args.t
tpr = args.s
ndx = args.ndx
dt = args.dt
co = args.cutoff
nc = args.num_coil

os.mkdir("./temp_clusters")
time_array = np.arange(args.b, args.e+1, args.dt)
output_file = open(args.output, 'w')
output_file.write("#Time(ps),cluster_size=0_coil,cluster_size=1_coil,cluster_size=2_coil, etc.")

for time in time_array:
    TIME = int(time)
    frame_coils = []
    gmx_clustsize = f"mpirun -np 1 gmx_mpi clustsize -f {xtc} -s {tpr} -n {ndx} -hc ./temp_clusters/hist_{TIME}.xvg -b {TIME} -e {TIME} -pbc yes -cut {co}"
    subprocess.run(gmx_clustsize.split())
    os.remove("avclust.xvg")
    os.remove("csize.xpm")
    os.remove("csizew.xpm")
    os.remove("maxclust.ndx")
    os.remove("maxclust.xvg")
    os.remove("nclust.xvg")
    os.remove("temp.xvg")


    convert_xvg = f"python3 {args.path} -f ./temp_clusters/hist_{TIME}.xvg -o ./temp_clusters/hist_{TIME}.csv"
    subprocess.call(convert_xvg.split())

    output_file.write(f"{TIME}")
    with open(f"./temp_clusters/hist_{TIME}.csv") as f:
        for line in f:
            if line[0] == "#":
                continue
            elif line[0] == "@":
                continue
            else:
                line_split = line.rstrip("\n").split(",")
                # only add cluster size values that correspond to whole coils
                
                if float(line_split[0]) % args.n == 0:
                    coil_oligo = int(line_split[0]) / 10.
                    frac = (float(line_split[1]) * coil_oligo) / float(nc)
                    output_file.write(f",{frac:.6f}")
                    
                else:
                    pass

    output_file.write("\n")

output_file.close()


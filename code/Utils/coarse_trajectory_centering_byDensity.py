"""
Name:               coarse_trajectory_centering_byDensity.py
Author:             Dominique A Ramirez
Date:               2024 01 22

[provide description of script here]


Updates:
"""

import subprocess
import os
import argparse
import numpy as np
import mdtraj as md



def read_density(FILE):
    zcoord_array = np.zeros(args.sl, dtype=float)
    density_array = np.zeros(args.sl, dtype=float)

    with open(FILE) as f:
        for line in f:
            line_count = 0
            if line[0] == "#":
                continue
            elif line[0] == "@":
                continue
            else:
                sp = line.lstrip(" ").rstrip("\n").split(" ")
                data = []
                for i in sp:
                    if len(i) == 0:
                        continue
                    else:
                        data.append(float(i))
                zcoord_array[line_count] = data[0]
                density_array[line_count] = data[1]
                line_count += 1
                data.clear()
    return zcoord_array, density_array


"""
Set up the argument parser and parse the arguments
"""
parser = argparse.ArgumentParser(description="Tool for coarsely 'centering' a LLPS cluster in a simulation box."
                                 " This works by finding the highest density in a coarse (i.e. low slice no.) "
                                 "density profile, translating the system so the coarse max position is now "
                                 "the center of the box, then repeating for all frames.")
parser.add_argument('-t', help="Input trajectory file for centering. This MUST be modified using pbc nojump "
                    "and pbc whole prior to centering.", 
                    required=True, type=str
)
parser.add_argument('-s', 
        help="A GROMACS .tpr run file for the associated trajectory. PROVIDE EXTENSION.", 
        required=True, type=str
)
parser.add_argument('-boxZ',
        help="The z-coordinate corresponding to the center of the box, in the z-direction.",
        required=True, type=float
)
parser.add_argument('-sl',
        help="The number of slices to cut up the density profile. This will be the same value "
        "passed to `gmx density -sl`. Please only provide an integer",
        required=True, type=int
)
parser.add_argument('-start', 
        help="[tu] The time at which to START the analysis.", 
        required=True, type=int)
parser.add_argument('-end', 
        help="[tu] The time at which to STOP the analysis.", 
        required=True, type=int)
parser.add_argument('-dt', 
        help="The intervals at which to do the analysis. ", 
        required=True, type=int, default=1)
parser.add_argument('-tu', 
        help="Time unit for the start and stop times. Must the same time unit for both times. "
        "Accepted values are 'ps' , 'ns' , 'us' ",
        required=True, type=str)
parser.add_argument('-mpi',
        help="Pass this flag to turn on the mpi settings for all the gmx functions",
        action="store_true",
        default=False)
        

args = parser.parse_args()

TRAJ = args.t
TPR = args.s

SLICE = args.sl
BOXZ = args.boxZ

START_FRAME = args.start
STOP_FRAME = args.end
FRAME_ITER = args.dt
TIME_UNIT = args.tu


# Define the TIME_FACTOR, which is important for multiplying things below
# in other words, what do you need to multiply the frames by to get
# to ps?
if TIME_UNIT == "us":
    TIME_FACTOR = 1000000
elif TIME_UNIT == "ns":
    TIME_FACTOR = 1000
elif TIME_UNIT == "ps":
    TIME_FACTOR = 1
else:
    print("The Time Unit (`-tu`) option you selected is incompatible with accepted values.")
    print("Please try again and choose only an acceptable value for `-tu`")
    print("(your provided value: {args.tu})")
    exit(1)


# Handle what kind of function to call, gmx or gmx_mpi
if args.mpi:
    FUNC = "mpirun -np 1 gmx_mpi"
else:
    FUNC = "gmx"


# This is just a silly print statement to make the user see what parameters are chosen
# more verbosity is never a bad thing...
print()
print("Performing coarse centering on files: ")
print(f"Trajectory: {TRAJ} | .tpr {TPR}")
print(f"Start time: {START_FRAME} {TIME_UNIT} | End time: {STOP_FRAME} {TIME_UNIT}")
print(f"Centering every {FRAME_ITER} frame")
print(f"Box dimensions: {BOXZ}")
print()




# loop through every DESIRED frame
for FRAME in range(START_FRAME, STOP_FRAME+FRAME_ITER, FRAME_ITER):
    # Step 1 -- calculate the density profile for the desired frame
    # and save to the special density trajectory
    density = (f"echo 0 | {FUNC} density -f {TRAJ} -s {TPR} "
               f"-b {FRAME} -e {FRAME} "
               f"-dens number -sl {SLICE} "
               f"-o ./density_untranslated/density_{FRAME}{TIME_UNIT}.xvg")
    subprocess.call(density, shell=True)

    
    # step 2 -- read the density file and find the zcoordinate location of the max
    # density
    zpositions, densities = read_density(f"./density_untranslated/density_{FRAME}{TIME_UNIT}.xvg")
    print(zpositions)
    print(densities)
    print(np.where(densities == np.max(densities))[0][0])
    max_index = np.where(densities == np.max(densities))[0][0]
    trans_z_by = BOXZ - zpositions[max_index]

    print(f"Slice with the largest density: {zpositions[max_index]}")
    print(f"Translating frame {FRAME}{TIME_UNIT} by {trans_z_by}")


    # step 3 -- translate the ENTIRE system by the translation coordinates and save out using -pbc mol!
    translated_frame = (f"echo 0 | {FUNC} trjconv "
            f"-f {TRAJ} " 
            f"-s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
            f"-trans 0.0 0.0 {trans_z_by} "
            f"-pbc mol "
            f"-o ./coarse_trans_traj/frame_{FRAME}{TIME_UNIT}_coarse_transl.xtc")
    subprocess.call(translated_frame, shell=True)



 
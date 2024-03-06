"""
Name:               coarse_trajectory_centering_byDensity.py
Author:             Dominique A Ramirez
Date:               2024 01 22

This script performs what I'm calling as a "coarse trajectory centering" where a LLPS slab 
simulation system is centered by translating the region of the highest density. The region
of highest density is determined by a "coarse density profile" where a simulation box is
first divided into slices of large width (approx. 15-25 nm slices). The approximate 
z-dimensions corresponding to the highest density is then centered in the simulation box
by translating the entire system.

IMPORTANTLY -- the updated version 


Updates:
2024 03 04  -- updated to do coarse analysis ONLY on a system with the largest cluster.
    Also changed one of the directory names!
"""

import subprocess
import os
import argparse
import numpy as np
#import mdtraj as md



def write_csv(FILENAME):
    out_file = open(f"./coarse_density_profile/{FILENAME}.csv", 'w')
    with open(f"./coarse_density_profile/{FILENAME}.xvg", "r") as f:
        for line in f:
            if line[0] == "#":
                out_file.write(line)
            elif line[0] == "@":
                out_file.write(line)
            else:
                sp = line.lstrip(" ").rstrip("\n").split(" ")
                for i in range(len(sp)):
                    if len(sp[i]) == 0:
                        continue
                    elif i != (len(sp)-1):
                        out_file.write(f"{sp[i]},")
                    else:
                        out_file.write(f"{sp[i]}")
                out_file.write("\n")
    out_file.close()
    return None


def read_density(FILE):
    zcoord_array = np.zeros(args.sl, dtype=float)
    density_array = np.zeros(args.sl, dtype=float)

    with open(FILE) as f:
        line_count = 0
        for line in f:
            if line[0] == "#":
                continue
            elif line[0] == "@":
                continue
            else:
                line_split = line.rstrip("\n").split(",")
                zcoord_array[line_count] = float(line_split[0])
                density_array[line_count] = float(line_split[1])
                line_count += 1
    return zcoord_array, density_array


def remove_files():
    try:
        os.remove("./cluster_analysis/avclust.xvg")
    except:
        pass
    try:
        os.remove("./cluster_analysis/csize.xpm")
    except:
        pass
    try:
        os.remove("./cluster_analysis/csizew.xpm")
    except:
        pass
    try:
        os.remove("./cluster_analysis/histo-clust.xvg")
    except:
        pass
    try:
        os.remove("./cluster_analysis/maxclust.xvg")
    except:
        pass
    try:
        os.remove("./cluster_analysis/nclust.xvg")
    except:
        pass
    try:
        os.remove("./cluster_analysis/temp.xvg")
    except:
        pass
    return None

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
parser.add_argument('-n',
        help="File corresponding to the default index file (.ndx) that comes from GROMACS. "
        "This file MUST exist in ./index_files/ in reference to where the trajectory files exist. "
        "Please include the extension.",
        default="standard.ndx",
        type=str)
        

args = parser.parse_args()

TRAJ = args.t
TPR = args.s
NDX = args.n

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



"""
Critical section!
Since I'm capturing output from the function calls, if there is an error
I won't see it well.

An easy place for this script to fail is when the necessary directories are missing.
Check each necessary directory one by one and make sure it exists
"""
if os.path.isdir("index_files"):
    pass
else:
    os.mkdir("index_files")

if os.path.isdir("cluster_analysis"):
    pass
else:
    os.mkdir("cluster_analysis")

if os.path.isdir("coarse_density_profile"):
    pass
else:
    os.mkdir("coarse_density_profile")

if os.path.isdir("coarse_trans_traj"):
    pass
else:
    os.mkdir("coarse_trans_traj")


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
    # Step 1 -- calculate the clustsize to determine the biggest cluster, if it exists,
    # and save the cluster index file of the largest cluster
    if os.path.isfile(f"./index_files/max_{FRAME}{TIME_UNIT}.ndx"):
        # this is when I've already done the cluster analysis and make the index files 
        # (they don't need to be remade obviously)
        pass
    else:
        os.chdir("./cluster_analysis/")
        clustsize = (f"{FUNC} clustsize -f ../{TRAJ} -s ../{TPR} "
                    f"-mcn ../index_files/max_{FRAME}{TIME_UNIT}.ndx " 
                    f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -mol -cut 0.9 -pbc")
        try:
            # if there is a cluster that is < N_all_molecules, then this will work
            # and produce some output from the clustsize analysis
            subprocess.run(clustsize.split(), check=True, capture_output=True)
        except:
            # Updated 01.22.24
            # The time this will fail is if ALMOST EVERY molecule is in a single cluster, i.e.
            # when all molecules are dispersed and not touching each other, which
            # gmx clustsize doesn't like. 
            # BUT, sometimes it will give an error and still make an index file if there's, like,
            # 99% of molecules in a single cluster then 1% in a smaller cluster
            # so I'm changing the exception checking
            print(f"Frame time {FRAME}{TIME_UNIT} produced an error with gmx clustsize ")
            print("Continuing with analysis and will check if an index file was still created.")
            
        # clean up the cluster_analysis files regardless of what happens
        os.chdir("../")
        remove_files()

    
    # Check to verify if gmx clustsize produced the expected max_{FRAME}{TIME_UNIT}.ndx index
    # file in the index_files directory. Sometimes, even if clustsize errors out, the index
    # file will still be created and I need to handle that.
    # notice that I'm in the head directory at this point, just as I should be
    if os.path.isfile(f"./index_files/max_{FRAME}{TIME_UNIT}.ndx"):
        # Make an index file specific for the largest cluster for use down the road
        if os.path.isfile(f"./index_files/all_{FRAME}{TIME_UNIT}.ndx"):
            pass
        else:
            os.chdir("./index_files/")
            concat = f"cat {NDX} max_{FRAME}{TIME_UNIT}.ndx >> all_{FRAME}{TIME_UNIT}.ndx"
            subprocess.run(concat, shell=True, capture_output=True)
            os.chdir("../")

        # Now, I have to set the CORRECT VARIABLE corresponding to the largest cluster
        MAXNDX = f"max_{FRAME}{TIME_UNIT}.ndx"

        # Now, run the density analysis on ONLY the largest cluster and save the output as normal
        density = (f"echo 10 | {FUNC} density -f {TRAJ} -s {TPR} "
                f"-n ./index_files/{MAXNDX} "
                f"-b {FRAME} -e {FRAME} "
                f"-dens number -sl {SLICE} "
                f"-o ./coarse_density_profile/density_{FRAME}{TIME_UNIT}.xvg")
        subprocess.run(density, shell=True, capture_output=True)
        # convert the file to .csv for easy reading
        write_csv(f"density_{FRAME}{TIME_UNIT}")

    
    else:
        # the max_{FRAME}{TIME_UNIT}.ndx index file was not created. This means the 
        # entire system is in a cluster and I can run the density analysis on the entire system
        density = (f"echo 0 | {FUNC} density -f {TRAJ} -s {TPR} "
                f"-b {FRAME} -e {FRAME} "
                f"-dens number -sl {SLICE} "
                f"-o ./coarse_density_profile/density_{FRAME}{TIME_UNIT}.xvg")
        subprocess.call(density, shell=True)
        # convert the file to .csv for easy reading
        write_csv(f"density_{FRAME}{TIME_UNIT}")
    

    # Regardless of which option is chose, a .csv file is written in ./coarse_density_profile/  
    # Read the density file and find the zcoordinate location of the max density
    zpositions, densities = read_density(f"./coarse_density_profile/density_{FRAME}{TIME_UNIT}.csv")
    print(f"Z-position slices:")
    print(zpositions)
    print(f"Densities of z-position slices:")
    print(densities)
    print(np.where(densities == np.max(densities))[0][0])
    max_index = np.where(densities == np.max(densities))[0][0]
    trans_z_by = BOXZ - zpositions[max_index]


    # A nice little debug statement
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
    subprocess.run(translated_frame, 
                   shell=True,
                   capture_output=True)

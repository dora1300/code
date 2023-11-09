"""
Name:               coordinate_based_trajectory_centering_mpi.py
Author:             Dominique A Ramirez
Date:               2023 09 20

This script performs centering on a trajectory so that the slab, if one exists, is at the center
of the box in all specified frames.

This will not use mdtraj so I don't have to load the trajectory into memory.



Updates:
2023 10 20          Updated to include an argument parser for general use!

2023 10 22          This has updated commands to be able to run on a supercomputer, if needed!

2023 11 09          Updated to remove the hard coded boundaries for verifying densities.
"""

import subprocess
import os
import numpy as np
import argparse

def read_density(FILENAME, low_point, high_point):
    slices = []; densities = []

    with open(FILENAME, 'r') as f:
        for line in f:
            if line[0] == "@" or line[0] == "#":
                continue
            else:
                sp = line.lstrip(" ").rstrip("\n").split(" ")
                for i in range(len(sp)):
                    if len(sp[i]) == 0:
                        continue
                    elif i != (len(sp)-1):
                        slices.append(float(sp[i]))
                    else:
                        densities.append(float(sp[i]))
   
    exterior = []; interior = []
    for i, dens in enumerate(densities):
        if i < low_point:
            exterior.append(dens)
        elif i >= high_point:
            exterior.append(dens)
        else:
            interior.append(dens)
    
    ext_max = np.max(np.array(exterior))
    int_max  = np.max(np.array(interior))
    if (ext_max / int_max) > 0.75 and (ext_max / int_max) < 1.25:
        box_center = "rect"
    elif (int_max / ext_max) >= 1.25:
        box_center = "rect"
    else:
        box_center = "zero"
     
    return box_center

def check_density(FILENAME, low_point, high_point):
    slices = []; densities = []

    with open(FILENAME, 'r') as f:
        for line in f:
            if line[0] == "@" or line[0] == "#":
                continue
            else:
                sp = line.lstrip(" ").rstrip("\n").split(" ")
                for i in range(len(sp)):
                    if len(sp[i]) == 0:
                        continue
                    elif i != (len(sp)-1):
                        slices.append(float(sp[i]))
                    else:
                        densities.append(float(sp[i]))
   
    exterior = []; interior = []
    for i, dens in enumerate(densities):
        # more strict check for what constitutes the edge of the box
        if i < low_point-5:
            exterior.append(dens)
        elif i >= high_point+5:
            exterior.append(dens)
        else:
            interior.append(dens)
    
    ext_max = np.max(np.array(exterior))
    int_max  = np.max(np.array(interior))

    # Here is the return statement. Return 0 == the interior has the highest density
    # Return 1 == the exterior has the highest density.
    if (int_max / ext_max) > 0.90 and (int_max / ext_max) < 1.10:
        # interior has higher density, but barely (within acceptable error)
        # strict error tolerance
        return 0
    elif (int_max / ext_max) >= 1.10:
        return 0
    else:
        return 1

def remove_files():
    os.remove("./cluster_analysis/avclust.xvg")
    os.remove("./cluster_analysis/csize.xpm")
    os.remove("./cluster_analysis/csizew.xpm")
    os.remove("./cluster_analysis/histo-clust.xvg")
    os.remove("./cluster_analysis/maxclust.xvg")
    os.remove("./cluster_analysis/nclust.xvg")
    os.remove("./cluster_analysis/temp.xvg")
    return None


"""
Set up the argument parser and parse the arguments
"""
parser = argparse.ArgumentParser(description="Tool for performing coordinate-based trajectory centering, which centers "
                                 "the trajectory on the center-of-mass of the largest cluster in any given frame. "
                                 "This writes *MANY* single frame .xtc files and then requires concatenating the "
                                 "many files. The finished product should be a trajectory where the largest cluster "
                                 "is in the middle of the box.")
parser.add_argument('-t', help="Input trajectory file for centering. PBC and other standard traj corrections "
                    "do not need to be made. PROVIDE EXTENSION.", required=True, type=str)
parser.add_argument('-s', 
        help="A GROMACS .tpr run file for the associated trajectory. PROVIDE EXTENSION.", 
        required=True, type=str)
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
parser.add_argument("-extLow",
        help="What is the *lower* value of the density profile slices corresponding to the "
        "dilute regime? e.g. in a box with 100 slices in the density, this value could be `15`.",
        required=True, type=int)
parser.add_argument("-extHigh",
        help="What is the *higher* value of the density profile slices corresponding to the "
        "dilute regime? e.g. in a box with 100 slices in the density, this value could be `15`.",
        required=True, type=int)
parser.add_argument("-sl", 
        help="How many slices do you want to divide the density profile into? This parameter "
        "will be passed into the `gmx_density -sl` function",
        required=True, type=int)


args = parser.parse_args()

TRAJ = args.t
TPR = args.s

START_FRAME = args.start
STOP_FRAME = args.end
FRAME_ITER = args.dt
TIME_UNIT = args.tu

LOW = args.extLow
HIGH = args.extHigh

SL = args.sl


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



"""
Here is where the centering takes place!!
The general procedure is as follows:
* for each frame in the trajectory*
    1. Calculate the index of the LARGEST cluster of the frame with `gmx clustsize`
    2. Prepare the primary index file with the cluster indices
    3. Calculate the density to determine what type of box to use
    4. Center the trajectory on the largest cluster using index file from (2)
    5. Verify that the centering worked! and that the largest 
"""
# loop through every DESIRED frame
for FRAME in range(START_FRAME, STOP_FRAME+FRAME_ITER, FRAME_ITER):
    print("######################      STEP 1")
    # Step 1 -- calculate the clustsize to determine the biggest cluster, if it exists
    os.chdir("./cluster_analysis/")
    clustsize = (f"mpirun -1 np 1 gmx_mpi clustsize -f ../{TRAJ} -s ../{TPR} "
                 f"-mcn ../index_files/max_{FRAME}{TIME_UNIT}.ndx " 
                 f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -mol -cut 0.9 -pbc")
    try:
        # if there is a cluster that is < N_all_molecules, then this will work
        # and produce some output from the clustsize analysis
        subprocess.run(clustsize.split(), check=True)
    except:
        # The time this will fail is if EVERY molecule is in a single cluster, which
        # gmx clustsize doesn't like. But this also means I can just use Group 1 to do
        # the centering, since it means that all protein is in a single cluster and is
        # already in a group to manipulate
        print(f"Frame time {FRAME}{TIME_UNIT} failed to calculate the cluster size. Moving on in the analysis.")
        os.chdir("../")

        print("######################      STEP 1-a")
        # Step 3 - calculate the density and determine what type of box center option to use
        density = (f"echo 1 | mpirun -np 1 gmx_mpi density -f {TRAJ} -s {TPR} "
                   f"-b {FRAME*TIME_FACTOR} -e {FRAME*TIME_FACTOR} "
                   f"-sl {SL} -dens number -o dens_{FRAME}.xvg")
        subprocess.call(density, shell=True)
        BOX_CENTER = read_density(f"dens_{FRAME}.xvg", LOW, HIGH)
        print(f"Initial selected box center spec: {BOX_CENTER}")
        os.remove(f"dens_{FRAME}.xvg")
        
        print("######################      STEP 1-b")
        # Step 4 - do the trajectory correction ON THE PROTEIN GROUP 1
        conv = (f"echo 1 1 | mpirun -np 1 gmx_mpi trjconv -f {TRAJ} -s {TPR} -b {FRAME} "
                f"-e {FRAME} -tu {TIME_UNIT} -pbc mol -center "
                f"-boxcenter {BOX_CENTER} -n ./index_files/standard.ndx -o ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
        subprocess.call(conv, shell=True)
        
        print("######################      STEP 1-c")
        # Step 5 - this is the quality check. This is where I see if the frame was correctly centered
        density_verify = (f"echo 1 | mpirun -np 1 gmx_mpi density "
                          f"-f ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc -s {TPR}  "
                          f"-b {FRAME*TIME_FACTOR} -e {FRAME*TIME_FACTOR} -sl {SL} "
                          f"-dens number -o dens_{FRAME}_verify.xvg")
        subprocess.call(density_verify, shell=True)
        verify_results = check_density(f"dens_{FRAME}_verify.xvg", LOW, HIGH)
        if verify_results == 0:
            # this means the interior has the highest density, which means the slab has hopefully been placed
            # in the center of the simulation box
            # this means do nothing! The initial attempt was correct
            os.remove(f"dens_{FRAME}_verify.xvg")
            continue
        else:
            # this means the initial centering attempt was incorrect, and I have to redo the centering.
            # first, figure out what the old box center parameter was and switch it up
            if BOX_CENTER == 'rect':
                NEW_CENTER = 'zero'
            else:
                NEW_CENTER = 'rect'
            # next, delete the old converted trajectory file
            os.remove(f"./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
            # now, redo the centering
            new_conv = (f"echo 1 1 | mpirun -np 1 gmx_mpi trjconv -f {TRAJ} -s {TPR} "
                        f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -pbc mol -center "
                        f"-boxcenter {NEW_CENTER} -n ./index_files/standard.ndx "
                        f"-o ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
            subprocess.call(new_conv, shell=True)
            continue



    # proceed with the rest of the code if the clustsize DOESNT throw an error
    os.chdir("../")
    remove_files()



    print("######################      STEP 2")
    # Step 2 - make the appropriate index files for the given frame
    os.chdir("./index_files/")
    concat = f"cat standard.ndx max_{FRAME}{TIME_UNIT}.ndx >> all_{FRAME}{TIME_UNIT}.ndx"
    subprocess.call(concat, shell=True)
    os.chdir("../")
 


    print("######################      STEP 3")
    # Step 3 - calculate the density and determine what type of box center option to use
    density = (f"echo 1 | mpirun -np 1 gmx_mpi density -f {TRAJ} -s {TPR} "
               f"-b {FRAME*TIME_FACTOR} -e {FRAME*TIME_FACTOR} "
               f"-sl {SL} -dens number -o dens_{FRAME}.xvg")
    subprocess.call(density, shell=True)
    BOX_CENTER = read_density(f"dens_{FRAME}.xvg", LOW, HIGH)
    print(f"Initial selected box center spec: {BOX_CENTER}")
    os.remove(f"dens_{FRAME}.xvg")


    
    print("######################      STEP 4")
    # Step 4 - do the trajectory correction
    conv = (f"echo 10 1 | mpirun -np 1 gmx_mpi trjconv -f {TRAJ} -s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -pbc mol -center "
            f"-boxcenter {BOX_CENTER} -n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
            f"-o ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
    subprocess.call(conv, shell=True)


    print("######################      STEP 5")
    # Step 5 - this is the quality check. This is where I see if the frame was correctly centered
    density_verify = (f"echo 1 | mpirun -np 1 gmx_mpi density "
                      f"-f ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc -s {TPR} "
                      f"-b {FRAME*TIME_FACTOR} -e {FRAME*TIME_FACTOR} "
                      f"-sl {SL} -dens number -o dens_{FRAME}_verify.xvg")
    subprocess.call(density_verify, shell=True)
    verify_results = check_density(f"dens_{FRAME}_verify.xvg", LOW, HIGH)
    if verify_results == 0:
        # this means the interior has the highest density, which means the slab has hopefully been placed
        # in the center of the simulation box
        # this means do nothing! The initial attempt was correct
        os.remove(f"dens_{FRAME}_verify.xvg")
        continue
    else:
        # this means the initial centering attempt was incorrect, and I have to redo the centering.
        # first, figure out what the old box center parameter was and switch it up
        if BOX_CENTER == 'rect':
            NEW_CENTER = 'zero'
        else:
            NEW_CENTER = 'rect'
        # next, delete the old converted trajectory file
        os.remove(f"./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
        # now, redo the centering
        new_conv = (f"echo 10 1 | mpirun -np 1 gmx_mpi trjconv -f {TRAJ} -s {TPR} "
                    f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -pbc mol -center "
                    f"-boxcenter {NEW_CENTER} -n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
                    f"-o ./traj_files/frame_{FRAME}{TIME_UNIT}.xtc")
        subprocess.call(new_conv, shell=True)




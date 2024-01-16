"""
Name:               trajectory_centering_by_COM.py
Author:             Dominique A Ramirez
Date:               2024 01 16

This script performs centering on a trajectory so that the slab, if one exists, is at the center
of the box in all specified frames.

This differs from previous attempts. Now I am performing centering by means of the COM of the slab.
The general workflow looks something like:
for each frame...
- Find the largest cluster of the simulation and extract out the index files for those molecules
- Save an .xtc file of ONLY the largest cluster, absolutely nothing else
- Calculate the COM coordinates of the largest cluster using mdtraj.compute_center_of_mass()
- Calculate the difference between the COM and the box center
    - transx = box_x - COM_x
    - transy = box_y - COM_y
    - transz = box_z - COM_z
- Make a new .xtc files of THE ENTIRE SYSTEM but translate all atoms by the difference using 
    `gmx trjconv -trans transx transy transz`
- Once finished, compile all of the trajectory files into a single .xtc


Updates:
"""

import subprocess
import os
import numpy as np
import argparse
import mdtraj as md



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
                    "CANNOT BE MADE. Original .xtc direct from simulation output ONLY. PROVIDE EXTENSION.", 
                    required=True, type=str
)
parser.add_argument('-s', 
        help="A GROMACS .tpr run file for the associated trajectory. PROVIDE EXTENSION.", 
        required=True, type=str
)
parser.add_argument('-p',
        help="A 'topology' file used for mdtraj. This is the .gro file DIRECTLY FROM simulation output. "
        "This is only used for connectivity for mdtraj. PROVIDE EXTENSION",
        required=True, type=str
)
parser.add_argument('-boxX',
        help="The x-coordinate corresponding to the center of the box, in the x-direction.",
        required=True, type=float
)
parser.add_argument('-boxY',
        help="The y-coordinate corresponding to the center of the box, in the y-direction.",
        required=True, type=float   
)
parser.add_argument('-boxZ',
        help="The z-coordinate corresponding to the center of the box, in the z-direction.",
        required=True, type=float
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
        

args = parser.parse_args()

TRAJ = args.t
TPR = args.s
TOP = args.p

BOXX = args.boxX
BOXY = args.boxY
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



"""
The general workflow looks something like:
for each frame...
- Find the largest cluster of the simulation and extract out the index files for those molecules
- Save an .xtc file of ONLY the largest cluster, absolutely nothing else
- Calculate the COM coordinates of the largest cluster using mdtraj.compute_center_of_mass()
- Calculate the difference between the COM and the box center
    - transx = box_x - COM_x
    - transy = box_y - COM_y
    - transz = box_z - COM_z
- Make a new .xtc files of THE ENTIRE SYSTEM but translate all atoms by the difference using 
    `gmx trjconv -trans transx transy transz`
- Once finished, compile all of the trajectory files into a single .xtc
"""


# loop through every DESIRED frame
for FRAME in range(START_FRAME, STOP_FRAME+FRAME_ITER, FRAME_ITER):
    print("######################      STEP 1")
    # Step 1 -- calculate the clustsize to determine the biggest cluster, if it exists,
    # and save the cluster index file of the largest cluster
    os.chdir("./cluster_analysis/")
    clustsize = (f"gmx clustsize -f ../{TRAJ} -s ../{TPR} "
                 f"-mcn ../index_files/max_{FRAME}{TIME_UNIT}.ndx " 
                 f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} -mol -cut 0.9 -pbc")
    try:
        # if there is a cluster that is < N_all_molecules, then this will work
        # and produce some output from the clustsize analysis
        subprocess.run(clustsize.split(), check=True)
    except:
        # The time this will fail is if EVERY molecule is in a single cluster, i.e.
        # when all molecules are dispersed and not touching each other, which
        # gmx clustsize doesn't like. 
        # So instead, I will still do a translation, but calculate the COM of the entire SYSTEM
        print(f"Frame time {FRAME}{TIME_UNIT} failed to calculate the cluster size. Moving on in the analysis.")
        os.chdir("../")
        remove_files()

        # Step 2a - No making a specific index file. Just save the frame!
        no_cluster_frame = (f"echo 1 | gmx trjconv -f {TRAJ} -s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
            f"-o ./clust_none_traj/frame_{FRAME}{TIME_UNIT}_clust_none.xtc")
        subprocess.call(no_cluster_frame, shell=True)

        # Step 3a - calculate the COM of the no cluster frame
        traj_no_clust = md.load(f"./clust_none_traj/frame_{FRAME}{TIME_UNIT}_clust_none.xtc",
                                top=TOP)
        no_clust_COM = md.compute_center_of_mass(traj_no_clust)

        # Step 4a - calculate the difference between COM and the actual box center
        transX = BOXX - no_clust_COM[0][0]
        transY = BOXY - no_clust_COM[0][1]
        transZ = BOXZ - no_clust_COM[0][2]

        # Step 5 -- make a new .xtc for just the given frame but TRANSLATE the ENTIRE SYSTEM
        # by the difference calculated in step 4
        # again, NO CENTERING
        translated_no_cluster_frame = (f"echo 0 | gmx trjconv -f {TRAJ} -s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
            f"-trans {transX} {transY} {transZ} "
            f"-o ./translated_trajs/frame_{FRAME}{TIME_UNIT}_transl.xtc")
        subprocess.call(translated_no_cluster_frame, shell=True)
        continue


    # proceed with the rest of the code if the clustsize DOESNT throw an error
    os.chdir("../")
    remove_files()


    # Step 2 -- save an .xtc file of ONLY THE CLUSTER, using the output from clustsize
    # to make an index file corresponding to only the largest cluster
    print("######################      STEP 2")
    os.chdir("./index_files/")
    concat = f"cat standard.ndx max_{FRAME}{TIME_UNIT}.ndx >> all_{FRAME}{TIME_UNIT}.ndx"
    subprocess.call(concat, shell=True)
    os.chdir("../")
    # now save out the .xtc file corresponding to only the cluster. ABSOLUTELY NOOOOOOOOOO
    # CENTERING
    only_cluster_frame = (f"echo 10 | gmx trjconv -f {TRAJ} -s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
            f"-n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
            f"-o ./clust_only_traj/frame_{FRAME}{TIME_UNIT}_clust_only.xtc")
    subprocess.call(only_cluster_frame, shell=True)
    

    # Step 3 -- calculate the COM on the largest cluster using mdtraj
    traj_clust_frame = md.load(f"./clust_only_traj/frame_{FRAME}{TIME_UNIT}_clust_only.xtc", 
                               top=TOP)
    clust_only_COM = md.compute_center_of_mass(traj_clust_frame)

    
    # Step 4 -- calculate the difference between COM and the actual box center
    transX = BOXX - clust_only_COM[0][0]
    transY = BOXY - clust_only_COM[0][1]
    transZ = BOXZ - clust_only_COM[0][2]


    # Step 5 -- make a new .xtc for just the given frame but TRANSLATE the ENTIRE SYSTEM
    # by the difference calculated in step 4
    # again, NO CENTERING
    translated_frame = (f"echo 0 | gmx trjconv -f {TRAJ} -s {TPR} "
            f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
            f"-trans {transX} {transY} {transZ} "
            f"-o ./translated_trajs/frame_{FRAME}{TIME_UNIT}_transl.xtc")
    subprocess.call(translated_frame, shell=True)

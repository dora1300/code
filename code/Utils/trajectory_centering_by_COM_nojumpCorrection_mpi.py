"""
Name:               trajectory_centering_by_COM_nojumpCorrection.py
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
import argparse
import mdtraj as md





"""
Set up the argument parser and parse the arguments
"""
parser = argparse.ArgumentParser(description="Tool for performing coordinate-based trajectory centering, which centers "
                                 "the trajectory on the center-of-mass of the largest cluster in any given frame. "
                                 "This writes *MANY* single frame .xtc files and then requires concatenating the "
                                 "many files. The finished product should be a trajectory where the largest cluster "
                                 "is in the middle of the box.")
parser.add_argument('-t', help="Input trajectory file for centering. This MUST be modified using pbc nojump "
                    "and pbc whole prior to centering.", 
                    required=True, type=str
)
parser.add_argument('-s', 
        help="A GROMACS .tpr run file for the associated trajectory. PROVIDE EXTENSION.", 
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
parser.add_argument('-mpi',
        help="Pass this flag to turn on the mpi settings for all the gmx functions",
        action="store_true",
        default=False)
parser.add_argument('-silent',
        help="Suppress output from GROMACS commands. Use with **caution**! If the centering fails, "
        "you won't know where or why.",
        action="store_true",
        default=False)
        

args = parser.parse_args()

TRAJ = args.t
TPR = args.s

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
if os.path.isdir("stepA_wholesys_mod"):
    pass
else:
    os.mkdir("stepA_wholesys_mod")

if os.path.isdir("translated_trajs"):
    pass
else:
    os.mkdir("translated_trajs")

if os.path.isdir("stepB_cluster_mod"):
    pass
else:
    os.mkdir("stepB_cluster_mod")

if os.path.isdir("stepA_wholesys_mod_noclust"):
    pass
else:
    os.mkdir("stepA_wholesys_mod_noclust")

if os.path.isdir("only_trans_cluster_trajs"):
    pass
else:
    os.mkdir("only_trans_cluster_trajs")

if os.path.isdir("index_files"):
    pass
else:
    print("The necessay index files do not exist. Cancelling script")
    exit(1)


# This is just a silly print statement to make the user see what parameters are chosen
# more verbosity is never a bad thing...
print()
print("Performing COM-based centering on files: ")
print(f"Trajectory: {TRAJ} | .tpr {TPR}")
print(f"Start time: {START_FRAME} {TIME_UNIT} | End time: {STOP_FRAME} {TIME_UNIT}")
print(f"Centering every {FRAME_ITER} frame")
print(f"Box dimensions: x = {BOXX} ; y = {BOXY} ; {BOXZ}")
print()

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
    # If there exists a unique index file for the frame that contains the max cluster,
    # then that means I will perform the COM based centering on the largest cluster
    # RECALL THIS SCRIPT MUST BE RUN AFTER THE COARSE DENSITY CENTERING WHICH MAKES THIS
    # INDEX FILE
    if os.path.isfile(f"./index_files/all_{FRAME}{TIME_UNIT}.ndx"):
        # Export the entire system with modifications, starting from the nojump_whole traj
        # this will be important for later
        # this step DOES use centering
        wholesys_export = (f"echo 0 | {FUNC} trjconv -f {TRAJ} -s {TPR} "
                    f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
                    f"-pbc mol -ur compact "
                    f"-n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
                    f"-o ./stepA_wholesys_mod/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc")
        if args.silent:
            subprocess.run(wholesys_export, shell=True, capture_output=True)
        else:
            subprocess.run(wholesys_export, shell=True)

        # export ONLY the largest cluster from the same frame, which will be important
        # for calculating the COM of the cluster
        # this step ALSO uses centering
        # removed the "-center -boxcenter tric" part of the trajectory conversion
        #   the xtc file is required as the trajectory for the md.load step
        cluster_export_xtc = (f"echo 10 10 | {FUNC} trjconv -f {TRAJ} -s {TPR} "
                    f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
                    f"-pbc mol -ur compact "
                    f"-n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
                    f"-o ./stepB_cluster_mod/frame_{FRAME}{TIME_UNIT}_cluster.xtc")
        if args.silent:
            subprocess.run(cluster_export_xtc, shell=True, capture_output=True)
        else:
            subprocess.run(cluster_export_xtc, shell=True)

        #   the gro file is required as the topology for the .xtc file in mdtraj
        cluster_export_gro = (f"echo 10 10 | {FUNC} trjconv -f {TRAJ} -s {TPR} "
                    f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
                    f"-pbc mol -ur compact "
                    f"-n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
                    f"-o ./stepB_cluster_mod/frame_{FRAME}{TIME_UNIT}_cluster.gro")
        if args.silent:
            subprocess.run(cluster_export_gro, shell=True, capture_output=True)
        else:
            subprocess.run(cluster_export_gro, shell=True)

        # calculate the COM of the cluster and determine the translational moves
        clust_TOP = f"./stepB_cluster_mod/frame_{FRAME}{TIME_UNIT}_cluster.gro"
        traj_clust_frame = md.load(f"./stepB_cluster_mod/frame_{FRAME}{TIME_UNIT}_cluster.xtc", 
                                top=clust_TOP)
        clust_only_COM = md.compute_center_of_mass(traj_clust_frame)

        transX = BOXX - clust_only_COM[0][0]
        transY = BOXY - clust_only_COM[0][1]
        transZ = BOXZ - clust_only_COM[0][2]

        print(f"Center of Mass of largest cluster for frame: {FRAME}{TIME_UNIT}")
        print(f"{clust_only_COM[0]}")
        print(f"Translation moves for frame: {FRAME}{TIME_UNIT}")
        print(f"X = {transX}")
        print(f"Y = {transY}")
        print(f"Z = {transZ}")

        
        # translate the modified whole system from above (see step A) using the translation
        # steps calculated directly above
        # AND perform PBC corrections on this
        translated_frame = (f"echo 0 | {FUNC} trjconv "
                f"-f ./stepA_wholesys_mod/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc " 
                f"-s {TPR} "
                f"-trans {transX} {transY} {transZ} "
                f"-pbc mol "
                f"-o ./translated_trajs/frame_{FRAME}{TIME_UNIT}_transl.xtc")
        if args.silent:
            subprocess.run(translated_frame, shell=True, capture_output=True)
        else:
            subprocess.run(translated_frame, shell=True)


        ## DEBUG
        # this is a debug step, to see if my procedure is even working
        # I am going to translate only the largest cluster to make a cluster only trajectory
        # to see if the PBC things I'm doing is even correct
        translated_cluster = (f"echo 10 | {FUNC} trjconv "
                f"-f ./stepA_wholesys_mod/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc " 
                f"-s {TPR} "
                f"-trans {transX} {transY} {transZ} "
                f"-n ./index_files/all_{FRAME}{TIME_UNIT}.ndx "
                f"-pbc mol "
                f"-o ./only_trans_cluster_trajs/frame_{FRAME}{TIME_UNIT}_transl_cluster.gro")
        if args.silent:
            subprocess.run(translated_cluster, shell=True, capture_output=True)
        else:
            subprocess.run(translated_cluster, shell=True)
        continue


    else:
        # Export the entire system but still do centering and what not to be consistent
        # with the other steps
        # and since theres separate cluster file, I will save out xtc and gro from this
        # step and calculate the COM on it
        # removed the "-center -boxcenter tric" part of the trajectory conversion
        wholesys_noclust_xtc = (f"echo 0 | {FUNC} trjconv -f {TRAJ} -s {TPR} "
                f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
                f"-pbc mol -ur compact "
                f"-o ./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc")
        subprocess.run(wholesys_noclust_xtc, shell=True)
        wholesys_noclust_gro = (f"echo 0 | {FUNC} trjconv -f {TRAJ} -s {TPR} "
                f"-b {FRAME} -e {FRAME} -tu {TIME_UNIT} "
                f"-pbc mol -ur compact "
                f"-o ./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.gro")
        subprocess.run(wholesys_noclust_gro, shell=True)


        # Step 3a - calculate the COM of the no cluster frame and determine the translational moves
        noclust_TOP = f"./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.gro"
        traj_no_clust = md.load(f"./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc",
                                top=noclust_TOP)
        no_clust_COM = md.compute_center_of_mass(traj_no_clust)


        print(f"Center of Mass of largest cluster for frame: {FRAME}{TIME_UNIT}")
        print(f"{clust_only_COM[0]}")
        print(f"Translation moves for frame: {FRAME}{TIME_UNIT}")
        print(f"X = {transX}")
        print(f"Y = {transY}")
        print(f"Z = {transZ}")

        transX = BOXX - no_clust_COM[0][0]
        transY = BOXY - no_clust_COM[0][1]
        transZ = BOXZ - no_clust_COM[0][2]


        # now, using the translational moves, translate the .xtc file generated above 
        # and perform pbc corrections on it!
        translated_noclust_frame = (f"echo 0 | {FUNC} trjconv "
            f"-f ./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc " 
            f"-s {TPR} "
            f"-trans {transX} {transY} {transZ} "
            f"-pbc mol "
            f"-o ./translated_trajs/frame_{FRAME}{TIME_UNIT}_transl.xtc")
        subprocess.run(translated_noclust_frame, shell=True, capture_output=True)

        """
        2025 02 26 Leave the following code commented out for now
        """
        # ## DEBUG
        # # this is a debug step, to see if my procedure is even working
        # # I am going to translate only the largest cluster to make a cluster only trajectory
        # # to see if the PBC things I'm doing is even correct
        # translated_noclust_frame_cluster = (f"echo 0 | {FUNC} trjconv "
        #     f"-f ./stepA_wholesys_mod_noclust/frame_{FRAME}{TIME_UNIT}_whole_sys.xtc " 
        #     f"-s {TPR} "
        #     f"-trans {transX} {transY} {transZ} "
        #     f"-pbc mol "
        #     f"-o ./only_trans_cluster_trajs/frame_{FRAME}{TIME_UNIT}_transl_cluster.gro")
        # subprocess.run(translated_noclust_frame_cluster, shell=True, capture_output=True)


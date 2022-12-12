"""
@Title:         dbscan_clustering.py
@Author:        Mando Ramirez
@Date:          20220307

@Description:   This code uses the DBSCAN clustering method (provided in scipy) to cluster a trajectory and determine
the medoid structure. This general code can handle all types of trajectories. It is modeled and adapted from the
cg_openmm DBSCAN clustering method implemented by Lenny and Chris of the MR Shirts lab.

@Updates:
20220328    Added a feature to control the frame stride, as needed.
"""

import os
import mdtraj as md
import dbscan_functions as dbscan
import argparse


## Set up the argument parser!
parser = argparse.ArgumentParser(description="DBSCAN clustering method to find medoid structures from the provided "
                                             "trajectory. PLEASE NOTE: this will find the medoid for the ENTIRE system"
                                             " in the trajectory. Not a subset. If you need a subset, you will need"
                                             " to write your own code to handle that.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-fs", help="frame stride i.e. No. every other frames to analyze, default=1",
                        default=1, type=int)
parser.add_argument("-od", help="output directory location. Please do not include any extraneous"
                                " dots or dashes", default="cluster_output")
parser.add_argument("-id", help="id for labeling purposes")


args = parser.parse_args()
trajectory= args.t
topology = args.p
frameStride = args.fs


# Load the trajectory
traj = md.load(trajectory, top=topology)


## Run the clustering!
medoids, cluster_size, cluster_rmsd, n_noise, silhouette_avg, labels = \
    dbscan.get_cluster_medoid_positions_DBSCAN(
        traj, min_samples=10, eps=0.50,
        frame_start=0, frame_stride=frameStride, frame_end=traj.n_frames,
        output_cluster_traj=True, filter=True,
        output_dir=args.od
    )

savefile = open(f"./{args.od}/cluster_results{args.id}.csv", 'w')
savefile.write(f"Cluster size,{cluster_size}\n")
savefile.write(f"Cluster RMSD,{cluster_rmsd}\n")
savefile.write(f"Noise members,{n_noise}\n")
savefile.write(f"Silhouette avgs,{silhouette_avg}\n")
savefile.write(f"Labels,{labels}\n")

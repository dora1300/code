"""
@Title:             dbscan_functions.py

@Authors:           Dominique Ramirez, Chris Walker, Lenny Fobe
@Date:              20220105

@Description:       This file contains the methods to perform the DBSCAN clustering method on a provided trajectory.
This file only contains the code necessary for executing the DBSCAN clustering, and not where the analysis actually
occurs.

@Acknowledgment:    Most of the code is adapted/copied from
(https://github.com/shirtsgroup/analyze_foldamers/blob/master/analyze_foldamers/ensembles/cluster.py)
which was written by Lenny Fobe and primarily by Chris Walker. I am thankful for both of them for making their code
available for me to adapt.
"""
import os
import mdtraj as md
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, silhouette_samples
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import minimize

"""
"""

def get_rmsd_matrix(trajectory, frame_start, frame_stride, frame_end):
    """
    Internal function for reading trajectory files and computing the RMSD.
    This code was adapted to handle trajectory objects that are not cgmodels, as well as only handling one single
    trajectory instead of several.
    :param trajectory:
    :param frame_start:
    :param frame_stride:
    :param frame_end:

    :return:
    """
    # load the trajectory with mdtraj
    # shorten the trajectory if need be, following the start end and stride variables
    traj_all = trajectory[frame_start:frame_end:frame_stride]

    # Align all the structures with the first frame as the reference
    for i in range(1,traj_all.n_frames):
        md.Trajectory.superpose(traj_all[i],traj_all[0])
        # This rewrites to traj_all

    # compute the pairwise RMSD matrix
    distances = np.empty((traj_all.n_frames, traj_all.n_frames))
    for i in range(traj_all.n_frames):
        distances[i] = md.rmsd(traj_all, traj_all, i)

    return distances, traj_all


def filter_distances(distances, filter_ratio=0.05):
    """
    Function for filtering out data points with few neighbors within a cutoff radius.
    This code was copied exactly from the analyze_foldamers code from Lenny and Chris

    :param distances: square matrix of pairwise RMSD distances
    :type distances: 2d numpy array

    :param filter_ratio: desired fraction of data remaining after neighborhood radius filtering
    :type filter_ratio: float

    :returns:
       - distances_filtered (2d numpy array) - distance matrix of data points satisfying filter parameters
       - neighbors_dense (1d numpy array) - indices of the original dataset which satisfy filter parameters
       - filter_ratio (float) - fraction of data remaining after filtering
    """

    filter_ratio_target = filter_ratio

    def get_filter_ratio(x0):
        # Function to minimize
        cutoff_radius = x0[0]
        density_cutoff = x0[1]

        neighbors = np.zeros((len(distances[:, 0])))

        for i in range(len(distances[:, 0])):
            neighbors[i] = (distances[:, i] <= cutoff_radius).sum() - 1
            # Excludes the self data point

        neighbors_dense = np.argwhere(neighbors >= density_cutoff)[:, 0]
        filter_ratio = len(neighbors_dense) / len(neighbors)

        return (filter_ratio - filter_ratio_target) ** 2

    # Optimize cutoff_radius, density_cutoff parameters to get desired filter ratio
    # A value of 0.05 is reasonable for rmsd distances, 75 is reasonable for torsion n-dimensional euclidean distances
    x0 = [np.mean(distances) / 2, 5]

    results = minimize(get_filter_ratio, x0, method='Nelder-Mead')

    cutoff_radius = results.x[0]
    density_cutoff = results.x[1]

    # Apply filtering parameters:

    neighbors = np.zeros((len(distances[:, 0])))

    for i in range(len(distances[:, 0])):
        neighbors[i] = (distances[:, i] <= cutoff_radius).sum() - 1
        # Excludes the self data point

    neighbors_dense = np.argwhere(neighbors >= density_cutoff)[:, 0]

    # Need to select 1 dimension at a time:
    distances_filtered = distances[:, neighbors_dense]
    distances_filtered = distances_filtered[neighbors_dense, :]

    filter_ratio_actual = len(neighbors_dense) / len(neighbors)

    print(f"filtered {(1 - filter_ratio_actual) * 100} % of data using:")
    print(f"cutoff radius = {cutoff_radius}")
    print(f"number neighbors cutoff: {density_cutoff}")

    return distances_filtered, neighbors_dense, filter_ratio_actual

def write_clusters_to_file(labels, traj_all, output_dir, output_format):
    """"""
    # Write ouput directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    clusters = np.unique(labels)

    for k in clusters:
        cluster_indices  = np.argwhere(labels == k)
        if k == -1:
            k = "noise"
        file_name = str(f"{output_dir}/clustermembers_{k}.{output_format}")
        cluster_traj = traj_all.slice(cluster_indices.reshape(-1))
        cluster_traj.save(file_name)

def make_silhouette_plot(
        cluster_fit, silhouette_sample_values, silhouette_avg,
        n_clusters, cluster_rmsd, cluster_sizes, plotfile):
    """Internal function for creating silhouette plot"""

    fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))

    y_lower = 10

    for k in range(n_clusters):
        kth_cluster_sil_values = silhouette_sample_values[cluster_fit.labels_ == k]

        kth_cluster_sil_values.sort()

        y_upper = y_lower + cluster_sizes[k]

        color = cm.nipy_spectral(float(k) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper), 0,
            kth_cluster_sil_values,
            facecolor=color, edgecolor=color, alpha=0.7)

        ax1.text(-0.05, y_lower + 0.5 * cluster_sizes[k], str(k))

        y_lower = y_upper + 10

    ax1.set_xlabel("Silhouette coefficient")
    ax1.set_ylabel("Cluster label")
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear y ticks/labels
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    xlim = ax1.get_xlim()
    ylim = ax1.get_ylim()

    for k in range(n_clusters):
        ax1.text(xlim[0] + 0.75 * (xlim[1] - xlim[0]), (0.9 - (k / 25)) * (ylim[1] - ylim[0]),
                 f"Cluster {k} RMSD: {cluster_rmsd[k]:.3f}",
                 {'fontsize': 8},
                 horizontalalignment='left')

    plt.tight_layout()
    plt.savefig(plotfile)

    return

def get_cluster_medoid_positions_DBSCAN(
    trajectory, min_samples=5, eps=0.5,
    frame_start=0, frame_stride=1, frame_end=-1, output_format="pdb",
    output_dir="cluster_output", output_cluster_traj = False, plot_silhouette=True,
    plot_rmsd_hist=True, filter=True, filter_ratio=0.25,
    core_points_only=True
):
    # check to see if output directory exists, and make if not
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # perform the pairwise RMSD matrix calculation
    distances, traj_all = get_rmsd_matrix(trajectory, frame_start, frame_stride, frame_end)

    # filter the distances to remove objects that are not close to each other
    if filter:
        distances, dense_indices, filter_ratio_actual = filter_distances(distances, filter_ratio)
        traj_all = traj_all[dense_indices]

    if plot_rmsd_hist:
        # plot the rmsd histogram
        distances_row = np.reshape(distances, (distances.shape[0] * distances.shape[1], 1))

        # Remove the diagonal 0 elements:
        distances_row = distances_row[distances_row != 0]

        figure = plt.figure()
        n_out, bin_edges_out, patch = plt.hist(
            distances_row, bins=1000, density=True)
        plt.xlabel('rmsd')
        plt.ylabel('probability density')
        plt.savefig(f'{output_dir}/distances_rmsd_hist.png')
        plt.close()

    # Cluster with sklearn DBSCAN
    dbscan = DBSCAN(min_samples=min_samples, eps=eps, metric='precomputed').fit(distances)
    # The produces a cluster labels from 0 to n_clusters-1, and assigns -1 to noise points
    # the precomputed metric means that I've precomputed sparse neighborhoods (thanks to the filtering)
    # which means DBSCAN doesn't have to do it and speeds up the computation

    # Get labels
    labels = dbscan.labels_

    # Get core sample indices:
    core_sample_indices = dbscan.core_sample_indices_

    # Number of clusters:
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    # Number of noise points:
    n_noise = list(labels).count(-1)

    # Get indices of frames in each cluster:
    cluster_indices = {}
    cluster_indices_core = {}
    cluster_sizes = []
    cluster_sizes_core = []

    for k in range(n_clusters):
        cluster_indices[k] = np.argwhere(labels == k)[:, 0]
        cluster_indices_core[k] = []
        for elem in cluster_indices[k]:
            if elem in core_sample_indices:
                cluster_indices_core[k].append(elem)
        cluster_sizes.append(len(cluster_indices[k]))
        cluster_sizes_core.append(len(cluster_indices_core[k]))

    # Get indices of frames classified as noise:
    noise_indices = np.argwhere(labels == -1)[:, 0]

    # Find the structure closest to each center (medoid):
    # OPTICS/DBSCAN does not have a built-in function to transform to cluster-distance space,
    # as the centroids of the clusters are not physically meaningful in general. However, as
    # RMSD between structures is our only clustering feature, the cluster centers (regions of
    # high density) will likely be representative structures of each cluster.

    # Following the protocol outlined in MDTraj example:
    # http://mdtraj.org/1.9.3/examples/centroids.html

    # Create distance matrices within each cluster:
    distances_k = {}

    if core_points_only:
        # this, I believe, finds the medoids considering only core samples from the clusters
        for k in range(n_clusters):
            distances_k[k] = np.zeros((cluster_sizes_core[k], cluster_sizes_core[k]))
            for i in range(cluster_sizes_core[k]):
                for j in range(cluster_sizes_core[k]):
                    distances_k[k][i, j] = distances[cluster_indices_core[k][i], cluster_indices_core[k][j]]

        # Compute medoid based on similarity scores:
        medoid_index = []  # Global index
        intra_cluster_medoid_index = []  # Index within cluster
        for k in range(n_clusters):
            intra_cluster_medoid_index.append(
                np.exp(-distances_k[k] / distances_k[k].std()).sum(axis=1).argmax()
            )
            # Here we need to use the global sample index to find the medoid structure:
            medoid_index.append(cluster_indices_core[k][intra_cluster_medoid_index[k]])

    else:
        # this part, I believe, finds the medoids but also takes into consideration the location of noise samples
        # in addition to the core samples.
        for k in range(n_clusters):
            distances_k[k] = np.zeros((cluster_sizes[k], cluster_sizes[k]))
            for i in range(cluster_sizes[k]):
                for j in range(cluster_sizes[k]):
                    distances_k[k][i, j] = distances[cluster_indices[k][i], cluster_indices[k][j]]

        # Compute medoid based on similarity scores:
        medoid_index = []  # Global index
        intra_cluster_medoid_index = []  # Index within cluster
        for k in range(n_clusters):
            intra_cluster_medoid_index.append(
                np.exp(-distances_k[k] / distances_k[k].std()).sum(axis=1).argmax()
            )
            # Here we need to use the global sample index to find the medoid structure:
            medoid_index.append(cluster_indices[k][intra_cluster_medoid_index[k]])

    medoid_xyz = np.zeros([n_clusters, traj_all.n_atoms, 3])
    for k in range(n_clusters):
        medoid_xyz[k, :, :] = traj_all[medoid_index[k]].xyz[0]


    # write the resulting medoids as their own PDB files
    number_medoids = medoid_xyz.shape[0]
    for k in range(number_medoids):
        positions = medoid_xyz[k]
        file_name = str(f"{output_dir}/medoid_{k}.{output_format}")

        # temp_traj = md.load(topology=traj.topology)
        # top = temp_traj.topology

        medoid_traj = md.Trajectory(xyz=positions, topology=trajectory.topology)
        medoid_traj.save_pdb(file_name)

    # 2022 01 06 --  the save mechanism worked! Yay!
    # the original code explicitly specified units of nanometer from the simtk package (provided by openmm code)
    # but I do not have that installed right now. Since the output should be in nanometers since that's the base
    # unit that mdtraj operates in, I will omit the explicit nanometer call for the time being

    if output_cluster_traj:
        write_clusters_to_file(labels, traj_all, output_dir, output_format)

    # Compute intra-cluster rmsd of samples to medoid based on structure rmsd
    cluster_rmsd = np.zeros(n_clusters)

    for k in range(n_clusters):
        cluster_rmsd[k] = np.sqrt(((distances_k[k][intra_cluster_medoid_index[k]]**2).sum())/len(cluster_indices[k]))

    # Get silhouette scores
    try:
        silhouette_sample_values = silhouette_samples(distances, labels)
        silhouette_avg = np.mean(silhouette_sample_values[labels != -1])

        if plot_silhouette:
            # Plot silhouette analysis
            plotfile = f"{output_dir}/silhouette_dbscan_min_sample_{min_samples}_eps_{eps}.png"

            make_silhouette_plot(
                dbscan, silhouette_sample_values, silhouette_avg,
                n_clusters, cluster_rmsd, cluster_sizes, plotfile
            )
    except ValueError:
        print(
            "There are either no clusters, or no noise points identified. Try adjusting DBSCAN min_samples, eps parameters.")
        silhouette_avg = None


    return medoid_xyz, cluster_sizes, cluster_rmsd, n_noise, silhouette_avg, labels


if __name__ == "__main__":
    t = "../sim_1/production_125/traj_comp.xtc"
    top = "../sim_1/production_125/coil_md_125.pdb"

    traj = md.load(t, top=top)

    meds = get_cluster_medoid_positions_DBSCAN(
        traj, min_samples=10, eps=0.125,
        frame_start=0, frame_stride=1, frame_end=traj.n_frames, output_cluster_traj=True,
        filter=True
    )
    print(meds)

"""
@Title:             rmsd_dimer_analysis.py
@Name:              Mando A Ramirez
@Date:              2022 07 14

@Description:       This script calculates the RMSD specifically for coil dimer simulations because this handles
parallel and anti-parallel symmetry in my coil models.
    This code builds off of rmsd_analysis.py and extensively uses the workflow provided by Chris Walker in the
analyze_foldamers github repo: https://github.com/shirtsgroup/analyze_foldamers/blob/master/analyze_foldamers/ensembles/cluster.py#L833
Thank you Chris.

THIS ONLY CALCULATES THE 1D RMSD, NOT PAIRWISE.

@Updates:
"""
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description="RMSD analysis of coil dimers, specifically to handle parallel/"
                                             "antiparallel invariance.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="topology file i.e. PDB file with extension", required=True)
parser.add_argument("-reft", help="reference structure to do RMSD against. NEEDS to be a single structure file,"
                                  " not another trajectory.", required=True)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file NO extension", default="output_rmsd_symmetry")
parser.add_argument("-plot_timeseries", help="Invoke to turn on time series specific RMSD data. Only "
                                             "meaningful for 1D RMSDs", action="store_true")
parser.add_argument("-plot_distribution", help="Invoke to plot distributions of RMSD values for the trajectory. "
                                             "Only meaningful for 1D RMSDs", action="store_true")

args = parser.parse_args()

trajectory= args.t
Topology = args.p
reference = args.reft
sim_code = args.id
output = args.o + ".csv"

REF = md.load(reference)        # reference structure is the provided reference in this case
traj = md.load(trajectory, top=Topology)    # sample trajectory to compute RMSD for


# Step 1 -- sort the trajectory and separate parallel coils from antiparallel coils. This is important to prevent
# double counting in later steps.
# Also turn off periodic images so that the correct distances are calculated and periodic images are not included
parallel_indices = []
antiparallel_indices = []
for i in range(traj.n_frames):
    if md.compute_distances(traj[i], [[0, 32]], periodic=False) < md.compute_distances(traj[i], [[0, 63]], periodic=False):
        parallel_indices.append(i)
    elif md.compute_distances(traj[i], [[0, 32]], periodic=False) > md.compute_distances(traj[i], [[0, 63]], periodic=False):
        antiparallel_indices.append(i)
    else:
        parallel_indices.append(i)

# Now make the two separate trajectories for parallel and antiparallel coils
parallel_traj = traj[parallel_indices[0]]
for i in parallel_indices[1:]:
    parallel_traj = parallel_traj.join(traj[i])

if len(antiparallel_indices) != 0:
    antiparallel_traj = traj[antiparallel_indices[0]]
    for i in antiparallel_indices[1:]:
        antiparallel_traj = antiparallel_traj.join(traj[i])
else:
    pass


# Superpose the parallel trajectory and do the regular parallel RMSD analysis
for i in range(parallel_traj.n_frames):
    md.Trajectory.superpose(parallel_traj[i], REF)
parallel_rmsd = md.rmsd(parallel_traj, REF)


if len(antiparallel_indices) != 0:
    # Now get ready to do the antiparallel analysis by making the inverted indices trajectory first
    xyz_rev = np.empty((antiparallel_traj.n_frames, antiparallel_traj.n_atoms, 3))

    for i in range(antiparallel_traj.n_frames):
        j_f = -1
        for j in range(0, 32):
            xyz_rev[i][j][:] = antiparallel_traj.xyz[i][j][:]
        for j in range(32, 64):
            xyz_rev[i][j][:] = antiparallel_traj.xyz[i][j_f][:]
            j_f -= 1

    # Now make the reversed indices trajectory object, superpose against the reference structure, and then do the analysis
    rev_ap_traj = md.Trajectory(xyz=xyz_rev, topology=antiparallel_traj.topology)
    for i in range(rev_ap_traj.n_frames):
        md.Trajectory.superpose(rev_ap_traj[i], REF)
    antiparallel_rmsd = md.rmsd(rev_ap_traj, REF)
else:
    pass


# Now combine the RMSD analyses!
if len(antiparallel_indices) != 0:
    combined = np.concatenate((parallel_rmsd, antiparallel_rmsd))
else:
    combined = np.copy(parallel_rmsd)


# Now do the plotting
if args.plot_timeseries:
    fig, ax = plt.subplots()
    ax.plot(traj.time, combined, 'b-', linewidth=0.35, label="1D RMSD")
    plt.xlabel("Time (ps)")
    plt.ylabel("Distance (nm)")
    plt.ylim(0, np.max(combined)+0.1)
    plt.title(f"RMSD !Symmetry! time series of {sim_code} trajectory in reference to \n{reference}")
    plt.savefig(f"{sim_code}_rmsd_symmetry_timeseries.png", dpi=300)
    plt.close("all")
if args.plot_distribution:
    fig, ax = plt.subplots()
    ax.hist(combined, density=True, bins="sqrt", color="blue", histtype="step")
    plt.xlabel("Distance (nm)")
    plt.ylabel("Frequency (counts)")
    #plt.xlim(0, np.max(combined)+0.5)
    plt.xlim(0, 2.75)
    plt.title(f"RMSD !Symmetry! distribution of {sim_code} trajectory \nin reference to {reference}")
    plt.savefig(f"{sim_code}_rmsd_symmetry_distribution.png", dpi=300)
    plt.close("all")

# Now its time to save the RMSD data to a file for later use, if required. As of 20211103 I am not saving out the 2D
#   RMSD since that will be a large file and will be kinda crazy to output. And I'm not sure I'd need that data for
#   any additional analysis.
with open(output, 'w') as outfile:
    outfile.write("Frame no., Time (ps), RMSD (nm)\n")
    for i in range(len(combined)):
        outfile.write(f"{i},{traj.time[i]},{combined[i]:.4f}\n")

avg_rmsd = np.average(combined)
std_rmsd = np.std(combined)
with open(args.o + "_short.csv", 'w') as shortout:
    shortout.write("Avg RMSD,Stddev RMSD\n")
    shortout.write(f"{avg_rmsd:.3f},{std_rmsd:.3f}\n")


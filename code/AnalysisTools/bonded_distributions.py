"""
@Title:             bonded_distributions.py
@Name:              Mando A Ramirez
@Date:              2022 05 03

@Description:       This script calculates the distributions of bonds, angles, and torsions for any system provided. It
makes no distinction between linkers, coils, or anything else. As of right now, it treats everything globally and will
plot everything on a single plot -- so be careful with its application.

@Updates:
20220509    - Changed the binning size of the bond angle histogram

"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description="Bonded Distributions Plotting -- plots the distributions of bond lengths,"
                                             " angles, and torsions of C-alpha coarse-grained systems. This analyzes "
                                             "*everything* in aggregate, so be careful when you interpret the "
                                             "distributions!!\n"
                                             "Also please note that currently all the models in the simulation need to"
                                             " be the same length.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual models in the simulation", required=True, type=int)
parser.add_argument("-l", help="number of beads in each individual model", required=True, type=int)
parser.add_argument("-f", help="stride of trajectory to analyze", default=1, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)

# print warning message to remove responsibility for poor user interpretation
print()
print("Just in case you didn't get the memo, this plots *AGGREGATE* distributions. Be careful in your interpretation "
      "and analysis and make sure you take this into account.")
print()

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
lmodel = args.l
nframes = args.f
sim_code = args.id

traj = md.load(trajectory, top=topology)



# Do the analysis for bond lengths
pairs = []
for n in range(1, nmodels + 1):
    for j in range(1, lmodel):
        i = j + (lmodel * (n - 1))
        pairs.append([i-1, i])
"""
The structure of these distance arrays will be:
[ [frame1: pair1 pair2 pair3 ...],
  [frame2: pair1 pair2 pair3 ...],
  [frame3: pair1 pair2 pair3 ...],
  ...
]
and so this is a 2D array
"""
pairs_dist = md.compute_distances(traj[::nframes], pairs)
pairs_avg = np.average(np.reshape(pairs_dist, -1))



# Do the analysis for Bond angles
total_bondangles = np.array([])

for n in range(1, nmodels + 1):
    for j in range(1, (lmodel - 1)):
        i = j + (lmodel * (n - 1))

        ang_coords = [i - 1, i, i + 1]
        # these are the angles for the given coordinates for the ENTIRE trajectory
        angle_deg = md.compute_angles(traj[::nframes], [ang_coords, ang_coords])[:, 0] * (180 / np.pi)

        total_bondangles = np.concatenate((total_bondangles, angle_deg), axis=0)
# The final array, total_bondangles, will be massive! but that's okay
average_bondangles = np.average(total_bondangles)



# Begin the analysis for torsions
total_torsions = np.array([])

for n in range(1, nmodels + 1):
    for j in range(1, (lmodel - 2)):
        i = j + (lmodel * (n - 1))

        dih_coords = [i - 1, i, i + 1, i + 2]
        # these are the angles for the given coordinates for the ENTIRE trajectory
        dihedral_deg = md.compute_dihedrals(traj[::nframes], [dih_coords, dih_coords])[:, 0] * (180 / np.pi)

        total_torsions = np.concatenate((total_torsions, dihedral_deg), axis=0)
# The final array, total_torsions, will be massive! but that's okay
average_torsions = np.average(total_torsions)



# Now move onto the plotting
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(15, 20))

ax1.hist(np.reshape(pairs_dist, -1), density=True, bins=150, color="tab:orange", alpha=0.5)
ax1.axvline(x=pairs_avg, linestyle="--", color="black", label=f"Avg: {pairs_avg:.4f}")
ax1.set_xlim(0.30, 0.45)
ax1.set_xlabel("Bond length distance (nm)")
ax1.set_ylabel("Frequency (normalized counts)")
ax1.set_title("Bond length distribution")
ax1.legend()

ax2.hist(total_bondangles, density=True, bins=100, color="royalblue", alpha=0.5)
ax2.axvline(x=average_bondangles, color="black", linestyle="--", label=f"Avg: {average_bondangles:.2f}")
ax2.set_xlim(-180, 180)
ax2.set_xlabel(r"Pseudo-angle $\theta$ (deg)")
ax2.set_ylabel("Frequency (normalized counts)")
ax2.set_title("Bond angle distribution")
ax2.legend()

ax3.hist(total_torsions, density=True, bins=100, color="darkorchid", alpha=0.5)
ax3.axvline(x=average_torsions, color="black", linestyle="--", label=f"Avg: {average_torsions:.2f}")
ax3.set_xlim(-180, 180)
ax3.set_xlabel(r"Pseudo-torsion $\alpha$ (deg)")
ax3.set_ylabel("Frequency (normalized counts)")
ax3.set_title("Torsion angle distribution")
ax3.legend()

plt.savefig(f"{sim_code}_distributions.png", dpi=600)

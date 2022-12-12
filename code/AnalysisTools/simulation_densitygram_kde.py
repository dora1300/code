"""
@Name:          Dominique Ramirez
@Date:          2021 11 29

@Title:         simulation_densitygram_kde.py

@Description:   This script is an extension and application from DAR3-7a Step 5.
This has been modified for use in DAR3-8. It reads the pickled Tozzini data from DAR3-7 and
then plots them along with the alpha/theta angles from a provided simulation.
It will also plot densitygrams and KDE plots for the distributions.

@Updates:
20211215 - updated to be moved from DAR3-7 into its own analysis script within the primary Code location for my work.
"""
import matplotlib.pyplot as plt
import numpy as np
import math as m
import pickle
import mdtraj as md
import pandas as pd
import seaborn as sns
import argparse

"""
Argument Parser set up to take in trajectory and topology information.
"""
parser = argparse.ArgumentParser(
    description="Tool to plot density-grams and KDE plots for provided simulation trajectory. Also plots the reference Tozzini landscape.")

parser.add_argument("-t", help="Trajectory file (.xtc) plus path and extension")
parser.add_argument("-p", help="Topology file (.pdb) plus path and extension")
parser.add_argument("-o", help="Output DIRECTORY where to save the goods")
parser.add_argument("-n", help="Number of coil models in the simulation trajectory", type=int)
parser.add_argument("-l", help="Length of the coil model in the simulation. All coil models must be the same length",
                    type=int)
parser.add_argument("-avg",
                    help="Switch to calculate the mean and standard deviation of the Tozzini and simulated angles distribution.",
                    type=bool, default=False)
parser.add_argument("-plot", help="Switch to plot the scatter plot of Tozzini and simulated angles.", type=bool,
                    default=False)
parser.add_argument("-kde", help="Switch to plot the KDE plots of Tozzini and simulated angles.", type=bool,
                    default=False)

args = parser.parse_args()

"""
Open the pickled lists that were generated in DAR3-7a. These are considered static in the 
~/Code/AnalysisTools/Tozzini_reference_pickles/ folder. If anything changes in DAR3-7a, these referenced pickle files 
must be updated accordingly. DAR3-7a is considered the master file for these data.
"""
with open("./Tozzini_reference_pickles/tozzini_converted_angles/total_helix_alpha.pkl", 'rb') as f:
    helix_alpha = pickle.load(f)
with open("./Tozzini_reference_pickles/tozzini_converted_angles/total_helix_theta.pkl", 'rb') as f:
    helix_theta = pickle.load(f)
with open("./Tozzini_reference_pickles/tozzini_converted_angles/total_nonhelix_alpha.pkl", 'rb') as f:
    nonhelix_alpha = pickle.load(f)
with open("./Tozzini_reference_pickles/tozzini_converted_angles/total_nonhelix_theta.pkl", 'rb') as f:
    nonhelix_theta = pickle.load(f)

"""
Open the provided trajectory and analyze the angles within.
"""
trajectory = args.t
topology = args.p
output_dir = args.o
nmodels = args.n
lhelix = args.l
traj = md.load(trajectory, top=topology)
total_sim_alphas = []
total_sim_thetas = []

# !!!! This part of the code analyzes my data!
# This code was taken from coil_rama_analysis.py and repurposed for my use. See the original script for a description
# of how the code works.
for n in range(1, nmodels + 1):
    for j in range(1, (lhelix - 2)):
        i = j + (lhelix * (n - 1))
        dih_coords = [i - 1, i, i + 1, i + 2]
        ang_coords = [i - 1, i, i + 1]
        # Compute the dihedrals for the given atom coordinates, grab the 0 column (0th or 1st, doesn't matter) and convert
        #    to degrees. Do the same for angles. These both produce numpy arrays
        dihedral_deg = md.compute_dihedrals(traj, [dih_coords, dih_coords])[:, 0] * (180 / m.pi)
        angle_deg = md.compute_angles(traj, [ang_coords, ang_coords])[:, 0] * (180 / m.pi)
        for k in range(len(dihedral_deg)):
            total_sim_alphas.append(dihedral_deg[k])
            total_sim_thetas.append(angle_deg[k])
    # there is a potential that these lists might become very, very large in memory. Hopefully truncating the simulation
    # will help mitigate potential memory issues
print(f"Total simulation alpha angles: {len(total_sim_alphas)}")
print(f"Total simulation theta angles: {len(total_sim_thetas)}")

"""
Section to control the calculating of distribution averages and printing them. Controlled with args.avg
"""
if args.avg:
    print("Variances and comparisons of reference helix distribution to simulated helix distribution.")
    avg_helix_ref_alpha = np.average(np.array(helix_alpha))
    std_helix_ref_alpha = np.std(np.array(helix_alpha))
    avg_helix_ref_theta = np.average(np.array(helix_theta))
    std_helix_ref_theta = np.std(np.array(helix_theta))

    avg_sim_alpha = np.average(np.array(total_sim_alphas))
    std_sim_alpha = np.std(np.array(total_sim_alphas))
    avg_sim_theta = np.average(np.array(total_sim_thetas))
    std_sim_theta = np.std(np.array(total_sim_thetas))

    print(f"Reference alpha: {avg_helix_ref_alpha:.2f} +/- {std_helix_ref_alpha:.2f}")
    print(f"Simulation alpha: {avg_sim_alpha:.2f} +/- {std_sim_alpha:.2f}")
    print(f"Reference theta: {avg_helix_ref_theta:.2f} +/- {std_helix_ref_theta:.2f}")
    print(f"Simulation theta: {avg_sim_theta:.2f} +/- {std_sim_theta:.2f}")
else:
    pass

# ### This section is designed to create the reference Ramachandran conversion map so that I can compare my data to what
# # is allowable
# switch = False
# if switch:
#     phi = np.arange(-180, 181, 0.5)
#     psi = np.arange(-180, 181, 0.5)
#     full_list = []
#
#     for i in range(len(phi)):
#         for k in range(len(psi)):
#             full_list.append([phi[i], psi[k]])
#
#     converted_alphas = []
#     converted_thetas = []
#
#     for i in range(len(full_list)):
#         phi_r = full_list[i][0] * (m.pi / 180)
#         psi_r = full_list[i][1] * (m.pi / 180)
#         t_theta = tc.theta(tc.t_r, tc.y1_r, tc.y2_r, phi_r, psi_r)
#         t_alpha = tc.alpha_complex(tc.t_r, tc.y1_r, tc.y2_r, phi_r, psi_r)
#         converted_alphas.append(t_alpha) ; converted_thetas.append(t_theta)
# else:
#     pass


"""
This is for plotting the scatter plot of Tozzini angles and simulation angles. This is controlled by switch args.plot
"""
if args.plot:
    # Plot the angles in Tozzini space, and then save!
    fig, ax = plt.subplots(figsize=(8, 8))
    # plot the nonhelix pairs first, those will be on bottom, then plot the simulated angles to cover the nonhelix, then
    # plot the helix angles to show how much of my simulated helix angles are covered by the Tozzini Top8000 analysis

    ax.plot(np.array(nonhelix_alpha), np.array(nonhelix_theta), markeredgecolor='blue', markerfacecolor='blue',
            linestyle='', marker='x', markersize='0.10', label='Non-helix pseudoangles (Top8000 dataset)',
            fillstyle='full')  # Non-helix Tozzini
    ax.plot(np.array(total_sim_alphas), np.array(total_sim_thetas), markeredgecolor='limegreen',
            markerfacecolor="limegreen",
            linestyle='', marker='x', markersize='0.10', label='Simulation Trajectory angles',
            fillstyle='full')  # Helix Tozzini
    ax.plot(np.array(helix_alpha), np.array(helix_theta), markeredgecolor='red', markerfacecolor='red',
            linestyle='', marker='x', markersize='0.10', label='Helix pseudoangles (Top8000 dataset)',
            fillstyle='full')  # Simulated angles

    plt.xlim(-180, 180);
    plt.ylim(40, 170)
    plt.xlabel(r"$\alpha$ (deg)");
    plt.ylabel(r"$\theta$ (deg)")
    ax.legend(markerscale=20.0, fontsize="small")
    plt.title("Plot of Tozzini angles from Top8000 dataset and angles from coil simulation")
    # plt.show()
    plt.savefig(f"{output_dir}/simulation_angles_scatterplot.png", dpi=600)
    plt.close("all")
else:
    pass

"""
This is for plotting the kernel density maps for the Tozzini angles and simulated angles. Controlled by args.kde
"""
if args.kde:
    fig, ax = plt.subplots()
    helix_df = pd.DataFrame(list(zip(helix_alpha, helix_theta)))
    nonhelix_df = pd.DataFrame(list(zip(nonhelix_alpha, nonhelix_theta)))
    simulation_df = pd.DataFrame(list(zip(total_sim_alphas, total_sim_thetas)))
    helix_df.columns = ['alpha', 'theta']
    nonhelix_df.columns = ['alpha', 'theta']
    simulation_df.columns = ['alpha', 'theta']
    sns.set_style("white")

    sns.kdeplot(x=nonhelix_df.alpha, y=nonhelix_df.theta, cmap='Blues', bw_adjust=0.35,
                label='Non-helix pseudoangles (Top8000 dataset)')
    sns.kdeplot(x=simulation_df.alpha, y=helix_df.theta, bw_adjust=0.35, cmap='Greens',
                label='Simulation trajectory angles')
    sns.kdeplot(x=helix_df.alpha, y=helix_df.theta, bw_adjust=0.35, cmap='Reds',
                label='Helix pseudoangles (Top8000 dataset)')

    plt.xlim(-180, 180);
    plt.ylim(40, 170)
    plt.legend()
    ax.set_xlabel(r"$\alpha$ (deg)");
    ax.set_ylabel(r"$\theta$ (deg)")
    plt.savefig(f"{output_dir}/simulation_angles_KDE.png", dpi=600)
    plt.close("all")
else:
    pass


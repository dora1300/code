"""
@Title:         KDE_plotter.py
@Author:        Mando Ramirez
@Date:          20211215

@Description:   This code implements (1) a KDE generator to make probability distribution function estimates for a
given set of angles, (2) makes a density plot for the generated KDE, and (3) calculates the Kullbeck-Leibler
divergence for the provided simulated angles in comparison to the Tozzini reference angles.

To speed up some of the calculations, I pre-calculated the KDE and probabilites for the helix Tozzini angles, so I
will load those pickles into the file for analysis so I don't have to analyze those every single time.

@Updates:
20211217 - added the Bhattacharya distance calculation, and fixed the KL divergence to take in the Tozzini angles as
the set X to calculate the divergence on.

20211220 - added a switch for plotting but set the default to True. Also added the Jensen-Shannon distance to my metric
calculations and modified the KL divergence to do an additional calculation.

20220112 - changed the code to ONLY make a plot of the densities of each distribution.
This is so I can make graphs without having to take all the time to do divergence calculations.

20220204 - changed the name, used to be titled "kl_divergence_plot_ONLY.py" but now it is more descriptive.
"""
import math as m
import numpy as np
import scipy.stats as st
import scipy.special as sp
import scipy.spatial.distance as spd
import pickle
import argparse
import mdtraj as md
import matplotlib.pyplot as plt

def bc_metric(p_probs, q_probs):
    """

    :param p_probs:
    :param q_probs:
    :return:
    """
    BC = 0
    for i in range(len(p_probs)):
        BC += np.sqrt(p_probs[i] * q_probs[i])
    BC_dist = -1 * m.log(BC)
    Hellinger_dist = m.sqrt((1 - BC))
    return BC, BC_dist, Hellinger_dist

def kl_divergence(p_probs, q_probs):
    """

    :param p_probs:
    :param q_probs:
    :return:
    """
    kl = st.entropy(pk=p_probs, qk=q_probs)
    kl_convex = sp.kl_div(p_probs, q_probs, out=None)
    return kl, np.sum(kl_convex)


"""
Set up of argument parser, and parsing of arguments
"""
parser = argparse.ArgumentParser(description="Kullbeck-Leibler Divergence analysis. This script determines the KDE"
                                             " for a distribution of angles from both the Tozzini reference set and"
                                             " a provided simulation and compares the distributions to see if they"
                                             " are similar.")
parser.add_argument("-t", help="trajectory file with extension (.xtc/.pdb)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-f", help="analyze every f'th frame in the trajectory", default=10, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_KLdivergence")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="Boolean switch to control plotting. Default is True.", action="store_true")

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially

"""
Open the pickled lists that were generated in DAR3-7a. These are considered static in the 
~/Code/AnalysisTools/Tozzini_reference_pickles/ folder. If anything changes in DAR3-7a, these referenced pickle files 
must be updated accordingly. DAR3-7a is considered the master file for these data.
"""
maclocation = "/Users/mramirez/Code/AnalysisTools/Tozzini_reference_pickles"
linuxlocation = "/home/mando/Code/AnalysisTools/Tozzini_reference_pickles"
with open(f"{linuxlocation}/total_helix_alpha.pkl", 'rb') as f:
    helix_alpha = pickle.load(f)
with open(f"{linuxlocation}/total_helix_theta.pkl", 'rb') as f:
    helix_theta = pickle.load(f)
with open(f"{linuxlocation}/total_nonhelix_alpha.pkl", 'rb') as f:
    nonhelix_alpha = pickle.load(f)
with open(f"{linuxlocation}/total_nonhelix_theta.pkl", 'rb') as f:
    nonhelix_theta = pickle.load(f)
# with open(f"{location}/reference_helix_kernel.pkl", 'rb') as f:
#     refkernel = pickle.load(f)
# with open(f"{location}/reference_helix_probabilities.pkl", 'rb') as f:
#     REF = pickle.load(f)


"""
Load the trajectory and analyze the angles within
"""
traj = md.load(trajectory, top=topology)
total_sim_alphas = []
total_sim_thetas = []

# !!!! This part of the code analyzes my data!
# This code was taken from simulation_densitygram_kde.py and repurposed for my use.
# See the original script for a description of how the code works.
for n in range(1, nmodels + 1):
    for j in range(1, (lhelix - 2)):
        i = j + (lhelix * (n - 1))
        dih_coords = [i - 1, i, i + 1, i + 2]
        ang_coords = [i - 1, i, i + 1]
        # Compute the dihedrals for the given atom coordinates, grab the 0 column (0th or 1st, doesn't matter) and convert
        #    to degrees. Do the same for angles. These both produce numpy arrays
        # This part is also set up to slice the trajectory based on the argument args.f, which means analyze every
        #    f'th frame. The slicing notations means exactly that.
        dihedral_deg = md.compute_dihedrals(traj[::args.f], [dih_coords, dih_coords])[:, 0] * (180 / m.pi)
        angle_deg = md.compute_angles(traj[::args.f], [ang_coords, ang_coords])[:, 0] * (180 / m.pi)
        for k in range(len(dihedral_deg)):
            total_sim_alphas.append(dihedral_deg[k])
            total_sim_thetas.append(angle_deg[k])
    # there is a potential that these lists might become very, very large in memory. Hopefully truncating the simulation
    # will help mitigate potential memory issues



"""
Generate the kernel density estimates for the Tozzini reference angles and the simulated angles
"""
# Make the simulated angles and reference angles into a (2, N) sized array
# EVERYTHING will be in the order of (alpha, theta) == (x, y)
ref_values = np.vstack([np.array(helix_alpha), np.array(helix_theta)])
sim_values = np.vstack([np.array(total_sim_alphas), np.array(total_sim_thetas)])

# Estimate the kernel using the scipy.stats.gaussian_kde() class described here:
# <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html>
# I am using the default "Scott" method for bin width
refkernel = st.gaussian_kde(ref_values)
simkernel = st.gaussian_kde(sim_values)



"""
Generate a plot of the resulting KDEs as density grams for comparison to each other. This behavior is now controlled by 
the user but it would be wise to include it.
"""
if args.plot:
    amin = np.array(total_sim_alphas).min() ; amax = np.array(total_sim_alphas).max()
    tmin = np.array(total_sim_thetas).min() ; tmax = np.array(total_sim_thetas).max()

    refamin = np.array(helix_alpha).min() ; refamax = np.array(helix_alpha).max()
    reftmin = np.array(helix_theta).min() ; reftmax = np.array(helix_theta).max()

    simA, simT = np.mgrid[amin:amax:100j, tmin:tmax:100j]
    refA, refT = np.mgrid[refamin:refamax:100j, reftmin:reftmax:100j]

    simPositions = np.vstack([simA.ravel(), simT.ravel()])
    refPositions = np.vstack([refA.ravel(), refT.ravel()])

    SIM = np.reshape(simkernel.pdf(simPositions).T, simA.shape)
    REF = np.reshape(refkernel.pdf(refPositions).T, refA.shape)

    fig, ax = plt.subplots()

    ax.imshow(np.rot90(SIM), cmap="Blues", extent=[amin, amax, tmin, tmax], alpha=0.50)
    cset1 = ax.contour(simA, simT, SIM, colors='blue', linewidths=0.75)

    ax.imshow(np.rot90(REF), cmap="Reds", extent=[refamin, refamax, reftmin, reftmax], alpha=0.75)
    # it is important to use the CORRECT mins and maxes or else the actual image will be off from where it should be
    cset2 = ax.contour(refA, refT, REF, colors='black', linewidths=0.75)

    ax.clabel(cset1, inline=False, fontsize=2)
    ax.clabel(cset2, inline=False, fontsize=2)
    ax.set_xlim(25, 90); ax.set_ylim(40, 120)
    plt.xlabel(r"$\alpha$ (deg)"); plt.ylabel(r"$\theta$ (deg)")
    plt.title(f"{sim_code} comparison of KDEs")
    plt.savefig(f"{sim_code}_KDEcountours.png", dpi=600)

exit()

"""
Calculate the Kullback-Leibler divergence (and other metrics). 
Remember, P(x) is the reference PDF i.e. Tozzini angles. Q(x) is the tested PDF i.e. simulated angles
I am making an important choice that the values of x from the set X are taken from the reference Tozzini angles dataset.
This is because I am wanting to know how much the simulated distribution differs from the reference data.

In the below sections, alpha = from the set of Tozzini reference angles dataset
                       gamma = from the simulated angles dataset
"""
sim_alpha_probs = simkernel.pdf(ref_values)
ref_alpha_probs = refkernel.pdf(ref_values)

sim_gamma_probs = simkernel.pdf(sim_values)
ref_gamma_probs = refkernel.pdf(sim_values)

# This part is important for normalizing the probabilities, because remember a set of probabilities should add to 1 (duh)
sim_alpha_probs = 1.0*sim_alpha_probs / np.sum(sim_alpha_probs, keepdims=True)
ref_alpha_probs = 1.0*ref_alpha_probs / np.sum(ref_alpha_probs, keepdims=True)

sim_gamma_probs = 1.0*sim_gamma_probs / np.sum(sim_gamma_probs, keepdims=True)
ref_gamma_probs = 1.0*ref_gamma_probs / np.sum(ref_gamma_probs, keepdims=True)

#Bhattacharya distance calculation
bc, bc_dist, h_dist = bc_metric(ref_alpha_probs, sim_alpha_probs)

# KL calculation, using Tozzini as the set of X and p = Ref, q = Sim
KL_alpha, KL_waste = kl_divergence(ref_alpha_probs, sim_alpha_probs)
KL_gamma, KL_waste = kl_divergence(ref_gamma_probs, sim_gamma_probs)

# And now the Jensen-Shannon distance (sqrt of divergence)
JS = np.sum(spd.jensenshannon(ref_alpha_probs, sim_alpha_probs))


outfile = open(f"{args.o}.csv", 'w')
outfile.write("BC metric,BC distance,Hellinger distance,KL on set of Tozzini,KL on set of simulated,JS\n")
outstring = f"{bc:.4f},{bc_dist:.4f},{h_dist:.4f},{KL_alpha:.4f},{KL_gamma:.4f},{JS:.4f}\n"
outfile.write(outstring)
outfile.close()

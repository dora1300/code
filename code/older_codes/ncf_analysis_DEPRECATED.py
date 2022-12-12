"""
@Title:         ncf_analysis_DEPRECATED.py
@Author:        Mando Ramirez
@Date:          20211209

@Description:   This code implements a native contact fraction analysis for coil helices to calculate the fraction of a
simulation when a coil helix is making or not making native contacts. Native contacts are defined as residues in the 'a'
and 'd' positions existing at appropriate distances from each other.

@Updates:
20220204 - As of Feb 4 2022 this code is decomissioned because it's not correct and is just for record keeping sake.
"""
import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import numpy as np


"""
Set up of argument parser, and parsing of arguments
"""
parser = argparse.ArgumentParser(description="Native contact fraction analysis. This calculates the fraction of a-d"
                                             " contacts in a coil trajectory in every frame that are 'native'.")
parser.add_argument("-t", help="trajectory file with extension (.xtc/.pdb)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-n", help="number of individual coil models in the simulation", required=True, type=int)
parser.add_argument("-f", help="analyze every f'th frame in the trajectory", default=1, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-o", help="name of output file, NO extension", default="output_NCF")
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="switch to control plotting", default=False, type=bool)

args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodels = args.n        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially


# Set up and open the output files
full_output = open(f"{args.o}_full.csv", 'w')
short_output = open(f"{args.o}_short.csv", 'w')
output_string = "Frame no"
for num in range(nmodels):
    output_string += f",M{num+1} AD avg(nm),M{num+1} AD stdev,M{num+1} DA avg(nm),M{num+1} DA stdev"
output_string += "\n"
full_output.write(output_string)
short_output.write("Model,AD-dist (nm),AD-dist stddev, DA-dist (nm), DA-dist stddev\n")


# Open the trajectory and load it
traj = md.load(trajectory, top=topology)



"""
Create the list of A-D atoms pairs and D-A atom pairs for the distance calculation. 
I will create a these lists but only for the first model, then will scale the atom pair indices for each additional model
during the main distance calculation section. *IMPORTANT* I am assuming each model is the same size/length which as of
20211209 has been the only situation I'm handling.
REMEMBER: MDTRAJ 0-INDEXES THE ATOMS IN THE TRAJECTORY. Helix length is passed by the user as 1-INDEXED
"""
# This section creates the list of a and d indices within the first model.
a_atoms = []
d_atoms = []
atom_count = 0
while atom_count < lhelix:
    a_atoms.append(atom_count)
    atom_count += 7
atom_count = 3
while atom_count < lhelix:
    d_atoms.append(atom_count)
    atom_count+= 7

# This section creates the list of atom pair indices that will be needed for the distance calculation. Only contains the
# atom pair indices for the first model.
ad_pairs = []
da_pairs = []
for i in range(len(a_atoms)):
    ad_pairs.append([a_atoms[i], d_atoms[i]])
for i in range(len(d_atoms)):
    try:
        da_pairs.append([d_atoms[i],a_atoms[i+1]])
    except:
        pass



"""
Define the global constants of native contact distances
"""
AD_DIST = 0.505         # (nm) constant, this is the distance going from 'a'-'d'
DA_DIST = 0.63          # (nm) constant, this in the distance in (nm) going from 'd'-'a'

AD_NUM = len(ad_pairs)  # This is the total NUMBER of 'a'-'d' native contacts that are possible
DA_NUM = len(da_pairs)  # This is the total NUMBER of 'd'-'a' native contacts that are possible

TOL = 0.05              # (nm) this is the tolerance factor i.e. the amount that a given a-d or d-a distance can deviate from
                        # the constants above and still count as as native contact


"""
Main distance calculation section
"""
# These four lists are for storing each model's distances, for every frame analyzed, so I can output ALL the data in
# a full_output file
total_avg_ad_dist = []; total_avg_da_dist = []; total_std_ad_dist = []; total_std_da_dist = []
# And these four lists are for storing the global avgs and stddevs for each model. There should only be one number for each
# list for each model i.e. these are 1-dimensional lists and the indices refer to the models in order
model_ad_avgs = []; model_da_avgs = []; model_ad_stds = []; model_da_stds = []


for n in range(nmodels):        # Just for clarification, nmodels is 1-INDEXED, and n is 0-INDEXED
    ad_avgs = []; ad_stds = []; da_avgs = []; da_stds = []
    # ^^ These lists will contain the AVERAGE and STDEV for each type of distance pair from each frame, for each analyzed
    # model i.e. indices are trajectory steps, and values are either Avg or Stdev of the corresponding distances. These
    # get reset for each new model.
    # first, make sure you scale the atom pair indices for the correct model
    if n == 0:
        pass
    else:
        for i in range(len(ad_pairs)):
            ad_pairs[i][0] = ad_pairs[i][0] + (n*lhelix)
            ad_pairs[i][1] = ad_pairs[i][1] + (n*lhelix)
        for i in range(len(da_pairs)):
            da_pairs[i][0] = da_pairs[i][0] + (n*lhelix)
            da_pairs[i][1] = da_pairs[i][1] + (n*lhelix)

    # Now, perform the distance calculations. I will analyze the distances frame by frame, average the total distances
    # in a given frame, and save the average and stddev of those distances from the fth frame into a separate list
    for frame in range(0, traj.n_frames, args.f):
        ad_dist = md.compute_distances(traj[frame], np.array(ad_pairs))
        da_dist = md.compute_distances(traj[frame], np.array(da_pairs))
        ad_avgs.append(np.mean(ad_dist)) ; ad_stds.append(np.std(ad_dist))
        da_avgs.append(np.mean(da_dist)) ; da_stds.append(np.std(da_dist))

    # Plot the individual model Avgs and Stdevs for the AD and DA distances
    if args.plot:
        fig, ax = plt.subplots(figsize=(10,5))
        ax.errorbar(np.arange(0, traj.n_frames, args.f), np.array(ad_avgs), yerr=np.array(ad_stds),
                    color="royalblue", linestyle="-", ecolor="deepskyblue", capsize=None,
                    linewidth=0.75, elinewidth=0.5, label="A-D pairs")
        ax.errorbar(np.arange(0, traj.n_frames, args.f), np.array(da_avgs), yerr=np.array(da_stds),
                    color="maroon", linestyle="-", ecolor="lightcoral", capsize=None,
                    linewidth=0.75, elinewidth=0.5, label="D-A pairs")

        ax.axhline(y=0.505, color="black", linestyle="--", linewidth=1.5, label="A-D reference")
        ax.axhline(y=np.mean(np.array(ad_avgs)), linestyle="--", linewidth=1.5, color="grey",
                   label=f"A-D avg: {np.mean(np.array(ad_avgs)):.3f}")

        ax.axhline(y=0.63, color="darkorange", linestyle="--", linewidth=1.5, label="D-A reference")
        ax.axhline(y=np.mean(np.array(da_avgs)), linestyle="--", linewidth=1.5, color="gold",
                   label=f"D-A avg: {np.mean(np.array(da_avgs)):.3f}")

        plt.xlabel("Trajectory frame")
        plt.ylabel("Distance (nm)")
        plt.legend()
        plt.title(f"Model {n+1} A-D and D-A distance averages and standard deviations")
        plt.savefig(f"{sim_code}_model{n+1}_NCF_distanceplot.png",dpi=600)

    # organize the data into the different datasets
    total_avg_ad_dist.append(ad_avgs)
    total_avg_da_dist.append(da_avgs)
    total_std_ad_dist.append(ad_stds)
    total_std_da_dist.append(da_stds)
    model_ad_avgs.append(np.mean(ad_avgs))
    model_da_avgs.append(np.mean(da_avgs))
    model_ad_stds.append(np.mean(ad_stds))
    model_da_stds.append(np.mean(da_stds))


# Now finish things up by writing stuff to files!
for i in range(nmodels):
    short_output.write(f"Model {n+1},{model_ad_avgs[i]:.3f},{model_ad_stds[i]:.3f},"
                       f"{model_da_avgs[i]:.3f},{model_da_stds[i]:.3f}\n")

for i in range(len(total_avg_ad_dist[0])):
    frame = i * args.f
    outwrite = f"{frame},"
    for j in range(len(total_avg_ad_dist)):
        outwrite += f"{total_avg_ad_dist[j][i]:.3f},{total_std_ad_dist[j][i]:.3f}," \
                    f"{total_avg_da_dist[j][i]:.3f},{total_std_da_dist[j][i]:.3f}\n"
    full_output.write(outwrite)


# Don't forget to close the files!
full_output.close()
short_output.close()
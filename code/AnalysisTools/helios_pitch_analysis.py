"""
@Title:         helios_pitch_analysis.py
@Author:        Mando Ramirez
@Date:          20211202

@Description:   This code is a modification of Chris Walker's code from the analyze_folders GitHub that uses the helios
package to calculate the pitch of a coil/helix in a molecular simulation trajectory.

@Citations:     The original paper that describes the Helios implementation: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00721
                Kevin Hauser's github: https://github.com/kehauser/heliosv1

@Updates:
"""

import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import numpy as np
import subprocess

"""
Function for writing the helios script that controls the helios run
"""
def write_helios(filename, outputname, helix_length):
    """
    This function writes and saves a bash run file for the Helios pitch analyzer program.
    :param filename: (string) the name of the PDB structure file to be analyzed
    :param helix_length: (int) the number of atoms in the PDB structure file "filename
    :return:
        - hfile_name (string) - the name of the helios run file that can be used in a subprocess thread
    """
    hfile_name = "helios_run_file.sh"
    hfile = open(f"{hfile_name}", 'w')
    # This code came from Chris Walker's analyze foldamers
    hfile.write("#!/bin/bash\n")
    hfile.write("\n")
    hfile.write("cat > input << EOF\n")
    hfile.write(f"inputhelix {filename}\n")
    hfile.write(f"helixout_name {outputname}\n")
    hfile.write("coord_type 1\n")
    hfile.write("num_grid 360\n")
    hfile.write(f"natoms {helix_length}\n")
    hfile.write("nframes 1\n")
    hfile.write("grid_phi_beg 0\n")
    hfile.write("grid_phi_end 180\n")
    hfile.write("grid_theta_beg 0\n")
    hfile.write("grid_theta_end 180\n")
    hfile.write("helix_atom_names CA\n")
    hfile.write("print_to_plot 1\n")
    hfile.write("EOF\n")
    # Thanks Chris
    hfile.write("/Users/mramirez/GitHub/heliosv1/bin/helios.o input\n")     # This part is hardcoded. Make sure this
    # is correct before using

    hfile.close()
    return hfile_name

def run_helios(helios_input):
    run_command = f"bash {helios_input}"
    subprocess.run(run_command.split(" "))
    return None

def process_helios(heliosfile):
    lineindex = 1
    with open(heliosfile) as f:
        for line in f:
            if lineindex == 43:
                line43 = line.split()
                radius = float(line43[3])
                pitch = float(line43[4])
            else:
                pass
            lineindex += 1
    return radius, pitch


"""
Set up of argument parser, and parsing of arguments
"""
parser = argparse.ArgumentParser(description="Coil pitch analysis using the Helios package.")
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-nmodel", help="the SINGLE model that you want analyzed i.e. if there are "
                                    "two models and you want model 1 analyzed, type -nmodel 1", required=True, type=int)
parser.add_argument("-f", help="analyze every f'th frame. Default=100. You will be prompted if you select something"
                               " smaller", default=100, type=int)
parser.add_argument("-id", help="moniker for the simulation being analyzed", required=True)
parser.add_argument("-l", help="the length (no. atoms) in the helices. Each helix must be the same length.",
                    required=True, type=int)
parser.add_argument("-plot", help="switch for plotting. Set to True to plot everything, or False to not", default=False,
                    type=bool)
parser.add_argument("-print", help="switch to choose to ONLY print to stdout, no saving files!", default=False, type=bool)


args = parser.parse_args()

trajectory= args.t
topology = args.p
nmodel = args.nmodel        # this is 1-indexed, essentially
sim_code = args.id
lhelix = args.l         # this is also 1-indexed, essentially
frames = args.f

if args.print and args.plot:
    # The Plot switch will override the print switch, since I am going to print regardless.
    args.print = False

# Make sure the user knows what they're doing in the event of smaller f
if frames < 100:
    print(f"You have selected to analyze every {frames}'th frame of the trajectory.")
    print(f"Due to the heavy read/write nature of this code, this could result in high memory usage and long completion times.")
    answer = input("Are you sure you want to continue? (Y/N)   ")
    if answer == "Y" or answer == "y":
        pass
    else:
        print("Exiting the analysis. Please select a different f.")



"""
Load the trajectory and prepare a new trajectory containing only the atoms of interest
"""
traj = md.load(trajectory, top=topology)

# Create a sliced topology that contains only the information for the model that is selected to be analyzed. This is
# necessary for making a selection of the trajectory that only contains the desired model.
atom_indices = []
for i in range(lhelix):
    index = i + ((nmodel-1)*lhelix)
    atom_indices.append(index)

newtraj = traj.atom_slice(atom_indices)

"""
Main computation code
"""
outputpdb = "output.pdb"
total_radii = []
total_pitch = []

# This part loops through the specified frames of the trajectory, saves the trajectory as a PDB (which should only be
# the helix), then analyze that PDB with helios. The function process_helios extracts the useful information. Save
# these numbers in the lists shown above
for i in range(0, newtraj.n_frames, frames):
    newtraj[i].save(outputpdb)
    helios_name = write_helios(outputpdb, "helios.out", lhelix)
    run_helios(helios_name)
    r, p = process_helios("helios.out")
    total_pitch.append(p) ; total_radii.append(r)

avg_pitch = np.mean(np.array(total_pitch))
std_pitch = np.std(np.array(total_pitch))
avg_radius = np.mean(np.array(total_radii))
std_radius = np.std(np.array(total_radii))

# 2021 12 10 - There is an interesting behavior of Helios. I suspect it is very sensitive to the curvature of the helix
# which produces some very, very large pitch values when it is very clear that the helix is not actually that large or
# unfolding, which the number would suggest.
excluded_pitch = []
for i in range(len(total_pitch)):
    if total_pitch[i] > 15:
        pass
    else:
        excluded_pitch.append(total_pitch[i])

if args.print or args.plot:
    print("Full data (no exclusions)")
    print(f"Average helix pitch: {avg_pitch:.4f}")
    print(f"Pitch std dev: {std_pitch:.4f}")
    print(f"Average helix radius: {avg_radius:.4f}")
    print(f"Radius std dev: {std_radius:.4f}")
    print()
    print("Pitch data but without massive values")
    print(f"Average helix pitch: {np.mean(np.array(excluded_pitch)):.4f}")
    print(f"Pitch std dev: {np.std(np.array(excluded_pitch)):.4f}")

if args.plot:
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2)
    ax1.plot(np.arange(0, newtraj.n_frames, frames), np.array(total_pitch), color="darkviolet",
            linestyle="-", linewidth=1.0)
    ax2.plot(np.arange(0, newtraj.n_frames, frames), np.array(total_pitch), color="darkviolet",
             linestyle="-", linewidth=1.0)
    ax2.set_ylim(0, 15)
    ax2.axhline(avg_pitch, color="black", linestyle="--", label=f"Full avg: {avg_pitch:.3f}")
    ax2.axhline(np.mean(np.array(excluded_pitch)), color="red", linestyle="--",
                label=f"Reduced avg: {np.mean(np.array(excluded_pitch)):.3f}")
    ax2.legend()
    fig.supxlabel("Frame no.")
    fig.supylabel("Pitch (Angstroms)")
    fig.suptitle(f"Model {nmodel} helix pitch analysis plot")
    plt.savefig(f"{sim_code}_pitchPlot.png", dpi=600)

"""
Open the save files and write to them
"""
fullout = open(f"{sim_code}_heliosAnalysis_full.csv", 'w')
shortout = open(f"{sim_code}_heliosAnalysis_short.csv", 'w')

shortout.write("Model,Full avg pitch (ang),Full pitch StdDev,Excluded avg pitch,Excluded pitch StdDev,"
               "Avg radius (ang),Radius StdDev\n")
shortout.write(f"{nmodel},{avg_pitch:.4f},{std_pitch:.4f},{np.mean(np.array(excluded_pitch)):.4f},"
               f"{np.std(np.array(excluded_pitch))},{avg_radius:.4f},{std_radius:.4f}\n")

fullout.write("Frame,Pitch (ang),Radius (ang)\n")
for i in range(len(total_pitch)):
    frame = i * frames
    fullout.write(f"{frame},{total_pitch[i]:.4f},{total_radii[i]:.4f}\n")

shortout.close()
fullout.close()
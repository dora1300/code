"""
@Author:            Dominique A Ramirez
@Date:              2022 12 22
@Title:             plot_singlemolecule_Rg.py

@Description:       This script is a general plotting function to plot the Rg analysis from any
simulation OF JUST ONE MOLECULE. This is not set up for multiple molecules.

"""

import argparse
import matplotlib.pyplot as plt
import numpy as np


# Set up arg parse and get the args
parser = argparse.ArgumentParser(description="Plotting tool for plotting the Rg of just ONE "
    "molecule, using the 'gmx gyrate' function to generate the data.")
parser.add_argument('-f', help="Input .csv file for analysis. Usually this comes after 'gmx gyrate'"
    ". INCLUDE EXTENSION!", required=True)
parser.add_argument('-acf', help="Input .csv file for the autocorrelation function analysis,"
    " which also comes from 'gmx gyrate'. Include extension", default=None)
parser.add_argument('-o', help="Name of the output png file.", default="rg_plot.png")
parser.add_argument('-o2', help="Name of the output ACF png file.", default="rg__acf_plot.png")
parser.add_argument('-color', help="Specify the color to use for the plot.", default="blue")
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()


# Open the file and parse the data
Rg = []
Time = []
with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            sp = line.lstrip(" ").rstrip("\n").split(",")
            Time.append(float(sp[0].rstrip()))
            Rg.append(float(sp[1].rstrip()))

# plot the Rg data
fig, ax1 = plt.subplots()
ax1.plot(Time, Rg, color=args.color, linewidth=1, marker="")
ax1.set_ylabel(r"Rg (nm)")
ax1.set_ylim(0, np.amax(np.array(Rg))+0.1)
plt.grid(color="black", linestyle=":", alpha=0.25)
plt.xlabel("Simulation time (ps)")
plt.title(args.title)
plt.savefig(args.o, dpi=300)
plt.close()


# Then plot the ACF data if it's provided
if args.acf is None:
    exit(0)
    
else:
    ACF = []
    Tau = []
    
    with open(args.acf, "r") as f:
        for line in f:
            if line[0] == "#" or line[0] == "@":
                continue
            else:
                sp = line.lstrip(" ").rstrip("\n").split(",")
                Tau.append(float(sp[0].rstrip()))
                ACF.append(float(sp[1].rstrip()))
                fig, ax1 = plt.subplots()
            
    ax1.plot(Tau, ACF, color=args.color, linewidth=1, marker="")
    ax1.set_ylabel(r"Autocorrelation")
    ax1.set_ylim(0, 1)
    plt.grid(color="black", linestyle=":", alpha=0.25)
    plt.xlabel("$\tau$ (frames)")
    plt.title(f"{args.title} - ACF")
    plt.savefig(args.o2, dpi=300)
    plt.close()





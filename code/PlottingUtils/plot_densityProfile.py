"""
@Author:            Dominique A Ramirez
@Date:              2022 06 29
@Title:             plot_densityProfile.py

@Description:       This script is a general plotting script for plotting density profiles for slab simulations.

@Updates:
2022 07 26 - added a color argument

2022 08 04 - fixed the ylimit because I'm only using number density and that can only max out at 1
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plotting tool for plotting density profiles for individual simulations!")
parser.add_argument('-f', help="Input .csv file for analysis. Usually this comes after 'gmx density'. INCLUDE"
                               " EXTENSION!", required=True)
parser.add_argument('-o', help="Name of the output png file.", default="densityprofile_plot.png")
parser.add_argument('-color', help="Specify the color to use for the plot.", default="blue")
parser.add_argument('-slice', help="The number of slices used in the analysis. Use the same number as for the '-sl' "
                                   "flag passed into 'gmx density'", required=True, type=int)
parser.add_argument('-center', help="Provide flag to mark the center of the box as '0'. Otherwise, "
                                    "z-axis will be listed using the coordinate slice from the "
                                    "analysis.", action="store_true", default=False)
parser.add_argument("-bl", help="The box length in nm. This must be provided if center is also provided.", type=int)
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()

if args.center and args.bl is None:
    print("If you choose to center the graph, you must provide the box length. Exiting now")
    exit()

Slices = []
Density = []

with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            sp = line.lstrip(" ").rstrip("\n").split(",")
            Slices.append(float(sp[0].rstrip()))
            Density.append(float(sp[1].rstrip()))

if args.center:
    arSlice = np.linspace(-args.bl/2, args.bl/2, args.slice)
else:
    arSlice = np.array(Slices)

fig, ax1 = plt.subplots()
ax1.plot(arSlice, np.array(Density), color=args.color, label = "Density",
         linewidth=1)
ax1.set_ylabel(r"$\rho$ (number density)")
ax1.set_ylim(0, 1)
#ax1.legend()
plt.grid(color="black", linestyle=":", alpha=0.5)
plt.xlabel("z-coordinate (nm)")
plt.title(args.title)
plt.savefig(args.o, dpi=300)
plt.close()

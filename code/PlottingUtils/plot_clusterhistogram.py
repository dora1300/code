"""
@Author:            Dominique A Ramirez
@Date:              2022 09 14
@Title:             plot_clusterhistogram.py

@Description:       This script is a general script to plot cluster distributions from single simulation runs. The
analysis typically comes from gmx clustsize

@Updates:

"""

import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plotting tool for plotting cluster size distributions from "
                                             "single simulations!")
parser.add_argument('-f', help="Input .csv file for analysis. Usually this comes after 'gmx clustsize'. INCLUDE"
                               " EXTENSION!", required=True)
parser.add_argument('-o', help="Name of the output png file.", default="clusthisto_plot.png")
parser.add_argument('-color', help="Specify the color to use for the plot.", default="blue")
parser.add_argument('-nmodels', help="The number of models in the system.", required=True, type=int)
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()


# parse through the file
x_axis = []
probabilities = []
with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            line_split = line.rstrip("\n").split(",")
            x_axis.append(int(line_split[0]))
            probabilities.append((int(line_split[0]) * float(line_split[1])) / args.nmodels)

fig, ax = plt.subplots()
ax.plot(x_axis, probabilities, color=args.color, linestyle="-", linewidth=2)
plt.xlim(0, args.nmodels+1)
plt.ylim(0, 1)
plt.yticks(ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.xticks(np.arange(0, args.nmodels+1, 2))
plt.title(args.title)
plt.grid(color="black", linestyle=":", alpha=0.25)
plt.xlabel("Cluster size")
plt.ylabel("Probability")

plt.savefig(args.o, dpi=300)

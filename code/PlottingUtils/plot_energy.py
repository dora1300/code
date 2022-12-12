"""
@Author:            Dominique A Ramirez
@Date:              2022 06 27
@Title:             plot_energy.py

@Description:       This script is a general plotting script for plotting the internal energies
for the system, using a csv file as input.

@Updates:
2022 06 29 - Added a title argparse feature to add to the plots if I desire
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plotting tool for plotting Potential, Kinetic, and Total energy"
                                             " for a simulation!")
parser.add_argument('-f', help="Input .csv file for analysis. U, K, T_e must be present, in addition to time.",
                    required=True)
parser.add_argument('-o', help="Name of the output png file.", default="energies_plot.png")
parser.add_argument('-colorU', help="Specify the color to use for the Potential energy.", default="blue")
parser.add_argument('-colorK', help="Specify the color to use for the Kinetic energy.", default="red")
parser.add_argument('-colorT', help="Specify the color to use for the Total energy.", default="black")
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()

Time = []
Potential = []
Kinetic = []
TotalE = []

with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            sp = line.lstrip(" ").rstrip("\n").split(",")
            Time.append(float(sp[0].rstrip()))
            Potential.append(float(sp[1].rstrip()))
            Kinetic.append(float(sp[2].rstrip()))
            TotalE.append(float(sp[3].rstrip()))

arTime = np.array(Time)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(15, 10))
ax1.plot(arTime, np.array(Potential), color=args.colorU, label="Potential energy",
         linewidth=0.5)
ax1.set_ylabel("kJ/mol")
ax1.legend()
ax1.grid(color="black", linestyle="-", alpha=0.25)
ax2.plot(arTime, np.array(Kinetic), color=args.colorK, label="Kinetic energy",
         linewidth=0.5)
ax2.legend()
ax2.set_ylabel("kJ/mol")
ax2.grid(color="black", linestyle="-", alpha=0.25)
ax3.plot(arTime, np.array(TotalE), color=args.colorT, label="Total energy",
          linewidth=0.5)
ax3.set_ylabel("kJ/mol")
ax3.grid(color="black", linestyle="-", alpha=0.25)
ax3.legend()
plt.xlabel("Time (ps)")
plt.suptitle(args.title)
plt.savefig(args.o, dpi=300)
plt.close()

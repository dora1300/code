"""
@Author:            Dominique A Ramirez
@Date:              2023 07 10
@Title:             plot_energy.py

@Description:       This script is a general plotting script for plotting the temperature 
of a simulation!

@Updates:

"""
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plotting tool for plotting the "
                                    "temperature of a simulation!")
parser.add_argument('-f', help="Input .csv file for analysis. Only temperature allowed in file.",
                    required=True)
parser.add_argument('-o', help="Name of the output png file.", default="temperature_plot.png")
parser.add_argument('-color', help="Specify the color to use for the plot", 
        default="forestgreen")
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()

Time = []
Temperature = []


with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            sp = line.lstrip(" ").rstrip("\n").split(",")
            Time.append(float(sp[0].rstrip()))
            Temperature.append(float(sp[1].rstrip()))


arTime = np.array(Time)

fig, ax1 = plt.subplots()
ax1.plot(arTime, np.array(Temperature), color=args.color, label="Temperature",
         linewidth=0.5)
plt.ylabel("Temperature (K)")
plt.grid(color="black", linestyle=":", alpha=0.35)
plt.xlabel("Time (ps)")
plt.title(args.title)
plt.tight_layout()
plt.savefig(args.o, dpi=300)
plt.close()

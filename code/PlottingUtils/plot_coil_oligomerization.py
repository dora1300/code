"""
@Author:            Dominique A Ramirez
@Date:              2022 09 12
@Title:             plot_coil_oligomerization.py

@Description:       This script is a companion script to the coil_oligomerization_analysis.py script. This allows the
user to plot the analyzed data as histograms and not have to re-run the analysis.

@Updates:

"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Plotting tool for plotting histograms of coil oligomers following"
                                             " the coil_oligomerization_analysis.py analysis.")
parser.add_argument('-f', help="Input .csv file for analysis INCLUDE"
                               " EXTENSION!", required=True)
parser.add_argument('-o', help="Name of the output png file.", default="oligomer_histo_plot.png")
parser.add_argument('-nmodels', help="The number of models in the system.", required=True, type=int)
parser.add_argument("-title", help="A title for your plot, if you so desire")

args = parser.parse_args()

# Read the input .csv file into a pandas workbook
df = pd.read_csv(args.f)

fig, ax = plt.subplots()
ax.hist(np.array(df.loc[2]), density=True, label="Dimers", color="royalblue", alpha=0.35)
ax.hist(np.array(df.loc[3]), density=True, label="Trimers", color="orangered", alpha=0.35)
ax.hist(np.array(df.loc[4]), density=True, label="Tetramers", color="limegreen", alpha=0.35)
ax.hist(np.array(df.loc[5]), density=True, label="Higher order", color="darkorchid", alpha=0.35)

plt.ylim(0, 1)
plt.xlim(1, args.nmodels+1)
plt.show()
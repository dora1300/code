"""
@Title:             _plotter_make_summary_heatmap.py
@Name:              Mando A Ramirez
@Date:              2025 10 30

@Description:       The purpose of this script is to plot a summary heatmap following the completion
and analysis of a set of validation simulations.

@Notes:             This code was enhanced with Microsoft Copilot (provided through CU)
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import math

parser = argparse.ArgumentParser(description="A tool to make a summary heatmap for a completed set of validation simulations")

parser.add_argument("-codenames", help="the list of protein codenames from the validation set that you want to make a "
                "heat map for.", nargs="+", type=str)
parser.add_argument("-analysis_path", help="The location (absolute path!) of the combined analysis file for the validation "
                    "in question. Don't include the final slash.", type=str)
parser.add_argument("-output_name", help="A name to add to the final heatmap", type=str)

args = parser.parse_args()

codenames = args.codenames
ANPATH = args.analysis_path
OUTPUTNAME = args.output_name

# Load arrays
total_array = np.empty(len(codenames), dtype=object)
for codei, code in enumerate(codenames):
    total_array[codei] = np.loadtxt(f"{ANPATH}/analysis_{code}/heatmap_{code}_correct_and_total_multimer.csv")

# Reshape combined data
redesigned_array = np.reshape(np.concatenate(total_array), (len(codenames), 7)).T

# Figure setup
fig_width = (len(codenames) / 10) * 7
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(fig_width, 6))

# First heatmap
im1 = ax[0].imshow(redesigned_array[3:5, :], cmap="coolwarm", interpolation="none", origin="lower", vmin=0.0, vmax=1.0)
ax[0].set_xticks(np.arange(len(codenames)))
ax[0].set_xticklabels(codenames, rotation=25, fontsize=8)
ax[0].set_yticks([0, 1])
ax[0].set_yticklabels(["Correct\nmultimer", "incorrect\nmultimer"], rotation=30)
ax[0].set_title("fraction time bound in some multimer", fontsize=10)
cbar1 = fig.colorbar(im1, ax=ax[0], shrink=0.75)
cbar1.set_label("fraction of frames")

# Second heatmap
im2 = ax[1].imshow(redesigned_array[5:7, :], cmap="coolwarm", interpolation="none", origin="lower", vmin=0.0, vmax=1.0)
ax[1].set_xticks(np.arange(len(codenames)))
ax[1].set_xticklabels(codenames, rotation=25, fontsize=8)
ax[1].set_yticks([0, 1])
ax[1].set_yticklabels(["correct ETE\ncoil 1", "correct ETE\ncoil 2"], rotation=30)
ax[1].set_title("fraction time coils have correct ETE", fontsize=10)
cbar2 = fig.colorbar(im2, ax=ax[1], shrink=0.75)
cbar2.set_label("fraction of frames")

# Third heatmap
im3 = ax[2].imshow(redesigned_array[0:3, :], cmap="coolwarm", interpolation="none", origin="lower", vmin=0.0, vmax=1.0)
ax[2].set_xticks(np.arange(len(codenames)))
ax[2].set_xticklabels(codenames, rotation=25, fontsize=8)
ax[2].set_yticks([0, 1, 2])
ax[2].set_yticklabels(["correctly\nbound", "incorrectly\nbound", "correct+\nincorrect"], rotation=30)
ax[2].set_title("fraction time bound X fraction correct ETE")
cbar3 = fig.colorbar(im3, ax=ax[2], shrink=0.75)
cbar3.set_label("fraction of time\nsimulation behaves")

# Save figure
plt.tight_layout()
plt.savefig(f"heatmap_{OUTPUTNAME}_summary.png", dpi=600)
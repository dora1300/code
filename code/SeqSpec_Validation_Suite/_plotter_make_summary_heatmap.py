"""
@Title:             _plotter_make_summary_heatmap.py
@Name:              Mando A Ramirez
@Date:              2025 10 30

@Description:       The purpose of this script is to plot a summary heatmap following the completion
and analysis of a set of validation simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="A tool to make a summary heatmap for a completed set of validation simulations")

parser.add_argument("-codenames", help="the list of protein codenames from the validation set that you want to make a " \
                "heat map for.", nargs="+", type=str)
parser.add_argument("-analysis_path", help="The location (absolute path!) of the combined analysis file for the validation " \
                    "in question. Don't include the final slash.", type=str)
parser.add_argument("-output_name", help="A name to add to the final heatmap", type=str)

args = parser.parse_args()

codenames = args.codenames
ANPATH = args.analysis_path
OUTPUTNAME = args.output_name

total_array = np.empty(len(codenames), dtype=object)

for codei, code in enumerate(codenames):
    total_array[codei] = np.loadtxt(f"{ANPATH}/analysis_{code}/heatmap_{code}_correct_and_total_multimer.csv")


figsize_scalar = len(codenames) / 6

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=((6 * figsize_scalar), 7))

im1 = ax[0].imshow(total_array.T[3:5,:], cmap="coolwarm", interpolation="none", origin="lower",
            vmin=0.0, vmax=1.0)
ax[0].set_xticks(np.arange(0, len(codenames)), [])
ax[0].set_yticks([0, 1], ["Correct\nmultimer", "incorrect\nmultimer"], rotation=30)
ax[0].set_title(f"are the coils correctly or\nincorrectly bound? (fraction time bound)")
cbar1 = fig.colorbar(im1, shrink=0.5)
cbar1.set_label("fraction of frames")

im2 = ax[1].imshow(total_array.T[5:7,:], cmap="coolwarm", interpolation="none", origin="lower",
            vmin=0.0, vmax=1.0)
ax[1].set_xticks(np.arange(0, len(codenames)), [])
ax[1].set_yticks([0, 1], [f"correct ETE\ncoil 1", f"correct ETE\ncoil 2"], rotation=30)
ax[1].set_title("do the coils have the\nright shape? (fraction time correct ETE)")
cbar2 = fig.colorbar(im2, shrink=0.5)
cbar2.set_label("fraction of frames")

im3 = ax[2].imshow(total_array.T[0:3, :], cmap="coolwarm", interpolation="none", origin="lower",
            vmin=0.0, vmax=1.0)
ax[2].set_xticks(np.arange(0, len(codenames)), codenames, rotation=45)
ax[2].set_yticks([0, 1, 2], ["correctly\nbound", "incorrectly\nbound", "correct+\nincorrect"], rotation=30)
ax[2].set_title("fraction time bound X fraction correct ETE")

cbar3 = fig.colorbar(im3, shrink=0.5)
cbar3.set_label("fraction of time simulation behaves")
plt.tight_layout()
plt.savefig(f"heatmap_{OUTPUTNAME}_summary.png", dpi=600)

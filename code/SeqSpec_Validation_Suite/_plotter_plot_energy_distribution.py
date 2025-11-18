"""
@Title:             _plotter_plot_energy_distribution.py
@Name:              Mando A Ramirez
@Date:              2025 11 18

@Description:       This script reads in a .csv file containing LJ energy in the second column
and time steps in the first column, which comes from converting an .xvg file after gmx energy.
That's it. Additional functionality maybe added later

@Notes:             
"""


import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12      #deviation form usual but I'm making small plots on purpose
import numpy as np
import argparse
import helperscript_level1_multimerization_dictionary as mult_dict


parser = argparse.ArgumentParser(description="A tool to make energy distributions for " \
                                 "Sequence Specific simulations.")

parser.add_argument("-energy_file", help="[csv] The file containing the LJ energy. " \
                    "Give the full path!!", type=str)

parser.add_argument("-protein_codename", 
                    help="The validation code in question. Useful for making a title.", type=str)

parser.add_argument("-output_name", help="A name to add to the final plot.", type=str)


args = parser.parse_args()
ENERGYFILE = args.energy_file
CODENAME = args.protein_codename
OUTPUTNAME = args.output_name

# I will need the total number of beads!
protein_dictionary = mult_dict[CODENAME]
total_system_beads = np.sum(np.array(protein_dictionary[3]))


energy_file = np.loadtxt(f"{ENERGYFILE}", 
                         delimiter=",", dtype=float, skiprows=24)
total_lj_energy = energy_file[:, 1].flatten()

avg_total_lj_energy_by_bead = total_lj_energy / total_system_beads


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(5, 4))
axs[0].hist(total_lj_energy, density=True, bins="sqrt", color="goldenrod", alpha=0.5)
axs[0].set_title(f"total LJ energy", fontstyle="italic")
axs[0].set_ylabel("density (counts)")

axs[1].hist(avg_total_lj_energy_by_bead, density=True, bins="sqrt", color="goldenrod",
            label=f"Avg: {np.average(avg_total_lj_energy_by_bead):.4f}")
axs[1].vlines(x=np.average(avg_total_lj_energy_by_bead), ymin=0, ymax=3, color="black", linewidth=1.5)
axs[1].set_xlabel("LJ energy (kJ/mol)")
axs[1].set_ylabel("density (counts)")
axs[1].set_title("avg LJ energy per individual bead", fontstyle="italic")
axs[1].legend()

plt.tight_layout()
plt.savefig(f"{OUTPUTNAME}.png", dpi=300)
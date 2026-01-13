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
plt.rcParams['font.size'] = 10      #deviation form usual but I'm making small plots on purpose
import numpy as np
import argparse
import helperscript_level1_multimerization_dictionary as mult_dict


parser = argparse.ArgumentParser(description="A tool to make energy distributions for " \
                                 "Sequence Specific simulations.")

parser.add_argument("-energy_file", help="[csv] The file containing the LJ energy. " \
                    "Give the full path!!", type=str)

parser.add_argument("-protein_codename", 
                    help="The validation code in question. Useful for making a title.", type=str)

parser.add_argument("-energy_dist_description", help="A brief description of the type of " \
                    "energy distribution being " \
                    "plotted. default = 'energy distribution'", type=str, default="energy distribution")

parser.add_argument("-num_beads_normalize", help="The number of beads you wish to divide the total energy by "
                    "to get the avg per-bead energy. If None, then the number of beads comes from dictionary "
                    "and is the total number of beads in the simulation.", default=None, type=int)

parser.add_argument("-normalization_mode", help="Instead of providing the number of beads to normalize to, " \
                    "you can also provide a mode of normalization. " \
                    "`n_multimer_beads` is the number of beads in the multimerization list. " \
                    "`n_total` is the same as not specifying the total number of beads. " \
                    "`n_rest_of_syst` is the difference between 'n_total' and 'n_multimer_beads'. " \
                    "Whatever is specified will supersede the number of beads flag.", default=None,
                    type=str)

parser.add_argument("-distribution_color", help="[default: 'goldenrod'] A custom color to give your energy "\
                    "distribution. It has to be a matplotlib color.", type=str, default='goldenrod')

parser.add_argument("-output_name", help="A name to add to the final plot.", type=str)


args = parser.parse_args()
ENERGYFILE = args.energy_file
CODENAME = args.protein_codename
OUTPUTNAME = args.output_name

# I will need the total number of beads!
if args.normalization_mode == "n_multimer_beads":
    protein_dictionary = mult_dict.validation_dictionary[CODENAME]
    total_system_beads = len(np.array(protein_dictionary[6]).flatten())
elif args.normalization_mode == "n_total":
    protein_dictionary = mult_dict.validation_dictionary[CODENAME]
    total_system_beads = np.sum(np.array(protein_dictionary[3]))
elif args.normalization_mode == "n_rest_of_syst":
    protein_dictionary = mult_dict.validation_dictionary[CODENAME]
    total_system_beads = np.sum(np.array(protein_dictionary[3])) - len(np.array(protein_dictionary[6]).flatten())
else:
    if args.num_beads_normalize is None:
        protein_dictionary = mult_dict.validation_dictionary[CODENAME]
        total_system_beads = np.sum(np.array(protein_dictionary[3]))
    else:
        total_system_beads = args.num_beads_normalize


energy_file = np.loadtxt(f"{ENERGYFILE}", 
                         delimiter=",", dtype=float, skiprows=24)
total_lj_energy = energy_file[:, 1].flatten()

avg_total_lj_energy_by_bead = total_lj_energy / total_system_beads


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(5, 5))
axs[0].hist(total_lj_energy, density=True, bins="sqrt", color=args.distribution_color, alpha=0.5)
axs[0].set_title(f"total {args.energy_dist_description}", fontstyle="italic")
axs[0].set_ylabel("density (counts)")

axs[1].hist(avg_total_lj_energy_by_bead, density=True, bins="sqrt", color=args.distribution_color,
            label=f"Avg: {np.average(avg_total_lj_energy_by_bead):.4f} per {total_system_beads} beads")
axs[1].vlines(x=np.average(avg_total_lj_energy_by_bead), ymin=0, ymax=3, color="black", linewidth=1.5)
axs[1].set_xlabel(f"{args.energy_dist_description} energy (kJ/mol)")
axs[1].set_ylabel("density (counts)")
axs[1].set_title("avg energy per individual bead", fontstyle="italic")
axs[1].legend(fontsize=10)

fig.suptitle(f"{args.energy_dist_description}:\n{CODENAME}")

plt.tight_layout()
plt.savefig(f"{OUTPUTNAME}.png", dpi=300)
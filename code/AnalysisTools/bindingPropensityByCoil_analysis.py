"""
@Title:         bindingPropensityByCoil_analysis.py
@Date:          2023 05 11
@Author:        Dominique A Ramirez

@Description:   This script uses the output from the **contact_analysis.py** script. It
specifically uses the percentage interaction time data for each coil to determine binding
propensity/interaction by coil, across all coils. This is a very simple bit of code but I'm
making it a separate script for workflow reasons. Maybe I'll throw it in to the contact analysis
when I'm done with it.
"""
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="This script analyzed binding propensity per coil in a protein, which"
                                             " averages the percentage coil interactions from *contact_analysis* "
                                             "across the same coil ID. Thus, this produces an average binding "
                                             "propensity for each of the coils in a protein.")
parser.add_argument("-data", help="Percentage coil interactions data file "
                                  "from *contact_analysis.py*. This DOES NOT use the contacts file "
                                  "that contact_analysis or selfVsOther tools use.", required=True)
#parser.add_argument("-ncoils", help="The total number of coils in the simulation", required=True,
#                    type=int)
#parser.add_argument("-nprots", help="The total number of individual proteins in the simulation", required=True,
#                    type=int)
parser.add_argument("-cpp", help="Not C++ lol! The number of coils *per protein*. This analysis"
                                 " assumes you're working with homopolymers only. Support for heteropolymers is not yet"
                                 " supported", type=int, required=True)
parser.add_argument("-output", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)

args = parser.parse_args()

RESULTS = args.data
COILpp = args.cpp
OUTPUT = args.output


binding_propensity_by_coil = np.empty((COILpp,), dtype=object)

for i, obj in enumerate(binding_propensity_by_coil):
    binding_propensity_by_coil[i] = []



with open(RESULTS, 'r') as f:
    header = 1
    for line in f:
        if header == 1:
            header += 1
            continue
        else:
            line_split = line.rstrip("\n").split(",")
            binding_propensity_by_coil[int(line_split[0]) % COILpp].append(float(line_split[1]))



average_binding_propensity = []
stdev_binding_propensity = []
for i, vals in enumerate(binding_propensity_by_coil):
    average_binding_propensity.append(np.average(vals))
    stdev_binding_propensity.append(np.std(vals))



error_config = {'elinewidth': 3,
                'capthick': 1.5,
                'capsize': 7.5}
plt.rcParams['font.size'] = 14
fig, ax = plt.subplots()
ax.bar(np.arange(1, COILpp + 1, 1), average_binding_propensity, 0.8, yerr=stdev_binding_propensity,
       color="grey", ecolor="black", edgecolor="black",
       error_kw=error_config)

plt.ylim(0, 110)
plt.xticks(np.arange(0, COILpp + 1, 1))
plt.xlabel("Protein coil ID")
plt.ylabel("Percentage of simulation time")
plt.tight_layout()
plt.savefig(f"{OUTPUT}_bindingPropensity_byCoil.png", dpi=600)

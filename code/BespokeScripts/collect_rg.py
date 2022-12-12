"""
This code goes through all the different analysis folders and collects the information from the Rg for output and
plotting
"""

import matplotlib.pyplot as plt
import numpy as np

runlist = []            # this contains the numbers associated with the temperature of each simulation
avglist = []
stdlist = []
auclist = []

outputfile = open("collated_Rg.csv", 'w')
outputfile.write("Run No., Avg Rg, StDev Rg, AUC\n")

for i in [200, 230, 250, 280, 293, 298, 310, 350, 400]:
    with open(f"./t_{i}/analysis/t_{i}_rg_analysis.txt") as f:
        line_count = 1
        for line in f:
            if line_count == 2:
                split = line.split("\t")
                try:
                    avg = float(split[2])
                except:
                    avg = None
            elif line_count == 3:
                split = line.split("\t")
                try:
                    std = float(split[2])
                except:
                    std = None
            elif line_count == 4:
                split = line.split("\t")
                try:
                    auc = float(split[1])
                except:
                    auc = None
            line_count += 1

        runlist.append(i) ; avglist.append(avg) ; stdlist.append(std) ; auclist.append(auc)
        outputfile.write(f"{i}, {round(avg,4)}, {round(std,4)}, {round(auc,4)}\n")
outputfile.close()




"""
Plot the Avg Rg of the runs
"""
xlabels = ['200', '230', '250', '280', '293', '298', '310', '350', '400']

fig,ax = plt.subplots(figsize=(12,6))
ax.errorbar(np.arange(len(runlist)), np.array(avglist), yerr=np.array(stdlist), ecolor='k', elinewidth=1, marker='.', markersize=5,
            linestyle="", markerfacecolor='k', markeredgecolor='k', capsize=2)
plt.xlabel("Simulation Temperature")
plt.ylabel("Radius of Gyration (nm)")
plt.title("Average radius of gyration for each run in DAR3-1f rep1")
plt.xticks(np.arange(len(runlist)), xlabels)
plt.savefig("avg_rg_DAR3-1f_rep1_plot.png", dpi=300)



"""
Plot the Rg AUC of the runs
"""
fig,ax = plt.subplots(figsize=(8,4))
ax.errorbar(np.arange(len(runlist)), np.array(auclist), marker='.', markersize=5,
            linestyle="", markerfacecolor='k', markeredgecolor='k')
plt.xlabel("Simulation Temperature")
plt.ylabel("Rg AUC (nm*ps)")
plt.title("Rg AUC for each run in DAR3-1f rep1")
plt.xticks(np.arange(len(runlist)), xlabels)
plt.savefig("auc_rg_DAR3-1f_rep1_plot.png", dpi=300)

"""
This code goes through all the different analysis folders and collects the information from the Rama
analysis for plotting and other things, if needed.
"""

import matplotlib.pyplot as plt
import numpy as np

coil1_minMax_list = []              # will be a list of lists, since I'm storing min and max
coil2_minMax_list = []
coil1_avg_list = []
coil2_avg_list = []
coil1_std_list = []
coil2_std_list = []

# outputfile = open("collated_rama.csv", 'w')
# outputfile.write("Run No., Avg Rg, StDev Rg, AUC\n")
h1_count = 0
h2_count = 0
for i in [200, 230, 250, 280, 293, 298, 310, 350, 400]:
    coil_list = [[], []]
    with open(f"./t_{i}/analysis/t_{i}_rama_analysis.txt") as f:
        line_count = 1
        for line in f:
            if line_count >= 4 and line_count <= 35:
                split = line.split("\t")
                # print(split)
                try:
                    coil_list[0].append(float(split[4]))
                except:
                    pass
            elif line_count == 36:
                split = line.split(" ")
                # print(split)
                try:
                    coil1_minMax_list.append([float(split[3])])
                except:
                    pass
            elif line_count == 37:
                split = line.split(" ")
                # print(split)
                try:
                    coil1_minMax_list[h1_count].append(float(split[3]))
                    h1_count += 1
                except:
                    pass
            elif line_count >= 43 and line_count <= 74:
                split = line.split("\t")
                # print(split)
                try:
                    coil_list[1].append(float(split[4]))
                except:
                    pass
            elif line_count == 75:
                split = line.split(" ")
                # print(split)
                try:
                    coil2_minMax_list.append([float(split[3])])
                except:
                    pass
            elif line_count == 76:
                split = line.split(" ")
                # print(split)
                try:
                    coil2_minMax_list[h2_count].append(float(split[3]))
                    h2_count += 1
                except:
                    pass
            line_count += 1

        coil1_avg_list.append(round(np.average(np.array(coil_list[0])),5))
        coil2_avg_list.append(round(np.average(np.array(coil_list[1])),5))
        coil1_std_list.append(round(np.std(np.array(coil_list[0])),5))
        coil2_std_list.append(round(np.std(np.array(coil_list[1])),5))


# """
# Plot the averages of each coil for each simulation PLUS/MINUS std dev
# """
# xlabels = ['200', '230', '250', '280', '293', '298', '310', '350', '400']
#
# fig, ax = plt.subplots(figsize=(15,8))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels - width/2, np.array(coil1_avg_list), width, yerr=np.array(coil1_std_list),
#             ecolor='k', capsize=3, color='darkorange', label="Coil 1 +/- SD")
# ax.bar(sim_labels + width/2, np.array(coil2_avg_list), width, yerr=np.array(coil2_std_list),
#            ecolor='k', capsize=3, color='royalblue', label="Coil 2 +/- SD")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average fraction of time each simulation coil pair\nhas a helix dihedral/pair conformation")
# plt.legend(loc="lower right")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.savefig("average_rama_bothCoils_stddev.png", dpi=300)
#
# """
# Plot the averages for each coil, but separated into coil1 and then coil2 for easier viewing
# """
# fig, ax = plt.subplots(figsize=(12,6))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels - width/2, np.array(coil1_avg_list), width, yerr=np.array(coil1_std_list),
#             ecolor='k', capsize=3, color='darkorange', label="Coil 1 +/- SD")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average fraction of time Coil 1\nhas a helix dihedral/pair conformation")
# plt.legend(loc="lower right")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.savefig("average_rama_coil1alone_stddev.png", dpi=300)
#
# fig, ax = plt.subplots(figsize=(12,6))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels + width/2, np.array(coil2_avg_list), width, yerr=np.array(coil2_std_list),
#            ecolor='k', capsize=3, color='royalblue', label="Coil 2 +/- SD")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average fraction of time Coil 1\nhas a helix dihedral/pair conformation")
# plt.legend(loc="lower right")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.savefig("average_rama_coil2alone_stddev.png", dpi=300)
#
# plt.close("all")
#
#
# """
# Now plot the same averages but with the MIN and MAX for each coil
# """
# coil1_minmax = [[], []]
# coil2_minmax = [[], []]
#
# for i in range(len(coil1_minMax_list)):
#     coil1_minmax[0].append(coil1_avg_list[i] - coil1_minMax_list[i][0])
#     coil1_minmax[1].append(coil1_minMax_list[i][1] - coil1_avg_list[i])
#     coil2_minmax[0].append(coil2_avg_list[i] - coil2_minMax_list[i][0])
#     coil2_minmax[1].append(coil2_minMax_list[i][1] - coil2_avg_list[i])
#
# #   Plot both coils in one plot
# fig, ax = plt.subplots(figsize=(15,8))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels - width/2, np.array(coil1_avg_list), width, yerr=np.array(coil1_minmax),
#             ecolor='k', capsize=3, color='darkorange', label="Coil 1 +/- Max/Min")
# ax.bar(sim_labels + width/2, np.array(coil2_avg_list), width, yerr=np.array(coil1_minmax),
#            ecolor='k', capsize=3, color='royalblue', label="Coil 2 +/- Max/Min")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average Fraction of time each simulation coil pair\nhas a helix dihedral/pair conformation +/- Max/Min")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.legend(loc="lower right")
# plt.savefig("average_rama_bothCoils_minMax.png", dpi=300)
#
# """
# Plot the averages for each coil, but separated into coil1 and then coil2 for easier viewing
# """
# fig, ax = plt.subplots(figsize=(12,6))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels - width/2, np.array(coil1_avg_list), width, yerr=np.array(coil1_minmax),
#             ecolor='k', capsize=3, color='darkorange', label="Coil 1 +/- Max/Min")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average Fraction of time Coil 1\nhas a helix dihedral/pair conformation +/- Max/Min")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.legend(loc="lower right")
# plt.savefig("average_rama_coil1alone_minMax.png", dpi=300)
#
# fig, ax = plt.subplots(figsize=(12,6))
# sim_labels = np.arange(1, len(xlabels)+1)
# width = 0.3
# ax.bar(sim_labels + width/2, np.array(coil2_avg_list), width, yerr=np.array(coil1_minmax),
#            ecolor='k', capsize=3, color='royalblue', label="Coil 2 +/- Max/Min")
# plt.xticks(sim_labels, xlabels)
# plt.xlabel("Simulation temperature")
# plt.ylabel("Fraction of simulation time")
# plt.title("Average Fraction of time Coil 2\nhas a helix dihedral/pair conformation +/- Max/Min")
# plt.axhline(y=0.5, color='k', linestyle='--')
# plt.legend(loc="lower right")
# plt.savefig("average_rama_coil2alone_minMax.png", dpi=300)
#
# plt.close("all")


"""
But all this is meaningless unless I write out a file!
"""
outfile = open("collated_rama.csv", 'w')
outfile.write("Sim No.,Coil1_avg,coil1_stdDev,coil1_max,coil1_min,Coil2_avg,coil2_stdDev,coil2_max,coil2_min\n")
for i in range(len(xlabels)):
    outfile.write(f"{xlabels[i]},{coil1_avg_list[i]},{coil1_std_list[i]},{coil1_minMax_list[i][1]},{coil1_minMax_list[i][0]},"
                  f"{coil2_avg_list[i]},{coil2_std_list[i]},{coil2_minMax_list[i][1]},{coil2_minMax_list[i][0]}\n")

outfile.close()
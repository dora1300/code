"""
This code goes through all the different analysis folders and collects the information from the rod length analysis
for combined output
"""

import matplotlib.pyplot as plt
import numpy as np


outputfile = open("collated_rodlength.csv", 'w')
outputfile.write("Run_id, h1_prop, h2_prop\n")

helix1_raw = []
helix2_raw = []
temp_str = []

for i in [200, 230, 250, 280, 293, 298, 310, 350, 400]:
    with open(f"./t_{i}/analysis/output_rodanalysis_proportions.csv") as f:
        line_count = 1
        for line in f:
            if line_count == 1:
                pass
            elif line_count == 2:
                split = line.rstrip("\n").split(",")
                helix1_raw.append(int(split[3]))
            else:
                split = line.rstrip("\n").split(",")
                helix2_raw.append(int(split[3]))
            line_count += 1
        temp_str.append(str(i))

for i in range(len(temp_str)):
    outputfile.write(f"{temp_str[i]},{helix1_raw[i]},{helix2_raw[i]}\n")

outputfile.close()

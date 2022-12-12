"""
@Author:            Dominique A Ramirez
@Date:              2022 12 12
@Title:             density_equilibration_metric_alpine.py

@Description:       This script calculates a density equilibration metric for slab simulations. It will perform 'gmx
density' on a trajectory for however many frames I want, extract the relevant densities and average them for a 'high
density' and 'low density' region, then plot over time. This will be useful for determining if I have actually reached
a stable density equilibrium in the simulation and will help in the decision of what is or is not phase separated.

This is specifically for the use on Alpine from CU RC. All the appropriate gmx modules are updated for the supercomputer.

@Updates:

"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os


parser = argparse.ArgumentParser()
parser.add_argument('-t', help="(.xtc) Trajectory file to perform analysis on. Make sure this has been PBC corrected and "
                               "verify a slab exists in the center before proceeding!")
parser.add_argument('-s', help="(.tpr) Gromacs run-file, used as the topology for the gmx density command")
parser.add_argument('-dt', help="[ps] Analyze every dt ps frame in the trajectory. Default is every frame.", default=1,
                    type=int)
parser.add_argument('-b', help="[ps] Starting simulation time to analyze density", default=0, type=int)
parser.add_argument('-e', help="[ps] Ending simulation time to analyze density", required=True, type=int)
parser.add_argument('-sl', help="[nm] The size of the box along the Z-axis. The entirety of the Z-box "
                                "is used for analysis so provide the full length", required=True,
                    type=int)
parser.add_argument('-high', help="The z-region that corresponds to the high density region. Please only"
                                  "provide 2 values", nargs='+', required=True, type=float)
parser.add_argument('-low', help='The z-regions that correspond to the low density region. Please only provide 4 '
                                 'values.', nargs='+', required=True, type=float)
parser.add_argument('-image', help="Name for the output graph of density vs time. Please include extension!",
                    default="output.png")
parser.add_argument('-title', help="ID to include in the density vs time plot.", default=None)


args = parser.parse_args()

xtc = args.t
tpr = args.s
dt = args.dt

high_z_lim = args.high; low_z_lim = args.low
assert len(high_z_lim) == 2, "Region for high density is not 2 values. Try again!"
assert len(low_z_lim) == 4, "Region for low density is not 4 values. Try again!"


os.mkdir("./temp_dens")

time_array = np.arange(args.b, args.e, args.dt)
high_dens_arr = []
low_dens_arr = []

for time in time_array:
    TIME = int(time)
    high_region = []; low_region = []
    gmx_density = f"mpirun -np 1 gmx_mpi density -f {xtc} -s {tpr} -o ./temp_dens/frame{TIME}.xvg -sl {args.sl} -dens number -b {TIME} -e {TIME}"
    dens_input = subprocess.Popen(('echo', '1'), stdout=subprocess.PIPE)
    subprocess.check_output(gmx_density.split(), stdin=dens_input.stdout)
    dens_input.wait()

    convert_xvg = f"python3 /projects/dora1300/code/code/Utils/xvg2csv.py -f ./temp_dens/frame{TIME}.xvg -o ./temp_dens/frame{TIME}.csv"
    subprocess.call(convert_xvg.split())

    with open(f"./temp_dens/frame{TIME}.csv") as f:
        for line in f:
            if line[0] == "#":
                continue
            elif line[0] == "@":
                continue
            else:
                line_split = line.rstrip("\n").split(",")

                if float(line_split[0]) >= high_z_lim[0] and float(line_split[0]) <= high_z_lim[1]:
                    high_region.append(float(line_split[1]))
                elif float(line_split[0]) >= low_z_lim[0] and float(line_split[0]) <= low_z_lim[1]:
                    low_region.append(float(line_split[1]))
                elif float(line_split[0]) >= low_z_lim[2] and float(line_split[0]) <= low_z_lim[3]:
                    low_region.append(float(line_split[1]))
                else:
                    pass

    high_dens_arr.append(np.average(high_region))
    low_dens_arr.append(np.average(low_region))



fig, ax = plt.subplots()
ax.plot(time_array/1000000, high_dens_arr, marker="", linestyle="-", color="crimson", label="High density region")
ax.plot(time_array/1000000, low_dens_arr, marker="", linestyle="-", color="darkblue", label="Low density region")
plt.ylim(0, 1)
plt.xlabel(r"Simulation time ($\mu$s)")
plt.ylabel(r"$\rho$ (number density)")
plt.title(f"Density equilibration plot: {args.title}")
plt.savefig(args.image, dpi=300)



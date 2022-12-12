"""
@Title:			kineticEnergy_validation.py
@Author:		Mando Ramirez
@Date:			20220509

@Description:	This script handles the kinetic energy equipartition assessment
using the physical validation package from Pascal, and was inspired by a workflow 
from Wei-Tse Hsu.
    This tool must take .trr files because that is where kinetic data is stored.

@Updates:
2022 06 23 - corrected the flag for the equipartition option

2022 06 24 - fixed some pathological behavior of flags not trigger the desired result
"""

import physical_validation
import argparse
import numpy as np
import time

# set up the arg parser like always!
parser = argparse.ArgumentParser(description="Kinetic energy validation script. This performs BOTH the Boltzmann "
                                             "Distribution validation, and the equipartition validation. You can select"
                                             " only one or the other with the appropriate flags. Default behavior is to"
                                             " do both.\n"
                                             "WARNING: Please note that this can ONLY be run for stable simulations"
                                             " where the temperature does not explode. Doing so will cause the code "
                                             "to never finish and it will eat up all your RAM.")

parser.add_argument("-t1", help="trajectory file with extension (.trr preferred)", required=True)
parser.add_argument("-p1", help="topology file", required=True)
parser.add_argument("-m1", help="mdp file (use the mdpout.mdp file)", required=True)
parser.add_argument("-e1", help="energy file", required=True)
parser.add_argument("-boltzmann", help="Run the Boltzmann distribution validation", action="store_true",
                        default=False)
parser.add_argument("-equipartition", help="Run the equipartition analysis validation", action="store_true",
                        default="False")
parser.add_argument("-strict", help="Switch to control if the equipartition analysis is done in strict, or not strict"
                                    " mode. Default is not strict", action="store_true", default=False)
parser.add_argument("-filename", help="Provide moniker to append to default filenames, otherwise none is used",
                        default="")

args = parser.parse_args()



# Set up the time analysis, since these take a loooong time
time_start = time.time()



# Load in the trajectory!
physparser = physical_validation.data.GromacsParser()

sim = physparser.get_simulation_data(
    trr=args.t1,
    top=args.p1,
    mdp=args.m1,
    edr=args.e1
)

#If both flags are provided, then do both tests:
if args.boltzmann and args.equipartition:
    if args.strict:
        ked = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=True,
            filename=f"KE_distribution_strict{args.filename}"
        )
        print("Kinetic energy distribution - strict")
        print(ked)
        print()
        keequi_strict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=True,
            filename=f"KE_equipartition_strict{args.filename}"
        )
        print("Equipartition analysis - strict")
        print(keequi_strict)
        print()
    else:       # not strict analysis
        ked_ns = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=False,
            filename=f"KE_distribution_notstrict{args.filename}"
        )
        print("Kinetic energy distribution - not strict")
        print(ked_ns)
        print()
        keequi_notstrict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=False,
            filename=f"KE_equipartition_notstrict{args.filename}"
        )
        print("Equipartition analysis - not strict")
        print(keequi_notstrict)
        print()
    time_finish = time.time()
    print(f"Time elapsed: {time_finish - time_start:.4f} seconds")
elif args.boltzmann and not args.equipartition:
    if args.strict:
        ked = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=True,
            filename=f"KE_distribution_strict{args.filename}"
        )
        print("Kinetic energy distribution - strict")
        print(ked)
        print()
    else:
        ked_ns = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=False,
            filename=f"KE_distribution_notstrict{args.filename}"
        )
        print("Kinetic energy distribution - not strict")
        print(ked_ns)
        print()
        time_finish = time.time()
        print(f"Time elapsed: {time_finish - time_start:.4f} seconds")
elif args.equipartition and not args.boltzmann:
    if args.strict:
        keequi_strict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=True,
            filename=f"KE_equipartition_strict{args.filename}"
        )
        print("Equipartition analysis - strict")
        print(keequi_strict)
        print()
    else:
        keequi_notstrict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=False,
            filename=f"KE_equipartition_notstrict{args.filename}"
        )
        print("Equipartition analysis - not strict")
        print(keequi_notstrict)
        print()
        time_finish = time.time()
        print(f"Time elapsed: {time_finish - time_start:.4f} seconds")
else:
    if args.strict:
        ked = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=True,
            filename=f"KE_distribution_strict{args.filename}"
        )
        print("Kinetic energy distribution - strict")
        print(ked)
        print()
        keequi_strict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=True,
            filename=f"KE_equipartition_strict{args.filename}"
        )
        print("Equipartition analysis - strict")
        print(keequi_strict)
        print()
    else:       # not strict analysis
        ked_ns = physical_validation.kinetic_energy.distribution(
            data=sim,
            strict=False,
            filename=f"KE_distribution_notstrict{args.filename}"
        )
        print("Kinetic energy distribution - not strict")
        print(ked_ns)
        print()
        keequi_notstrict = physical_validation.kinetic_energy.equipartition(
            data=sim,
            strict=False,
            filename=f"KE_equipartition_notstrict{args.filename}"
        )
        print("Equipartition analysis - not strict")
        print(keequi_notstrict)
        print()
    time_finish = time.time()
    print(f"Time elapsed: {time_finish - time_start:.4f} seconds")

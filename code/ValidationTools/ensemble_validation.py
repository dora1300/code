"""
@Title:			ensemble_validation.py
@Author:		Mando Ramirez
@Date:			20220509

@Description:	This script handles the ensemble dynamics validation using the physical validation
package from Pascal, and was inspired by a workflow from Wei-Tse Hsu.
    I believe this ensemble validation tool is designed for .trr files but can handle .xtc if that
is all that is provided.

@Updates:
2022 06 23 - Added a feature to turn on/off passing uncorrelated data to the analysis engine,
which will allow me to turn off/on auto-equilibration detection.
"""

import physical_validation
import argparse
import numpy as np
import time

# set up the arg parser like always!
parser = argparse.ArgumentParser(description="Ensemble validation of simulations. 2 simulations, identical in all"
                                             " set ups except for a change in one state variable, is required.\n"
                                             "This script can also be used to tell you how to change the state variable"
                                             " of the second simulation for proper analysis.\n"
                                             "This analysis directs all output directly to stdout, it DOES NOT save a"
                                             " separate file. Make sure to handle this appropriately.")

parser.add_argument("-check", help="provide this flag to only run the interval estimator for the ensemble validation"
                                   " analysis. Only information for one trajectory is necessary, obviously.",
                    action="store_true")
parser.add_argument("-t1", help="trajectory file with extension (.trr preferred) - 1st simulation", required=True)
parser.add_argument("-p1", help="topology file - 1st simulation", required=True)
parser.add_argument("-m1", help="mdp file (use the mdpout.mdp file) - 1st simulation", required=True)
parser.add_argument("-e1", help="energy file - 1st simulation", required=True)
parser.add_argument("-t2", help="trajectory file with extension (.trr preferred) - 2nd simulation")
parser.add_argument("-p2", help="topology file - 2nd simulation")
parser.add_argument("-m2", help="mdp file (use the mdpout.mdp file) - 2nd simulation")
parser.add_argument("-e2", help="energy file - 2nd simulation")
parser.add_argument("-uncorr", help="provide flag to turn off auto-equilibration feature", action="store_true", 
                     default=False)
args = parser.parse_args()

# Begin the timing. This is important because these analyses can take a looooong time
time_start = time.time()

# load the first trajectory simulation --> this must happen regardless if the check happens
physparser = physical_validation.data.GromacsParser()

sim1 = physparser.get_simulation_data(
    trr=args.t1,
    top=args.p1,
    mdp=args.m1,
    edr=args.e1
)

# Handle the switch case for interval estimation
if args.check:
    print(f"Doing interval estimation for trajectory: {args.t1}")
    interval = physical_validation.ensemble.estimate_interval(
        data=sim1,
    )
    print(interval)
    print()
    time_finish = time.time()
    print(f"Time elapsed: {time_finish-time_start:.4f} seconds")
    exit(0)
else:
    # now this is the case for doing the actual ensemble validation
    if args.t2 is not None and args.p2 is not None and args.m2 is not None and args.e2 is not None:
        pass
    else:
        print("Files for second trajectory not provided, exiting now.")
        exit(1)
    print("Doing ensemble validation of following trajectories:")
    print(f"{args.t1}  and  {args.t2}")

    sim2 = physparser.get_simulation_data(
        trr=args.t2,
        top=args.p2,
        mdp=args.m2,
        edr=args.e2
    )

    ensem = physical_validation.ensemble.check(
        data_sim_one=sim1,
        data_sim_two=sim2,
        screen=False,
        filename="ensemble_validation_plot",
        verbosity=2,
        data_is_uncorrelated=args.uncorr
    )

    time_finish = time.time()
    print(f"Time elapsed: {time_finish-time_start:.4f} seconds")

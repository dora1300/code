"""
@Title:             mismatch_coil_interaction_frequency.py
@Name:              Mando A Ramirez
@Date:              2024 10 29

@Description:       This script is used for the nonspecific mismatch coil simulations.
This calculates the frequency that each type of coil interacts with all other types of coils.
Consider the following mismatched 5-coil trimer protein:
        5(i) -- 5(ii) -- 5(iii) -- 3(iv) -- 3(v)

This analysis will generate a matrix and show the frequency that coil i interacts with all other
coils i--v. Of course this compiles both intra and inter chain interactions and really doesn't 
discriminate between the two.

This is useful to help understand how the mismatch coils are interacting. The hypothesis
is that 5--5 and 3--3 are more frequent, stable, and likely than 5--3. There's only one way
to find out, though.
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import time



"""
Set up the argparser!
"""
parser = argparse.ArgumentParser(description="This script calculates one type of a Frustration Metric: "
                                 "The number of coils that are (1) unbound; (2) bound in inter-chain "
                                 "interactions; and (3) bound in intra-chain interactions.")
parser.add_argument("-contacts_data", help="Contacts list data file multimerization analysis -- list of coil multimer"
                                  " contacts, each line marks an ANALYZED frame (.csv)."
                                  " Please provide path if not in working directory.", required=True)
parser.add_argument("-coils_per_protein", help="The number of coil segments per protein. CURRENTLY, this analysis only"
                    " supports simulations that use proteins of the same type, i.e. all proteins are C6L5. "
                    "Mixed coils-per-protein types will be added later.", required=True, type=int)
parser.add_argument("-total_coils", help="The total number of coils in the simulation", required=True,
                    type=int)
# parser.add_argument("-ft", help="[ps] Frame time, i.e. the amount of time that each frame in the Contacts Data File "
#                                 "is worth in simulation. Integers only. Only necessary for making the final graph",
#                             type=int, default=None)
parser.add_argument("-output_name", help="A common name to give to output files. Default = [output]",
                    default="output", type=str)


args = parser.parse_args()
CONTACTS_FILE = args.contacts_data
CPP = args.coils_per_protein
TOT_COILS = args.total_coils
output_name = args.output_name
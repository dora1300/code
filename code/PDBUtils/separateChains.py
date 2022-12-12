
import prody as pd
import argparse
import os.path

"""
Name of script: separateChains.py
Author: Dominique Ramirez
Date: 2021 07 30
Version: 1.0

This script contains a function for separating chains from a *single* into different pdb files, one for each chain.
"""

def separateChains(pdb, user_chain, outPrefix):
    """
    This function extracts chains from a PDB file and then saves them as separate PDB files.
    :param pdb: a valid PDB file for parsing with ProDy
    :param user_chain: a list of chains to be extracted from the file,
    :param outPrefix: a string for prefixing output files, if desired
    :return: Nothing, but it does save the chains as different pdbs
    """
    protein = pd.parsePDB(pdb)
    hv = protein.getHierView()
    if len(user_chain) == 0:
        chains = list(hv)
        for i in range(len(chains)):
            file_name = outPrefix+"_"+"chain_"+str((i+1))+".pdb"
            pd.writePDB(file_name, chains[i])
    else:
        for i in range(len(user_chain)):
            file_name = outPrefix+"_"+"chain_"+user_chain[i]
            pd.writePDB(file_name, hv[user_chain[i]])
    return None

# -- Set up the argument parser
parser = argparse.ArgumentParser("PDB Utils - Separate Chains.\n"
                                 "Separate chains from a given PDB file and save them! Automated or user provided chains.")
parser.add_argument("pdbFile", action="store", type=str, help="PDB file to be parsed. Must include extension. Give full path if not in directory!")
parser.add_argument("-u", "--user", type=list, default=[], help="(Optional) User provided chains to extract, if known. Please provide as a list.")
parser.add_argument("-op","--outputPrefix", type=str, default="", help="(Optional) Prefix for output files, if desired.")

args = parser.parse_args()

# -- Parse the args and make sure the PDB exists.
if args.pdbFile is None:
    print("No input file given, exiting now.")
    exit(0)
else:
    if os.path.exists(args.pdbFile):
        input_pdb = args.pdbFile
    else:
        print("Input file does not exist! Please try again.")
        exit(0)

user_provided_chains = args.user
prefix = args.outputPrefix

separateChains(input_pdb, user_provided_chains, prefix)
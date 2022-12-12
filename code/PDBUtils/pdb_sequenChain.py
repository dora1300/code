"""
Name of script: pdb_sequenChain.py
Author: Dominique Ramirez
Date: 2021 10 12
Version: 1.0

This script changes the chain identifiers of multiple structures in a PDB file, but requires that each one is separated
by a TER statement. This does not handle MODEL/ENDMDL commands.
"""

import argparse
import string

parser = argparse.ArgumentParser("PDB Utils - Sequential Chain Modifier. Changes every chain ID to be sequential "
                                 "in order of appearance.")
parser.add_argument("file", help="The PDB file that needs chain re-sequentiation.")

args = parser.parse_args()

inputfile_name = args.file.split(".")
outputfile = open(f"{inputfile_name[0]}_seqCh.pdb", 'w')

chain_counter = 0
with open(args.file, 'r') as f:
    for line in f:
        line_split = line.split(" ")

        if line[0:4] != "ATOM" and line[0:3] != "TER":
            # This is mostly for the start of the file. Keep header stuff consistent.
            outputfile.write(line)
        elif line[0:4] == "ATOM":
            # Now we're dealing with actual atoms, time to start changing the chain ID.
            new_line = line[0:21]+ string.ascii_uppercase[chain_counter] + line[22:]
            outputfile.write(new_line)
        elif line[0:3] == "TER":
            # This is my swtich to increment the chain counter.
            outputfile.write(line)
            chain_counter += 1
        else:
            # in case there are other commands in there that I am not accounting for, I will just output them and
            # prevent myself from doing anything strange to the output.
            outputfile.write(line)

outputfile.close()
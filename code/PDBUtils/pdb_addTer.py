"""
Name of script: pdb_addTer.py
Author: Dominique Ramirez
Date: 2021 10 12
Version: 1.0

This script is useful for adding TER statements into PDB files, especially after packmol since it fails to do so.
"""

import argparse

parser = argparse.ArgumentParser("PDB Utils - Termination Adder. Adds TER commands in the right place!")
parser.add_argument("file", help="The PDB file that needs to have TER statements added into it.")

args = parser.parse_args()

inputfile_name = args.file.split(".")
outputfile = open(f"{inputfile_name[0]}_ter.pdb", 'w')

ter_counter = 0
with open(args.file, 'r') as f:
    for line in f:
        line_split = line.split(" ")
        if line_split[0] != "ATOM":
            # this is for header information and other commands, put those directly into the output file and then skip
            # onto next line
            outputfile.write(line)
            continue

        if int(line[23:27]) == 1:
            # check to see if we're at the start of a model. This checks the residue index of a given model
            if ter_counter == 0:
                # if this is true, then we're at the very start of the file and I should just write to the outputfile
                outputfile.write(line)
                ter_counter += 1
                continue
            else:
                # If this is true, then we have models already in the file, and a termination is needed.
                outputfile.write("TER\n")
                outputfile.write(line)
                ter_counter += 1
                continue
        else:
            # in the case where the index is not 1, that means we're in the middle of the structure file! Output these
            # coordinates.
            outputfile.write(line)
        if line_split[0] == "END":
            # we've reached the end of the file, and by default the final TER won't be added. So, I have to add it
            # now
            outputfile.write("TER\n")
            outputfile.write(line)

outputfile.close()



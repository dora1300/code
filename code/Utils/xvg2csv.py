"""
@Title:         xvg2csv.py
@Author:        Mando Ramirez
@Date:          20220405

@Description: This script takes a XVG file from gromacs analysis and converts into a CSV file, which is handy.
One of the issues with XVG files is that the number of spaces that separates entries in the file are often not consistent,
which means that scripts to parse xvg files need to be crafted to deal with that. This is one of way of dealing with that.

Spaces between data in rows are replaced by a singular comma for the csv file, which is substantially easier to handle
than the spaces provided by XVG.
"""

import argparse

# Set up the argparser and then parse the arguments!
parser = argparse.ArgumentParser(description="Converts your .xvg files into .csv files!")
parser.add_argument("-f", help="Input .xvg file. (provide path if necessary)", required=True)
parser.add_argument("-o", help="Name of output .csv file.", default="output_conversion.csv")

args = parser.parse_args()

in_file = args.f
out_file = open(args.o, "w")

# Parse through the input XVG file and then write to the outfile as the input file is read
with open(args.f, "r") as f:
    for line in f:
        if line[0] == "#":
            out_file.write(line)
        elif line[0] == "@":
            out_file.write(line)
        else:
            sp = line.lstrip(" ").rstrip("\n").split(" ")
            for i in range(len(sp)):
                if len(sp[i]) == 0:
                    continue
                elif i != (len(sp)-1):
                    out_file.write(f"{sp[i]},")
                else:
                    out_file.write(f"{sp[i]}")
            out_file.write("\n")

out_file.close()
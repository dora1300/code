
import prody as pd
import argparse
import os.path

"""
Name of script: translateChains.py
Author: Dominique Ramirez
Date: 2021 07 30
Version: 1.0

This script the functionality to translate the coordinates of a PDB by a set linear vector for x, y, or z. This 
corresponds to a geometric translation.
"""

def translateCoords(input_file, x, y, z, chain_id):
    """

    :param input_file:
    :param x:
    :param y:
    :param z:
    :param chain_id:
    :return:
    """
    # PDB = preparsed PDB file from prody. This must happen outside the function
    PDB = pd.parsePDB(input_file)
    hv = PDB.getHierView()
    chains = list(hv)
    chain = chains[chain_id-1]
    for i in range(1, len(chain)+1):
        # this loops through the indices for the residues within the selected chain. Every loop will select a single
        # residue, one after the next
        res_coords = chain.getResidue(i).getCoords()
        for j in range(len(res_coords)):
            # this loops through the array of coordinates for each of the atoms in the selected residue
            res_coords[j][0] = res_coords[j][0] + x     # translate the X coords by value x
            res_coords[j][1] = res_coords[j][1] + y     # translate the Y coords by value y
            res_coords[j][2] = res_coords[j][2] + z     # translate the Z coords by value z
            # at the end of this loop, it will have changed the coords of all the atoms only of the given
            # residue in question
        chain.getResidue(i).setCoords(res_coords)     # save the new coords to the chain info

    new_prot = chain
    if len(chains) > 1:
        for i in range(len(chains)):
            if i == chain_id-1:
                continue
            else:
                new_prot += chains[i]
    pd.writePDB("translated_"+input_file, new_prot)
    return None

# -- Set up the argument parser
parser = argparse.ArgumentParser("PDB Utils - Translate Coords v1.0.\n"
                                 "Translates an entire chain by provided (x,y,z) scalars. Currently only translates one"
                                 " chain at a time. Default is the first chain")
parser.add_argument("pdbFile", action="store", type=str, help="PDB file to be parsed. Must include extension. Give full path if not in directory!")

parser.add_argument("-x", type=float, default=0., help="The amount to scale X by.")
parser.add_argument("-y", type=float, default=0., help="The amount to scale Y by.")
parser.add_argument("-z", type=float, default=0., help="The amount to scale Z by.")
parser.add_argument("-c", "--chainid", type=int, default=1, help="The chain to translate, given as the 1-based index.")

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

translateCoords(input_pdb, args.x, args.y, args.z, args.chainid)
print("PDB File: "+input_pdb+" was successfully translated by coordinates given. "
                             "Output file: "+"translated_"+input_pdb)
print()
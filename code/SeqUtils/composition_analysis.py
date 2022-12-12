"""
Name: composition_analysis.py
Date: 2021 08 26
Author: DA Ramirez

Description: This script contains several functions useful for composition analysis of protein sequences especially
as it relates to metrics useful for intrinsically disordered proteins and those that phase separate.

This uses the CIDER package provided by the Pappu lab. Citation:
Holehouse AS, Das RK, Ahad JN, Richardson MOG, and Pappu RV (2017) CIDER: Resources to Analyze Sequence-Ensemble Relationships of Intrinsically Disordered Proteins. Biophys J 112:16–21.

"""

import argparse
from localcider.sequenceParameters import SequenceParameters

# - Universal script constants/lists
nonpolar = ["A", "G", "I", "L", "P", "V", "M"]
aromatic = ["F", "W", "Y"]
polar_acidic = ["S", "C", "T", "N", "Q"]
polar_neutral = ["S", "C", "T", "N", "Q", "H"]
positive_acidic = ["R", "H", "K"]
positive_neutral = ["R", "K"]
negative = ["D", "E"]

# - Function definitions
def composition(sequence, acidic):
    total_count = 0
    seq_length = len(sequence)
    nonpolar_count = 0
    aromatic_count = 0
    positive_count = 0
    polar_count = 0
    negative_count = 0


    print("Non-Polar residues")
    for i in range(len(nonpolar)):
        print("--- Residue: " + nonpolar[i])
        res_count = sequence.count(nonpolar[i])
        total_count += res_count
        nonpolar_count += res_count
        print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))

    print()
    print("Aromatic residues")
    for i in range(len(aromatic)):
        print("--- Residue: " + aromatic[i])
        res_count = sequence.count(aromatic[i])
        total_count += res_count
        aromatic_count += res_count
        print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))

    print()
    if acidic == "acidic":
        print("Polar residues")
        for i in range(len(polar_acidic)):
            print("--- Residue: " + polar_acidic[i])
            res_count = sequence.count(polar_acidic[i])
            total_count += res_count
            polar_count += res_count
            print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))
    else:
        print("Polar residues")
        for i in range(len(polar_neutral)):
            print("--- Residue: " + polar_neutral[i])
            res_count = sequence.count(polar_neutral[i])
            total_count += res_count
            polar_count += res_count
            print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))

    print()
    if acidic == "acidic":
        print("Positively ionizable residues / Fraction Positive Charge")
        for i in range(len(positive_acidic)):
            print("--- Residue: " + positive_acidic[i])
            res_count = sequence.count(positive_acidic[i])
            total_count += res_count
            positive_count += res_count
            print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))
    else:
        print("Positively ionizable residues / Fraction Positive Charge")
        for i in range(len(positive_neutral)):
            print("--- Residue: " + positive_neutral[i])
            res_count = sequence.count(positive_neutral[i])
            total_count += res_count
            positive_count += res_count
            print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))

    print()
    print("Negatively ionizable residues / Fraction Negative Charged")
    for i in range(len(negative)):
        print("--- Residue: " + negative[i])
        res_count = sequence.count(negative[i])
        total_count += res_count
        negative_count += res_count
        print("Percentage of Sequence: " + str((round(res_count / seq_length * 100, 2))))

    print()
    print(f"Total amino acids counted: {total_count} / {seq_length}")
    print(f"Total nonpolar amino acids: {nonpolar_count}")
    print(f"Total polar amino acids: {polar_count}")
    print(f"Total aromatic amino acids: {aromatic_count}")
    print(f"Total positively charged amino acids: {positive_count}")
    print(f"Total negatively charged amino aicds: {negative_count}")
    print()

def FCR(sequence, acidity):
    """
    Calculates the Fraction of Charged Residues.
    :param sequence:
    :param acidity:
    :return:
    """
    seq_length = float(len(sequence))
    fcr_count = 0
    if acidity == "acidic":
        for i in range(len(positive_acidic)):
            fcr_count += sequence.count(positive_acidic[i])
    else:
        for i in range(len(positive_neutral)):
            fcr_count += sequence.count(positive_neutral[i])
    for i in range(len(negative)):
        fcr_count += sequence.count(negative[i])
    fcr = round(fcr_count / seq_length, 4)
    return fcr

def NCPR(sequence, acidity):
    """
    Calculates the Net Charge Per Residue
    :param sequence:
    :param acidity:
    :return:
    """
    seq_length = float(len(sequence))
    ncpr_positive = 0
    ncpr_negative = 0
    if acidity == "acidic":
        for i in range(len(positive_acidic)):
            ncpr_positive += sequence.count(positive_acidic[i])
    else:
        for i in range(len(positive_neutral)):
            ncpr_positive += sequence.count(positive_neutral[i])

    for i in range(len(negative)):
        ncpr_negative += sequence.count(negative[i])

    pos_frac = ncpr_positive / seq_length
    neg_frac = ncpr_negative / seq_length
    NCPR = round(pos_frac - neg_frac, 4)
    return NCPR

def meanNetCharge(sequence, acidity):
    seq_length = float(len(sequence))
    charged_res = 0

    if acidity == "acidic":
        for i in range(len(positive_acidic)):
            charged_res += sequence.count(positive_acidic[i])
    else:
        for i in range(len(positive_neutral)):
            charged_res += sequence.count(positive_neutral[i])

    for i in range(len(negative)):
        charged_res += sequence.count(negative[i])

    MNC = round(charged_res / seq_length, 4)
    return MNC




# - Main Script section, and the code that does stuff
parser = argparse.ArgumentParser(
    "Set of tools for compositional sequence analysis. - I am Ubik. Before the universe was, I am.")

parser.add_argument("sequence", type=str, help="Positional argument for the sequence file, in .fasta format. "
                                               "Include the path if not in the current directory. Only handles files"
                                               " with ONE sequence at a time.")
parser.add_argument('--list', nargs='+', type=int, help="A switch to define which analysis method to run. Input is a list"
                                              "and can include multiple analyses if provided correctly:\n"
                                             "1 - Compositional Analysis;\n"
                                             "2 - FCR/NCPR;\n"
                                             "3 - Pappu Kappa value;\n"
                                             "4 - SCD (Ghosh); \n"
                                             "5 - Uversky Hydropathy", required=True)
parser.add_argument("-pH", type=str, help="A string to convey whether the conditions are acidic (<6) or neutral "
                                              "(>6). Use either 'acidic' or 'neutral', no other strings. Default is"
                                          " neutral.", default="neutral")
parser.add_argument("-print", type=bool, help="Do you want to print the sequence to stdout? Type '-print True'!",
                    default=False)

args = parser.parse_args()

input_file = args.sequence
fxn_switches = args.list
fxn_switches.sort()
ph = args.pH
printer = args.print

# -- parse the input file and read the sequence, and check the validity of the acidic parameter
sequence = ""
try:
    with open(input_file, 'r') as file:
        for line in file:
            if line[0] == ">":
                seq_ID = line[1:]
            else:
                seq_split = line.rstrip("\n").split(" ")
                for i in range(len(seq_split)):
                    sequence += seq_split[i]
except FileNotFoundError:
    print("Sequence file could not be found using the given path and file name. Please try again.")
    exit(0)
if ph == "acidic" or ph == "neutral":
    pass
else:
    print("The pH argument passed in cannot be interpreted. Please try again")
    exit(0)

# -- prepare the sequence for analysis by making a Pappu class
SeqOb = SequenceParameters(sequence)

# -- the function switch board, where we decide which analyses to run
if printer:
    print("--- Sequence ---")
    print(sequence)
for switch in fxn_switches:
    if switch == 1:
        print()
        print("--- COMPOSITIONAL ANALYSIS ---")
        composition(sequence, ph)
    elif switch == 2:
        print()
        print("--- FCR and NCPR ---")
        pappu_fcr = SeqOb.get_FCR()
        pappu_ncpr = SeqOb.get_NCPR()
        frac_charge = FCR(sequence, ph)
        ncpr_frac = NCPR(sequence, ph)
        print(f"My FCR: {frac_charge}")
        print(f"Pappu FCR: {round(pappu_fcr,4)}")
        print(f"NCPR: {ncpr_frac}")
        print(f"Pappu NCPR: {round(pappu_ncpr,4)}")
        if frac_charge < 0.300:
            print("Sequence is a weak polyampholyte (FCR < 0.300)")
        else:
            print("Sequence is a strong polyampholyte (FCR >= 0.300)")
    elif switch == 3:
        print()
        print("--- KAPPA VALUE ---")
        print(round(SeqOb.get_kappa(),4))
    elif switch == 4:
        print()
        print("--- Ghosh SCD ---")
        print(round(SeqOb.get_SCD(),4))
    elif switch == 5:
        print()
        print("--- UVERSKY HYDROPATHY SCORE ---")
        print(round(SeqOb.get_uversky_hydropathy(),4))
        print(f"Mean Net Charge: {meanNetCharge(sequence, ph)}")
    else:
        print(f"Invalid method selection - {switch} - skipping input.")


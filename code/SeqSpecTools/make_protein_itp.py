"""
@Title:             make_protein_itp.py
@Name:              Mando A Ramirez
@Date:              2025 02 03

@Description:       This script generates the BIG KAHUNA -- the ITP file for the individual
protein sequence that I will simulate using the sequence specific interaction terms.

Though this itp shouldn't change much in between interaction parameters, there are many places
for intervention in this ITP should I want to change strength of force constants and assignment
of secondary structure.

AND THAT -- assignment of seconday structure -- is the crux of this code. I have to know something
about which parts of a protein are likely to be helical/coiled-coil and which are not, because I have
to know how to assign the torsional parameters.

This script DOES NOT make predictions about secondary structure. Do that somewhere else. This 
script just takes some prediction and translates that into a GROMACS-ready .itp file.
"""

import argparse


""" CONSTANTS """
# The parameters listed below are taken from the CC-LLPS simulation framework [Ramirez et al. 2024].
DAMP = 0.01                         # From DAR3-25p, used for damping linker parameters

# MASS = 109.0              # this isn't used. Use amino acid specific masses!
BOND_LENGTH = 0.385
BOND_FORCE = 50000.000

COIL_P14_C6 = 0.069130          # 1-4 pairs C6 (attractive) term
COIL_P14_C12 = 5.97404E-04      # 1-4 pairs C12 (repulsive) term
COIL_P15_C6 = 0.249383          # 1-5 pairs C6 (attractive) term
COIL_P15_C12 = 7.773996E-03     # 1-5 pairs C6 (attractive) term
LINK_P14_C6 = 0.069130 * DAMP          # 1-4 pairs C6 (attractive) term FOR LINKERS
LINK_P14_C12 = 5.97404E-04 * DAMP      # 1-4 pairs C12 (repulsive) term FOR LINKERS
LINK_P15_C6 = 0.249383 * DAMP          # 1-5 pairs C6 (attractive) term FOR LINKERS
LINK_P15_C12 = 7.773996E-03 * DAMP     # 1-5 pairs C6 (attractive) term FOR LINKERS

THETA_HELIX = 93
THETA_HELIX_FORCE = 50.0
ALPHA1_HELIX = -130.0
ALPHA1_HELIX_FORCE = 20
ALPHA3_HELIX = -30.0
ALPHA3_HELIX_FORCE = 20

THETA_LINK_FORCE = 50 * DAMP
ALPHA1_LINK_FORCE = 20 * DAMP
ALPHA3_LINK_FORCE = 20 * DAMP


aa_mass_dict = {
    "Ala": 71.03711,
    "Arg": 156.10111,
    "Asn": 114.04293,
    "Asp": 115.02694,
    "Cys": 103.00919,
    "Glu": 129.04259,
    "Gly": 57.02146,
    "His": 137.05891,
    "Ile": 113.08406,
    "Leu": 113.08406,
    "Lys": 128.09496,
    "Met": 131.04049,
    "Phe": 147.06841,
    "Pro": 97.05276,
    "Ser": 87.03203,
    "Thr": 101.04768,
    "Trp": 186.07941,
    "Tyr": 163.06333,
    "Val": 99.06841
}

aas_list = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys',
            'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro',
            'Ser', 'Thr', 'Trp', 'Tyr', 'Val']

aas_one_letter = ['A', 'R', 'N', 'D', 'C',
                  'Q', 'E', 'G', 'H', 'I',
                  'L', 'K', 'M', 'F', 'P',
                  'S', 'T', 'W', 'Y', 'V']




"""     FUNCTIONS       """

def parse_deepcoil_prediction():
    return None

def parse_s4pred_prediction():
    return None

def parse_psipred_prediction():
    return None
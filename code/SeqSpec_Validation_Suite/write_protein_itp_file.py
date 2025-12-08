"""
@Title:             write_protein_itp_file.py
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
import importonly_itp_helper_functions as helper


""" CONSTANTS  and  DATA CONTAINERS """
# The parameters listed below are taken from the CC-LLPS simulation framework [Ramirez et al. 2024].
DAMP = 0.01                         # From DAR3-25p, used for damping linker parameters

# MASS = 109.0              # this isn't used. Use amino acid specific masses!
BOND_LENGTH = 0.385
BOND_FORCE = 50000.000

COIL_P14_C6 = 0.069130          # 1-4 pairs C6 (attractive) term
COIL_P14_C12 = 5.97404E-04      # 1-4 pairs C12 (repulsive) term
COIL_P15_C6 = 0.249383          # 1-5 pairs C6 (attractive) term
COIL_P15_C12 = 7.773996E-03     # 1-5 pairs C6 (attractive) term
# LINK_P14_C6 = 0.069130 * DAMP          # 1-4 pairs C6 (attractive) term FOR LINKERS
# LINK_P14_C12 = 5.97404E-04 * DAMP      # 1-4 pairs C12 (repulsive) term FOR LINKERS
# LINK_P15_C6 = 0.249383 * DAMP          # 1-5 pairs C6 (attractive) term FOR LINKERS
# LINK_P15_C12 = 7.773996E-03 * DAMP     # 1-5 pairs C6 (attractive) term FOR LINKERS

THETA_HELIX = 93                # theta corresponds to the pseudo-angle between three c-alpha CG beads
THETA_HELIX_FORCE = 50.0
ALPHA1_HELIX = -130.0           # alpha corresponds to the pseudo-torsion between four c-alpha CG beads
ALPHA1_HELIX_FORCE = 20
ALPHA3_HELIX = -30.0
ALPHA3_HELIX_FORCE = 20

THETA_LINK_FORCE = 50 * DAMP
ALPHA1_LINK_FORCE = 20 * DAMP
ALPHA3_LINK_FORCE = 20 * DAMP


# Masses taken from [https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html]
# units in amu
aa_mass_dict = {
    "Ala": 71.03711,
    "Arg": 156.10111,
    "Asn": 114.04293,
    "Asp": 115.02694,
    "Cys": 103.00919,
    "Glu": 129.04259,
    "Gln": 128.05858,
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

# Important! The order of the 'aas_list' is identical to the order
# of 'aas_one_letter'! For easy finding of amino acids!
aas_list = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys',
            'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro',
            'Ser', 'Thr', 'Trp', 'Tyr', 'Val']

aas_one_letter = ['A', 'R', 'N', 'D', 'C',
                  'Q', 'E', 'G', 'H', 'I',
                  'L', 'K', 'M', 'F', 'P',
                  'S', 'T', 'W', 'Y', 'V']



"""  FUNCTIONS FOR WRITING THE ITP  """
"""
I am choosing to put the writing into a function instead of having it free-flowing
in this script. I suppose that makes this more pythonic

This is the order that the .itp file will take
        * openening text *
        [ moleculetype ]
        [ atoms ]
        [ bonds ]
        [ pairs ]
        [ exclusions ]
        [ angles ]
        [ dihedrals ]

"""

def fasta_parser(FASTAFILE):
    """
    Description:    Parses the given protein sequence and returns the sequence in string format
    Arguments:
        :FASTAFILE:
    Returns:        A string of amino acids (one letter code)
    """
    with open(FASTAFILE, 'r') as fasta:
        for line in fasta:
            if line[0] == ">":
                sequence_name = line[1:].rstrip("\n")
            elif line[0] == "#":
                # I'm allowing comments in my fasta file, which is not standard format
                # but I might want to make notes
                continue
            elif line[0] != ">" and line[0] != "#":
                sequence = line.rstrip("\n")
            else:
                break

        return sequence_name, sequence
    

def backbone_file_parser(BACKBONEFILE):
    """
    Description:
        Parses the file that contains the MAPPED per residue backbone scalar values
        for the protein sequence in question
    Arguments:
        BACKBONEFILE   :   the file to parse [.csv]
    Returns:
        a list that contains the backbone scalar values, for each residue
    """
    backbone_scalar_values = []
    with open(BACKBONEFILE, 'r') as backfile:
        for line in backfile:
            split_line = line.rstrip("\n").split(",")
            if len(split_line) < 3:
                # this if statement captures what might happen at the end of the file, if there's a blank
                # new line (which there certainly will be)
                # that constitutes the end of the file and then it's time to return
                break
            else:
                try:
                    backbone_scalar_values.append(float(split_line[2]))
                except:
                    backbone_scalar_values.append(float('nan'))
    
    return backbone_scalar_values


def write_itp_text(protein_name,
                   protein_sequence,
                   molecule_name,
                   backbone_scalars,
                   itp_filename,
                   nrexcl_val,
                   list_of_pairs):
    """
    Description:
    Arguments:
        protein_name
        protein_sequence
        molecule_name
        backbone_scalars
        itp_filename
        value_of_nrexcl
    Returns:        Nothing. This writes all the generated text to a file.
    """

    # First things first -- get the length of the sequence.
    seq_length = len(protein_sequence)


    # Second things second
    opener_text = f""";
; Title: {itp_filename}
; Individual topology file (.itp) for a coiled-coil protein
;  (named {protein_name}, which was included in the .fasta)
;  made for simulation with a Sequence Specific Force Field.
; This file generated by make_protein_itp.py.
; Verify the name of this file is properly included in the overall .top
;  before starting a simulation.
;
"""

    moleculetype_text = f"""
    
[ moleculetype ]
;name          nrexcl
{molecule_name}     {nrexcl_val}
"""

    atoms_text = helper.generate_atom_text(protein_sequence, 
                                           aas_one_letter,
                                           aas_list,
                                           aa_mass_dict)


    bonds_text = helper.generate_bond_text(seq_length,
                                           BOND_LENGTH,
                                           BOND_FORCE)


    pairs_text = f"""

[ pairs ]
; the value of func, V and W are not included because they are read from pairtypes.itp
;i    j    func    V (attrac)    W (repul)
"""
    if list_of_pairs is None:
        pairs_text += f"; no pairs added\n"
    else:
        for pair_num in list_of_pairs:
            pairs_text += f";    1-{pair_num} pairs\n"
            for pair_start in range(1, (seq_length-(pair_num-1))+1):
                # the above range seems complicated. Here's why it's this way:
                # seq_length for the length of the sequence. +1 because end of range is exclusive.
                # -(pair_num-1) because then this puts the correct upper bound on the sequence indices that I
                # can supply to the pairs section. e.g. if pair_num = 4, then the final pair that I can have
                # is (seq_len-4, seq_len)
                pairs_text += f"{pair_start} {pair_start+(pair_num-1)}\n"
    pairs_text += f"\n"

# an update as of 2025 12 06 -- Removing the generation of og_pair types because 
# these are not sequence specific and aren't relevant. However, I'm keeping it commented out
# just for posterity.
#     if generate_og_pairs:
#         pairs_text = f"""
#
# [ pairs ]
# ;OG 1-4 & 1-5 pairs from the CC-LLPS simulation framework
# ;i    j    func    C6(attrac)    C12(repul)
# """
#         pairs_14_text = helper.generate_pairs_text(seq_length,
#                                                 "4",
#                                                 backbone_scalars,
#                                                 COIL_P14_C6,
#                                                 COIL_P14_C12,
#                                                 combining_method="mean")
#         pairs_15_text = helper.generate_pairs_text(seq_length,
#                                                 "5",
#                                                 backbone_scalars,
#                                                 COIL_P15_C6,
#                                                 COIL_P15_C12,
#                                                 combining_method="mean")
#         pairs_text += pairs_14_text
#         pairs_text += pairs_15_text



    exclusions_text = f"""

[ exclusions ]
; no exclusions for a sequence specific protein!
"""
   

    angles_text = f"""

[ angles ]
;i    j    k    func   theta0 (deg)    Ktheta (kJ/mol/rad^2)
"""
    
    angles_text += helper.generate_angles_text(seq_length, THETA_HELIX, THETA_HELIX_FORCE,
                                               backbone_scalars, combining_method="mean")



    dihedrals_text = f"""

[ dihedrals ]
;i   j   k   l   func   phi0(deg)   kb (kJ/mol)  mult.
"""
    dihedrals_text += helper.generate_torsions_text(seq_length, ALPHA1_HELIX, ALPHA1_HELIX_FORCE, 1, 
                                                    backbone_scalars, combining_method="mean")

    dihedrals_text += helper.generate_torsions_text(seq_length, ALPHA3_HELIX, ALPHA3_HELIX_FORCE, 3, 
                                                    backbone_scalars, combining_method="mean")


    # Now that I've made everything, write that bitch out!
    writer_file = open(itp_filename, 'w')
    writer_file.write(opener_text)
    writer_file.write(moleculetype_text)
    writer_file.write(atoms_text)
    writer_file.write(bonds_text)
    writer_file.write(pairs_text)
    writer_file.write(exclusions_text)
    writer_file.write(angles_text)
    writer_file.write(dihedrals_text)
    writer_file.close()
    return None












if __name__ == '__main__':
    # run this code if I'm running this script by itself!
    """
        Argument Parser Set up
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-fastafile', help="File name (and path, if appropriate) of the protein sequence you wish to "
                        "generate topology files for, in FASTA format!", type=str, required=True)
    parser.add_argument('-molecule_type_name', help="The name to give your molecule! i.e. [ moleculetype ]",
                        type=str, default='coil-model')
    parser.add_argument('-backbone_file', help="The file containing the **mapped** backbone scalar values on a per residue "
                        "basis, for the protein sequence in question. This file must have 3 columns to be valid!",
                        type=str, required=True)
    parser.add_argument('-output_file', help="Name that you wish to give to the .itp file that you've just made!"
                        " Include extension!", type=str)

    args = parser.parse_args()
    FASTA = args.fastafile
    MOLTYPE = args.molecule_type_name
    BACKBONE = args.backbone_file
    OUTNAME = args.output_file



    """
        PROTEIN ITP MAKER
            This is where the itp file will be assembled!!
    """
    # it would behoove me to run some tests before I get too carried away with things, to make sure I have my
    # files set up as they should be

    # Verification -- make sure the sequence I pull from .fasta has the same number of amino acids as numbers
    # that I extract from the backbone file
    input_protein_name, input_protein_sequence = fasta_parser(FASTA)
    input_protein_backbones = backbone_file_parser(BACKBONE)
    if len(input_protein_sequence) != len(input_protein_backbones):
        raise ValueError(f"The length of the sequence you gave in the .fasta file ({len(input_protein_sequence)} aa) " \
                         f"does not match the number of mapped backbone scalar values ({len(input_protein_backbones)}) " \
                         "given in the backbone file. Please check your files and try again.")


    # Now that I've done my validation, I can go about writing my itp file!
    write_itp_text(input_protein_name, input_protein_sequence, MOLTYPE, input_protein_backbones, OUTNAME)

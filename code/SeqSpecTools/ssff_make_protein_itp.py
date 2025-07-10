"""
@Title:             ssff_make_protein_itp.py
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
import import_itp_helper_functions as helper


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
    with open(FASTAFILE) as fasta:
        line_counter = 1
        for line in fasta:
            if line_counter == 1:
                sequence_name = line[1:].rstrip("\n")
            elif line_counter == 2:
                sequence = line.rstrip("\n")
            else:
                break
            line_counter += 1

        return sequence_name, sequence

def write_itp_text(fasta_file, 
                   molecule_name,
                   structure_tool, 
                   structure_data, 
                   itp_filename, 
                   use_secondary_structure_info, 
                   force_adjuster,
                   index_end,
                   index_start=1):
    """
    Description:
    Arguments:
    Returns:        Nothing. This writes all the generated text to a file.
    """

    # First things first -- get the sequence and name of the protein!
    # one letter code sequence only
    protein_name, protein_sequence = fasta_parser(fasta_file)
    seq_length = len(protein_sequence)


    # Second things second
    opener_text = f""";
; Title: {itp_filename}
; Individual topology file (.itp) for a coiled-coil protein
; made for simulation with a Sequence Specific Force Field.
; This file generated by make_protein_itp.py
; This file MUST BE MANUALLY INCLUDED IN THE OVERALL .top FILE!
;
"""

    moleculetype_text = f"""
[ moleculetype ]
;name          nrexcl
{molecule_name}     4
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
;i    j    func    C6(attrac)    C12(repul)
"""
    pairs_14_text = helper.generate_pairs_text(seq_length,
                                               "4",
                                               COIL_P14_C6,
                                               COIL_P14_C12)
    pairs_15_text = helper.generate_pairs_text(seq_length,
                                               "5",
                                               COIL_P15_C6,
                                               COIL_P15_C12)
    pairs_text += pairs_14_text
    pairs_text += pairs_15_text


    exclusions_text = f"""
[ exclusions ]
"""
    
    # at this point I need the information from secondary structure prediction tools
    if structure_tool == "deepcoil":
        pseudo_angles, pseudo_torsions = helper.parse_deepcoil_prediction(structure_data)
    
    if use_secondary_structure_info:
        angles_text = helper.generate_angles_text(pseudo_angles, THETA_HELIX, THETA_HELIX_FORCE, predictions=True)
    else:
        THETA_MOD_FORCE = THETA_HELIX_FORCE * force_adjuster
        pseudo_angles = helper.basic_angle_list_generator(index_start, index_end)
        pseudo_torsions = helper.basic_torsion_list_generator(index_start, index_end)
            # yes I know I'm only handling angles in this section but I need to define the torsions now in this
            # else block so that I can properly use the torsions in the next section
        angles_text = helper.generate_angles_text(pseudo_angles, THETA_HELIX, THETA_MOD_FORCE, predictions=False)


    dihedrals_text = f"""
[ dihedrals ]
;i   j   k   l   func   phi0(deg)   kb (kJ/mol)  mult.
"""
    if use_secondary_structure_info:
        torsion_mult1 = helper.generate_torsions_text(pseudo_torsions, ALPHA1_HELIX, ALPHA1_HELIX_FORCE, 1, predictions=True)
        torsion_mult3 = helper.generate_torsions_text(pseudo_torsions, ALPHA3_HELIX, ALPHA3_HELIX_FORCE, 3, predictions=True)
    else:
        ALPHA1_MOD_FORCE = ALPHA1_HELIX_FORCE * force_adjuster
        ALPHA3_MOD_FORCE = ALPHA3_HELIX_FORCE * force_adjuster
        torsion_mult1 = helper.generate_torsions_text(pseudo_torsions, ALPHA1_HELIX, ALPHA1_MOD_FORCE, 1, predictions=False)
        torsion_mult3 = helper.generate_torsions_text(pseudo_torsions, ALPHA3_HELIX, ALPHA3_MOD_FORCE, 3, predictions=False)
    dihedrals_text += torsion_mult1
    dihedrals_text += torsion_mult3


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
                        "generate topology files for, in FASTA format!", type=str)
    parser.add_argument('-molecule_type_name', help="The name to give your molecule! i.e. [ moleculetype ]",
                        type=str, default='coil-model')
    parser.add_argument('-no_structure_prediction', help="A switch to turn off the feature that scales force constants "
                        "by the secondary structure prediction mechanisms. Passing this flag sets the variable to FALSE. "
                        "Default value is True.",
                        action="store_false", default=True)
    parser.add_argument('-number_of_protein_beads', help="The number of total CG beads/atoms that will be in the "
                        "protein represented by the .itp file. This is necessary if you pass no_structure_prediction. "
                        "I'm making it required so that I never forget it.",
                        required=True, type=int)
    parser.add_argument('-structure_prediction', help="Please provide the name of the tool used to make secondary structure "
                        "predictions. Choose out of these options: 'deepcoil', 'psipred', 's4pred'", type=str)
    parser.add_argument('-sp_data', help="File name of domain/secondary structure prediction output. Current"
                        " supported tools include DeepCoil. Please provide prediction tool name in the "
                        "other flag. Also, include full path if file is not in directory.", type=str)
    parser.add_argument('-output_file', help="Name that you wish to give to the .itp file that you've just made!"
                        " Include extension!", type=str)
    parser.add_argument("-force_scale_factor", help="[Default None] A factor by which to universally scale the force constant "
                        "values for ANGLES and TORSIONS **ONLY**. No other force constants are changed. "
                        "The default/reference force constants for angles and torsions are the same as used for the default "
                        "coil segments in the Ramirez et al 2024 BPJ paper. "
                        "Each force constant is scaled multiplicatively by the scale factor, so please "
                        "provide fractional scale factors. e.g. 1 = no change, 0.90 = 90%% of regular "
                        "FC, 1.25 means 25%% stronger than reference FC etc.", default=None,
                        type=float)

    args = parser.parse_args()
    FASTA = args.fastafile
    MOLTYPE = args.molecule_type_name
    doIUsePrediction = args.no_structure_prediction
    NUMBEADS = args.number_of_protein_beads
    PREDICTION = args.structure_prediction
    STRFILE = args.sp_data
    OUTNAME = args.output_file
    FSCALAR = args.force_scale_factor



    """
        PROTEIN ITP MAKER
            This is where the itp file will be assembled!!
    """
    write_itp_text(FASTA, MOLTYPE, PREDICTION, STRFILE, OUTNAME, doIUsePrediction, FSCALAR, NUMBEADS)

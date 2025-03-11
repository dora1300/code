"""
@Title:             import_itp_helper_functions.py
@Name:              Mando A Ramirez
@Date:              2025 03 05

@Description:       This code holds the helper scripts for the make_protein_itp.py file.
Basically all the functions that are used by make_protein_itp.py live in here. The reason is for clarity.
"""
import numpy as np

def generate_atom_text(PROT_SEQUENCE, AAS_ONE_LETTER, AAS_LIST, AA_MASS):
    """
    Description:
        [ atoms ]
        ;nr    type    resnr    residue    atom    cgnr    charge    mass
    Arguments:
    Returns:        [string] all the text for the [ atoms ] directive that will go into the 
    itp file.
    """

    output_text = f"""
[ atoms ]
; to keep everything consistent, 'type' and 'residue' are identical
; though 'type' is what matters when referencing the nonbonded.itp file.
;nr  type  resnr residue atom  cgnr  charge   mass
"""
    
    for aa_posi, aa in enumerate(PROT_SEQUENCE):
        AA = AAS_LIST[AAS_ONE_LETTER.index(aa)] #convert one letter to 3 letter
        line_text = f" {aa_posi+1}   {AA}   {aa_posi+1}     {AA}     CA     {aa_posi+1}    0.000    {AA_MASS[AA]}\n"
        output_text += line_text
    
    return output_text


def generate_bond_text(sequence_length, bond_length, bond_force):
    """
    Description:
        [ bonds ]
        ;i    j    func   r0 (nm)        kb (kJ/mol/nm^2)
    Arguments:
    Returns:       
    """

    output_text = f"""
[ bonds ]
;i    j    func  r0 (nm)   kb (kJ/mol/nm^2)
"""

    for aa_posi in range(1, sequence_length):
        line_text = f" {aa_posi}    {aa_posi+1}    1     {bond_length}     {bond_force}\n"
        output_text += line_text

    return output_text


def generate_pairs_text(sequence_length, pair_type, pair_attractive, pair_repulsive):
    """
    Description:
        [ pairs ]
        ;i    j    func    C6(attrac)    C12(repul)
    Arguments:
    Returns:       
    """
    pair_type_int = int(pair_type)
    increment = pair_type_int - 1

    output_text =f""";       1-{pair_type_int} pairs
"""
    for i in range(1, sequence_length-(increment-1)):
        line_text = f" {i}    {i+increment}     1      {pair_attractive}       {pair_repulsive}\n"
        output_text += line_text

    return output_text


def generate_angles_text(pseudoangles, equil_value, angle_force):
    """
    Description:
        [ angles ]
        ;i    j    k    func   theta0 (deg)    Ktheta (kJ/mol/rad^2)
    Arguments:
    Returns:       
    """
    
    output_text = f"""
[ angles ]
;i    j    k    func   theta0 (deg)    Ktheta (kJ/mol/rad^2)
"""
    for angle_set in pseudoangles:
        force_adjustment = angle_force * angle_set[3]
        line_text = f" {angle_set[0]}    {angle_set[1]}    {angle_set[2]}    1      {equil_value}           {force_adjustment:.4f}\n"
        output_text += line_text

    return output_text


def generate_torsions_text(pseudotorsions, torsion_value, torsion_force, multiplicity):
    """
    Description:
        [ dihedrals ]
        ;i   j   k   l   func   phi0(deg)   kb (kJ/mol)  mult.
    Arguments:
    Returns:       
    """

    output_text = f""";     multiplicity = {multiplicity}
"""
    
    for torsion_set in pseudotorsions:
        force_adjustment = torsion_force * torsion_set[4]
        line_text = f" {torsion_set[0]}   {torsion_set[1]}   {torsion_set[2]}   {torsion_set[3]}    1     {torsion_value}      {force_adjustment:.4f}        {multiplicity}\n"
        output_text += line_text

    return output_text










"""
    Functions that parse structure prediction output files and return 
    angles and torsions with their associated window averaged score values.
"""
def parse_deepcoil_prediction(FILENAME):
    """
    Description:

        Important! A score of 0.000 doesn't make sense (I don't want torsion force constants to be
        0, just weak!)! But at the same time, a score of 0.005 isn't that different from 0.010 if
        I'm using the score to scale the force constant.
        So I'm converting every score <0.01 to 0.010.
        This corresponds to using a dampening of 0.01 for linkers in my old framework, i.e. this 
        corresponds to a primarily disordered state.

        Update 2025 03 11 --> all scores are now getting "binned" to their nearest tenth place.
        All values < 0.05 become 0.01. Values >= 0.05 and < 0.1 turn into 0.1.
        This should avoid having values like 0.25 which seem to be causing numerical stability issues.
        This code doesn't actually get updated, but I'm putting a comment here because this function
        *IS* affected.
        
    Arguments:
    Returns:
        :pseudoangles:      [list]
        :pseudotorsions:    [list]
        Important! This only returns the pseudoangles and the pseudotorsions!
    """
    deepcoil_data_file = np.loadtxt(FILENAME,
                                    delimiter="\t",
                                    skiprows=1, 
                                    dtype=str).T
    len_sequence = len(deepcoil_data_file[0])

    # Define the containers to hold the data I care about
    aa_sequence = np.zeros(len_sequence, dtype=str)
    deepcoil_scores = np.zeros(len_sequence, dtype=float)
    index_array = np.arange(1, len_sequence+1, step=1)
    
    # Loop through the file and sort amino acids and the raw_cc scores
    # So I'm converting all 0s to 0.010
    for INDEX, AA in enumerate(deepcoil_data_file[0]):
        aa_sequence[INDEX] = AA
        if float(deepcoil_data_file[2][INDEX]) < 0.010:
            deepcoil_scores[INDEX] = 0.001
        else:
            deepcoil_scores[INDEX] = float(deepcoil_data_file[2][INDEX])

    pseudoangles = deepcoil_make_pseudoangles(index_array, deepcoil_scores, len_sequence)
    pseudotorsions = deepcoil_make_pseudotorsions(index_array, deepcoil_scores, len_sequence)

    return pseudoangles, pseudotorsions


def deepcoil_make_pseudoangles(arr_of_aa_indicies, arr_of_dc_scores, num_of_aas):
    """
    Description:
        Update 2025 03 11 --> all scores are now getting "binned" to their nearest tenth place.
        All values < 0.05 become 0.01. Values >= 0.05 and < 0.1 turn into 0.1.
        This should avoid having values like 0.25 which seem to be causing numerical stability issues.
        This update takes place in the "average_dc_score(...)" line since that incorporates rounding.
    Arguments:
    Returns:
    """
    list_of_pseudoangles = []

    for I in range(0, num_of_aas-2):
        average_dc_score = np.round(np.average([arr_of_dc_scores[I],
                                       arr_of_dc_scores[I+1],
                                       arr_of_dc_scores[I+2]]),
                                       1)
        if average_dc_score == 0.0:
            average_dc_score = 0.01     # this is to be consistent with the disordered linker state
        list_of_pseudoangles.append([arr_of_aa_indicies[I],
                                     arr_of_aa_indicies[I+1],
                                     arr_of_aa_indicies[I+2],
                                     average_dc_score])
            
    return list_of_pseudoangles


def deepcoil_make_pseudotorsions(arr_of_aa_indicies, arr_of_dc_scores, num_of_aas):
    """
    Description:
        Update 2025 03 11 --> all scores are now getting "binned" to their nearest tenth place.
        All values < 0.05 become 0.01. Values >= 0.05 and < 0.1 turn into 0.1.
        This should avoid having values like 0.25 which seem to be causing numerical stability issues.
        This update takes place in the "average_dc_score(...)" line since that incorporates rounding.
    Arguments:

    Returns:        Returns an array of window averaged DeepCoil prediction scores
    for all the possible torsions indicated by indicies/num_of_aas. 
    """
    list_of_pseudotorsions = []

    for I in range(0, num_of_aas-3):
        average_dc_score = np.round(np.average([arr_of_dc_scores[I],
                                       arr_of_dc_scores[I+1],
                                       arr_of_dc_scores[I+2],
                                       arr_of_dc_scores[I+3]]),
                                       1)
        if average_dc_score == 0.0:
            average_dc_score = 0.1
        list_of_pseudotorsions.append([arr_of_aa_indicies[I],
                                     arr_of_aa_indicies[I+1],
                                     arr_of_aa_indicies[I+2],
                                     arr_of_aa_indicies[I+3],
                                     average_dc_score])
        
    return list_of_pseudotorsions


def scale_deepcoil_values():
    return None


def parse_s4pred_prediction():
    return None

def parse_psipred_prediction():
    return None
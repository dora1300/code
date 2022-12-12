"""
@Title:             topology_assistant_specific.py
@Name:              Mando A Ramirez
@Date:              2022 07 19

@Description:       This script writes the topology for any type of coil model I want, with coils and linkers specified,
to assist me in the making of files for simulation. It will have helper scripts to make general .top files and .itp
files for individual models. **IMPORTANT** == this code writes for a restricted interaction scheme, meaning that coils
can only interact with the same type of coil segment on other molecules. What I mean is if I have a 2 coil model, then
it is structured like:
                                    coil1 -- linker -- coil2
Coil1 cannot interact with Coil2 either intra- or inter-molecularly. Coil1 can only interact with coil1 on OTHER models.
I think this best represents what should be seen in real coil proteins.

This code goes back and forth between using 1 indexing and 0 indexing, which ever is most convenient for the application.
Be mindful of this! I will try to denote what indexing is being used.

This code builds off of the other topology assistants.

@Updates:

2022 10 14 - updated file name to reflect the specific interaction scheme supported

"""

import argparse

""" Constants """
# Since I have now found parameters to use in the model, I will make those parameters constant
DAMP = 0.01                         # From DAR3-25p, used for damping linker parameters

MASS = 109.0
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




""" Function Definitions """
def write_topology(top_filename, aa_c6, aa_c12, coilNumber):
    """
    This writes the .top file for a coil model.
    :param top_filename: the user provided filename for the topology
    :param aa_c6: the value of the attractive part of the A-A vdW interaction
    :param aa_c12: the value of the repulsive part of the A-A vdW interaction
    :param coilNumber: an int corresponding to the number of coil segments in the model
    :return: None (this just saves the file)
    """
    output = open(top_filename, 'w')
    top = f"""; Topology file for a {coilNumber} coil model
; .itp files must be MANUALLY input
; updated to include new VdW parameter format
; Also includes updated B-B repulsion interaction set at 2kJ/mol
; but only for the nonbond_param
;       RESTRICTED INTERACTION SCHEME

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               1               yes              1.0     0

[ atomtypes ]
;name  mass     charge   ptype  V(C6, attra)  W(C12, repul)\n"""

    for CN in range(coilNumber):
        at_line = f" A{CN+1}     109.0    0.000   A        6.9132E-03    5.97404E-05\n"
        top += at_line
    top += " B      109.0    0.000   A        0             5.97404E-05\n"

    top += """
[ nonbond_params ]
;i     j     func    V(C6, attra)    W(C12, repul)\n"""
    # section for attractiveness between similar coil segments
    for CN in range(coilNumber):
        np_line = f" A{CN+1}     A{CN+1}     1        {aa_c6}          {aa_c12}\n"
        top += np_line
    # section for the repulsion between un-similar coil segments
    for CN in range(coilNumber):
        next_coil = CN + 2
        while next_coil <= coilNumber:
            np_repul_coil = f" A{CN+1}     A{next_coil}     1        0            5.97404E-04\n"
            top += np_repul_coil
            next_coil += 1
    # section for the repulsion to all B beads
    for CN in range(coilNumber):
        np_repul_line = f" A{CN+1}     B      1        0            5.97404E-04\n"
        top += np_repul_line
    top += " B      B      1        0            5.97404E-04\n"

    top += """
; Include the individual .itps for different molecules
; This is set MANUALLY - DONT FORGET!
#include

[ system ]
; Name
Coil-model-simulation

[ molecules ]
;moleculetype-name    # molecules
; fill in manually here

"""

    output.write(top)
    output.close()
    return None

def write_itp(itp_filename, Length, segments_list, positions_list, coilNumber):
    """

    This is the order that the .itp file will take
        [ moleculetype ]
        [ atoms ]
        [ bonds ]
        [ pairs ]
        [ exclusions ]
        [ angles ]
        [ dihedrals ]
    :return: None
    """
    A_BEADS = []
    A_BEADS_BY_COIL = []

    opener = f""";
; Individual topology file (.itp) for a {coilNumber}-coil model
; Title: {itp_filename}
; WITH PAIRS -- parameters in C6 and C12 format!
; This file is generated automatically by "coil_topology_assistant.py"
; This file MUST BE MANUALLY INCLUDED IN THE OVERALL .top FILE!
;           RESTRICTED INTERACTION SCHEME

[ moleculetype ]
;name          nrexcl
coil_model     4

"""

    # next up do the atoms
    atom_information = f"""

[ atoms ]
;nr  type  resnr  residue  atom  cgnr  charge  mass
"""
    coil_counter = 0
    for index in range(len(segments_list)): # this is 0-indexed
        if segments_list[index] == "coil":
            coil_counter += 1
            temp_string, temp_a_posits, temp_d_posits = atom_helper(
                positions_list[index][0],
                positions_list[index][1],
                segment_type="coil",
                heptad_start=positions_list[index][2],
                coil_id = coil_counter
            )   # all these are 0-indexed too
            atom_information += temp_string
            coilN_beads = temp_a_posits; coilN_beads.extend(temp_d_posits)
            coilN_beads.sort()
            A_BEADS_BY_COIL.append(coilN_beads)
            A_BEADS.extend(temp_a_posits)           # all these positions are 1-indexed
            A_BEADS.extend(temp_d_posits)           # all these positions are 1-indexed
        else:
            temp_string = atom_helper(
                positions_list[index][0],
                positions_list[index][1]
            )   # all these are 0-indexed too
            atom_information += temp_string

    # now write the bond information
    bond_information = f"""

[ bonds ]
;i    j    func   r0 (nm)        kb (kJ/mol/nm^2)
"""
    for posb in range(1, Length):    # 1-indexed
        bond_information += f" {posb}    {posb+1}    1        {BOND_LENGTH}        {BOND_FORCE}\n"


    # Now write out ALL the pairs
    pair_information = f"""

[ pairs ]
;i   j   func     C6(attrac)    C12(repul)
;    1-4 pairs
"""
    for index1 in range(len(segments_list)): # this is 0-indexed
        if segments_list[index1] == "coil": # but the reference to positions_list values are 1-indexed
            temp_string = pairs_helper(
                positions_list[index1][0],
                positions_list[index1][1],
                segment_type="coil",
                pair_type="1-4"
            )
            pair_information += temp_string
        else:
            temp_string = pairs_helper(
                positions_list[index1][0],
                positions_list[index1][1],
                segment_type="linker",
                pair_type="1-4"
            )
            pair_information += temp_string

    pair_information += ";    1-5 pairs\n"

    for index2 in range(len(segments_list)): # this is 0-indexed
        if segments_list[index2] == "coil": # but the reference to positions_list values are 1-indexed
            temp_string = pairs_helper(
                positions_list[index2][0],
                positions_list[index2][1],
                segment_type="coil",
                pair_type="1-5"
            )
            pair_information += temp_string
        else:
            temp_string = pairs_helper(
                positions_list[index2][0],
                positions_list[index2][1],
                segment_type="linker",
                pair_type="1-5"
            )
            pair_information += temp_string


    # now handle the Exclusions
    exclusions_information = f"""

[ exclusions ]
; This coil interacts through a **SPECIFIC INTERACTION SCHEME**
; this model achieves specific interactions through using
; different sticky coil beads to keep the coil segments separated
; therefore, intra-model coil interactions will be removed through atom-type
; and nonbond_param specification
; Exclude only A-A interactions within the same individual coil segments

"""
    for coil_no, beads in enumerate(A_BEADS_BY_COIL):
        exclusions_information += f"; Exclusions for coil segment: {coil_no+1}\n"
        for ai in range(len(beads)):  # this is handled by 0-indexed
            exclusion_string = f""
            if ai != len(beads)-1:
                for tai in range(ai, len(beads)):
                    exclusion_string += f"{beads[tai]} "
                exclusion_string += "\n"
                exclusions_information += exclusion_string
            else:
                pass
        exclusions_information += "\n"


    # Now time to handle the angles
    angle_information = f"""
[ angles ]
;i    j    k    func   theta0 (deg)    Ktheta (kJ/mol/rad^2)
"""
    for i3 in range(len(segments_list)):    # this is 0-indexed
        if segments_list[i3] == "coil":     # the values in positions_list are 1-indexed
            coil_start = positions_list[i3][0]
            coil_stop = positions_list[i3][1]
            temp_string = angle_helper(
                coil_start,
                coil_stop,
                segment_type="coil"
            )
            angle_information += temp_string
            # now I need to handle the manual angles (2 need to be specified)
            if i3 == (len(segments_list)-1):
                continue
            else:
                angle_information += f" {coil_stop - 1}    {coil_stop}    {coil_stop + 1}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^'_CL\n"
                angle_information += f" {coil_stop}    {coil_stop + 1}    {coil_stop + 2}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^''_CL\n"
        else:
            linker_start = positions_list[i3][0]
            linker_stop = positions_list[i3][1]
            temp_string = angle_helper(
                linker_start,
                linker_stop,
                segment_type="linker"
            )
            angle_information += temp_string
            # now I need to handle the manual angles (2 need to be specified)
            if i3 == (len(segments_list)-1):    # stop if I'm at the end of the model obviously
                continue
            else:
                angle_information += f" {linker_stop - 1}    {linker_stop}    {linker_stop + 1}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^'_LC\n"
                angle_information += f" {linker_stop}    {linker_stop + 1}    {linker_stop + 2}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^''_LC\n"


    # now it's time to handle the dihedrals!
    dihedral_information = f"""

[ dihedrals ]
;i   j   k   l   func   phi0(deg)   kb (kJ/mol)    mult.
;    multiplicity = 1
"""
    for i4 in range(len(segments_list)):
        if segments_list[i4] == "coil":
            coil_d_start = positions_list[i4][0]
            coil_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                coil_d_start,
                coil_d_stop,
                1,
                segment_type="coil",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {coil_d_stop - 2}   {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^0_C\n"
                dihedral_information += f" {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^'_CL\n"
                dihedral_information += f" {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}   {coil_d_stop + 3}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^''_CL\n"
        else:
            link_d_start = positions_list[i4][0]
            link_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                link_d_start,
                link_d_stop,
                1,
                segment_type="linker",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {link_d_stop - 2}   {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^0_C\n"
                dihedral_information += f" {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^'_CL\n"
                dihedral_information += f" {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}   {link_d_stop + 3}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^''_CL\n"

    dihedral_information += f";    multiplicity = 3\n"

    for i4 in range(len(segments_list)):
        if segments_list[i4] == "coil":
            coil_d_start = positions_list[i4][0]
            coil_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                coil_d_start,
                coil_d_stop,
                3,
                segment_type="coil",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {coil_d_stop - 2}   {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^0_C\n"
                dihedral_information += f" {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^'_CL\n"
                dihedral_information += f" {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}   {coil_d_stop + 3}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^''_CL\n"
        else:
            link_d_start = positions_list[i4][0]
            link_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                link_d_start,
                link_d_stop,
                3,
                segment_type="linker",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {link_d_stop - 2}   {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^0_C\n"
                dihedral_information += f" {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^'_CL\n"
                dihedral_information += f" {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}   {link_d_stop + 3}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^''_CL\n"


    # NOW -- BIG FINALE -- write all the parts to a file
    with open(itp_filename, 'w') as f:
        f.write(opener)
        f.write(atom_information)
        f.write(bond_information)
        f.write(pair_information)
        f.write(exclusions_information)
        f.write(angle_information)
        f.write(dihedral_information)

    return None


def atom_helper(aStart, aStop, segment_type="linker", heptad_start=0, coil_id=0):  # arguments provided as 1-indexed
    """

    :param aStart: the start (1-indexed) index of the segment in question
    :param aStop: the end (1-indexed) index of the segment in question
    :param segment_type: the switch for either linker or coil segment, behavior changes based on segment
    :param heptad_start: the location of the FIRST A bead in the coil segment, only has meaning for segment_type="coil"
    :param coil_id: the identifying number to assign to each A bead, important for the restricted interaction scheme
    :return: a complete string of every atom in the segment provided, bounded inclusive by aStart and aStop.
    Also return the list a_posits and d_posits for the given segment, if a coil, which includes the indices of the A and
    D position A-beads of the segment (1-indexed)
    """
    atom_string = ""
    a_posits = []           # stores A position indices by 1-indexing
    d_posits = []           # stores D position indices by 1-indexing
    if segment_type == "linker":    # control for making linker information
        for posl in range(aStart, aStop+1):  # 1-indexing
            atom_string += f" {posl}   B    {posl}    GLY     CA    {posl}    0.000    {MASS}\n"
        return atom_string
    else:                           # control for making coil information
        a_posits = [heptad_start]
        while len(a_posits) < 200:
            new_pos = a_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                a_posits.append(new_pos)
            else:
                break
        d_posits = [heptad_start+3]
        while len(d_posits) < 200:
            new_pos = d_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                d_posits.append(new_pos)
            else:
                break

        for posc in range(aStart, aStop+1):        # uses 1-indexing
            if posc in a_posits:
                atom_string += f" {posc}   A{coil_id}    {posc}    ILE     CA    {posc}    0.000    {MASS}\n"
            elif posc in d_posits:
                atom_string += f" {posc}   A{coil_id}    {posc}    ILE     CA    {posc}    0.000    {MASS}\n"
            else:
                atom_string += f" {posc}   B     {posc}    ALA     CA    {posc}    0.000    {MASS}\n"

        return atom_string, a_posits, d_posits

def pairs_helper(pStart, pStop, segment_type="linker", pair_type="1-4"):
    """

    :param pStart:
    :param pStop:
    :param segment_type:
    :param pair_type:
    :return:
    """
    pairs_string = ""
    pairs_list = []
    # This part handles if the pairs are 1-4 pairs or not
    if pair_type == "1-4":
        for i in range(pStart, pStop-2):        # this is 1-indexed!
            pairs_list.append([i, i+3])
        if segment_type == "coil":
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {COIL_P14_C6}     {COIL_P14_C12}\n"
        else:
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {LINK_P14_C6:.10f}     {LINK_P14_C12:.10f}\n"
        return pairs_string
    elif pair_type == "1-5":
        for i in range(pStart, pStop-3):        # this is 1-indexed!
            pairs_list.append([i, i+4])
        if segment_type == "coil":
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {COIL_P15_C6}     {COIL_P15_C12}\n"
        else:
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {LINK_P15_C6:.10f}     {LINK_P15_C12:.10f}\n"
        return pairs_string
    else:
        print("The pair_type specifer doesn't match any expected pairs for the coil model. There must be a problem.")
        return None

def angle_helper(anStart, anStop, segment_type="linker"):
    """

    :param anStart:
    :param anStop:
    :param segment_type:
    :return:
    """
    angle_string = ""
    angle_list = []
    for i in range(anStart, anStop-1):
        angle_list.append([i, i+1, i+2])
    if segment_type == "coil":
        for j in range(len(angle_list)):
            angle_string += f" {angle_list[j][0]}    {angle_list[j][1]}    {angle_list[j][2]}     1       {THETA_HELIX}        {THETA_HELIX_FORCE}\n"
        return angle_string
    else:
        for j in range(len(angle_list)):
            angle_string += f" {angle_list[j][0]}    {angle_list[j][1]}    {angle_list[j][2]}     1       {THETA_HELIX}        {THETA_LINK_FORCE}\n"
        return angle_string

def dihedral_helper(dStart, dStop, multiplicity, segment_type="linker"):
    """

    :param dStart:
    :param dStop:
    :param multiplicity:
    :param segment_type:
    :return:
    """
    dihedral_string = ""
    dihedral_list = []
    for i in range(dStart, dStop-2):        # 1-indexed
        dihedral_list.append([i, i+1, i+2, i+3])
    if segment_type == "coil":
        if multiplicity == 1:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA1_HELIX}        {ALPHA1_HELIX_FORCE}      {multiplicity}\n"
            return dihedral_string
        if multiplicity == 3:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA3_HELIX}        {ALPHA3_HELIX_FORCE}      {multiplicity}\n"
            return dihedral_string
    else:   # linker segments
        if multiplicity == 1:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      {multiplicity}\n"
            return dihedral_string
        if multiplicity == 3:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      {multiplicity}\n"
            return dihedral_string

if __name__ == "__main__":
    # I don't anticipate needing to use the functions in this script elsewhere, so the
    # name==main style is a little overkill, but I'm including it here for good pythonic practice

    # Set up the argument parser!
    parser = argparse.ArgumentParser(description="Topology Writer Tool -- Use this to make topology files for coil "
                                                 "models, easy as pie! Cheap as Ubik, too!")
    parser.add_argument("-df", help="Definition file: this contains the information on how to build the model, with"
                                    " alternating coil and linker segments.")
    parser.add_argument("-n", help="The number of beads in the coil model.", type=int, required=True)
    parser.add_argument("-top", help="Activate to turn on .top file writer.", action="store_true")
    parser.add_argument("-top_filename", help="File name for the .top file. Please provide extension.",
                        default="top.top", type=str)
    parser.add_argument("-aaAttract", help="The C6 attractive value for A-A bead interactions. Value is handled "
                                          "as string", type=str)
    parser.add_argument("-aaRepul", help="The C12 repulsive value for A-A bead interactions. Value is handled "
                                          "as string", type=str)
    parser.add_argument("-itp", help="Activate to turn on .itp file writer.",
                        action="store_true")
    parser.add_argument("-itp_filename", help="File name for the .itp file. Please provide extension.",
                        default="itp.itp")

    # parse the arguments
    args = parser.parse_args()

    # set up the list of segments, in order, and their corresponding indices, in order
    # Segments contains the names of each of the segments as they appear in order in the provided file
    # Positions contains the start/stop (and heptad) indices for the segments that correspond to the same index in the
    # Segments list
    Segments = []; Positions = []
    Coil_Number = 0
    with open(args.df, "r") as file:
        for line in file:
            line_split = line.rstrip("\n").split(",")
            Segments.append(str(line_split[0]))
            if str(line_split[0]) == "coil":   # uses 0-indexing
                Positions.append([int(line_split[1]), int(line_split[2]), int(line_split[3])])
                Coil_Number += 1
            else:
                Positions.append([int(line_split[1]), int(line_split[2])])

    if args.top:
        if args.aaAttract is None or args.aaRepul is None:
            print("A-A interaction energies are not provided. Please do so. Exiting now.")
            exit(1)
        else:
            write_topology(args.top_filename, args.aaAttract, args.aaRepul, Coil_Number)
    if args.itp:
        write_itp(args.itp_filename, args.n, Segments, Positions, Coil_Number)
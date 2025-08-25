"""
@Title:             ssff_model_maker.py
@Name:              Mando A Ramirez
@Date:              2025 03 06

@Description:       This script is the main model maker for my coil models -- this produces the 3D atomistic and CA
coarse-grained starting structures for each model. This is essentially the same script I use for my first 
framework, but now repurposed for the SSFF framework.

@Updates:
2025 03 07      -   I am adding a new feature to convert and existing PDB structure into CA representation so I 
don't have to open a PDB structure, pull out the CAs, and then resave. Hopefully this will save me time in the
long run.

"""

import math as m
import argparse
import numpy as np
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import prody as pd


"""
Tozzini converter functions
"""
#GLOBAL ANGLE CONSTANTS -- provided in DEGREES but then converted into radians
tau = 111
omega = 180
y1 = 20.7
y2 = 14.7            # y2 should not be 20.7, causes the pseudoangles to be too low
t_r = tau * (m.pi / 180)
o_r = omega * (m.pi / 180)
y1_r = y1 * (m.pi / 180)
y2_r = y2 * (m.pi / 180)

def THETA(T, Y1, Y2, Phi1, Psi1):
    """
    The structure of this equation is complex. I will shorten it into this so that I can keep better
    track of it
        cos(theta) = cos(T)*A - B + sin(T)*C
    where,
        A = [cos(Y1)cos(Y2) - sin(Y1)sin(Y2)cos(Psi1)cos(Phi1)]
        B = sin(Y1)sin(Y2)sin(Psi1)sin(Phi1)
        c = [cos(Phi1)sin(Y1)cos(Y2) + cos(Psi1)cos(Y1)sin(Y2)]

    All variables are provided in RADIANS.
    The answer is returned in DEGREES
    """
    A = m.cos(Y1)*m.cos(Y2) - m.sin(Y1)*m.sin(Y2)*m.cos(Psi1)*m.cos(Phi1)
    B = m.sin(Y1)*m.sin(Y2)*m.sin(Psi1)*m.sin(Phi1)
    C = m.cos(Phi1)*m.sin(Y1)*m.cos(Y2) + m.cos(Psi1)*m.cos(Y1)*m.sin(Y2)

    cos_theta = m.cos(T)*A - B + m.sin(T)*C
    theta = m.acos(cos_theta) * (180 / m.pi)
    return theta


def ALPHA(Y1, Y2, Phi1, Psi1, Phi2, Psi2):
    """
    All variables are provided as RADIANS. It is lastly converted into degrees after the computation
    step.
    """
    a = (m.pi + Psi1 + Phi2 + Y1*m.sin(Psi2) + Y2*m.sin(Phi1)) * (180 / m.pi)
    return a

def ALPHA_COMPLEX(T, Y1, Y2, Phi, Psi):
    A = 0.25 * Y1**2 * m.sin(2*Psi)
    B = 0.25 * Y2**2 * m.sin(2*Phi)
    C = Y1 * Y2 * m.sin(Phi + Psi)
    D = Y1 * (T - (m.pi / 2)) * m.sin(Psi)
    E = Y2 * (T - (m.pi / 2)) * m.sin(Phi)

    ac_r = Phi + Psi + m.pi + (Y1 * m.sin(Psi)) + (Y2 * m.sin(Phi)) + A + B + C - D - E
    ac = ac_r * (180 / m.pi)    # answer is returned in DEGREES

    if ac > 180:
        ac = ac - 360
    else:
        pass
    return ac


"""
Peptide Builder Functions
"""
def fasta_parser(FASTA_FILE):
    with open(FASTA_FILE, 'r') as sequence_file:
        line_counter = 1
        for line in sequence_file:
            if line_counter == 1:
                pass
            elif line_counter == 2:
                protein_sequence = line.rstrip("\n")
            else:
                break
            line_counter += 1
    return protein_sequence

def build_peptide(PeptideName, Sequence):
    """
    This function makes the starting structure of a given coil model. It only operates in a single directory,
    with no functionality to specify different directory output, so be careful where you use it.
    :param PeptideName: the name given to the built model, user input
    :param Sequence: the sequence that corresponds to the model to be built. This is parsed by the main code.
    :return: None -- all the necessary files are handled and saved. Nothing is returned.
    """
    aa_name = f"{PeptideName}.pdb"
    ca_name = f"{PeptideName}_ca.pdb"
    output_file = f"{PeptideName}_aa_dihedrals.csv"

    # Start by setting up the structure and defining the phi/psi angles beyond the defaults. Leave the defaults for the
    # bond angles and bond lengths
    geo = Geometry.geometry(Sequence[0])
    # These angles were from DAR3-15
    geo.phi = -64.088
    geo.psi_im1 = -41.260
    linker = PeptideBuilder.initialize_res(geo)
    for res in Sequence[1:]:
        linker = PeptideBuilder.add_residue(linker, res, geo.phi, geo.psi_im1)

    out_linker = Bio.PDB.PDBIO()
    out_linker.set_structure(linker)
    out_linker.save(aa_name)

    # Now this part of the code will reopen this fully atomistic structure, read the torsion angles of it and convert into
    # Tozzini space, then save the C-alpha coarse grained model. It will also save the torsions and the pseudotorsions into
    # and output file
    linker_pdb = pd.parsePDB(aa_name)
    hv = pd.HierView(linker_pdb)
    outfile = open(output_file, 'w')
    outfile.write("Residue No.,Residue Name,phi (deg),psi (deg),alpha (deg),theta (deg)\n")
    #    Calculate the torsion angles and the pseudo angles and then put into the output file
    #    Keep in mind that ProDy 1-INDEXES the residues
    for i in range(1, len(Sequence) + 1):
        if i == 1 or i == len(Sequence):
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()}\n")
        else:
            phi = pd.calcPhi(hv.getResidue('A', i))
            psi = pd.calcPsi(hv.getResidue('A', i))
            phi_r = phi * (np.pi / 180)
            psi_r = psi * (np.pi / 180)
            theta = THETA(t_r, y1_r, y2_r, phi_r, psi_r)
            alpha = ALPHA_COMPLEX(t_r, y1_r, y2_r, phi_r, psi_r)
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()},{phi:.2f},{psi:.2f},{alpha:.2f},{theta:.2f}\n")
    outfile.close()
    # Lastly, save the structure in CA form
    linker_ca = linker_pdb.select("name CA")
    pd.writePDB(ca_name, linker_ca)
    return None


def convert_pdb(pdb_filename, Sequence, PeptideName):
    # ATTENTION -- this is a complete copy from the code above! Because it's already there.
    # Now this part of the code will reopen this fully atomistic structure, read the torsion angles of it and convert into
    # Tozzini space, then save the C-alpha coarse grained model. It will also save the torsions and the pseudotorsions into
    # and output file
    ca_name = f"{PeptideName}_ca.pdb"
    output_file = f"{PeptideName}_aa_dihedrals.csv"

    linker_pdb = pd.parsePDB(pdb_filename)
    hv = pd.HierView(linker_pdb)
    outfile = open(output_file, 'w')
    outfile.write("Residue No.,Residue Name,phi (deg),psi (deg),alpha (deg),theta (deg)\n")
    #    Calculate the torsion angles and the pseudo angles and then put into the output file
    #    Keep in mind that ProDy 1-INDEXES the residues
    for i in range(1, len(Sequence) + 1):
        if i == 1 or i == len(Sequence):
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()}\n")
        else:
            phi = pd.calcPhi(hv.getResidue('A', i))
            psi = pd.calcPsi(hv.getResidue('A', i))
            phi_r = phi * (np.pi / 180)
            psi_r = psi * (np.pi / 180)
            theta = THETA(t_r, y1_r, y2_r, phi_r, psi_r)
            alpha = ALPHA_COMPLEX(t_r, y1_r, y2_r, phi_r, psi_r)
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()},{phi:.2f},{psi:.2f},{alpha:.2f},{theta:.2f}\n")
    outfile.close()
    # Lastly, save the structure in CA form
    linker_ca = linker_pdb.select("name CA")
    pd.writePDB(ca_name, linker_ca)
    return None





"""
Run the code here!
"""
if __name__ == "__main__":
    # I don't anticipate needing to use the functions in this script elsewhere, so the
    # name==main style is a little overkill, but I'm including it here for good pythonic practice

    # Set up the argument parser!
    parser = argparse.ArgumentParser(description="Model Maker Tool -- use this to turn your desired CC protein into"
                                                 " a starting PDB structure! Easy, and corresponds with the Topology"
                                                 " Assistant too. Spreads on easy, like Ubik!")
    parser.add_argument("-fasta", help="The fasta file containing the sequence of your protein of interest!")
    parser.add_argument("-conversion", help="Switch. Pass this flag to turn on conversion only.",
                        action="store_true")
    parser.add_argument("-convert_pdb", help="The file name (plus extension) (plus path) of a PDB file you "
                        "wish to convert into C-alpha coarse graining through Tozzini calculations.", 
                        default=None)
    parser.add_argument("-n", help="The number of beads in the CC protein. This is useful for sanity checking.",
                        type=int)
    parser.add_argument("-name", help="The name you'd like to give the CC protein. No "
                                      "extensions, please!", type=str, default="coil_model")
    

    # parse the arguments
    args = parser.parse_args()
    fasta_input = args.fasta
    protein_length = args.n
    protein_name = args.name
    pdb_file = args.convert_pdb

    prot_fasta_sequence = fasta_parser(fasta_input)

    if args.conversion:
        # this step will only happen if I turn on conversion, since all I have to do is convert the existing file!
        convert_pdb(pdb_file, prot_fasta_sequence, protein_name)
    else:
        # otherwise, make a brand new structure file for the provided sequence.
        assert len(prot_fasta_sequence) == args.n, f"The sequence detected in the fasta file does not match the number of beads "\
                                            f"that are expected."
        build_peptide(protein_name, pprot_fasta_sequence)

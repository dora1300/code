"""
@Author:            Dominique Ramirez
@Date:              2021 11 12

@Name:              dihedral_topology_selector.py

@Description:       This script implements the solution of the gmx cosine dihedral potential so that, for a given
dihedral PHI, I can find the reference PHI_0 that I need to put into the topology
"""

import argparse
import math as m

def angle_calc(phi, n):
    """
    The function to calculate the *reference angle* used in the GMX cosine dihedral potential.
    :param phi: IN DEGREES - the *minimum angle* where the potential will be at a minimum
    :param n: the multiplicity of the cosine potential
    :return:
    """
    phi_0 = (n*phi) - (m.acos(-1) * (180 / m.pi))
    return phi_0


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Simple script to calculate the Phi_naught used in the Gromacs"
                                                 " cosine potential. The angle given is the angle where the potential "
                                                 "will be a minimum.")
    parser.add_argument("-p", help="Angle (termed phi) where the cosine potential is at a minimum.", type=float)
    parser.add_argument("-n", help="Multiplicity of the cosine potential", type=int)

    args = parser.parse_args()

    print(f"Provided minimum angle: {args.p}")
    print(f"Multiplicity of potential: {args.n}")
    print(f"Reference angle of potential: {angle_calc(args.p, args.n)}")

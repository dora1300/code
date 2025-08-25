"""
@Title:             write_atomtypes_for_top.py
@Name:              Mando A Ramirez
@Date:              2025 03 06

@Description:       This script generates the file that will take the place
of the [atomtypes] directive in the primary .top file.
"""

import argparse         # this is for when I test the script solo
import numpy as np

""" GLOBAL DEFINITIONS AND OTHER USEFUL THINGS """
aas_list = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys',
            'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro',
            'Ser', 'Thr', 'Trp', 'Tyr', 'Val']

aas_one_letter = ['A', 'R', 'N', 'D', 'C',
                  'Q', 'E', 'G', 'H', 'I',
                  'L', 'K', 'M', 'F', 'P',
                  'S', 'T', 'W', 'Y', 'V']

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


"""  THE FUNCTION TO WRITE THE ITP FILES  """

def write_atomtypes_itp(output_file, attractive_term=6.9132E-03, repulsive_term=5.97404E-05):
    """
    Description:
        [ atomtypes ]
        ;name   mass        charge   ptype  V(C6, attra)  W(C12, repul)
    Arguments:
    Returns:        
    """
    output_text = f"""[ atomtypes ]
;name   mass        charge   ptype  V(C6, attra)    W(C12, repul)
"""
    for aa_index, aa in enumerate(aas_list):
        line_text = f" {aas_one_letter[aa_index]}    {aa_mass_dict[aa]}    0.000    A      {attractive_term:.5e}       {repulsive_term:.5e}\n"
        output_text += line_text

    writer_file = open(output_file, 'w')
    writer_file.write(output_text)
    writer_file.close()

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A short little script to write the atomtypes.itp file to include with your SSFF topology!")
    parser.add_argument("-itp_filename")
    parser.add_argument("-attractive_term", default=None)
    parser.add_argument("-repulsive_term", default=None)

    args = parser.parse_args()
    FILENAME = args.itp_filename
    ATTR = args.attractive_term
    REPUL = args.repulsive_term

    if (ATTR is None) and (REPUL is None):
        write_atomtypes_itp(FILENAME)
    elif (ATTR is not None) and (REPUL is None):
        write_atomtypes_itp(FILENAME, attractive_term=ATTR)
    elif (ATTR is None) and (REPUL is not None):
        write_atomtypes_itp(FILENAME, repulsive_term=REPUL)
    else:
        write_atomtypes_itp(FILENAME, attractive_term=ATTR, repulsive_term=REPUL)
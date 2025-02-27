"""
@Title:             make_nonbonded_parameters_itp.py
@Name:              Mando A Ramirez
@Date:              2025 02 01

@Description:       This script generates the big kahuna -- the nonbonded_parameters.itp
file that is used with the sequence specific simulation framework and associated topologies.

This file is generated from a 20x20 matrix of epsilsons and also possible sigmas. That is up
to the user.
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

dignon_sigmas = [0.504, 0.656, 0.568, 0.558, 0.548,
                 0.602, 0.592, 0.450, 0.608, 0.618,
                 0.618, 0.636, 0.618, 0.636, 0.556,
                 0.518, 0.562, 0.678, 0.646, 0.586]

# units of dignon_sigmas are in nm
# the Dignon sigmas come from his PLOS One 2018 paper describing his
# coarse-grained sequence specific model. They are in alphabetical order.
# They also seen to be based on the Kim and Hummer (2008) sigma values

""" FUNCTION DEFINITIONS """

def write_nonbonded_itp(epsilon_file_name, sigma_table):
    """
    Description:
        This is the main function for generating the text that will be written
        to the final .itp file.
    
    Arguments:

    Returns:
    """

    # load the table of epsilons. This is mandatory! I can't make the .itp file
    # without it!
    epsilon_table = np.loadtxt(epsilon_file_name, delimiter=',')    # might need to add encoding thing

    # Now the magic happens. It's file to make the file contents
    nonbond_params_text = f"""[ nonbond_params ] 
;i    j     func    V(sigma)   W(epsilon)\n"""

    for I in range(len(aas_list)):
        for J in range(I, len(aas_list)):
            new_string = f"{aas_list[I]}   {aas_list[J]}   1       {sigma_table[I][J]:.5f}    {epsilon_table[I][J]:.5f}\n"
            nonbond_params_text += new_string

    return nonbond_params_text



def load_calculated_sigmas(SIGMA_BY_MATRIX=True, SIGMA_MATRIX_FILE=None,
                           SIGMA_BY_AA=False, SIGMA_AA_NAME=None,
                           DIGNON=False, combining_rule="lorentz"):
    """
    Description:

    Arguments:
    
    Returns:
        A 20x20 numpy array containing the pairwise amino acid
        sigma values for use in the LJ functional in GROMACS.
        Please recognize this is a SYMMETRIC array BUT!!!!!!!!!!!!!!!
        if SIGMA_BY_AA is true or DIGNON is true then the array will only
        be filled out above the diagnoal! This is fine, and will cause no
        problems just be aware of it!
    """
    if SIGMA_BY_MATRIX is True:
        # this simply reads a file. The file name is provided in argument
        # SIGMA_MATRIX_FILE. This stops the function after loading because
        # the sigmas have been made, game over.
        sigma_matrix = np.loadtxt(SIGMA_MATRIX_FILE, delimiter=",", dtype=float)
        return sigma_matrix
    

    # now, I need to generate the sigma matrix if I'm not reading from a file. 
    # Thus, I need to make a 20x20 numpy array first to store the values
    sigma_matrix_manual = np.zeros((20, 20), dtype=float)
    if combining_rule == "lorentz":
        CR = lambda s1, s2 : (s1 + s2) / 2.

    if SIGMA_BY_AA is True:
        # calculate a matrix of sigmas based on a list of sigmas for each amino acid.
        # these must be combined. The list of sigmas MUST be in alphabetical order as
        # specified above. The file MUST be .csv
        sigma_aa_list = np.loadtxt(SIGMA_AA_NAME, delimiter=',', dtype=float)
        for I in range(len(sigma_aa_list)):
            for J in range(I, len(sigma_aa_list)):
                sigma_matrix_manual[I][J] = CR(sigma_aa_list[I], sigma_aa_list[J])
        return sigma_matrix_manual
    
    if DIGNON is True:
        # calculate sigmas based on dignon then return. No need to go onto others
        # these must be combined. Use whatever is specified by the function arguments.
        for I in range(len(dignon_sigmas)):
            for J in range(I, len(dignon_sigmas)):
                sigma_matrix_manual[I][J] = CR(dignon_sigmas[I], dignon_sigmas[J])
        return sigma_matrix_manual
    


def save_itp_file_to_disk(nonbonded_file_text, itp_file_name):
    """
    Description:
        This is the main function for generating the text that will be written
        to the final .itp file.
    
    Arguments:

    Returns:
    """

    with open(itp_file_name, 'w') as OUTPUTFILE:
        OUTPUTFILE.write(nonbonded_file_text)
    
    return None




""" MAIN area for testing this script solo """
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-name_of_itp_file", help="The name that you'd like to give to the file! "
                        "Please include the extension and the path if appropriate!")
    
    parser.add_argument("-epsilon_file_name", help="[csv] the 20x20 matrix/CSV file containing the "
                        "epsilons for the itp file.")
    
    parser.add_argument("-sigma_by_matrix", help="Switch. Are you using a matrix of pre-calculated sigmas?",
                        default=False, action="store_true")
    
    parser.add_argument("-sigma_matrix_file", help="[csv] Name (and path) of the file containing the 20x20 "
                        "matrix of sigma values to use.")
    
    parser.add_argument("-sigma_by_aa", help="Switch. Are you providing a list of sigma values but only for "
                        "each individual amino acids?",
                        default=False, action="store_true")
    
    parser.add_argument("-sigma_aa_name", help="[csv] Name of the file which contains the sigma values for "
                        "each individual amino acid. Make sure each sigma gets its own line AND that the "
                        "sigmas are in alphabetical order corresponding to their amino acids! No labels please!")
    
    parser.add_argument("-sigma_by_dignon", help="Switch. Do you want to generate sigmas based on the values "
                        "used by Dignon et al. 2018 PLOS One?",
                        default=False, action="store_true")
    
    parser.add_argument("-combining_rule", help="The code-name of the combining rule to use for calculating "
                        "sigmas. Only 'lorentz' is supported right now.",
                        default="lorentz")

    args = parser.parse_args()

    # Step 1 -- calculate the sigmas based on the setting I chose
    array_of_sigmas = load_calculated_sigmas(args.sigma_by_matrix, args.sigma_matrix_file,
                                             args.sigma_by_aa, args.sigma_aa_name,
                                             args.sigma_by_dignon, args.combining_rule)
    
    # Step 2 -- generate the text that will make up the file contents
    itp_file_contents = write_nonbonded_itp(args.epsilon_file_name, array_of_sigmas)

    # Step 3 -- save the itp file to disk!
    save_itp_file_to_disk(itp_file_contents, args.name_of_itp_file)

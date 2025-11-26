"""
@Title:             write_seqspec_top_file.py
@Name:              Mando A Ramirez
@Date:              20250306

@Description:       This script generates a topology file for a sequence specific force field
simulation of coiled-coils. Hence the abbreviation "ssff" at the start of the code.

This is not the only code you can or should run to generate the things you need. This only 
generates a .top file

@Updates:
20250708            Originally part of stand-alone testing for sequence specific parameters.
Now this is going to be part of an automated workflow for setting up validation simulations for
testing sequence specific parameters.
This will allow the user to specify multiple molecule types as necessary and all the relevant 
copies and .itp files. 
This script should be generalizable for any type of set up that I want, since the protein .itp files
are handled elsewhere.
"""
import argparse         # this is for when I test the script solo




""" FUNCTION DEFINITIONS """
def write_text_top(atomtypes_name, 
                   nonbonded_name, 
                   itp_names:list, 
                   moleculetype_name:list, 
                   molecule_number:list,
                   system_name,
                   output_name,
                   genpairs_val,
                   fudgelj_val):
    text = f""";
;           TOPOLOGY FILE
; This is a topology file for a sequence specific simulation
; of coiled-coil protein(s)
;

[ defaults ]
;nbfunc     comb-rule      gen-pairs        fudgeLJ     fudgeQQ
1           2              {genpairs_val}              {fudgelj_val}         0.0

; this is where [ atomtypes ] directive would live, but is incorporated 
; through an include statement
#include "{atomtypes_name}"

; this is where the [ nonbond_params ] would live, and it is also
; incorporated with an include statement
#include "{nonbonded_name}"

"""
    
    # handle all of the itp files now
    text += f"""; this is where the .itp files for all molecules lives
; they must be handled with include statements
"""
    for itp_index, itp_file in enumerate(itp_names):
        text += f"""#include "{itp_file}"
"""
        

    # now it's time to set up the system directive
    text += f"""
[ system ]
{system_name}
"""
        
    # handle all the molecule types now
    text += f"""
[ molecules ]
; provide the moleculetype-name from each molecule itp and also the total number in the simulation
;moleculetype-name      # molecules
"""
    for index, molecule_name in enumerate(moleculetype_name):
        text += f"""{molecule_name}                 {molecule_number[index]}
"""

    writer_file = open(output_name, 'w')
    writer_file.write(text)
    writer_file.close()



""" MAIN area for running the code! """
if __name__ == '__main__':
    # run this code if I'm running this script by itself!
    """
        Argument Parser Set up
    """
    parser = argparse.ArgumentParser(description="The script to write a custom topology (.top) file for a sequence specific simulation")
    parser.add_argument('-atomtypes_name', help="The name of the atomtypes.itp file you want to include in this topology."
                        " Please include extension.", type=str)
    parser.add_argument('-nonbonded_name', help="The name of the nonbonded.itp file you want to include in this toplogy."
                        " Please include extension.", type=str)
    parser.add_argument('-itp_file_name', help="The name of the itp files that you want to be included in the top file. "
                        "Please include extension. The order of the itp files must match the order of the molecule names and the "
                        "order of the number of molecules!", required=True, nargs="+")
    parser.add_argument('-molecule_type_name', help="The name of the first molecule type that you want to be included in the top file. You "
                        "can add more manually later.", required=True, nargs="+")
    parser.add_argument('-number_molecule_type', help="The number of the molecule_type_name to include in your topology. This of course"
                        " can be changed manually later.", required=True, nargs="+")
    parser.add_argument('-system_name', help="This is the name for the simulation given under the [ system ] directive. Doesn't really "
                        "do anything so you can use default if you wish", default="coiled coil simulation")
    parser.add_argument('-output_file', help="Name that you wish to give to the .top file that you've just made!"
                        " Include extension!", type=str)

    args = parser.parse_args()
    ATOMTYPE_NAME = args.atomtypes_name
    NONBONDED_NAME = args.nonbonded_name
    ITP_NAME = args.itp_file_name
    MOLECULETYPE_NAME = args.molecule_type_name
    NUM_MOLTYPE = args.number_molecule_type
    SYS_NAME = args.system_name
    OUTPUT = args.output_file

    if len(MOLECULETYPE_NAME) != len(ITP_NAME):
        raise ValueError("The number of molecule_type names does not match the number of .itp files provided. Please review your arguments. "\
                         "and try again.")
    else:
        pass

    if len(MOLECULETYPE_NAME) != len(NUM_MOLTYPE):
        raise ValueError("The number of molecule_type names does not match the number of molecule_numbers provided. There is a mismatch. "\
                         "Please review your arguments and try again.")
    else:
        pass

    write_text_top(ATOMTYPE_NAME, NONBONDED_NAME, ITP_NAME, MOLECULETYPE_NAME, NUM_MOLTYPE, SYS_NAME, OUTPUT)
"""
@Title:             ssff_full_top_writer.py
@Name:              Mando A Ramirez
@Date:              20250306

@Description:       This script generates a topology file for a sequence specific force field
simulation of coiled-coils. Hence the abbreviation "ssff" at the start of the code.

This is not the only code you can or should run to generate the things you need. This only 
generates a .top file
"""
import argparse         # this is for when I test the script solo




""" FUNCTION DEFINITIONS """
def write_text_top(atomtypes_name, nonbonded_name, output_name):
    text = f""";
;           TOPOLOGY FILE
; This is a topology file for a sequence specific simulation
; of coiled-coil protein(s)
;
; made using 'ssff_full_top_writer.py', est. 20250306

[ defaults ]
;nbfunc     comb-rule      gen-pairs        fudgeLJ     fudgeQQ
1           2              yes              1.0         0.0

; this is where [ atomtypes ] directive would live, but is incorporated 
; through an include statement
#include "{atomtypes_name}"

; this is where the [ nonbond_params ] would live, and it is lso
; incorporated with an include statement
#include "{nonbonded_name}"

; at this point, you must include the .itps for all the unique molecules
; you are going to simulating! 
; this must be done manually
#include "INSERT_FILE_HERE.itp"

[ system ]
; feel free to rename
Sequence specific simulation

[ molecules ]
; provide the moleculetype-name from each molecule itp and also the total number in the simulation
;moleculetype-name      # molecules
EXAMPLE                 1
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
                        " Please include extension", type=str)
    parser.add_argument('-output_file', help="Name that you wish to give to the .top file that you've just made!"
                        " Include extension!", type=str)

    args = parser.parse_args()
    ATOMTYPE_NAME = args.atomtypes_name
    NONBONDED_NAME = args.nonbonded_name
    OUTPUT = args.output_file

    write_text_top(ATOMTYPE_NAME, NONBONDED_NAME, OUTPUT)
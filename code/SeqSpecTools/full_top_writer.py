"""
@Title:             top_writer.py
@Name:              Mando A Ramirez
@Date:              2024 01 16

@Description:       This script generates a topology file 
"""
import argparse         # this is for when I test the script solo




""" FUNCTION DEFINITIONS """
def write_text_top():
    text = f""";
;           TOPOLOGY FILE
; This is a topology file for a sequence specific simulation
; of coiled-coil protein(s)
;
; made using 'top_writer.py', est. 20240116

[ defaults ]
;nbfunc     comb-rule      gen-pairs        fudgeLJ     fudgeQQ
1           2              yes              1.0         0.0

; this is where [ atomtypes ] directive would live, but is incorporated 
; through an include statement
#include "{customName}_atomtypes.itp"

; this is where the [ nonbond_params ] would live, and it is lso
; incorporated with an include statement
#include "{customeName}_nonbonded.itp

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

# def WRITE_FULL_TOP():




""" MAIN area for testing this script solo """
if __name__ == '__main__':
    # do something
    pass
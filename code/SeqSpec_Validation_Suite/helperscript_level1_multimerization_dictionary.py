"""
Dictionary to help store important variables for each of the different proteins that can be simulated
in the validation suite for Level 1 of multimerization validation.

each entry takes the form of
codename : [["protein_names_that_go_into_.top", ...], 
            [numbers_of_each_type_of_protein,...],
            COM_pulling_number:str,
            [length_of_each_coil,...],
            [[n,term,torsion,indices,...]],
            [[c,term,torsion,indices,...]],
            [[distance,ref,indices,...]],
            [comments about the structure, for posterity]
            ]
"""

validation_dictionary = {

#
#       POSITIVE CONTROLS -- THESE SHOULD BIND
#
    "s1hs2h":[["s1h", "s2h"],
              [1, 1],
              "2",
              [21, 21],
              [[1, 8, 29, 22]],
              [[18, 11, 32, 39]],
              [[1, 22], [8, 29], [11, 32], [18, 39]]
    ],
    
    "p5fp6f":[["p5f", "p6f"],
              [1, 1],
              "2",
              [28, 28],
              [[1, 11, 39, 29]],
              [[25, 15, 43, 53]],
              [[1, 29], [11, 39], [15, 43], [25, 53]]

    ],

    "gcnh2gcnh2":[["gcnh2"],
                  [2],
                  "2",
                  [21, 21],
                  [[1, 8, 29, 22]],
                  [[18, 11, 32, 39]],
                  [[1, 22], [8, 29], [11, 32], [18, 39]]
    ],

    "aph2aph2":[["aph2"],
                [2],
                "2",
                [21, 21],
                [[1, 8, 32, 39]],
                [[18, 11, 32, 22]],
                [[1, 39], [8, 32], [11, 29], [18, 22]]
    ],

    "4dzm":[["4dzm"],
            [2],
            "2",
            [32, 32],
            [[2, 12, 44, 34]],
            [[26, 16, 48, 58]],
            [[1, 34], [12, 44], [16, 48], [26, 58]]
    ],

    "7q1r":[["7q1r"],
            [2],
            "2",
            [30, 30], 
            [[2, 13, 46, 57]],
            [[27, 16, 43, 32]],
            [[2, 57], [13, 46], [16, 43], [27, 32]]
    ],

    "4dzl" : [["4dzl"],
              [3],
              "3",
              [32, 32, 32],
              [[2, 12, 34, 44], [2, 12, 66, 76]],
              [[26, 16, 48, 58], [26, 16, 80, 90]],
              [[2, 34], [12, 44], [16, 48], [26, 58],
               [2, 66], [12, 76], [16, 80], [26, 90]]
    ],

    "gcn4" : [["gcn4"],
              [2],
              "2",
              [31, 31],
              [[1, 11, 42, 32]],
              [[25, 18, 49, 56]],
              [[1, 32], [8, 39], [11, 42], [15, 46], [18, 49], [25, 56]],
              ["This comes from PDB code 1ysa. This does not exist in the CC+ database for reasons I don't understand because GCN4 is considered a canonical coiled-coil (leucine zipper).",
               "delete the final ER from the sequence in 4dmd. The sequence will only be 31 amino acids"]

    ],

    "3n4x" : [["3n4x"],
              [2],
              "2",
              [47, 47],
              [[0, 14, 61, 47]],
              [[46, 32, 79, 93]],
              [[0,47], [4,51], [7,54], [11,58], [14,61], [18,65], 
               [21,68], [25,72], [28,75], [32,79], [35,82], 
               [39,86], [42,89], [46,93]],
              ["from the CC+ database, a longerish coil",
               "parallel orientation"]

    ],

#
#       NEGATIVE CONTROLS -- THESE SHOULD NOT BIND
#
    "gcnh2aph2":[["gcnh2", "aph2"],
                 [1, 1],
                 "2",
                 [21, 21],
                 [[1, 8, 29, 22]],
                 [[18, 11, 32, 39]],
                 [[1, 22], [8, 29], [11, 32], [18, 39]]
    ],

    "gcnh2s2h":[["gcnh2", "s2h"],
                [1, 1],
                "2",
                [21, 21],
                [[1, 8, 29, 22]],
                [[18, 11, 32, 39]],
                [[1, 22], [8, 29], [11, 32], [18, 39]]
    ],

    "p6f4dzm":[["p6f", "4dzm"],
               [1, 1],
               "2",
               [28, 32],
               [[1, 15, 44, 30]],
               [[25, 15, 44, 54]],
               [[1, 30], [11, 40], [15, 44], [25, 54]]
    ],

    "p5fgcn4":[["p5f", "gcn4"],
               [1, 1],
               "2",
               [28, 30],
               [[1, 11, 39, 29]],
               [[25, 15, 43, 53]],
               [[1, 29], [11, 39], [15, 43], [25, 53]]
    ],

    "4dzlgcn4":[["4dzl", "gcn4"],
                [1, 1],
                "2",
                [32, 30],
                [[2, 12, 42, 33]],
                [[25, 16, 47, 57]],
                [[1, 33], [12, 43], [16, 47], [26, 57]]
    ]
}


if __name__ == "__main__":
    print("Keys of valid protein codenames:")
    print(validation_dictionary.keys())
    print()
    print("Dictionary of possible protein codenames to run for Level 1 Multimerization validation:")
    print(validation_dictionary)
    print()
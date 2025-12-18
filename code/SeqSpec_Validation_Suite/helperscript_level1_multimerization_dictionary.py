"""
Dictionary to help store important variables for each of the different proteins that can be simulated
in the validation suite for Level 1 of multimerization validation.

each entry takes the form of
codename : [["protein_names_that_go_into_.top", ...], 
            [numbers_of_each_type_of_protein,...],
            COM_pulling_number:str,
            [length_of_each_coil,...],
            [[n,term,torsion,indices,...]], 0-indexed!!!!!
            [[c,term,torsion,indices,...]], 0-indexed!!!!!
            [[distance,ref,indices,...]],   0-indexed!!!!!
            [heptad length for each coil,...],
            ["positive_or_negative_control"]
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
              [[1, 22], [8, 29], [11, 32], [18, 39]],
              [3, 3],
              ["positive"],
              ["3 hetpad long designer coils from Ramsak et al."]
    ],
    
    "p5fp6f":[["p5f", "p6f"],
              [1, 1],
              "2",
              [28, 28],
              [[1, 11, 39, 29]],
              [[25, 15, 43, 53]],
              [[1, 29], [11, 39], [15, 43], [25, 53]],
              [4, 4],
              ["positive"],
               ["4 hetpad long designer coils from Ramsak et al."]

    ],

    "gcnh2gcnh2":[["gcnh2"],
                  [2],
                  "2",
                  [21, 21],
                  [[1, 8, 29, 22]],
                  [[18, 11, 32, 39]],
                  [[1, 22], [8, 29], [11, 32], [18, 39]],
                  [3, 3],
                  ["positive"],
                  ["3 hetpad long designer coils from Ramsak et al.",
                    "parallel coils"]
    ],

    "aph2aph2":[["aph2"],
                [2],
                "2",
                [21, 21],
                [[1, 8, 32, 39]],
                [[18, 11, 32, 22]],
                [[1, 39], [8, 32], [11, 29], [18, 22]],
                [3, 3],
                ["positive"],
                ["3 hetpad long designer coils from Ramsak et al.",
                    "antiparallel coils"]
    ],

    "4dzm":[["4dzm"],
            [2],
            "2",
            [32, 32],
            [[2, 12, 44, 34]],
            [[26, 16, 48, 58]],
            [[1, 34], [12, 44], [16, 48], [26, 58]],
            [4.5, 4.5],
            ["positive"]
    ],

    "7q1r":[["7q1r"],
            [2],
            "2",
            [30, 30], 
            [[2, 13, 46, 57]],
            [[27, 16, 43, 32]],
            [[2, 57], [13, 46], [16, 43], [27, 32]],
            [4, 4],
            ["positive"]
    ],

    "4dzl" : [["4dzl"],
              [3],
              "3",
              [32, 32, 32],
              [[2, 12, 34, 44], [2, 12, 66, 76]],
              [[26, 16, 48, 58], [26, 16, 80, 90]],
              [[2, 34], [12, 44], [16, 48], [26, 58],
               [2, 66], [12, 76], [16, 80], [26, 90]],
               [4.5, 4.5, 4.5],
               ["positive"]
    ],

    "gcn4-2" : [["gcn4-2"],
              [2],
              "2",
              [31, 31],
              [[1, 11, 42, 32]],
              [[25, 18, 49, 56]],
              [[1, 32], [8, 39], [11, 42], [15, 46], [18, 49], [25, 56]],
              [4, 4],
              ["positive"],
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
               [6.7, 6.7],
               ["positive"],
              ["from the CC+ database, a longerish coil",
               "parallel orientation. Csm1 from S. Cerevisiae"]

    ],

    "2fxm" : [["2fxm"],
              [2],
              "2",
              [106, 106], 
              [[0, 25, 131, 106]],
              [[105, 81, 187, 211]],
              [[0, 106], [4, 110], [7, 113], [11, 117], [14, 120], [18, 124], [21, 127], 
               [25, 131], [28, 134], [32, 138], [35, 141], [39, 145], [42, 148], [46, 152], 
               [49, 155], [53, 159], [56, 162], [60, 166], [63, 169], [67, 173], [70, 176], 
               [74, 180], [77, 183], [81, 187], [84, 190], [88, 194], [91, 197], [95, 201], 
               [98, 204], [102, 208], [105, 211]],
              [15.1, 15.1],
              ["positive"],
              ["from the CC+ database, a coil with KK and EE interactions along the interface"]
    ],

    "7d9r" : [["7d9r_chain1", "7d9r_chain2"],
              [1, 1],
              "2",
              [46, 46],
              [[0, 14, 60, 46]],
              [[45, 31, 77, 91]],
              [[0, 46], [3, 49], [7, 53], [10, 56], [14, 60], [17, 63], 
               [21, 67], [24, 70], [28, 74], [31, 77], [35, 81], 
               [38, 84], [42, 88], [45, 91]],
              [6.6, 6.6],
              ["positive"],
              ["from the CC+ database, a coil with KK and EE interactions along the interface"]    ],

#
#       NEGATIVE CONTROLS -- THESE SHOULD NOT BIND
#
    "gcnh2aph2":[["gcnh2", "aph2"],
                 [1, 1],
                 "2",
                 [21, 21],
                 [[1, 8, 29, 22]],
                 [[18, 11, 32, 39]],
                 [[1, 22], [8, 29], [11, 32], [18, 39]],
                 [3, 3],
                 ["negative"]
    ],

    "gcnh2s2h":[["gcnh2", "s2h"],
                [1, 1],
                "2",
                [21, 21],
                [[1, 8, 29, 22]],
                [[18, 11, 32, 39]],
                [[1, 22], [8, 29], [11, 32], [18, 39]],
                [3, 3],
                ["negative"]
    ],

    "p6f4dzm":[["p6f", "4dzm"],
               [1, 1],
               "2",
               [28, 32],
               [[1, 15, 44, 30]],
               [[25, 15, 44, 54]],
               [[1, 30], [11, 40], [15, 44], [25, 54]],
               [4, 4.5],
               ["negative"]
    ],

    "p5fgcn4":[["p5f", "gcn4"],
               [1, 1],
               "2",
               [28, 30],
               [[1, 11, 39, 29]],
               [[25, 15, 43, 53]],
               [[1, 29], [11, 39], [15, 43], [25, 53]],
               [4, 4.3],
               ["negative"]
    ],

    "4dzlgcn4":[["4dzl", "gcn4"],
                [1, 1],
                "2",
                [32, 30],
                [[2, 12, 42, 33]],
                [[25, 16, 47, 57]],
                [[1, 33], [12, 43], [16, 47], [26, 57]],
                [4.5, 4.3],
                ["negative"]
    ]
}


if __name__ == "__main__":
    print("Keys of valid protein codenames:")
    print(validation_dictionary.keys())
    print()
    print("Dictionary of possible protein codenames to run for Level 1 Multimerization validation:")
    print(validation_dictionary)
    print()
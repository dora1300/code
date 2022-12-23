"""
Name:           autocorrelate_analysis.py

Author:         Dominique A. Ramirez
Date:           2022 12 22

Description:    This script parses a .csv file and calculates the autocorrelation of the provided data. The data file must be formulated such that: Column 1 = x-axis units that are in the time-lag,
and Column 2 = the data that the ACF is calculated on. If there is any other data it is ignored.

"""

import argparse
import numpy as np
import os.path



""" Function Definitions """
def autocorrelate(array, start_lag, end_lag):
    """
    This function calculates the autocorrelations for different lag/shift phases of the given array
    :param array: (np.array) the array of data that will be used for autocorrelation analysis
    :param start_lag: (int) starting lag time for autocorrelation
    :param end_lag: (int) end lag time, i.e. number of cycles to progress through
    :return: (list) autocorrelation values
    """
    auto_results = []
    if start_lag == 0:
        auto_results.append(np.corrcoef(array, array)[0,1])
        for index in range(start_lag+1, end_lag+1):
            auto_results.append(np.corrcoef(array[:-index], array[index:])[0,1])
    else:
        for index in range(start_lag, end_lag+1):
            auto_results.append(np.corrcoef(array[:-index], array[index:])[0,1])
    return auto_results




""" ArgParser set up """
parser = argparse.ArgumentParser(description="Autocorrelation analysis of provided data (must be"
    " in .csv file format).")

# -- Positional arguments
parser.add_argument("-f", action="store", help="Name of file (in directory) or complete path to file, WITH extension (.csv)", required=True)
parser.add_argument("-o", action="store", help="Name of autocorrelation results output .csv file (include extension)", default="acf_output.csv")
parser.add_argument("-beg", help="The starting lag-tau index range for the autocorrelation"
    " analysis.", required=True, type=int)
parser.add_argument("-end", help="The max lag-tau index range.", required=True, type=int)
parser.add_argument("-u", help="Units of the x-axis in the series data (e.g. time, or frames)",
    required=True, type=str)
    
args = parser.parse_args()

message = """Please note! The starting and ending lag indicies you provided refer to the position INDEX of the data, and has nothing to do with the x-axis units of the series the data comes from.
E.g. if there are 20000 objects in the series data, you could specify a starting index = 0 and ending index = 10000 (half index is most common for acf). 
But, if each of these entries corresponds to 0.25 ns of simulation data (i.e. a total of 5 us simulation time), then you COULD NOT use an end value of 2.5E+06 because there are not that many entries in the series data.
Please make sure you understand this before proceeding.
"""
print(message)



""" Parse data and do the analysis """
# -- extract the data from the input file.
input_x = [] ; input_data = []
with open(args.f, 'r') as f:
    for line in f:
        if line[0] == "#" or line[0] == "@":
            continue
        else:
            line_split = line.lstrip("").rstrip("\n").split(",")
            input_x.append(float(line_split[0]))
            input_data.append(float(line_split[1]))
            
            

# -- now do the analysis
if args.end >= len(input_data):
    print("The provided max lag-index is out of range of the data provided.")
    print(f"Length of data array: {len(input_data)}")
    print(f"Max lag index range: {args.end}")
    exit(1)
acf_results = autocorrelate(np.array(input_data), args.beg, args.end)



# -- now save the data into a csv file.
with open(args.o, 'w') as orf:
    orf.write(f"#lag-tau (original units {args.u}), acf_data\n")
    for i in range(len(acf_results)):
        orf.write(f"{args.beg + i},{acf_results[i]}\n")

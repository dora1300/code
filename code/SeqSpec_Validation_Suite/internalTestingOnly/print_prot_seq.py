
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f")

args = parser.parse_args()

FASTAFILE = args.f

with open(FASTAFILE) as fasta:
    for line in fasta:
        if line[0] == ">":
            sequence_name = line[1:].rstrip("\n")
        elif line[0] == "#":
            # I'm allowing comments in my fasta file, which is not standard format
            # but I might want to make notes
            continue
        elif line[0] != ">" and line[0] != "#":
            sequence = line.rstrip("\n")
        else:
            break

for i in range(len(sequence)):
    print(f"{i+1},{sequence[i]},")
import random
import argparse

parser = argparse.ArgumentParser(
    description="Silly file! Not really useful but if you need a randomly generated sequence of amino acids "
    "then here you go."
)

parser.add_argument("-n", help="the number of amino acids to randomly generate", type=int, required=True)
args = parser.parse_args()

aas_one_letter = ['A', 'R', 'N', 'D', 'C',
                  'Q', 'E', 'G', 'H', 'I',
                  'L', 'K', 'M', 'F', 'P',
                  'S', 'T', 'W', 'Y', 'V']


random_sequence = ""

while len(random_sequence) < args.n:
    random_integer = random.randint(0, 19)
    random_sequence += aas_one_letter[random_integer]

print()
print(random_sequence)
print()
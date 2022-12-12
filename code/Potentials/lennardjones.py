"""
Simple implementation of the GROMACS standard LJ potential for studying/didactic purposes.

@Author: DA Ramirez
@Date: 20210916

Updated: 20211014 - added the 12-10 potential for comparison
"""

import numpy
import math
import matplotlib.pyplot as plt

def lj(epsilon, sigma, r):
	"""
	Implementation of the Type 2 nonbonded function for GROMACS.
	
	@epsilon	: the epsilon value for atoms 1 and 2. Can be the same or different (kJ / mol)
	@sigma		: the sigma value for atoms 1 and 2. Can be the same or different (nm)
		Both of these need to be combined according to a specific rule before using.
	@r			: the distance between atoms 1 and 2
	
	$return y	: an array or single value of the LJ potential
	"""
	y = 4*epsilon * ((sigma/r)**12 - (sigma/r)**6)
	return y
	
def ca_lj(A, C, r):
	"""
	Implementation of the SMOG 12-10 LJ potential form, which uses A and C as seen in the topology.
	
	@twelve		: the 5*r_ij**12 term corresponding to the repulsive form of the 12-10
	@ten		: the 6*r_ij**10 term corresponding to the attractive form.
	@r			: the distance between the two particles
	"""
	y = (C/r**12) - (A/r**10)
	return y


"""
Global Constants for the gromacs Lj potential
"""
E1 = 0.1
E2 = 0.1
S1 = 0.575625
S2 = 0.575625
combo_rule = 2


"""
Global constants for the 12-10 potential
"""
a = 0.570448629E-01
c = 0.187347694E-01


"""
Combination rule section
"""
if combo_rule == 1:
	print("Haven't implemented this type yet.")
	exit()
elif combo_rule == 2:
	# Lorentz-Berthelot combination rules
	sigma = 0.5 * (S1 + S2)			# arithmetic mean
	epsilon = math.sqrt(E1 * E2)	# geometric mean
else:
	# Geometric mean rule
	sigma = math.sqrt(S1 * S2)		# geometric mean
	epsilon = math.sqrt(E1 * E2)	# geometric mean

print(f"Combined sigma: 	{sigma} nm")
print(f"Combined epsilon:	{epsilon} (kJ/mol)")


"""
Calculation of LJ potentials
"""
switch = 2
r = numpy.arange(0.001, 5, 0.01)

# This is the gromacs potential
LJ = lj(epsilon, sigma, r) 
# This is the 12-10 potential
ca_LJ = ca_lj(a, c, r)


"""
Plotting section
"""

x_axis = numpy.array(r)
# x_axis = numpy.array(r) / sigma

fig, ax = plt.subplots()
ax.plot(x_axis, LJ, 'b-', label="Gromacs Potential")
# ax.plot(x_axis, ca_LJ, 'r', label="12-10 Potential")
plt.ylim(min(LJ)-1 , 4)

plt.xlabel("r, dist between atoms (nm)")
# plt.xlabel("r/sigma, dist between atoms (unitless)")

plt.ylabel("LJ Potential Energy (kJ/mol)")
plt.title(f"Lennard-Jones Potentials")
plt.legend()

plt.show()

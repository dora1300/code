"""
Simple implementation of the JE Shea dihedral potential for studying/didactic purposes.

@Author: DA Ramirez
@Date: 20211007

!!!!! Math.cosine must take *radians* NOT degrees!
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def dihedral(Phi, Phi_0, D, G):
	"""
	Implementation of the Shea standard dihedral potential.
	@Phi 	: the current angle of the dihedral (degrees)
	@Phi_0	: the reference angle of the dihedral (degrees)
	@D		: the first force constant term, kcal/mol
	@G		: the second force constant term, kcal/mol
	
	$returns	: potential, in kJ/mol
	"""
	if type(Phi) is float:
		D_radians = (3*Phi - Phi_0) * (math.pi / 180)
		G_radians = (Phi-Phi_0) * (math.pi / 180)
		potential = D*math.cos(D_radians) - G*math.cos(G_radians)
		return potential
	else:
		potential = []
		for i in range(len(Phi)):
			D_radians = (3*Phi[i] - Phi_0) * (math.pi / 180)
			G_radians = (Phi[i] - Phi_0) * (math.pi / 180)
			function = D*math.cos(D_radians) - G*math.cos(G_radians)
			potential.append(function)
		return potential


phi_0 = 180.
d = 0
g = -2


# phi = float(50.9)
phi = np.arange(0, 360, 0.1)

potentials = dihedral(phi, phi_0, d, g)

if type(phi) is float:
	print(f"{potentials:.3f} kcal/mol")
else:
	fig, ax = plt.subplots()
	ax.plot(phi, np.array(potentials), 'b-', linewidth=0.75)
	ax.set_ylim(min(potentials)-1, max(potentials)+0.1)
	ax.set_xlim(0, 360)
	ax.set_xlabel("Phi (degrees)")
	ax.set_ylabel("Dihedral Potential, V (kcal/mol)")
	plt.show()

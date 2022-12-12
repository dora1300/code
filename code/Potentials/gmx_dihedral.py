"""
Simple implementation of the GROMACS dihedral potential for studying/didactic purposes.

@Author: DA Ramirez
@Date: 20210928

!!!!! Math.cosine must take *radians* NOT degrees!
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def dihedral(Phi, Phi_0, K, N):
	"""
	Implementation of the Gromacs standard dihedral potential.
	@Phi 	: the current angle of the dihedral (degrees)
	@Phi_0	: the equilibrium angle of the dihedral (degrees)
	@K		: the force constant for the dihedral potential (kJ/mol)
	@N		: the multiplicity of the angle term
	
	$returns	: potential, in kJ/mol
	"""
	if type(Phi) is float:
		radians = (N*Phi - Phi_0) * (math.pi / 180)
		potential = K*(1 + math.cos(radians))
		return potential
	else:
		potential = []
		for i in range(len(Phi)):
			radians = (N*Phi[i] - Phi_0) * (math.pi / 180)
			function = K*(1 + math.cos(radians))
			potential.append(function)
		return potential


phi_0 = -130.0

# phi = float(52.0)
phi = np.arange(-360, 360, 0.1)

p1 = dihedral(phi, -130.0, 4, 1)
p2 = dihedral(phi, -30, 4, 3)
p3 = []
for i in range(len(p1)):
    p3.append((p1[i]+p2[i]))

if type(phi) is float:
	print(f"{potentials:.3f} kJ/mol")
else:
    fig, ax = plt.subplots()
    ax.plot(phi, np.array(p1), 'b-', linewidth=0.75,
        label="Multiplicity=1")
    ax.plot(phi, np.array(p2), '-r', linewidth=0.75,
        label="Multiplicity=3")
    ax.plot(phi, np.array(p3), 'g-', linewidth=0.75,
        label="Combined potentials")
    plt.axhline(y=10, color='tab:orange', linestyle="--",
        label="Reference max")
    ax.set_ylim(-0.1, max(p3)+0.1)
    ax.set_xlim(-360, 360)
    plt.legend()
    ax.set_xlabel("Phi (degrees)")
    ax.set_ylabel("Dihedral Potential, V (kJ/mol)")
    plt.show()

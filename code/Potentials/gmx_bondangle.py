"""
Implementation of the GROMACS harmonic bond angle potential.

@Author: DA Ramirez
@Date: 20211109

!!!!! Math.cosine must take *radians* NOT degrees!
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

def harmonicangle(Theta, Theta_0, K):
	"""
	Implementation of the Gromacs harmonic bond angle potential
	@Theta 		: the current bond angle(degrees)
	@Theta_0	: the reference bond angle (degrees)
	@K			: the force constant for the bond angle potential (kJ/mol/rad^2)
	
	$returns	: potential, in kJ/mol
	"""
	t_r = Theta * (m.pi / 180); t0_r = Theta_0 * (m.pi / 180)
	if type(Theta) is float:
		potential = 0.5 * K * (t_r - t0_r)**2
		return potential
	else:
		potential = []
		for i in range(len(Theta)):
			p = 0.5 * K * (t_r[i] - t0_r)**2
			potential.append(p)
		return potential


theta_0 = 92.4
k = 1.000

# phi = float(50.9)
theta = np.arange(0, 360, 0.1)

potentials = harmonicangle(theta, theta_0, k)

if type(theta) is float:
	print(f"{potentials:.3f} kJ/mol")
else:
	fig, ax = plt.subplots()
	ax.plot(theta, np.array(potentials), 'b-', linewidth=0.75)
	ax.set_ylim(-0.1, max(potentials)+0.1)
	ax.set_xlim(0, 360)
	ax.set_xlabel(r"C$\alpha$ bond angle $\theta$ (degrees)")
	ax.set_ylabel("Bond angle potential, V (kJ/mol)")
	ax.set_title("Bond angle potential for linker domain pseudo bond angles")
	plt.show()

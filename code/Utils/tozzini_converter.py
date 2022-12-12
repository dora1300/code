"""
@Name:			Mando A Ramirez
@Date:			2021 11 02

@Description:	This script converts traditional phi/psi dihedral angles into the Tozzini 2006
pseudodihedral/bond angle parameters.

@Citations:		Tozzini V, Rocchia W, and McCammon JA (2006) Mapping All-Atom Models onto One-Bead Coarse-Grained Models:  General Properties and Applications to a Minimal Polypeptide Model. J Chem Theory Comput 2:667–673, American Chemical Society.

@Notes:			Recall that math.cos/sin must take arguments as radians (and outputs in rads)!!

@Updates:
20211116 - added clarity to the functions and code

20211117 - corrected the alpha calculation so that the 360 offset can be calculated in the function itself, which is
essential if I want to convert angles correctly into Tozzini space.
"""

import math as m
import argparse

def theta(T, Y1, Y2, Phi1, Psi1):
	"""
	The structure of this equation is complex. I will shorten it into this so that I can keep better
	track of it
		cos(theta) = cos(T)*A - B + sin(T)*C
	where,
		A = [cos(Y1)cos(Y2) - sin(Y1)sin(Y2)cos(Psi1)cos(Phi1)]
		B = sin(Y1)sin(Y2)sin(Psi1)sin(Phi1)
		c = [cos(Phi1)sin(Y1)cos(Y2) + cos(Psi1)cos(Y1)sin(Y2)]

	All variables are provided in RADIANS. This function assumes (Phi,Psi)_i = (Phi,Psi)_i+1
	:param T: the value of Tau from the paper, but in RADIANS
	:param Y1: the value of y1 from the paper, but in RADIANS
	:param Y2: the value of y2 from the paper, but in RADIANS
	:param Phi1: the value of the Phi angle, in RADIANS
	:param Psi1: the value of the Psi angle, in RADIANS
	:return theta: the theta value in DEGREES
	"""
	A = m.cos(Y1)*m.cos(Y2) - m.sin(Y1)*m.sin(Y2)*m.cos(Psi1)*m.cos(Phi1)
	B = m.sin(Y1)*m.sin(Y2)*m.sin(Psi1)*m.sin(Phi1)
	C = m.cos(Phi1)*m.sin(Y1)*m.cos(Y2) + m.cos(Psi1)*m.cos(Y1)*m.sin(Y2)

	cos_theta = m.cos(T)*A - B + m.sin(T)*C
	theta = m.acos(cos_theta) * (180 / m.pi)
	return theta


def alpha(Y1, Y2, Phi1, Psi1, Phi2, Psi2):
	"""
	All variables are provided as RADIANS. It is lastly converted into degrees after the computation
	step.
	This is the simple expression of the alpha function, and makes the approximation that (Phi,Psi)_i = (Phi,Psi)_i+1

	This can produce an error of 3.7% compared to the analytical solution, and sometimes requires "optimized" values.
	:param Y1: the value of y1 from the paper, but in RADIANS
	:param Y2: the value of y2 from the paper, but in RADIANS
	:param Phi1: the value of the Phi i angle, in RADIANS
	:param Psi1: the value of the Psi i angle, in RADIANS
	:param Phi2: the value of the Phi i+1 angle, in RADIANS
	:param Psi2: the value of the Psi  i+1 angle, in RADIANS
	:return a: the "alpha" value in DEGREES
	"""
	a = (m.pi + Psi1 + Phi2 + Y1*m.sin(Psi2) + Y2*m.sin(Phi1)) * (180 / m.pi)
	if a > 180:
		a = a - 360
	return a

def alpha_complex(T, Y1, Y2, Phi, Psi):
	"""
	This is the complex expression for alpha from the paper, but this expression is the one from the supplement.
	The whole function is broken up into different parts (A - E) so that it is more manageable in the code.

	This can produce an error of 1.4% compared to the analytical solution, but does not require "optimizing" values
	:param T: the value of Tau from the paper, but in RADIANS
	:param Y1: the value of y1 from the paper, but in RADIANS
	:param Y2: the value of y2 from the paper, but in RADIANS
	:param Phi: the value of the Phi angle, in RADIANS
	:param Psi: the value of the Psi angle, in RADIANS
	:return ac: the "alpha complex" value converted into DEGREES
	"""
	A = 0.25 * Y1**2 * m.sin(2*Psi)
	B = 0.25 * Y2**2 * m.sin(2*Phi)
	C = Y1 * Y2 * m.sin(Phi + Psi)
	D = Y1 * (T - (m.pi / 2)) * m.sin(Psi)
	E = Y2 * (T - (m.pi / 2)) * m.sin(Phi)

	ac_r = Phi + Psi + m.pi + (Y1 * m.sin(Psi)) + (Y2 * m.sin(Phi)) + A + B + C - D - E
	ac = ac_r * (180 / m.pi)
	if ac > 180:
		ac = ac - 360
	return ac


"""
GLOBAL ANGLE CONSTANTS -- provided in DEGREES but then converted into radians. These are from the paper
"""
tau = 111
omega = 180
y1 = 20.7
y2 = 14.7			# y2 should not be 20.7, causes the pseudoangles to be too low

t_r = tau * (m.pi / 180)
o_r = omega * (m.pi / 180)
y1_r = y1 * (m.pi / 180)
y2_r = y2 * (m.pi / 180)


def main_func():
	"""
	Set up the argparser
	"""
	parser = argparse.ArgumentParser(description="Converter to transform Phi/Psi"
	    "angles into Tozzini Pseudo-angles")
	parser.add_argument("phi", help="the phi angle, in degrees", type=float)
	parser.add_argument("psi", help="the psi angle, in degrees", type=float)

	args = parser.parse_args()

	phi1, phi2 = args.phi, args.phi
	psi1, psi2 = args.psi, args.psi

	phi1_r = phi1 * (m.pi / 180)
	psi1_r = psi1 * (m.pi / 180)
	phi2_r = phi2 * (m.pi / 180)
	psi2_r = psi2 * (m.pi / 180)


	"""
	Calculate the angles and print!
	"""
	pseudoBond = round(theta(t_r, y1_r, y2_r, phi1_r, psi1_r), 4)
	pseudoDihedral = round(alpha(y2_r, y2_r, phi1_r, psi1_r, phi2_r, psi2_r), 4)
		# why is this? because the citation for this work says to make both sets of phi/psi
		# identical, and to use the equation with y1=y2 and y1=y2=14.7 reproduces the values the best
		# in relation to the paper
	pseudoDihedral_c = round(alpha_complex(t_r, y1_r, y2_r, phi1_r, psi1_r), 4)
		# I am including this because it adds a little extra precision without the y1=y2 assumption. Note, I am using
		# both of the separate y1 and y2 values.


	# if pseudoDihedral > 180:
	# 	pseudoDihedral = pseudoDihedral - 360
	# 	pseudoDihedral_c = pseudoDihedral_c - 360
	# else:
	# 	pass


	print(f"Phi i: {phi1}	|Psi i: {psi1}")
	print(f"Phi i+1: {phi2}	|Psi i+1: {psi2}")
	print(f"Theta: {pseudoBond:.0f}	|Alpha: {pseudoDihedral:.0f}	|Complex Alpha: {pseudoDihedral_c:.0f}")

	return None


if __name__ == "__main__":
	main_func()

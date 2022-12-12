import math as m
import numpy as np
import scipy.stats as st
import pickle
import argparse
import mdtraj as md
import matplotlib.pyplot as plt


print("opening the pickled angles")
location = "/Users/mramirez/Code/AnalysisTools/Tozzini_reference_pickles"
with open(f"{location}/total_helix_alpha.pkl", 'rb') as f:
    helix_alpha = pickle.load(f)
with open(f"{location}/total_helix_theta.pkl", 'rb') as f:
    helix_theta = pickle.load(f)
print("done!")
    
print("making the KDE for the pickled angles")
ref_values = np.vstack([np.array(helix_alpha), np.array(helix_theta)])
refkernel = st.gaussian_kde(ref_values)
print("Done!")

print("writing the kernel object to pickle file")
with open("reference_helix_kernel.pkl", 'wb') as output:
	pickle.dump(refkernel, output, pickle.HIGHEST_PROTOCOL)
print("done!")
	
print("making the list of positions to analyze with the KDE and writing it")
refamin = np.array(helix_alpha).min() ; refamax = np.array(helix_alpha).max()
reftmin = np.array(helix_theta).min() ; reftmax = np.array(helix_theta).max()
refA, refT = np.mgrid[refamin:refamax:1000j, reftmin:reftmax:1000j]
refPositions = np.vstack([refA.ravel(), refT.ravel()])
REF = np.reshape(refkernel.pdf(refPositions).T, refA.shape)

with open("reference_helix_probabilities.pkl", 'wb') as outp:
	pickle.dump(REF, outp, pickle.HIGHEST_PROTOCOL)
print("done!")

print("making the plot of the positions from the KDE")
fig, ax = plt.subplots()
ax.imshow(np.rot90(REF), cmap="Reds", extent=[amin, amax, tmin, tmax], alpha=0.5)
ax.set_xlim(-180, 180);
ax.set_ylim(40, 170)

plt.show()

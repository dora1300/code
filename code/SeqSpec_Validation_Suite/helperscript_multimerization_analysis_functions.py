import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
import mdtraj as md

# def analyze_orientation_torsions(traj, 
#                        coil_length_list,
#                        n_cis_torsion,
#                        c_cis_torsion,
#                        distance_reference_points):
def analyze_angles_and_dihedrals(codename,
                       traj, 
                       coil_length_list):

    # calculate the angle distributions
    # make all the secondary structure backbone angle indices
    # this is generalized so that I can handle multiple types of multimers
    pangle_list = []; ptorsion_list = []

    for len_i, coil_length in enumerate(coil_length_list): # loop over all the coil lengths in the list
        for i in range(coil_length-2):    # loop over the indices corresponding to the current coil_length
            pangle_list.append(
                [i+(len_i*coil_length_list[len_i-1]),
                 i+(len_i*coil_length_list[len_i-1])+1,
                 i+(len_i*coil_length_list[len_i-1])+2])
                # this is convoluted, but it's basically saying that I modify the indices to take into account 
                # the length of the previous coil so that the indices scale continuosly throughout the entire
                # range of atoms in the system

    for len_i, coil_length in enumerate(coil_length_list): # loop over all the coil lengths in the list
        for i in range(coil_length-3):    # loop over the indices corresponding to the current coil_length
            ptorsion_list.append(
                [i+(len_i*coil_length_list[len_i-1]),
                 i+(len_i*coil_length_list[len_i-1])+1,
                 i+(len_i*coil_length_list[len_i-1])+2,
                 i+(len_i*coil_length_list[len_i-1])+3])
                # this is convoluted, but it's basically saying that I modify the indices to take into account 
                # the length of the previous coil so that the indices scale continuosly throughout the entire
                # range of atoms in the system

    pseudoangles = md.compute_angles(traj, np.array(pangle_list), periodic=True) * 180 / np.pi
    avgAngle = np.average(pseudoangles.flatten()); stdAngle = np.std(pseudoangles.flatten())
    pseudotorsions = md.compute_dihedrals(traj, np.array(ptorsion_list), periodic=True) * 180 / np.pi
    avgTorsion = np.average(pseudotorsions.flatten()); stdTorsion = np.std(pseudotorsions.flatten())


    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
    ax1.hist(pseudoangles.flatten(), density=True, color="goldenrod", bins="sqrt", alpha=0.75,
        label=f"{avgAngle:.3f}+/-{stdAngle:.3f}")
    ax2.hist(pseudotorsions.flatten(), density=True, color="grey", bins="sqrt", alpha=0.75,
        label=f"{avgTorsion:.3f}+/-{stdTorsion:.3f}")
    ax1.set_xlabel("pseudoangle (angles)")
    ax1.set_ylabel("counts (density)")
    ax1.grid(color="black", alpha=0.35, linestyle=":")
    ax1.set_xlim(-185, 185)
    ax1.legend()
    ax2.set_xlabel("pseudotorsion (angles)")
    ax2.set_ylabel("counts (density)")
    ax2.grid(color="black", alpha=0.35, linestyle=":")
    ax2.set_xlim(-185, 185)
    ax2.legend()
    plt.tight_layout()
    plt.savefig(f"plot_{codename}_backbone_angles.png", dpi=300)
    plt.close()

    return None



def analyze_orientation_torsions(codename,
                       traj, 
                       n_cis_torsion,
                       c_cis_torsion):
    # Calculate the orientation torsions and plot!
    #   *don't forget to convert to angles*
    # Keep in mind that n_cis_torsion and c_cis_torsion can have multiple arrays which define 
    # torsions, this is applicable for trimers and greater multimers
    zero_time_n_torsions = md.compute_dihedrals(traj[0], np.array(n_cis_torsion), periodic=True) * 180 / np.pi
    zero_time_c_torsions =md.compute_dihedrals(traj[0], np.array(c_cis_torsion), periodic=True) * 180 / np.pi
    n_orient_torsions = md.compute_dihedrals(traj, np.array(n_cis_torsion), periodic=True) * 180 / np.pi
    c_orient_torsions = md.compute_dihedrals(traj, np.array(c_cis_torsion), periodic=True) * 180 / np.pi

    LEN_TOR = len(n_cis_torsion)
    fig, axs = plt.subplots(ncols=LEN_TOR, nrows=1, figsize=(5*LEN_TOR, 5))
    if LEN_TOR > 1:
        for I in range(len(n_cis_torsion)):
            axs[I].hist(n_orient_torsions[:,I], density=True, color="mediumturquoise", bins="sqrt", alpha=0.75)
            axs[I].set_xlabel("N-terminal torsion (angles)")
            axs[I].set_ylabel("counts (density)")
            axs[I].set_xlim(-185, 185)
            # axs[I].set_title(f"coil{I+1}--coil{I+2}")
            axs[I].grid(color="black", alpha=0.35, linestyle=":")
        if LEN_TOR > 1:
            axs[1].set_title(f"coil{1}--coil{3}")
        axs[0].set_title(f"coil{1}--coil{2}")
    else:
        axs.hist(n_orient_torsions, density=True, color="mediumturquoise", bins="sqrt", alpha=0.75)
        axs.set_xlabel("N-terminal torsion (angles)")
        axs.set_ylabel("counts (density)")
        axs.set_xlim(-185, 185)
        # axs[I].set_title(f"coil{I+1}--coil{I+2}")
        axs.grid(color="black", alpha=0.35, linestyle=":")
    fig.suptitle("N-terminal orientation torsions")
    plt.tight_layout()
    plt.savefig(f"plot_{codename}_Nterminal_orientation_torsions.png", dpi=300)
    plt.close()

    LEN_TOR = len(c_cis_torsion)
    fig, axs = plt.subplots(ncols=LEN_TOR, nrows=1, figsize=(5*LEN_TOR, 5))
    if LEN_TOR > 1:
        for I in range(len(c_cis_torsion)):
            axs[I].hist(c_orient_torsions[:,I], density=True, color="mediumturquoise", bins="sqrt", alpha=0.75)
            axs[I].set_xlabel("N-terminal torsion (angles)")
            axs[I].set_ylabel("counts (density)")
            axs[I].set_xlim(-185, 185)
            # axs[I].set_title(f"coil{I+1}--coil{I+2}")
            axs[I].grid(color="black", alpha=0.35, linestyle=":")
        if LEN_TOR > 1:
            axs[1].set_title(f"coil{1}--coil{3}")
        axs[0].set_title(f"coil{1}--coil{2}")
    else:
        axs.hist(c_orient_torsions, density=True, color="mediumturquoise", bins="sqrt", alpha=0.75)
        axs.set_xlabel("N-terminal torsion (angles)")
        axs.set_ylabel("counts (density)")
        axs.set_xlim(-185, 185)
        # axs[I].set_title(f"coil{I+1}--coil{I+2}")
        axs.grid(color="black", alpha=0.35, linestyle=":")
    fig.suptitle("C-terminal orientation torsions")
    plt.tight_layout()
    plt.savefig(f"plot_{codename}_Cterminal_orientation_torsions.png", dpi=300)
    plt.close()

    return None



def analyze_reference_point_distances(codename,
                       traj, 
                       distance_reference_points):
    # Calculate the distances of the reference pairs
    reference_distances = md.compute_distances(traj, 
                            np.array(distance_reference_points),
                            periodic=True)
    
    # keep in mind there can be multiple sets of reference points for multimers greater than a dimer
    # that correspond to the reference poitns between more than two coils
    # But there will always be a multiple of 4 reference points if this is done correctly
    num_unique_ref_points = int(len(distance_reference_points)/4)

    for i in range(num_unique_ref_points):
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, figsize=(14, 4))
        ax1.hist(reference_distances[:,(0+(i*4))], density=True, color="grey", alpha=0.75, bins="sqrt")
        ax1.set_title("termini 1")
        ax1.grid(color="black", alpha=0.35, linestyle=":")

        ax2.hist(reference_distances[:,(1+(i*4))], density=True, color="grey", alpha=0.75, bins="sqrt")
        ax2.set_title("midpoint 1")
        ax2.grid(color="black", alpha=0.35, linestyle=":")

        ax3.hist(reference_distances[:,(2+(i*4))], density=True, color="grey", alpha=0.75, bins="sqrt")
        ax3.set_title("midpoint 2")
        ax3.grid(color="black", alpha=0.35, linestyle=":")

        ax4.hist(reference_distances[:,(3+(i*4))], density=True, color="grey", alpha=0.75, bins="sqrt")
        ax4.set_title("termini 2")
        ax4.grid(color="black", alpha=0.35, linestyle=":")

        fig.supxlabel("distance (nm)")
        fig.supylabel("counts (density)")
        plt.grid(color="black", alpha=0.35, linestyle=":")
        fig.suptitle(f"distance of reference points between coil{i+1}--coil{1+2}")
        if num_unique_ref_points == 2:
            fig.suptitle(f"coil{1}--coil{3}")
        fig.suptitle(f"coil{1}--coil{2}")
        plt.tight_layout()
        plt.savefig(f"plot_{codename}_key_ref_distances_{i}.png", dpi=300)
        plt.close()

    return None
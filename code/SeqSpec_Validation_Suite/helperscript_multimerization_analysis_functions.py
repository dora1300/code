import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
import mdtraj as md


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
    n_orient_torsions = md.compute_dihedrals(traj, np.array(n_cis_torsion), periodic=True) * 180 / np.pi
    c_orient_torsions = md.compute_dihedrals(traj, np.array(c_cis_torsion), periodic=True) * 180 / np.pi

    zero_time_n_torsions = n_orient_torsions[0] # converts to degrees!
    zero_time_c_torsions = c_orient_torsions[0]


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



def helper_are_distances_a_real_multimer(dist_arr, cutoff, acceptable_frac_dists):
    """
    :args:
        dist_arr
        cutoff
        acceptable_frac_dists
    """
    acceptable = 0
    unacceptable = 0
    for distance in dist_arr:
        if distance <= cutoff:
            acceptable += 1
        else:
            unacceptable += 1
    
    if (acceptable / (acceptable + unacceptable)) >= acceptable_frac_dists:
        return True
    else:
        return False
    

def helper_are_all2all_distances_a_multimer(all2all_dists, cutoff, acceptable_frac_dists):
    beads_that_are_making_contacts = 0
    beads_not_making_contacts = 0
    for beadi_distances in all2all_dists:
        if ((beadi_distances <= cutoff)).any():
            beads_that_are_making_contacts += 1
        else:
            beads_not_making_contacts += 1
    
    if (beads_that_are_making_contacts / (beads_that_are_making_contacts + beads_not_making_contacts)) >= acceptable_frac_dists:
        return True
    else:
        return False



def analyze_reference_point_distances(codename,
                       traj, 
                       distance_reference_points,
                       core_cutoff:float=0.9,
                       acceptable_threshold:float=0.75,
                       mismatch_cutoff:float=1.2):
    """
    :args:
        individual_coil_lengths : (list)
        core_cutoff : [nm]
            the distance cut-off for determining if a given set of reference points is 
            participating in a multimer. Only counts for reference points that are not defined as the termini
        acceptable_threshold :

        com_cutoff : [nm]

    """

    # Calculate the distances of the reference pairs
    reference_distances = md.compute_distances(traj, 
                            np.array(distance_reference_points),
                            periodic=True)
        # this calculates the distances for all reference points for ALL the frames in the provided trajectory


    # for me to do this analysis correctly and robustly, I need to separate out the reference beads into
    # arrays so that all beads within the same coil stay together
    # rows = the number of pairs of beads used to define reference points between coils
    # columns = the number of coils (which is the length of the distance_reference_points array)
    # 2025 10 29 == this analysis is only suitable for dimers as of right now!
    ref_points_by_coil = np.zeros((len(distance_reference_points), 2), dtype=int)
    for bead_pairs_i, bead_pairs in enumerate(distance_reference_points):
        ref_points_by_coil[bead_pairs_i][0] = bead_pairs[0]
        ref_points_by_coil[bead_pairs_i][1] = bead_pairs[1]
    
    # I am going to store my multimer analysis into a numpy array
    # each row is a frame of the trajectory.
    # Column 1 = correct multimer; Column 2 = bound but incorrect multimer; Column 3 = not bound at all
    multimer_analysis_table = np.zeros((traj.n_frames, 3))

    for frame_i, frame_ref_dists in enumerate(reference_distances):
        # Do the reference distances suggest a proper multimer?
        if helper_are_distances_a_real_multimer(frame_ref_dists, core_cutoff, acceptable_threshold):
            multimer_analysis_table[frame_i, 0] = 1
            continue
        else:
            # calculate the distances between all the pairs. I gotta see if there is an interaction of any
            # type happening
            # This step gets the xyz coordinates for the reference points, which is important for getting the distances
            # in a vectorized way
            frame_coil1_ref_points = traj.xyz[frame_i, ref_points_by_coil[0], :]
            frame_coil2_ref_points = traj.xyz[frame_i, ref_points_by_coil[1], :]

            # this calculates the distances of all points to each other using the vectorized method
            coil1_to_coil2_ref_point_dists = np.linalg.norm(frame_coil1_ref_points - frame_coil2_ref_points[:, None], axis=-1).T
            #   STOP!!
            #   PLEASE REMEMBER!!
            #   Using the above method -- I need to take the transpose before I do any analysis!!
            # this checks to see if any of the coils are making a contact that could construed as a multimer.
            if helper_are_all2all_distances_a_multimer(coil1_to_coil2_ref_point_dists.T, mismatch_cutoff, acceptable_threshold):
                multimer_analysis_table[frame_i, 1] = 1
            else:
                multimer_analysis_table[frame_i, 2] = 1


    # at this point, save out the multimer analysis table so I can do things with it
    np.savetxt(f"{codename}_areCoilsInMultimer_by_frame.csv", multimer_analysis_table, delimiter=",", fmt="%i")

    # also make a plot of the fraction of the type of multimer that coils are in, based on the table above
    fraction_of_frames_in_multimer = np.sum(multimer_analysis_table.T[0]) / traj.n_frames
    fraction_of_frames_in_wrong_multimer = np.sum(multimer_analysis_table.T[1]) / traj.n_frames
    fraction_of_frames_no_multimer = np.sum(multimer_analysis_table.T[2]) / traj.n_frames

    fig, axs = plt.subplots(figsize=(5, 4))
    plt.bar(np.array([1]), np.array([fraction_of_frames_in_multimer]),
                color="goldenrod", width=0.65)
    plt.bar(np.array([2]), np.array([fraction_of_frames_in_wrong_multimer]),
                color="black", width=0.65)
    plt.bar(np.array([3]), np.array([fraction_of_frames_no_multimer]),
                color="grey", width=0.65),
    plt.xticks([1, 2, 3], [f"Correctly\nbound", f"Incorrectly\nbound", "Unbound"], rotation=50)
    plt.ylim(0, 1)
    plt.ylabel("fraction of simulation frames")
    plt.tight_layout()
    plt.savefig(f"plot_{codename}_areCoilsInMultimer.png", dpi=300)
    plt.close()

    
    # keep in mind there can be multiple sets of reference points for multimers greater than a dimer
    # that correspond to the reference poitns between more than two coils
    # But there will always be a multiple of 4 reference points if this is done correctly
    # TODO -- 2025 10 29, I need to fix this to accurately handle trimers etc. I need a new way of storing my
    # reference points
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

    return fraction_of_frames_in_multimer, fraction_of_frames_in_wrong_multimer




def analyze_ETE_distances(codename,
                       traj,
                       list_of_protein_lengths,
                       ETE_error_threshold:float=0.80):
    """
    :args:
        - list_of_protein_lengths : [list]
            a list of 1-indexed lengths (i.e. total beads) in the coils in the system
            The numbers are in order of the coils present in the system.
    """
    # first I will be generating the list of indices that correspond to the starts and ends
    # of the coils
    # don't forget = mdtraj = 0-indexed
    coil_starts_and_ends = []
    start_tracker = 0
    for length in list_of_protein_lengths:
        end_of_coil = start_tracker + (length - 1) # -1 to account for the 0-indexing
        coil_starts_and_ends.append([start_tracker, end_of_coil])
        start_tracker += length
    
    # now I can calculate the ETEs
    ETEs = md.compute_distances(traj, np.array(coil_starts_and_ends))

    # The ETEs themselves are interesting, but I also want some metric of how the ETEs change througout the simulation
    coil1_initial_ete = ETEs[0][0]; coil2_initial_ete = ETEs[0][1]
    #   this gives each frame's ETE as a fraction of the initial, which is useful for seeing the percent/fractional deviation
    #   from expected
    coil1_ete_fraction_of_init = ETEs[:, 0] / coil1_initial_ete
    coil2_ete_fraction_of_init = ETEs[:, 1] / coil2_initial_ete
    #   this turns the fractions into 1s and 0s corresponding to correct and incorrect on a frame basis
    coil1_fraction_frame_correct = np.where(coil1_ete_fraction_of_init < ETE_error_threshold, 0, 1)
    coil2_fraction_frame_correct = np.where(coil2_ete_fraction_of_init < ETE_error_threshold, 0, 1)

    # now save my data and make a plot
    np.savetxt(f"{codename}_ETEs.csv", ETEs, delimiter=",", fmt="%1.4f")
    np.savetxt(f"{codename}_ETEs_fraction_of_initial.csv", 
               np.array([coil1_ete_fraction_of_init, coil2_ete_fraction_of_init]), delimiter=",", fmt="%1.4f")
    np.savetxt(f"{codename}_ETEs_fraction_frame_correct.csv", 
               np.array([coil1_fraction_frame_correct, coil2_fraction_frame_correct]), delimiter=",", fmt="%1.4f")


    coil1_frac_correct = np.sum(coil1_fraction_frame_correct)/traj.n_frames
    coil2_frac_correct = np.sum(coil2_fraction_frame_correct)/traj.n_frames

    fig, ax = plt.subplots(figsize=(5, 3))
    plt.bar(np.array([1, 4]), np.array([coil1_frac_correct, 1-coil1_frac_correct]),
            label="Coil 1", color="goldenrod", width=0.75, align="center")
    plt.bar(np.array([2, 5]), np.array([coil2_frac_correct, 1-coil2_frac_correct]),
            label="Coil 2", color="grey", width=0.75, align="center")
    
    plt.xticks([1.5, 4.5], [f"correct ETE", "incorrect ETE"], rotation=30)
    plt.ylabel(f"fraction of \nsimulation frames")
    plt.tight_layout()
    plt.savefig(f"plot_{codename}_ETE_fraction_frame.png", dpi=300)
    plt.close()

    return coil1_frac_correct, coil2_frac_correct

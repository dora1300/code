"""
@Title:             _pulling_simulation_preparer.py
@Name:              Mando A Ramirez
@Date:              2025 07 08

@Description:       This script is the main "director" for setting up the files necessary to run a simulation
for the sequence specific validations. This file prepares **a single multimer simulation** for ***COM PULLING
SIMULATIONS***.

@Updates:
2025 08 21          This is a significant update! I am completely revamping this code so that this
            code actually *runs* the MD simulations and validation simulations

"""

import argparse
import os
import shutil
import subprocess
import helperscript_level1_multimerization_dictionary as mult_dict
import write_protein_itp_file as write_prot_itp
import write_nonbond_params_for_top as write_nonbonded
import write_pairtypes_for_top as write_pairtypes
import write_atomtypes_for_top as write_atomtypes
import write_seqspec_top_file as write_top
import importonly_write_com_pulling_mdp_files as mdp_writer
import helperscript_multimerization_analysis_functions as analyzer
import mdtraj as md
import numpy as np
# import numpy as np
# import matplotlib.pyplot as plt



if __name__ == "__main__":
    # for running this script solo, or it can be run from a separate python script if needed 
    # for even more automated flow.

    """
    X -- Handle the argument parsing to get everything set up.
    SORRY IN ADVANCE THIS SECTION IS THICCC
    """
    parser = argparse.ArgumentParser(description="Hello! This code is part of the Sequence Specific Force Field validation "
                                     "suite, in which you can set up a COM pulling simulation between coiled-coils "
                                     "in a box.")
    parser.add_argument("-override", action="store_true", default=False, 
                        help="I'm not describing this. If you know why and how to use then then feel free. "
                        "Otherwise, don't touch.")
    
    parser.add_argument("-rerun_analysis", help="[T/F] [default: False] switch. Toggle this if you want to just rerun the analysis "
                        "on a finished simulation. Should save some time if you need to reanalyze.",
                        action="store_true", default=False)

    parser.add_argument("-reconvert_trajectories", help="[T/F] [default: False] switch. Toggle this if you need to reconvert "
                        "the trajectories for any reason. This command will only make sense if you also pass "
                        "-rerun_analysis.",
                        action="store_true", default=False)

    #
    #   The basic general arguments
    #
    parser.add_argument("-simulation_path", help="OMIT FINAL SLASH -- The absolute path for where you wish to run the validation "
                        "simulation. This is the head directory where a *NEW* directory will be made for your "
                        "protein of interest.", 
                        type=str, required=True)
    
    parser.add_argument("-protein_codename", help="[Required always] The *codename* of the protein you wish to simulate, but this "
                        "does not necessarily correlate to an actual name that would go into a topology file.",
                        type=str)

    parser.add_argument("-run_on_alpine", help="[T/F] [default: False] Switch: is this being run locally on a laptop, or on Alpine?",
                        action="store_true", default=False)

    parser.add_argument("-simulation_T", help="[default: 298] The temperature at which to run the COM pulling production simulation.",
                        type=float, default=298.0)
    
    parser.add_argument("-system_name", help="[default: 'coiled-coil-system'] The GROMACS system name that goes "
                        "into the topology file. Can use default.", default="coiled-coil-system", type=str)
    
    parser.add_argument("-nt", help="[default: 1] The number of threads to assign to running the simulation.", 
                        default=1, type=str)
    

    #
    #   Specific arguments for the protein.itp file!
    #
    parser.add_argument("-nrexcl", help="[default: 3] The value to set for the `nrexcl` component of the protein.itp file. " \
                        "This is relevant for determining pairs and excluded nonbonded interactions.",
                        type=str, default="3")

    parser.add_argument("-use_og_pair_terms", help="[T/F] [default: False] switch. Activate this switch (!) to use the 1-4 AND "\
                        "1-5 pair terms that come from my OG CC-LLPS simulation framework. The OG pairs were used to help stabilize "\
                        "the backbone helix but can be removed as desired. This will of course overwrite nonbonded interactions "\
                        "for 1-4 and 1-5 regardless of what nrexcl is set to. **These pairs were developed with nrexcl=4 just FYI**",
                        action="store_true", default=False)

    #
    #   Specific arguments for specifying nonbonded information for nonbonded.itp
    #
    parser.add_argument("-name_of_nonbonded_file", help="[default: itp_nonbonded.itp] The name that you'd like to give to the nonbonded_params.itp. "
                        "Please include the extension and the path if appropriate!",
                        default="itp_nonbonded.itp")
    
    parser.add_argument("-epsilon_matrix_file", help="[csv] the 20x20 matrix/CSV file containing the "
                        "epsilons for the itp file. PLEASE GIVE THE ABSOLUTE PATH")
    
    parser.add_argument("-sigma_matrix_file", help="[csv] Name (and path) of the file containing the 20x20 "
                        "matrix of sigma values to use. PLEASE GIVE THE ABSOLUTE PATH")
    
    parser.add_argument("-combining_rule", help="[default: 'lorentz']The code-name of the combining rule to use for calculating "
                        "sigmas. Only 'lorentz' is supported right now.",
                        default="lorentz")
    
    parser.add_argument("-epsilon_scale_factor", help="[default: 1.0] A global factor by which to scale the epsilons provided in the "
                        "epsilon matrix. Each epsilon is scaled multiplicatively by the scale factor, so please "
                        "provide fractional scale factors. e.g. 1 = no change, 0.90 = 90%% of regular "
                        "epsilons, 1.25 means 25%% stronger than reference epsilons etc.", default=1.0,
                        type=float)
    
    parser.add_argument("-sigma_scale_factor", help="[default: 1.0] A global factor by which to scale the sigmas provided in the "
                        "sigma matrix. Each sigma is scaled multiplicatively by the scale factor, so please "
                        "provide fractional scale factors. e.g. 1 = no change, 0.90 = 90%% of regular "
                        "sigmas, 1.25 means 25%% stronger than provided sigmas etc.", default=1.0,
                        type=float)
    
    #  
    #   Specific arguments for specifying the atomtypes.itp information
    #
    parser.add_argument("-name_of_atomtypes_file", help="[default: itp_atomtypes.itp] The name that you'd like to "
                        "give to the atomtypes.itp file! "
                        "Please include the extension and the path if appropriate!",
                        default="itp_atomtypes.itp")
    
    parser.add_argument("-atomtypes_attr_term", help="[default: 6.9132E-03] The generic attractive term for the C6 term in the atomtypes.itp "
                        "file. This is mostly useless but must be provided in the topology.",
                        default=6.9132E-03, type=float)
    
    parser.add_argument("-atomtypes_repul_term", help="[default: 5.97404E-05] The generic repulsive term for the C12 term in the atomtypes.itp "
                        "file. This is mostly useless but must be provided in the topology.",
                        default=5.97404E-05, type=float)
    

    #  
    #   Specific arguments for specifying the pairtypes.itp information
    #
    parser.add_argument("-include_pairtypes_file", help="[T/F] [default: False] switch. Toggle this if you want to include a separate "\
                        "pairtypes.itp file in the topology.",
                        action="store_true", default=False)
    
    parser.add_argument("-name_of_pairtypes_file", help="[default: itp_pairtypes.itp] The name that you'd like to "
                        "give to the pairtypes.itp file! "
                        "Please include the extension and the path if appropriate!",
                        default="itp_pairtypes.itp")
    
    parser.add_argument("-epsilon_14_scale_factor", help="[default: 1.0] A factor by which to multiplicatively scale the epsilons "\
                        "used in the pairtypes.itp file, which will affect ONLY the 1-4 pairs if gen-pairs = yes", default=1.0,
                        type=float)
    
    parser.add_argument("-sigma_14_scale_factor", help="[default: 1.0] A factor by which to multiplicatively scale the sigmas "\
                        "used in the pairtypes.itp file, which will affect ONLY the 1-4 pairs if gen-pairs = yes", default=1.0,
                        type=float)


    #  
    #  Specific arguments for specifying the total.top information
    #        
    parser.add_argument('-name_of_topology_file', help="[default: top_cc.top] Name that you want for the "
                        "overall topology (.top) file.", 
                        type=str, default="top_cc.top")
    
    parser.add_argument("-genpairs", help="[default: 'yes'] The value to pass to the `gen-pairs` component of the .top file. " \
                        "Be mindful of how you use this if you choose to provide additional pairtypes.itp files.",
                        default="yes", type=str)
    
    parser.add_argument("-fudgelj", help="[default: 1.0] The numeric value to pass to `fudgeLJ` term in the general .top "\
                        "file. This will be used to scale generated 1-4 pairs when `gen-pairs`='yes'",
                        default=1.0, type=float)
    

    #
    #   This is an interesting one -- options for specifying different parameters in the .mdp files
    #       You shouldn't really need to change these values, 
    #       but they can be changed if there is a good reason to do so.
    #
    parser.add_argument('-pulling_coordinate_distance', help="[default: 2.0] The distance at which the COM "
                        "pulling starts and the potential well kicks in [nm]. ",
                        type=float, default=2.0)
    

    #
    #   Arguments for changing different analysis parameters. You don't have to use this but can if you want
    #
    parser.add_argument('-cutoff_for_correct_multimerization', help="[default: 0.9] [nm] The distance between reference beads " \
                        "to classify a correct multimer.",
                        type=float, default=0.9) 
    parser.add_argument('-cutoff_for_wrong_multimerization', help="[default: 0.1.2] [nm] The distance between reference beads " \
                        "to classify a multimer that is wrongly aligned but still a multimer.",
                        type=float, default=1.2) 
    parser.add_argument('-threshold_for_multimerization', help="[default: 0.75] The fraction of reference beads that " \
                        "need to be within the provided cut-offs for structure to be called a multimer. This can be for correct " \
                        "or in correct multimer arrangements.",
                        type=float, default=0.75)
    parser.add_argument('-threshold_for_anycontact', help="[default: 0.5] The fraction of reference beads that " \
                        "need to be within the provided cut-offs for structure to be classified as any type of contact. "\
                        "This threshold will only be used if a given structure is already determined to not be 'correct'.",
                        type=float, default=0.5)

    parser.add_argument('-ete_tolerance_threshold', help="[default: 0.8] The threshold at which a given coil's ETE relative to " \
                        "starting frame can still be considered a tolerable ETE, " \
                        "e.g. if ETE_frame <= (ETE_initial * ete_tolerance_threshold) then the coil does not have correct ETE.",
                        type=float, default=0.8)



    """
    X -- Parse my args and get everything set up that I need to!
    """
    # 1- parse the args
    args = parser.parse_args()
    rerun_switch = args.rerun_analysis

    # 2- set up my environment variable that will hold all the information relevant for the protein codename I want
    #   to simulate. This information comes from the dictionary in helperscript_level1_...
    ENV_INFO = mult_dict.validation_dictionary[args.protein_codename]
    list_of_protein_names = ENV_INFO[0]
    list_of_ncoils = ENV_INFO[1]

    # 3- set some other variables that will be important
    WORKDIR = args.simulation_path
    # NCOILS = args.number_of_coils
    TEMP = args.simulation_T

    # 4- set Global variables, non-arguments
    if args.run_on_alpine:
        CODEDIR = "/projects/dora1300/code/code/SeqSpec_Validation_Suite"
    else:
        CODEDIR = "/home/mando/code/code/SeqSpec_Validation_Suite"
        # this houses all the code and directories necessary 
        # specific for my linux distribution

    # 5- define the path of the working simulation directory specific to this protein
    # PROTNAME will be a combination of all the proteins that are getting simulated
    pro_sim_dir = f"{WORKDIR}/{args.protein_codename}"   # this is the main directory that holds all the simulation things
                    # for the given protein validation simulation

    # not False = True
    # not True = False
    # if not rerun_switch(False) == if True
    # if not rerun_switch(True) == False
    if not rerun_switch:    # basically, if not True
        """
        X -- Set up the directories and necessary workspace for running the simulation. This involves copying all
        the relevant files from their storage location.
        
        Some other copying will take place elsewhere. This covers the starting structure portion.
        """
        # I am going to include a safety feature to stop me from overwriting things if I make
        # a mistake and provdie a path that already exists:
        if os.path.exists(f"{pro_sim_dir}"):
            if args.override:
                pass
            else:
                raise FileExistsError("The provided simulation path already contains a simulation folder for the specified protein. "
                "Please provide a different path to prevent ruining existing data OR override this feature.")
        else:
                os.mkdir(f"{pro_sim_dir}")
                os.chdir(f"{pro_sim_dir}")
                os.mkdir("top"); os.mkdir("mdp"); os.mkdir("em"); os.mkdir("md")
                os.mkdir(f"analysis_{args.protein_codename}")
                os.mkdir("starting_structure")

        # updated 20250821 - simpler (I think) way of copying the starting box_structure and the index files
        os.chdir(f"{pro_sim_dir}/starting_structure")
        shutil.copy(f"{CODEDIR}/structures_validation_multimerization/level1/box_{args.protein_codename}.gro",
                    f"{pro_sim_dir}/starting_structure")
        shutil.copy(f"{CODEDIR}/structures_validation_multimerization/level1/index_{args.protein_codename}.ndx",
                    f"{pro_sim_dir}/starting_structure")        


        """
        X -- Generate the protein .itp and .top files for the provided protein
            This code does not do anything with secondary structure prediction, or converting
            secondary structure predictions into mapped backbone scalar values. That has to be done
            elsewhere and beforehand.
            This includes the atomtypes, nonbonded, and the actual protein .itp, as well as the total .top

            This is all done in the ./top folder
        """
        os.chdir(f"{pro_sim_dir}/top")

        #
        #   Copy the relevant (.fasta, backbone) files from the validation suite into this directory
        #
        for NAME in list_of_protein_names:
            shutil.copy(f"{CODEDIR}/structures_validation_multimerization/level1/fasta_{NAME}.fasta",
                        f"{pro_sim_dir}/top/")
            shutil.copy(f"{CODEDIR}/structures_validation_multimerization/level1/backbones_{NAME}.csv",
                        f"{pro_sim_dir}/top/")


        #
        #   Generate itp_nonbonded.itp
        #   this is the sequence specific parameters!!!
        #
        # Step 1 -- calculate the sigmas based on the setting I chose
            # 2025 08 24 -- I am turning off the ability to do anything other than generate
            # sigmas from a matrix file. But I'm not changing the underlying function. I'm just
            # default passing in parameters which only allow generating sigmas from a matrix file.
        array_of_sigmas = write_nonbonded.load_calculated_sigmas(True, args.sigma_matrix_file,
                                                False, None,
                                                False, args.combining_rule)
        
        # Step 2 -- generate the text that will make up the file contents of itp_nonbonded
        nonbond_file_contents = write_nonbonded.write_nonbonded_itp(args.epsilon_matrix_file, array_of_sigmas, 
                                                args.epsilon_scale_factor, args.sigma_scale_factor)

        # Step 3 -- save the itp file to disk!
        write_nonbonded.save_itp_file_to_disk(nonbond_file_contents, args.name_of_nonbonded_file)



        #
        #   Generate itp_paritypes.itp
        #   this is the sequence specific parameters used for handling 1-4 pairs!!
        #
        #   But only do this if I pass the argument specifying to do this.
        #
        if args.include_pairtypes_file:
            array_of_pairtype_sigmas = write_pairtypes.load_calculated_sigmas(True, args.sigma_matrix_file,
                                                    False, None,
                                                    False, args.combining_rule)
            
            pairtype_file_contents = write_pairtypes.write_pairtypes_itp(args.epsilon_matrix_file, array_of_pairtype_sigmas, 
                                                    args.epsilon_14_scale_factor, args.sigma_14_scale_factor)

            # Step 3 -- save the itp file to disk!
            write_pairtypes.save_itp_file_to_disk(pairtype_file_contents, args.name_of_pairtypes_file)



        #
        #   Generate itp_atomtypes.itp
        #
        write_atomtypes.write_atomtypes_itp(args.name_of_atomtypes_file, 
                                            attractive_term=args.atomtypes_attr_term, 
                                            repulsive_term=args.atomtypes_repul_term)
        

        #
        #   Generate protein.itp. This will have to be done for each of the proteins provided!
        #       .fasta and the .backbone file have to be copied into this directory
        #
        list_of_itp_files = []
            # this list might seem a little convoluted. I'm not giving the user the choice to name the itp files.
            # instead, they give me a list of protein names then I turn those names into corresponding .itp file names.
            # I need these names later, though, which is why I'm saving them
        for prot_index, protein_name in enumerate(list_of_protein_names):
            input_protein_name, input_protein_sequence = write_prot_itp.fasta_parser(f"fasta_{protein_name}.fasta")
            input_protein_backbones = write_prot_itp.backbone_file_parser(f"backbones_{protein_name}.csv")
            if len(input_protein_sequence) != len(input_protein_backbones):
                raise ValueError(f"The length of the sequence you gave in the .fasta file ({len(input_protein_sequence)} aa) " \
                                f"does not match the number of mapped backbone scalar values ({len(input_protein_backbones)}) " \
                                "given in the backbone file. Please check your files and try again.")
            write_prot_itp.write_itp_text(input_protein_name, input_protein_sequence, 
                                        protein_name, input_protein_backbones, 
                                        f"itp_{protein_name}.itp", args.nrexcl, args.use_og_pair_terms)
            list_of_itp_files.append(f"itp_{protein_name}.itp")


        # 
        #    Generate the protein .top file.
        #       This is where sequence specificity happens i.e. where I incorporate my epsilon and sigma parameters!
        #
        #     write_text_top(ATOMTYPE_NAME, NONBONDED_NAME, ITP_NAME, MOLECULETYPE_NAME, NUM_MOLTYPE, SYS_NAME, OUTPUT, GENPAIRS, FUDGELJ)
        write_top.write_text_top(args.name_of_atomtypes_file, args.name_of_nonbonded_file ,
                                list_of_itp_files, list_of_protein_names, list_of_ncoils,
                                args.system_name, args.name_of_topology_file, args.genpairs, args.fudgelj,
                                include_pairtypes=args.include_pairtypes_file, 
                                pairtypes_name=args.name_of_pairtypes_file)




        """
        X -- Generate .mdp files for the COM pulling simulations
            The EM mdp files can be universal, but the actual NVT MD simulations with the COM pulling will need
            to be handled carefully so that I get the correct number of pulling index groups
        """
        os.chdir(f"{pro_sim_dir}/mdp")
        mdp_writer.write_em_mdp_text()
        mdp_writer.write_md_pulling_mdp_text(TEMP, ENV_INFO[2], pull_init=args.pulling_coordinate_distance)
        


        """
        X -- Run energy minimization! (located in ../em)
        """
        # first up = EM
        os.chdir(f"{pro_sim_dir}/em")
        if args.run_on_alpine:
            em_grompp_text = f"mpirun -np 1 gmx_mpi grompp -f ../mdp/em.mdp " \
                f"-c ../starting_structure/box_{args.protein_codename}.gro " \
                f"-p ../top/{args.name_of_topology_file} -n ../starting_structure/index_{args.protein_codename}.ndx " \
                f"-o em_{args.protein_codename}.tpr"
        else:
            em_grompp_text = f"gmx grompp -f ../mdp/em.mdp " \
                f"-c ../starting_structure/box_{args.protein_codename}.gro " \
                f"-p ../top/{args.name_of_topology_file} -n ../starting_structure/index_{args.protein_codename}.ndx " \
                f"-o em_{args.protein_codename}.tpr"
        em_grompp_cmd = em_grompp_text.split()
        em_grompp_process = subprocess.Popen(em_grompp_cmd)
        em_grompp_process.wait()

        if args.run_on_alpine:
            em_md_text = f"mpirun -np 1 --mca opal_common_ucx_opal_mem_hooks 1 gmx_mpi mdrun " \
                f"-ntomp {args.nt} -s em_{args.protein_codename}.tpr -deffnm {args.protein_codename}_em"
        else:
            em_md_text = f"gmx mdrun " \
                f"-nt {args.nt} -pin on -s em_{args.protein_codename}.tpr -deffnm {args.protein_codename}_em"
        em_md_cmd = em_md_text.split()
        em_md_process = subprocess.Popen(em_md_cmd)
        em_md_process.wait()



        """
        X -- Run the md simulation now! (located in ../md)
        """
        # second step = MD
        os.chdir(f"{pro_sim_dir}/md")
        if args.run_on_alpine:
            md_grompp_text = f"mpirun -np 1 gmx_mpi grompp -f ../mdp/md_pulling.mdp " \
                f"-c ../em/{args.protein_codename}_em.gro " \
                f"-p ../top/{args.name_of_topology_file} -n ../starting_structure/index_{args.protein_codename}.ndx " \
                f"-o md_{args.protein_codename}.tpr"
        else:
            md_grompp_text = f"gmx grompp -f ../mdp/md_pulling.mdp " \
                f"-c ../em/{args.protein_codename}_em.gro " \
                f"-p ../top/{args.name_of_topology_file} -n ../starting_structure/index_{args.protein_codename}.ndx " \
                f"-o md_{args.protein_codename}.tpr"
        md_grompp_cmd = md_grompp_text.split()
        md_grompp_process = subprocess.Popen(md_grompp_cmd)
        md_grompp_process.wait()      

        if args.run_on_alpine:
            md_md_text = f"mpirun -np 1 --mca opal_common_ucx_opal_mem_hooks 1 gmx_mpi mdrun " \
                f"-ntomp {args.nt} -rdd 1.6 " \
                f"-s md_{args.protein_codename}.tpr -deffnm {args.protein_codename}_md"
        else:
            md_md_text = f"gmx mdrun " \
                f"-nt {args.nt} -pin on -rdd 1.6 " \
                f"-s md_{args.protein_codename}.tpr -deffnm {args.protein_codename}_md"  
        md_md_cmd = md_md_text.split()
        md_md_process = subprocess.Popen(md_md_cmd)
        md_md_process.wait()


        """
        X -- At this point, it's time to clean up the trajectories and run some analyses!! 
        (still in ../md)
        """
        input_selection = b"10 0\n"
        if args.run_on_alpine:
            convert_text = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o {args.protein_codename}_centered.xtc " \
                f"-n ../starting_structure/index_{args.protein_codename}.ndx"
            
            final_frame_txt = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o final_frame.pdb " \
                f"-b 500000 -e 500000 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
            start_frame_txt = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o starting_frame.pdb " \
                f"-b 0 -e 0 -n ../starting_structure/index_{args.protein_codename}.ndx"
        else:
            convert_text = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o {args.protein_codename}_centered.xtc " \
                f"-n ../starting_structure/index_{args.protein_codename}.ndx"
            
            final_frame_txt = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o final_frame.pdb " \
                f"-b 500000 -e 500000 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
            start_frame_txt = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o starting_frame.pdb " \
                f"-b 0 -e 0 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
        # This part converts the trajectory so that it's centered on the first coil 'coil-1'
        convert_cmd = convert_text.split()
        convert_process = subprocess.Popen(convert_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = convert_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())

        # this generates the final_frame.pdb for mdtraj stuff
        finalframe_cmd = final_frame_txt.split()
        finalframe_process = subprocess.Popen(finalframe_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = finalframe_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())

        # this generates the starting_frame.pdb because thats also useful
        startframe_cmd = start_frame_txt.split()
        startframe_process = subprocess.Popen(startframe_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = startframe_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())


    """
    X -- Now we need to do the analysis. Load the trajectories first and then pass the mdtraj 
    object into the different functions that I need
    (manually changes into ../md
        this seems unnecessary but this allows me to neatly rerun analyses)
    """
    os.chdir(f"{pro_sim_dir}/md")


    # first step is to reconvert trajectories if necessary. Just a simple flag to deal with in the 
    # arguments.
    if args.reconvert_trajectories:
        input_selection = b"10 0\n"
        if args.run_on_alpine:
            convert_text = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o {args.protein_codename}_centered.xtc " \
                f"-n ../starting_structure/index_{args.protein_codename}.ndx"
            
            final_frame_txt = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o final_frame.pdb " \
                f"-b 500000 -e 500000 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
            start_frame_txt = f"mpirun -np 1 gmx_mpi trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o starting_frame.pdb " \
                f"-b 0 -e 0 -n ../starting_structure/index_{args.protein_codename}.ndx"
        else:
            convert_text = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o {args.protein_codename}_centered.xtc " \
                f"-n ../starting_structure/index_{args.protein_codename}.ndx"
            
            final_frame_txt = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o final_frame.pdb " \
                f"-b 500000 -e 500000 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
            start_frame_txt = f"gmx trjconv -f {args.protein_codename}_md.xtc " \
                f"-s md_{args.protein_codename}.tpr -pbc whole -center -o starting_frame.pdb " \
                f"-b 0 -e 0 -n ../starting_structure/index_{args.protein_codename}.ndx"
            
        # This part converts the trajectory so that it's centered on the first coil 'coil-1'
        convert_cmd = convert_text.split()
        convert_process = subprocess.Popen(convert_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = convert_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())

        # this generates the final_frame.pdb for mdtraj stuff
        finalframe_cmd = final_frame_txt.split()
        finalframe_process = subprocess.Popen(finalframe_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = finalframe_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())

        # this generates the starting_frame.pdb because thats also useful
        startframe_cmd = start_frame_txt.split()
        startframe_process = subprocess.Popen(startframe_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = startframe_process.communicate(input=input_selection)
        print("STDOUT:", stdout.decode())
        print("STDERR:", stderr.decode())

    

    # Now I can handle the analysis

    simtraj = md.load(f"{args.protein_codename}_centered.xtc", top="final_frame.pdb")
    analyzer.analyze_angles_and_dihedrals(args.protein_codename,
                                          simtraj,
                                          ENV_INFO[3])
    shutil.copyfile(f"plot_{args.protein_codename}_backbone_angles.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_backbone_angles.png")
    
    analyzer.analyze_orientation_torsions(args.protein_codename,
                                          simtraj,
                                          ENV_INFO[4],
                                          ENV_INFO[5])
    shutil.copyfile(f"plot_{args.protein_codename}_Nterminal_orientation_torsions.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_Nterminal_orientation_torsions.png")
    shutil.copyfile(f"plot_{args.protein_codename}_Cterminal_orientation_torsions.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_Cterminal_orientation_torsions.png")

    fracframe_correct_multimer, fracframes_wrong_multimer = analyzer.analyze_reference_point_distances(args.protein_codename,
                                               simtraj,
                                               ENV_INFO[6],
                                               ENV_INFO[3],
                                               args.cutoff_for_correct_multimerization,
                                               args.cutoff_for_wrong_multimerization,
                                               args.threshold_for_multimerization,
                                               args.threshold_for_anycontact)
    shutil.copyfile(f"plot_{args.protein_codename}_key_ref_distances_0.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_key_ref_distances_0.png")
    shutil.copyfile(f"plot_{args.protein_codename}_areCoilsInMultimer.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_areCoilsInMultimer.png")
    shutil.copyfile(f"{args.protein_codename}_areCoilsInMultimer_by_frame.csv",
                    f"../analysis_{args.protein_codename}/{args.protein_codename}_areCoilsInMultimer_by_frame.csv")
    
    coil1_frac_correct_ETE, coil2_frac_correct_ETE = analyzer.analyze_ETE_distances(args.protein_codename,
                                   simtraj,
                                   ENV_INFO[3],
                                   args.ete_tolerance_threshold)
    shutil.copyfile(f"plot_{args.protein_codename}_ETE_fraction_frame.png",
                    f"../analysis_{args.protein_codename}/plot_{args.protein_codename}_ETE_fraction_frame.png")
    shutil.copyfile(f"{args.protein_codename}_ETEs_fraction_frame_correct.csv",
                    f"../analysis_{args.protein_codename}/{args.protein_codename}_ETEs_fraction_frame_correct.csv")
    

    # At this point I'm also adding some preparative work for the heat map information
    score_correct_multimer = fracframe_correct_multimer * ((coil1_frac_correct_ETE+coil2_frac_correct_ETE)/2)
    score_incorrect_multimer = fracframes_wrong_multimer * ((coil1_frac_correct_ETE+coil2_frac_correct_ETE)/2)
    score_total_multimer = (fracframe_correct_multimer + fracframes_wrong_multimer) * ((coil1_frac_correct_ETE+coil2_frac_correct_ETE)/2)

    np.savetxt(f"heatmap_{args.protein_codename}_correct_and_total_multimer.csv", 
               np.array([score_correct_multimer, score_incorrect_multimer, score_total_multimer,
                         fracframe_correct_multimer, fracframes_wrong_multimer,
                         coil1_frac_correct_ETE, coil2_frac_correct_ETE]),
               delimiter=",", fmt="%1.5f")
    shutil.copyfile(f"heatmap_{args.protein_codename}_correct_and_total_multimer.csv",
                    f"../analysis_{args.protein_codename}/heatmap_{args.protein_codename}_correct_and_total_multimer.csv")

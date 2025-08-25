"""
@Title:             importonly_write_bash_script.py
@Name:              Mando A Ramirez
@Date:              20250804


@Description:       This script writes the bash script for Alpine simulation to run the simulation of the desired
sequence specific validation simulation.
"""

def write_bash_script_text(script_file_name,
                           simulation_name,
                           starting_structure_file,
                           topology_file,
                           index_file,
                           sim_dir):
    bash_file_text = f"""#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --account=ucb342_asc3
#SBATCH --qos=normal
#SBATCH --constraint=ib
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --job-name={simulation_name}
#SBATCH --output=output_{simulation_name}_%j.txt

module purge
module load gcc/11.2.0
module load openmpi/4.1.1
export SLURM_EXPORT_ENV=ALL

#
#       Run EM
#
cd {sim_dir}/em
mpirun -np 1 gmx_mpi grompp -f ../mdp/em.mdp -c ../starting_structure/{starting_structure_file} \\
-p ../top/{topology_file} -n ../starting_structure/{index_file} \\
-o em_{simulation_name}.tpr

mpirun -np 1 --mca opal_common_ucx_opal_mem_hooks 1 gmx_mpi mdrun \\
-ntomp 1 -s em_{simulation_name}.tpr -deffnm {simulation_name}_em
cd ../


#
#       Run production MD
#
cd {sim_dir}/md
mpirun -np 1 gmx_mpi grompp -f ../mdp/md_pulling.mdp -c ../em/{simulation_name}_em.gro \\
-p ../top/{topology_file} -n ../starting_structure/{index_file} \\
-o md_{simulation_name}.tpr

mpirun -np 1 --mca opal_common_ucx_opal_mem_hooks 1 gmx_mpi mdrun \\
-ntomp 1 -rdd 1.6 \\
-s md_{simulation_name}.tpr -deffnm {simulation_name}_md


# Convert the trajectory while we're here
echo 10 0 | mpirun -np 1 gmx_mpi trjconv -f {simulation_name}_md.xtc \
-s md_{simulation_name}.tpr -pbc whole -center -o {simulation_name}_centered.xtc

echo 10 0 | mpirun -np 1 gmx_mpi trjconv -f {simulation_name}_md.xtc \
-s md_{simulation_name}.tpr -pbc whole -center \
-b 500000 -e 500000 -o final_frame.pdb
cd ../

"""
    
    with open(script_file_name, 'w') as output_file:
        output_file.write(bash_file_text)

    return None

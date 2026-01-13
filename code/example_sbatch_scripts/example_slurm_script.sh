#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --job-name=myExperiment
#SBATCH --output=output_my_experiment_%j.txt

# the %j in the output job name should be included because that assigns the
# slurm id to your output file, which can be really helpful in debugging problems
# associated with Alpine/slurm/hardware, if any pop up
#
# the flags that you can/should change are --ntasks and --time
# You should almost never need to change --nodes for any simulations we do.
#   1 node gives you access to 64 ntasks (i.e. cpus) for your simulations
#
# --ntasks should be set based on the number of cpus you need for your simulation
# which might be high for atomistic simulations of dimers
# --time can be set based on how much time you need for the simulation, but note that 
# you can only request 24 hours of time with --qos=normal. See Alpine documentation for
# other time request restraints



# purge and prep modules
module purge
# module load anaconda
# conda activate myEnvironment
module load gcc/11.2.0
module load openmpi/4.1.1
export OMPI_TMPDIR=/scratch/alpine/$USER
export SLURM_EXPORT_ENV=ALL

# The anaconda/conda lines are included for when you need to run python analysis in a sbatch file.
# You should be using conda environments when analyzing trajectories using something like mdtraj
# the gcc and openmpi lines are important to properly run the parallelized version of GROMACS I installed
# the two export lines are important to make sure openMPI runs correctly and doesn't end up crashing
# by running out of temporary storage space


SIMDIR="/projects/dora1300/pkgs/gromacs-2022-cpu-mpi/bin"
# this is the directory location of my 'gmx_mpi' executable you can use.
#   IMPORTANT! this version of gromacs can only use CPUs, and not GPUs, for simulation
# you can add this directory to your $PATH variable in your .bashrc if you feel comfortable
# doing that, which can make your life a little easier.



# an example set up for GROMACS md simulations. The backslashes allow you to use multiple lines for one single command
mpirun -np 1 ${SIMDIR}/gmx_mpi grompp \
-f parameters.mdp -p topol.top -c config.gro \
-n index.ndx -o md_system.tpr


# if you only want one core/cpu to run your simulation
mpirun -np 1 gmx_mpi --mca opal_common_ucx_opal_mem_hooks 1 ${SIMDIR}/gmx_mpi mdrun \
-ntomp 1 -rdd 1.6 -s md_system.tpr -deffnm md_output
#   OR
# if you want multiple cores/cpus to run your simulation --> make sure -ntomp <= --ntasks above!
mpirun -np 1 gmx_mpi --mca opal_common_ucx_opal_mem_hooks 1 ${SIMDIR}/gmx_mpi mdrun \
-ntomp 16 -rdd 1.6 -s md_system.tpr -deffnm md_output

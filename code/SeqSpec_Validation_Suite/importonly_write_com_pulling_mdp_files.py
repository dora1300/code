"""
@Title:             importonly_write_com_pulling_mdp_files.py
@Name:              Mando A Ramirez
@Date:              2025 07 08

@Description:       This script handles the automated writing of mdp files relevant for the COM
pulling validation MD simulations. 

This script DOES NOT SAVE THE ACTUAL MDP FILE TEXT!!!!
It only generates the text, then returns the text which can be used at a separate time.

"""
def write_em_mdp_text(em_mdp_file_name="em.mdp",
        em_tol = 50.0,
        em_step = 0.002,
        vbt = "1E-7"
):
    
    em_file_text = f""";
;    MDP file for Energy Minimization
;       part of Sequence-specific validation suite (est. 20250804)
;

; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = {em_tol}          ; Stop minimization when the maximum force < {em_tol} kJ/mol/nm
emstep      = {em_step}         ; Minimization step size (nm)
nsteps      = 500000        ; Maximum number of (minimization) steps to perform

nstlog                 = 1000            ; frequency to write to log file
nstenergy              = 1000            ; frequency to write to energy file
nstxout-compressed     = 1000            ; frequency to write coords to .xtc file


; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = Verlet
nstlist                  = 10
pbc                      = xyz
verlet-buffer-tolerance  = {vbt}

vdw-type                 = Cut-off
vdw-modifier             = Potential-shift
rvdw                     = 1.1

coulombtype              = Cut-off
rcoulomb                 = 1.1"""

    with open(em_mdp_file_name, 'w') as output_file:
        output_file.write(em_file_text)
    
    return None



def write_md_pulling_mdp_text(
    temp,
    pull_groups,
    md_mdp_file_name="md_pulling.mdp",
    vbt="1E-12",
    pull_init="2"
):
    """
    Description:
    Parameters:
        pull_groups     int; the number of pull groups (this function only supports 
                        2 and 3 right now)
    Returns:
        None. Writes a file with title "md_mdp_file_name" and then exits the function.
    """

    md_file_text = f""";
;    Production MD .mdp file
;       part of Sequence-specific validation suite (est. 20250804)
;
;   for COM pulling between multiple coil groups
;

; ---- Run control
integrator            = sd            ; stochastic dynamics integrator
                                      ; when used, tcoupl and nsttcouple are ignored
tinit                 = 0
dt                    = 0.010         ; units of ps
nsteps                = 20000000      ; units of (ps) ==> total time = 0.2 us == 200 ns
comm-mode             = Linear        ; remove the center of mass translational velocity
nstcomm               = 10            ; frequency for center of mass motion removal
ld-seed               = -1            ; -1 means pseudo-randomly chosen, picked by GROMACS


; ---- Output control
nstlog                = 10000         ; 5000 steps
nstcalcenergy         = 10
nstenergy             = 1000          ; 50000 steps, 10 ps value
nstxout-compressed    = 10000         ; 5000 steps, 100 ps value


; ---- Neighbor searching and short-range nonbonded interactions
cutoff-scheme            = Verlet
nstlist                  = 10        ; default value for Gromacs 2021
pbc                      = xyz
verlet-buffer-tolerance  = {vbt}     ; kJ mol^-1 ps^-1

; van der Waals
vdw-type                = Cut-off
vdw-modifier            = Potential-shift
rvdw                    = 2.0

; Electrostatics
coulombtype             = Cut-off
rcoulomb                = 2.0


; ---- Pressure controlling and coupling
pcoupl                  = no


; ---- Velocity generation
gen-vel                 = no         ; only useful for the md integrator


; ---- Temperature controlling and Langevin noise
tc-grps                 = System
tau-t                   = 5        
ref-t                   = {temp}      ; in units of (K)

"""

    if int(pull_groups) == 2:
        md_file_text += f"""
; ---- COM Pulling
pull                    = yes
pull-nstfout            = 100

pull-ngroups            = 2
pull-ncoords            = 1
pull-group1-name        = coil-1
pull-group2-name        = coil-2

pull-coord1-type        = flat-bottom
pull-coord1-geometry    = distance
pull-coord1-groups      = 1 2
pull-coord1-start       = no
pull-coord1-init        = {pull_init}           ; [nm]
pull-coord1-dim         = Y Y Y
pull-coord1-k           = 2000            ; [kJ/mol/nm^2]

"""
    else:
        md_file_text += f"""
; ---- COM Pulling
pull                    = yes
pull-nstfout            = 100

pull-ngroups            = 3
pull-ncoords            = 2
pull-group1-name        = coil-1
pull-group2-name        = coil-2
pull-group3-name        = coil-3

pull-coord1-type        = flat-bottom
pull-coord1-geometry    = distance
pull-coord1-groups      = 1 2
pull-coord1-start       = no
pull-coord1-init        = {pull_init}           ; [nm]
pull-coord1-dim         = Y Y Y
pull-coord1-k           = 2000          ; [kJ/mol/nm^2]

pull-coord2-type        = flat-bottom
pull-coord2-geometry    = distance
pull-coord2-groups      = 1 3
pull-coord2-start       = no
pull-coord2-init        = {pull_init}           ; [nm]
pull-coord2-dim         = Y Y Y
pull-coord2-k           = 2000          ; [kJ/mol/nm^2]

"""

    with open(md_mdp_file_name, 'w') as output_file:
        output_file.write(md_file_text)
    
    return None
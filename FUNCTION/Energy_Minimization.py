import subprocess
import os
import logging

def energy_min(minim_mdp_path, gro_name, top_name, output_name, gmx_path, number_of_run=1, max_warn=0):
    """
    Perform energy minimization using GROMACS tools (grompp and mdrun).

    Parameters:
    - minim_name: MDP file for energy minimization.
    - gro_name: Input structure file.
    - top_name: Topology file.
    - output_name: Output file name for the minimized structure.
    - number_of_run: Number of energy minimization steps (default is 1).
    - max_warn: 
    """

    log_file = "energy_minimization.out"

    # Run the first 
    grompp_cmd = f"{gmx_path} grompp -f {minim_mdp_path} -c {gro_name}.gro -p {top_name}.top -o {gro_name}_EM1.tpr -maxwarn 1"
    try:
        subprocess.run(grompp_cmd, shell=True, check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with 1st GROMPP!")
        return

    # Run the first 
    mdrun_cmd = f"{gmx_path} mdrun -s {gro_name}_EM1.tpr -c {gro_name}_EM1.gro -v"
    try:
        subprocess.run(mdrun_cmd, shell=True, check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with 1st MDRUN!")
        return

    # Further runs if number_of_run > 1
    if number_of_run > 1:
        start = 1
        for run in range(2, number_of_run + 1):
            # Run grompp for the next step
            grompp_cmd = f"{gmx_path} grompp -f {minim_mdp_path} -c {gro_name}_EM{start}.gro -p {top_name}.top -o {gro_name}_EM{run}.tpr -maxwarn {max_warn}"
            try:
                subprocess.run(grompp_cmd, shell=True, check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                print(f"Something went wrong with {run}st GROMPP!")
                return

            # Run mdrun for the next step
            mdrun_cmd = f"{gmx_path} mdrun -s {gro_name}_EM{run}.tpr -c {gro_name}_EM{run}.gro -v"
            try:
                subprocess.run(mdrun_cmd, shell=True, check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                print(f"Something went wrong with {run}st MDRUN!")
                return

            start += 1

    # Rename the final output
    final_gro_name = f"{gro_name}_EM{number_of_run}.gro"
    os.rename(final_gro_name, f"{output_name}_minim.gro")

    # Check if the minimization was successful
    with open(log_file, 'r') as f:
        minim_test = None
        for line in f:
            if "Norm of force" in line:
                minim_test = line.split('=')[1].strip()
                break
        
        if minim_test == "inf":
            print("\tSomething went wrong with the energy minimization.")
            return
    logging.info("Energy minimization completed successfully.")

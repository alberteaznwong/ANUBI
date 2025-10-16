import subprocess
import os
import time
import logging

def run_gromacs_command(command, error_message, pipe_file, output_file=None):
    """
    Runs a GROMACS command and logs the output to a file.
    
    command (str): The GROMACS command to execute.
    error_message (str): Error message to display if the command fails.
    pipe_file (str): File to write "exit" if the command fails.
    output_file (str): File to log command output.
    """
    try:
        if output_file:
            with open(output_file, "w") as out:
                subprocess.run(command, check=True, shell=True, stdout=out, stderr=subprocess.STDOUT)
        else:
            subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        print(f"Something went wrong: {error_message}")
        with open(pipe_file, "w") as f:
            f.write("exit")
        raise
'''
def make_new_minim_config_samd(input_structure_file, samd_mdp, output_gro, sequence, top_name="topol", pipe_file="samd_out.out"):
    """
    Performs SAMD using GROMACS commands, and logs outputs to files.
    
    input_structure_file (str): Path to the input structure file.
    samd_mdp (str): Path to the SAMD MDP file.
    output_gro (str): Path to the output GRO file.
    top_name (str): Topology file name without extension. Default is "topol".
    pipe_file (str): Pipe file for signaling errors. Default is "samd_out.out".
    """
    # SAMD 
    #print(f"{time.strftime('%H:%M:%S')} -- Running SAMD ")
    logging.info("Running Simulated Annealing MD.")
    grompp_samd_out = f"gromppSAMD_seq{sequence}.out"
    mdrun_samd_out = f"mdoutSAMD_seq{sequence}.out"

    # GROMPP step for SAMD
    samd_grompp_command = (
        f"gmx grompp -f {samd_mdp} -c {input_structure_file} -r {input_structure_file} "
        f"-p {top_name}.top -o system_SAMD.tpr -maxwarn 1"
    )
    run_gromacs_command(samd_grompp_command, "Something wrong on SAMD GROMPP", pipe_file, output_file=grompp_samd_out)

    # MDRUN step for SAMD
    samd_mdrun_command = (
        f"gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu "
        f"-s system_SAMD.tpr -c {output_gro}.gro -cpo state_SAMD.cpt -x traj_SAMD.xtc -e SAMD.edr -v"
    )
    run_gromacs_command(samd_mdrun_command, "Something wrong on SAMD MDRUN", pipe_file, output_file=mdrun_samd_out)

    # Energy checks for SAMD
    try:
        subprocess.run(f"echo 'Temperature\nPressure\nDensity\n\n' | gmx energy -f SAMD.edr -o press_SAMD.xvg", check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Something went wrong during the energy check")
        return
    try:
        subprocess.run(f"printf '1\n1\n' | gmx rms -s system_SAMD.tpr -f traj_SAMD.xtc -o rmsd_SAMD{sequence}.xvg -tu ps", check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Something went wrong during the RMSD check")
        return

    print("SAMD completed successfully!")

    # Copy output files to the specified folder for checking
    os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
    subprocess.run(f"cp ./grompp*.out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=True, shell=True)
'''
def run_md(md_mdp, tpr_file, trj_name, sequence, cycle_number, gmx_path, top_name = "topol", pipe_file="Prod_MD.out"):
    """
    Runs a full MD cycle including energy and RMSD checks for the first cycle.

    Parameters:
        grompp_path (str): Path to grompp executable.
        mdrun_path (str): Path to mdrun executable.
        energy_path (str): Path to gmx energy command.
        rms_path (str): Path to gmx rms command.
        md_engcomp_ff14sb_name (str): Name of MDP file for production MD.
        top_name (str): Topology file name.
        tpr_file (str): Name for output TPR file.
        sequence (str): Sequence identifier for files.
        cycle_number (int): Current cycle number.
        double_check_folder (str): Folder for storing double-check files.
        cycle_number_md_folder (str): Folder for current cycle's temporary files.

    Returns:
        None
    """
    # MD
    logging.info(f"Running MD Production cycyle{cycle_number}")
    grompp_md_out = f"gromppPROD_seq{sequence}.out"
    mdrun_md_out = f"mdoutPROD_seq{sequence}.out"

    # Run GROMPP
    grompp_command = (
        f"{gmx_path} grompp -f {md_mdp} -c system_equil.gro -r system_equil.gro "
        f"-p {top_name}.top -o {tpr_file}.tpr -maxwarn 1"
    )
    run_gromacs_command(grompp_command, "Error in GROMPP", pipe_file, output_file = grompp_md_out)

    # Run MDRUN
    mdrun_command = (
        # -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu
        f"{gmx_path} mdrun -s {tpr_file}.tpr -c system_Compl_MD.gro -x {trj_name}.xtc -e PROD.edr -v"
    )
    run_gromacs_command(mdrun_command, "Something wrong on MD MDRUN", pipe_file, output_file = mdrun_md_out)

    if cycle_number < 19:
        try:
            # Energy Check
            energy_command = (
                f"echo 'Temperature\nPressure\nDensity\n0\n' | {gmx_path} energy -f PROD.edr -o PROD{sequence}.xvg"
            )
            run_gromacs_command(
                command=energy_command,
                error_message="Something went wrong during the energy check",
                pipe_file=pipe_file
            )

            # RMSD Check
            rms_command = (
                f"printf '1\n1\n' | {gmx_path} rms -s {tpr_file}.tpr -f {trj_name}.xtc -o rmsd_PROD{sequence}.xvg "
                f"-a avgPROD{sequence}.xvg -tu ps"
            )
            run_gromacs_command(
                command=rms_command,
                error_message="Something went wrong during RMSD Check",
                pipe_file=pipe_file
            )
        
            os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
            subprocess.run(f"cp ./mdout.mdp DOUBLE_CHECK_FOLDER/mdoutPROD_seq{sequence}.mdp", check=True, shell=True)
            subprocess.run(f"cp ./grompp*.out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=True, shell=True)

        except Exception as e:
            print("Something wrong during Check for MD")
            raise
    
    logging.info(f"MD Production cycyle{cycle_number} Done")
# run_md(md_mdp_path,system_Compl_MD, traj_MD, 0, 1)

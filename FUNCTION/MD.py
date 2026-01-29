import subprocess
import os
import time
import logging


def run_gromacs_command(command, error_message, pipe_file, output_file=None, env=None):
    """Run a shell command (typically gmx) and record timing.

    - `command` is a shell command string.
    - `pipe_file` will be written with "exit" on failure to signal upstream.
    - `output_file` captures stdout/stderr when provided.
    - `env` is an optional dict of environment variables to use for the subprocess.
    """
    start = time.time()
    logging.info(f"Running command: {command}")
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    try:
        if output_file:
            with open(output_file, "w") as out:
                subprocess.run(command, check=True, shell=True, stdout=out, stderr=subprocess.STDOUT, env=run_env)
        else:
            subprocess.run(command, check=True, shell=True, env=run_env)
    except subprocess.CalledProcessError:
        logging.exception(error_message)
        with open(pipe_file, "w") as f:
            f.write("exit")
        raise
    finally:
        end = time.time()
        logging.info(f"Command finished (elapsed {end-start:.3f}s): {command}")


def run_md(md_mdp, tpr_file, trj_name, sequence, cycle_number, gmx_path, top_name="topol", pipe_file="Prod_MD.out", ntmpi=1, ntomp=None, gpu_flags=None, env=None):
    """Run production MD: grompp then mdrun with optional threading and GPU flags.

    Parameters mirror the ANUBI caller; additional knobs:
    - `ntmpi`, `ntomp` control `gmx mdrun` `-ntmpi`/`-ntomp` options.
    - `gpu_flags` is a string appended to the mdrun command (e.g. "-nb gpu -pme gpu ...").
    - `env` is an optional dict of extra environment variables to set for the subprocess.
    """
    logging.info(f"Running MD Production cycle {cycle_number}")
    grompp_md_out = f"gromppPROD_seq{sequence}.out"
    mdrun_md_out = f"mdoutPROD_seq{sequence}.out"

    # Prepare environment
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    # If ntomp is provided, set OMP_NUM_THREADS for the subprocess environment
    if ntomp is not None:
        run_env["OMP_NUM_THREADS"] = str(ntomp)
        # set a conservative MKL control if present
        run_env.setdefault("MKL_NUM_THREADS", "1")

    # Run GROMPP to produce the TPR
    grompp_command = (
        f"{gmx_path} grompp -f {md_mdp} -c system_equil.gro -r system_equil.gro -p {top_name}.top -o {tpr_file}.tpr -maxwarn 1"
    )
    run_gromacs_command(grompp_command, "Error in GROMPP", pipe_file, output_file=grompp_md_out, env=run_env)

    # Build mdrun command
    mdrun_base = f"{gmx_path} mdrun -s {tpr_file}.tpr -c system_Compl_MD.gro -x {trj_name}.xtc -e PROD.edr -v"
    mdrun_parts = [mdrun_base]
    if ntmpi:
        mdrun_parts.append(f"-ntmpi {ntmpi}")
    if ntomp:
        mdrun_parts.append(f"-ntomp {ntomp}")
    if gpu_flags:
        mdrun_parts.append(gpu_flags)

    mdrun_command = " ".join(mdrun_parts)
    run_gromacs_command(mdrun_command, "Something wrong on MD MDRUN", pipe_file, output_file=mdrun_md_out, env=run_env)

    # Post-checks for early cycles
    if cycle_number < 19:
        try:
            energy_command = f"echo 'Temperature\\nPressure\\nDensity\\n0\\n' | {gmx_path} energy -f PROD.edr -o PROD{sequence}.xvg"
            run_gromacs_command(energy_command, "Something went wrong during the energy check", pipe_file, env=run_env)

            rms_command = (
                f"printf '1\\n1\\n' | {gmx_path} rms -s {tpr_file}.tpr -f {trj_name}.xtc -o rmsd_PROD{sequence}.xvg -a avgPROD{sequence}.xvg -tu ps"
            )
            run_gromacs_command(rms_command, "Something went wrong during RMSD Check", pipe_file, env=run_env)

            os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
            subprocess.run(f"cp ./mdout.mdp DOUBLE_CHECK_FOLDER/mdoutPROD_seq{sequence}.mdp", check=False, shell=True)
            subprocess.run(f"cp ./grompp*.out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=False, shell=True)

        except Exception:
            logging.exception("Something wrong during Check for MD")
            raise

    logging.info(f"MD Production cycle {cycle_number} Done")
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
    def run_gromacs_command(command, error_message, pipe_file, output_file=None, env=None):

    # MDRUN step for SAMD
    samd_mdrun_command = (
        f"gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu "
        f"-s system_SAMD.tpr -c {output_gro}.gro -cpo state_SAMD.cpt -x traj_SAMD.xtc -e SAMD.edr -v"
    )
    run_gromacs_command(samd_mdrun_command, "Something wrong on SAMD MDRUN", pipe_file, output_file=mdrun_samd_out)

        start = time.time()
        logging.info(f"Running command: {command}")
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
        finally:
            end = time.time()
            logging.info(f"Command finished (elapsed {end-start:.1f}s): {command}")

    print("SAMD completed successfully!")
    def run_md(md_mdp, tpr_file, trj_name, sequence, cycle_number, gmx_path, top_name = "topol", pipe_file="Prod_MD.out", ntmpi=1, ntomp=None, gpu_flags=None):

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

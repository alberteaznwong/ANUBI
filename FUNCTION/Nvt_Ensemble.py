import subprocess
import os
import time
import logging

def run_gromacs_command(command, error_message, pipe_file, output_file=None):
    """
    command (str): The GROMACS command to execute.
    error_message (str): Error message to display if the command fails.
    pipe_file (str): 
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
        raise SystemExit(error_message)

def make_nvt(input_structure_file, nvt_mdp, output_gro, sequence, gmx_path, npt_mdp = "",top_name="topol", pipe_file="out.out"):
    """
    Performs NVT and NPT equilibration using GROMACS commands, and logs outputs to files.

    input_structure_file (str): input structure file.
    nvt_mdp (str): Path to the NVT MDP file.
    npt_mdp (str): Path to the NPT MDP file.
    output_gro (str): output GRO file.
    top_name (str): Topology file name without extension. Default is "topol".
    pipe_file (str): Pipe file for signaling errors. Default is "out.out".
    """
    # NVT
    logging.info("Running NVT MD for temperature equilibration.")
    #print(f"{time.strftime('%H:%M:%S')} -- Running NVT MD for temperature equilibration...")
    grompp_nvt_out = f"gromppNVT_seq{sequence}.out"
    mdrun_nvt_out = f"mdrun_NVT_MD_seq{sequence}.out"

    nvt_grompp_command = (
        f"{gmx_path} grompp -f {nvt_mdp} -c {input_structure_file} -r {input_structure_file} "
        f"-p {top_name}.top -o system_NVT_MD.tpr -maxwarn 1"
    )
    run_gromacs_command(nvt_grompp_command, "Something wrong on NVT GROMPP", pipe_file, output_file=grompp_nvt_out)
    # -nb gpu -pme gpu -bonded gpu -update gpu 
    # "-s system_NVT_MD.tpr -c system_NVT_MD.gro -cpo state_NVT_MD.cpt -e NVT.edr -v"
    nvt_mdrun_command = (
        f"{gmx_path} mdrun "
        f"-s system_NVT_MD.tpr -c {output_gro}.gro -cpo state_NVT_MD.cpt -e NVT.edr -v"
    )


    run_gromacs_command(nvt_mdrun_command, "Something wrong on NVT MDRUN", pipe_file, output_file=mdrun_nvt_out)
    """
    # NPT
    logging.info("Running NPT MD for pressure equilibration")
    #print(f"{time.strftime('%H:%M:%S')} -- Running NPT MD for pressure equilibration...")
    grompp_npt_out = f"gromppNPT_seq{sequence}.out"
    mdrun_npt_out = f"mdrun_NPT_MD_seq{sequence}.out"

    npt_grompp_command = (
        f"{gmx_path} grompp -f {npt_mdp} -c system_NVT_MD.gro -r system_NVT_MD.gro "
        f"-p {top_name}.top -o system_NPT_MD.tpr -t state_NVT_MD.cpt -maxwarn 2"
    )
    run_gromacs_command(npt_grompp_command, "Something wrong on NPT GROMPP", pipe_file, output_file=grompp_npt_out)

    npt_mdrun_command = (
        f"{gmx_path} mdrun "
        f"-s system_NPT_MD.tpr -c {output_gro}.gro -cpo state_NPT_MD.cpt -x traj_NPT_MD.xtc -e NPT.edr -v"
    )
    run_gromacs_command(npt_mdrun_command, "Something wrong on NPT MDRUN", pipe_file, output_file=mdrun_npt_out)

    # Energy checks
    try:
        # Correctly select Temperature for NVT energy file
        subprocess.run(f"echo 'Temperature\n\n' | {gmx_path} energy -f NVT.edr -o temp_NVT.xvg", check=True, shell=True)
    
        # Correctly select Pressure and DEnsity for NPT energy file
        subprocess.run(f"echo 'Pressure\nDensity\n\n' | {gmx_path} energy -f NPT.edr -o press_NPT.xvg", check=True, shell=True)
    
    except subprocess.CalledProcessError:
        print("Something went wrong on the energy check...")
        return
    
    logging.info("Equilibration completed successfully.")
    """
    # Copy output files to the specified folder for checking
    os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
    #subprocess.run(f"cp ./grompp*.out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=True, shell=True)
    subprocess.run(f"cp ./grompp*.out ./*edr DOUBLE_CHECK_FOLDER", check=True, shell=True)


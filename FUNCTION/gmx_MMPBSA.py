import subprocess
import os
import time
import logging
# ForceField="amber99sb-ildn"
# NP=32
def gmx_mmpbsa(cycle_number, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_md_folder, conf_name, root_name, 
                     force_field, FORCE_FIELD_PATH, mmpbsa_infile, MMPBSA_INFILE_PATH, np_value, number_of_frames):
    """Main function to run the gmxMMPBSA cycle"""
    try:
        # Print log for starting the cycle
        #print(f"{time.strftime('%H:%M:%S')} --running gmxMMPBSA_cycle{cycle_number}..")
        logging.info(f"Running gmxMMPBSA_cycle{cycle_number}")

        # Activate conda environment
        #subprocess.check_call(f"source {conda_activate_path}/activate {conda_gmxmmpbsa_name}", shell=True, executable = "/bin/bash")
        #subprocess.check_call(f"conda activate {conda_gmxmmpbsa_name}", shell = True, executable = "/bin/bash")
        #%conda activate {conda_gmxmmpbsa_name}

        # Search for the ForceField in GMXLIB

        # Prepare files for MMPBSA
        cycle_folder = os.path.join(cycle_number_md_folder, root_name)
        if not os.path.isdir(cycle_folder):
            raise RuntimeError(f"Cannot enter '{root_name}' folder")
        os.chdir(cycle_folder)

        subprocess.check_call(f"cp {MMPBSA_INFILE_PATH}/{mmpbsa_infile} {cycle_folder}/", shell=True, executable = "/bin/bash")

        # Calculate NP_used
        np_half = np_value // 2
        np_used = min(np_half, int(number_of_frames)) 
        print(f"NP_value={np_value} \t number_of_frames={number_of_frames} \t NP_used={np_used}")

        # Run gmxMMPBSA calculation
      
        gmxMMPBSA_env = conda_gmxmmpbsa_name

        command = f"mpirun -np {np_used} gmx_MMPBSA MPI -O -i {mmpbsa_infile} -cs {conf_name}_newGRO.tpr -ci index.ndx " \
          f"-cg 0 1 -ct {conf_name}_noPBC.xtc -cr ./{conf_name}_starting_protein.pdb -cp topol_protein.top " \
          f"-eo gmx_MMPBSA_plot.csv -nogui"
        
        
        
        #full_command = f"conda run -n {gmxMMPBSA_env} {command}"
        full_command = f"source {conda_activate_path} {gmxMMPBSA_env} && {command}"
     
        # output file
        output_file = "gmx_mmpbsa.out"

        with open(output_file, "w") as out_file:
            subprocess.check_call(full_command, shell=True, executable="/bin/bash", stdout=out_file, stderr=out_file)

    except subprocess.CalledProcessError as e:
        if "Segmentation fault" in str(e):
            print("Some error occurred with gmx_MMPBSA.. waiting 2min and try again")
            time.sleep(2 * 60)  # wait for 2 minutes
            subprocess.check_call(command, shell=True, executable = "/bin/bash")
        else:
            raise RuntimeError(f"gmx_MMPBSA failed: {e}")
    print(f"Finished gmx_MMPBSA on cycle{cycle_number}")
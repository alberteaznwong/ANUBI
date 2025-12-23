########################################## IMPORT MODULES ######################################
import os
import datetime
import logging
import shutil
import pandas as pd
import glob
import multiprocessing
import re
import subprocess
import numpy as np
import math
import random
import yaml
import argparse
import sys
from Bio.PDB import PDBParser
from FUNCTION import make_top_protein, fill_water_ions, energy_min, make_equi, run_md
from FUNCTION import files_gmxmmpbsa, gmx_mmpbsa, Data_Analysis_Pre, Data_Analysis_Cal, clean_for_each_cycle
from FUNCTION import Data_Analysis_Cal_child, peptide_mode, extract_lastframe_and_rename




############################################### FUNCTION DEFINATION #####################################
def load_config(config_file):
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)
    return config

def get_version(command):
    try:
        results = subprocess.run([command,"--version"],capture_output=True, text=True, check=True)
        return results.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def which_program(command):
    try:
        results = subprocess.run(["which",command],capture_output=True, text=True, check=True)
        return results.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def check_conda_env(env):
    try:
        check_env = subprocess.run(['conda', 'info', '--envs'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        if env in check_env.stdout:
            return True
        else:
            return False
    except subprocess.CalledProcessError:
        return False
'''
def create_output_directory():
    
    current_dir = os.getcwd()
    
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir_name = f"output_{timestamp}"
    output_dir_path = os.path.join(current_dir, output_dir_name)
    
    os.mkdir(output_dir_path)
    print(f"Created directory: {output_dir_path}")

    os.chdir(output_dir_path)
    return output_dir_path
'''

def create_output_directory():
    current_dir = os.getcwd()
    
    # Get the current date formatted as MMDD
    date_str = datetime.datetime.now().strftime('%m%d')
    base_name = f"output_{date_str}_"
    index = 1
    
    # Generate a new directory name that doesn't already exist
    while True:
        output_dir_name = f"{base_name}{index}"
        output_dir_path = os.path.join(current_dir, output_dir_name)
        
        if not os.path.exists(output_dir_path):
            os.mkdir(output_dir_path)
            print(f"Created directory: {output_dir_path}")
            break
        
        index += 1

    os.chdir(output_dir_path)
    return output_dir_path




def build_folders(current_dir, cycle_num):
    # Create folder for each cycle
    folders  ={}
    
    for cycle_n in range (1,cycle_num + 1):
        folder_name = f"cycle{cycle_n}_MD"
        folder_path = os.path.join(current_dir, folder_name)
        os.makedirs(folder_path, exist_ok = True)
        folders[f"cycle{cycle_n}_MD"] = folder_path

    folders["repository"] = os.path.join(current_dir,"REPOSITORY")
    folders["TEMP_FILES_FOLDER"] = os.path.join(current_dir,"TEMP_FILES_FOLDER")
    folders["REMOVED_FILES_FOLDER"] = os.path.join(current_dir,"REMOVED_FILES_FOLDER")
    folders["results"] = os.path.join(current_dir,"RESULTS")

    for folder in folders.values():
        os.makedirs(folder,exist_ok = True)
    '''
    header = [
    "#RUNnumber", "DeltaG(kcal/mol)", "Coul(kcal/mol)", "vdW(kcal/mol)",
    "PolSol(kcal/mol)", "NpoSol(kcal/mol)", "ScoreFunct", "ScoreFunct2",
    "Canonica_AVG", "MedianDG", "DeltaG_2s"]

    '''
    header = [
    "#RUNnumber", "DeltaG(kcal/mol)", "Coul(kcal/mol)", "vdW(kcal/mol)",
    "PolSol(kcal/mol)", "NpoSol(kcal/mol)", "ScoreFunct", "ScoreFunct2",
    "MedianDG", "DeltaG_2s"]
    df = pd.DataFrame(columns=header)
    results_file_path = os.path.join(folders["results"], "MoleculesResults.dat")
    df.to_csv(results_file_path, sep='\t', index=False, header=True)

    return folders

def add_ter_to_pdb(pdb_file_name):
    temp_file_name = f"{pdb_file_name}_temp"  

    with open(pdb_file_name, 'r') as f:
        lines = f.readlines()

    new_lines = []
    prev_chain_id = None  
    num_lines = len(lines)

    for i, line in enumerate(lines):
  
        if not line.startswith("ATOM"):
            new_lines.append(line)
            continue


        current_chain_id = line[21]  # get chain ID (the 22rd column)

        # add
        new_lines.append(line)


        if prev_chain_id is not None:

            if (i == num_lines - 1 or 
                (lines[i + 1].startswith("ATOM") and lines[i + 1][21] != current_chain_id) or 
                not lines[i + 1].startswith("ATOM")):
                new_lines.append("TER\n")

 
        prev_chain_id = current_chain_id


    with open(temp_file_name, 'w') as f:
        f.writelines(new_lines)

  
    os.rename(temp_file_name, pdb_file_name)

def replace_his_residues_flexible(input_pdb, output_pdb):
    with open(input_pdb, "r", encoding="utf-8") as infile, open(output_pdb, "w", encoding="utf-8") as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                match = re.match(r"(.{17})(HISE|HISD|HISP)(.*)", line)
                if match:
                    line = f"{match.group(1)}{'HIS':<4}{match.group(3)}"
            
            outfile.write(line.rstrip('\n') + '\n')


def MD_for_each_cycle(work_dir, cycle_number,sequence, md_mdp_path, tpr_file, trj_name, gmx_path):
    #cycle_number = 1
    #while cycle_number <= cycle_num:
        
    #cycle_MD_path = folders[f"cycle{cycle_number}_MD"] 
    os.chdir(work_dir)
    shutil.copy(os.path.join(folders["repository"], "system_equil.gro"), "./")
    shutil.copy(os.path.join(folders["repository"], "topol.top"), "./")

    for itp_file in glob.glob(os.path.join(folders["repository"], "*rotein_chain_*.itp")):
        shutil.copy(itp_file, "./")

    for itp_file in glob.glob(os.path.join(folders["repository"], "posres_*.itp")):
        shutil.copy(itp_file, "./")

    for cpt_file in glob.glob(os.path.join(folders["repository"], "*NPT*.cpt")):
        shutil.copy(cpt_file, "./")

    #make_new_minim_config_samd("system_equil.gro", samd_mdp_path, "system_Compl_MDstart", 0)
    #make_new_minim_config_samd(input_structure_file, samd_mdp_path, output_gro, sequence)
    #run_md(md_mdp_path,"system_Compl_MD", "traj_MD", 0, 1)
    run_md(md_mdp_path, tpr_file, trj_name, sequence, cycle_number, gmx_path)
    shutil.copy("system_Compl_MD.gro", f"LastFrame_cycle{cycle_number}.gro")
    #cycle_number += 1


def gmx_mmpbsa_for_each_cycle(work_dir, cycle_number,only_protein_md_mdp_path,temp_files_folder, FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, REMOVED_FILES_FOLDER, results_folder, repository_folder, current_conf_path):
    #cycle_number = 1
    #while cycle_number <= cycle_num:
    ConfName = f"cycle{cycle_number}"
    RootName = f"cycle{cycle_number}_BE"
    cycle_number_MD_FOLDER = folders[f"cycle{cycle_number}_MD"]
    # print
    print(f"Cycle Number: {cycle_number}")
    print(f"Configuration Name: {ConfName}")
    print(f"Root Name: {RootName}")
    print(f"MD Folder Path: {cycle_number_MD_FOLDER}")
    #os.chdir(cycle_number_MD_FOLDER )
    os.chdir(work_dir)
        
    repository_pdb_file = os.path.join(repository_folder, f"{protein_infile}.pdb")
    #startingFrameGMXPBSA="2000"
    # make files for gmx_mmpbsa
    # files_gmxmmpbsa(starting_gro_file, repository_pdb_file, trj_file, tpr_file, top_file, mdp_name, root_name, conf_name, vmd_function_folder, temp_files_folder)

    files_gmxmmpbsa("system_Compl_MD", repository_pdb_file, "traj_MD", "system_Compl_MD", "topol", only_protein_md_mdp_path, RootName, ConfName,temp_files_folder, cycle_number, startingFrameGMXPBSA, receptorFRAG, ABchains,gmx_path)
    # get number of frames
    try:
        with open("trj_check.out", "r") as file:
            number_of_frames = next(
                (line.split()[1] for line in file if line.startswith("Step")), None
            )
    except FileNotFoundError:
        print(f"Error: File trj_check.out not found.")
        number_of_frames = None
    #conda_activate_path="/home/bio/ls/bin"
 
    number_of_frames = int(number_of_frames)
 
    #conda_gmxmmpbsa_name="gmxMMPBSA"
    forcefield="amber99sb-ildn"
    
    mmpbsa_inFILE="mmpbsa_LinearPB_amber99SB_ILDN.in"
    
    np_value = config['run']['num_processors']
    #gmx_mmpbsa(1, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, 
    #             mmpbsa_inFILE, MMPBSA_INFILE_PATH , np_value, number_of_frames)
    gmx_mmpbsa(cycle_number, conda_actiavte_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, mmpbsa_inFILE, MMPBSA_INFILE_PATH, np_value, number_of_frames)
    # data analysis
    NUMframe = "all"
    Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe)
    Data_Analysis_Cal(cycle_number, results_folder)
    # clean and move files
    clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_FOLDER, RootName, REMOVED_FILES_FOLDER, current_conf_path)
    #cycle_number += 1
    #conf_name = f"cycle{cycle_number}"
    #root_name = f"cycle{cycle_number}_BE"
    #cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")

def run_cycle(cycle_number, cycle_num, md_args, gmx_args):
    """
    deal with gmx_MMPBSA for current cycle and MD for the next cycle
    """
    process = []

    # current cycle gmx_mmpbsa
    gmx_process = multiprocessing.Process(target=gmx_mmpbsa_for_each_cycle, args=(folders[f"cycle{cycle_number}_MD"],cycle_number, *gmx_args))
    process.append(gmx_process)
    gmx_process.start()

    # if we have next cycle，run MD for the next cycle 
    if cycle_number < cycle_num:
        next_md_process = multiprocessing.Process(target=MD_for_each_cycle, args=(folders[f"cycle{cycle_number+1}_MD"], cycle_number + 1, *md_args))
        process.append(next_md_process)
        next_md_process.start()

    # all processes finished
    for p in process:
        p.join()

def get_last_chain_residue_count(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    # Model 1
    model = structure[0]

    # Last chain
    chains = list(model.get_chains())
    last_chain = chains[-1]

    # residue number
    residue_count = sum(1 for res in last_chain if res.id[0] == ' ')
    
    #return last_chain.id, residue_count
    return residue_count


def get_peptide_residue_info(pdb_file):
    """
    Get residue positions from the last chain in a PDB file as space-separated string.
    
    Parameters:
        pdb_file (str): Path to PDB file
    
    Returns:
        str: Space-separated residue positions in format "resnum:chain" (e.g., "99:B 100:B 101:B")
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_file)
    model = structure[0]
    chains = list(model.get_chains())
    last_chain = chains[-1]
    chain_id = last_chain.id
    
    # Process residues (standard residues only)
    residues = [res for res in last_chain if res.id[0] == ' ']
    
    # Generate space-separated string
    residue_string = ' '.join(f"{res.id[1]}:{chain_id}" for res in residues)
    
    return residue_string

############################################# LOAD FILES ###########################################

# Command-line argument parsing
parser = argparse.ArgumentParser(description='ANUBI Pipeline')
parser.add_argument('-i', '--infile', type=str, required=True,
                   help='Input YAML configuration file')

args = parser.parse_args()

# Load configuration from the provided YAML file
config_file = args.infile
config = load_config(config_file )

#PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.getcwd()

DATA_DIR = os.path.join(PROJECT_ROOT, "DATA")
#VMD_DIR = os.path.join(PROJECT_ROOT, "VMD_FUNCTION")
FUNCTION_DIR = os.path.join(PROJECT_ROOT, "FUNCTION")
FORCE_FIELD_PATH = os.path.join(PROJECT_ROOT, "FORCE_FIELD")
MMPBSA_INFILE_PATH = os.path.join(PROJECT_ROOT, "gmx_mmpbsa_in")
# pdb file
#protein_infile = "HLA_BiAB_protein_50ns" 
#protein_infile = "mtbind"


#protein_file_path = os.path.join(DATA_DIR, f"{protein_infile}.pdb")

make_mutation_modeller_py = os.path.join(FUNCTION_DIR,"MakeNewMutant_Modeller.py") 
# MDP files
ions_mdp_file = "ions"
minim_mdp_file = "minim"
nvt_mdp_file = "NVT"
npt_mdp_file = "NPT"
#samd_mdp_file = "SAMD"
md_mdp_file = "MD_MMPBSA"
only_protein_md_mdp_file = "Protein_MD"

ions_mdp_path = os.path.join(DATA_DIR, f"{ions_mdp_file}.mdp")
minim_mdp_path = os.path.join(DATA_DIR, f"{minim_mdp_file}.mdp")
nvt_mdp_path = os.path.join(DATA_DIR, f"{nvt_mdp_file}.mdp")
npt_mdp_path = os.path.join(DATA_DIR, f"{npt_mdp_file}.mdp")
#samd_mdp_path = os.path.join(DATA_DIR, f"{samd_mdp_file}.mdp")
md_mdp_path = os.path.join(DATA_DIR, f"{md_mdp_file}.mdp")
only_protein_md_mdp_path = os.path.join(DATA_DIR, f"{only_protein_md_mdp_file}.mdp")
                                                     
ROOT_OUTPUT = create_output_directory()

# create a CSV form to get the results directly
final_csv_name = "MonteCarlo_OUTPUT.csv"
final_csv_columns = ["Mutant_Name", "Binding Energy", "delta_delta_G", "Status", "Sequence"]
final_csv_file = os.path.join(ROOT_OUTPUT, final_csv_name)

logging.basicConfig(
    filename = "OUTPUT.out",
    level = logging.INFO,
    format="%(asctime)s - %(levelname)s -%(message)s"
)

logging.info(f"ROOT FOLDER PATH: {ROOT_OUTPUT}")

####################################### CHECK SYSTEM #############################################
#config_file = os.path.join(PROJECT_ROOT, 'infile.yaml')
#config = load_config(config_file )
conda_actiavte_path = config['Basic_setting']['conda_activate_script_path']
#VMD_path = config['Basic_setting']['VMD_path']
gmx_path = config['Basic_setting']['GROMACS_executable_path']

conda_gmxmmpbsa_name = config['Basic_setting']['conda_gmx_MMPBSA_name']
conda_modeller_name = config['Basic_setting']['conda_Modeller_name']


# check FINAL_CSV
if not os.path.exists(final_csv_file):
    pd.DataFrame(columns=final_csv_columns).to_csv(final_csv_file, index=False)
    print(f"Created MonteCarlo_PROCESS_OUTPUT RESULTS file")
else:
    print(f"MonteCarlo_PROCESS_OUTPUT RESULTS file found")

# check conda
if not os.path.isfile(conda_actiavte_path):
    logging.error(f"ERROR: cannot find conda activate scrpt path as {conda_actiavte_path}")
else:
    logging.info(f"conda activate path --> {conda_actiavte_path} version: {get_version('conda')}")
    
# check python
python_version = get_version("python")
logging.info(f"Python --> {which_program('python')} version: {python_version}")

'''
# check VMD
if VMD_path:
    if os.path.isfile(VMD_path) and os.access(VMD_path, os.X_OK):
        logging.info(f"VMD path --> {VMD_path} ")
    else:
        logging.error(f"ERROR: cannot find VMD path as {VMD_path}")
else:
    logging.info(f"VMD path --> {which_program('vmd')} ")
'''

# check gromacs
if gmx_path:
    if os.path.isfile(gmx_path) and os.access(gmx_path, os.X_OK):
        logging.info(f"gmx path --> {gmx_path} ")
    else:
        logging.error(f"ERROR: cannot find gmx path as {gmx_path}")
else:
    logging.info(f"gmx path --> {which_program('gmx')} ")

# check gmx_mmpbsa

# check modeller
if check_conda_env(conda_gmxmmpbsa_name):
    logging.info(f"{conda_gmxmmpbsa_name} is installed")
else:
    logging.error(f"No {conda_gmxmmpbsa_name} founded, please check that")
'''
if check_conda_env(conda_modeller_name):
    logging.info(f"{conda_modeller_name} is installed")
else:
    logging.error(f"No {conda_modeller_name} founded, please check that")
'''
# check parameters
receptorFRAG = str(config['gmx_mmpbsa']['receptorFRAG'])
ABchains = str(config['gmx_mmpbsa']['ABchains'])
startingFrameGMXPBSA = config['gmx_mmpbsa']['startingFrameGMXPBSA']
pipeline_mode = config['input_files']['pipeline_mode']

min_peptide = config['run']['peptide_mode_length_min']
max_peptide = config['run']['peptide_mode_length_max']
frequency_peptide = config['run']['peptide_mode_frequency_control']
#protein_infile = config['input_files']['structure_infile_name']
protein_file_path = config['input_files']['structure_file_path']
if not os.path.exists(protein_file_path):
    raise ValueError(f"No path {protein_file_path}")
else:
    logging.info(f"The structure file is {protein_file_path}")
protein_infile = os.path.basename(protein_file_path)
protein_infile, _ =os.path.splitext(protein_infile)


max_mutant = config['modeller']['max_mutant']
cycle_num = config['modeller']['cycle_num'] # the run cycle numbers for each configuration  Default:10
MUTANT_signal = False
#Stored Average BE from the last configuration. - Default: no
#Stored_AVG= -92.8

#Stored BE standard deviation from the last configuration. - Default: no
#Stored_STD= 4.3
#Metropolis Temperature - Default: 2  1.5
#Metropolis_temp = 1
Metropolis_temp = config['run']['Metropolis_Temperature']
#Metropolis Temperature top limit - Default: 4
Metropolis_Temp_cap= 4

#Metropolis Temperature Used during the calculations. It could change. - Default: 
Eff_Metropolis_Temp= Metropolis_temp

#Number of consecutive discarded results. - Default: 0
Consecutive_DISCARD_Count= 4

MP_correction = False

attempts = 0

peptide_length = get_last_chain_residue_count(protein_file_path )

######################################### MAIN PROCESS ######################################


for sequence in range (0,max_mutant+1):
    try:
        os.chdir(ROOT_OUTPUT)
    except OSError:
        logging.error(f"Cannot enter {ROOT_OUTPUT} folder")
        exit()

    if sequence == 0:
        
        # create configuration folder
        configuration_path = os.path.join(ROOT_OUTPUT,"configuration")
        os.mkdir(configuration_path)
        current_path_store = configuration_path
        
        print(f"Create directory: {configuration_path}")
        os.chdir(configuration_path)
        logging.info(f"#### Begin with configuration{sequence} ####")
        logging.info(f"PATH : {configuration_path}")
    else:
        mutant_folder_path = os.path.join(ROOT_OUTPUT,f"Mutant{sequence}")
        os.mkdir(mutant_folder_path)
        current_path_store = mutant_folder_path
        
        print(f"Create directort: {mutant_folder_path}")
        os.chdir(mutant_folder_path)
        logging.info(f"#### Begin with Mutant{sequence} ####")
        logging.info(f"PATH : {mutant_folder_path}")

    if MUTANT_signal == True:

        new_mutant = True
        while new_mutant == True:
            attempts += 1
            pdb_file = os.path.join(ROOT_OUTPUT,f"{protein_infile}.pdb") #LastFRame_xxxx.pdb

            if pipeline_mode == 'peptide' and attempts % (frequency_peptide+1) == 0 and min_peptide <= peptide_length <= max_peptide:

                peptide_options = ["ADD_FIRST", "ADD_END", "REMOVE_LAST", "REMOVE_FIRST"]
                p_choice = random.choice(peptide_options)
                if p_choice in ["ADD_FIRST", "ADD_END"]:
                    peptide_length += 1
                else:
                    peptide_length -= 1
                
                logging.info(f"Making a mutation of {p_choice}")
                #print(f"The option for this configuration is {p_choice}")
                peptide_mode(pdb_file,1,p_choice,sequence)
                sequence_log = p_choice
                
            else:
                

                # pdb_file, res_position, chain, new_restype, res_pos_list,res_weight_files, new_restype_list, keep_hydration, output_name
                #pdb_file = os.path.join(ROOT_OUTPUT,f"{protein_infile}.pdb") #LastFRame_xxxx.pdb
                #res_position = None
                #chain = None
                #new_restype = None
                #res_pos_list = config['modeller']['res_pos_list']
                if pipeline_mode == 'peptide':
                    res_pos_list = get_peptide_residue_info(pdb_file)
                    logging.info(f"peptide sequence now:{res_pos_list}")
                else:
                    res_pos_list = config['modeller']['res_pos_list']

                #new_restype_list = ['LEU', 'VAL', 'ILE', 'MET', 'PHE', 'TYR', 'TRP','GLU', 'ASP','ARG', 'LYS','SER', 'THR', 'ASN', 'GLN', 'HIS']

                # NO CYS GLY PRO
                new_restype_list = ['ALA', 'ARG', 'ASN', 'ASP', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MET']
                output_name = f"Mutant{sequence}"

                logging.info("Making a new mutation.")
                #keep_hydration = False
                #make_new_mutation(pdb_file, res_position, chain, new_restype, res_pos_list,res_weight_files, new_restype_list, keep_hydration, output_name)
                command_mutant = (f"python {make_mutation_modeller_py} {pdb_file} -o ./Mutant{sequence} -rl {res_pos_list} -v")
                subprocess.run(command_mutant, shell =True, check = True)
                #attempts += 1
            
                sequence_mutant_log = "OUTPUT_MUTANT.log"
           
                with open(sequence_mutant_log, 'r') as file:
                    sequence_log = file.readline().strip()  
            
            new_mutant =False
        
        # Check that the sequence hasn't be tested already (self avoiding walk)
        #os.remove(pdb_file)
        protein_infile= f"Mutant{sequence}"

        protein_file_path= os.path.join(current_path_store, f"{protein_infile}.pdb")
        #destination_file = os.path.join(current_path_store, f"{protein_infile}_noH.pdb")
        #shutil.copy(protein_file_path, destination_file)


    Metropolis_flag= 0
    current_dir = os.getcwd()
    folders = build_folders(current_dir,cycle_num)

    # generating a topology and build box
    #make_top_protein(protein_file_path, "amber14sb", "tip3p", "system", "topol", gmx_path)
    make_top_protein(protein_file_path, "amber99sb-ildn", "tip3p", "system", "topol", gmx_path)

    # cp system.pdb {protein_infile}.pdb in current folder
    source = os.path.join(current_dir, "system.pdb")
    destination = os.path.join(current_dir, f"{protein_infile }.pdb")
    try:
        shutil.copy(source,destination)
    except Exception:
        print("Copy system.pdb failed.")

    add_ter_to_pdb(f"{protein_infile }.pdb")
    output_pdb = os.path.join(ROOT_OUTPUT, f"{protein_infile}.pdb")
    replace_his_residues_flexible(f"{protein_infile}.pdb",output_pdb)
    # Adding water and ions
    fill_water_ions("system", "topol", ions_mdp_path, gmx_path)
    # Energy Minimiization
    energy_min(minim_mdp_path, "system_ions", "topol", "system_compl",gmx_path)

    make_equi("system_compl_minim.gro", nvt_mdp_path, npt_mdp_path, "system_equil", sequence, gmx_path)
    #make_nvt("system_compl_minim.gro", nvt_mdp_path, "system_equil", sequence, gmx_path)
    # Move .cpt, .top, and .itp files to repository folder
    for file_pattern in [f"{current_dir}/*.cpt", f"{current_dir}/*.top", f"{current_dir}/*.itp"]:
        for file in glob.glob(file_pattern):
            shutil.move(file, folders["repository"])

    # Move specific files to repository folder
    shutil.move(f"{current_dir}/{protein_infile}.pdb", folders["repository"])
    shutil.move(f"{current_dir}/system_compl_minim.gro", folders["repository"])
    shutil.move(f"{current_dir}/system_equil.gro", folders["repository"])


    # Move temp* and *out files to removed files folder
    for file in glob.glob("./*temp*.*") + glob.glob("./*.temp") + glob.glob("./*out"):
        shutil.move(file, folders["REMOVED_FILES_FOLDER"])

    # Remove files with # in their name
    for file in glob.glob("./#*"):
        os.remove(file)

    md_args = (sequence, md_mdp_path, "system_Compl_MD", "traj_MD", f"{gmx_path}")
    gmx_args = (only_protein_md_mdp_path,folders["TEMP_FILES_FOLDER"], FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, folders["REMOVED_FILES_FOLDER"], folders["results"], folders["repository"], current_path_store)
    # 1st cycle MD
    MD_for_each_cycle(folders["cycle1_MD"],1, *md_args)

    # each cycle: gmx_mmpbsa and next MD
    for cycle_number in range(1, cycle_num + 1):
        run_cycle(cycle_number, cycle_num, md_args, gmx_args)

    
    last_cycle_MD_FOLDER = os.path.join(folders["repository"],f"cycle{cycle_number}_MD")
    last_cycle_gro = os.path.join(last_cycle_MD_FOLDER,f"LastFrame_cycle{cycle_number}.gro")
    shutil.copy(last_cycle_gro, os.path.join(folders["repository"],f"LastFrame_cycle{cycle_number}.gro"))
    logging.info(f"Making the starting PDB for the next Mutation from LastFrame_cycle{cycle_number}.gro")
    
    os.chdir(current_path_store)
   
    repository_pdb_file = os.path.join(folders["repository"], f"{protein_infile}.pdb")
    '''
    pathGRO = folders["repository"]
    fileNameGRO = f"LastFrame_cycle{cycle_number}"
    pathPDB = os.path.dirname(repository_pdb_file)
    pdb_name_with_extension = os.path.basename(repository_pdb_file) #xxxx.pdb
    pdb_name_without_extension = os.path.splitext(pdb_name_with_extension)[0] #xxxx
    fileNamePDB = pdb_name_without_extension
    FileNamePDB_OUT = f"LastFrame_cycle{cycle_number}"

    GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, FileNamePDB_OUT, VMD_DIR, folders["TEMP_FILES_FOLDER"], VMD_path)
    '''

    '''
    #gmx trjconv -s system_Compl_MD.tpr -f traj_MD.xtc -o last_frame.pdb -dump 9999 -pbc mol -ur compact
    # work in repository
    ref_pdb_extract = repository_pdb_file
    tpr_file_extract = os.path.join(last_cycle_MD_FOLDER,"system_Compl_MD.tpr")
    xtc_file_extract = os.path.join(last_cycle_MD_FOLDER,"traj_MD.xtc")
    time_ps_extract = "9999"
    tmp_pdb_extract = os.path.join(folders["repository"],"lastframe.pdb")
    out_pdb_extract = os.path.join(folders["repository"],f"LastFrame_cycle{cycle_number}.pdb")
    extract_lastframe_and_rename(gmx_path,ref_pdb_extract,tpr_file_extract,xtc_file_extract,time_ps_extract,tmp_pdb_extract,out_pdb_extract)
    '''
    cycle10_start_pdb = os.path.join(last_cycle_MD_FOLDER, f"cycle{cycle_number}_starting_protein.pdb")
  
    last_cycle_pdb = os.path.join(folders["repository"], f"LastFrame_cycle{cycle_number}.pdb")
    shutil.copy(cycle10_start_pdb,last_cycle_pdb)
    add_ter_to_pdb(last_cycle_pdb)        
    output_last_cycle_pdb = os.path.join(ROOT_OUTPUT, f"Mutant{sequence}_cycle{cycle_number}_LastFrameMD.pdb")
    replace_his_residues_flexible(last_cycle_pdb,output_last_cycle_pdb)

    protein_infile = f"Mutant{sequence}_cycle{cycle_number}_LastFrameMD"
    logging.info("Making the average of the cycles results.")
    
    all_cycle_data = "All_cycle_data.out"
    MoleculesResults_data = os.path.join(folders["results"], "MoleculesResults.dat")
    data_analysis_temp = "DataAnalysis_temp.csv"

    if os.path.exists(all_cycle_data):
        os.remove(all_cycle_data)

    flag_header = True

    with open(MoleculesResults_data, 'r') as infile, open(all_cycle_data, 'w') as outfile:
        for line in infile:
            if flag_header == True:
                # head row
                outfile.write(f"#{'configNum':<10} \t{line}")
                flag_header = False
            else:
                # data row
                outfile.write(f"{'avg':<10} \t{line}")
    # get the data from the second line
    df = pd.read_csv(all_cycle_data, sep='\t')  
    df_filtered = df.iloc[:, 1:]  # get the data from the second column

    # save to csv
    df_filtered.to_csv(data_analysis_temp, sep = '\t', index=False,header = False)
    #Data_Analysis_Signal = False
    Data_Analysis_Cal_child(data_analysis_temp, "AllData.temp", False)

    frame_count = 0
    # get AVG and STD to AllData.out
    with open("AllData.temp", 'r') as temp_file, open(all_cycle_data, 'a') as outfile:
        for line in temp_file:

            if line.startswith("#frame"):
                frame_count +=1
                if frame_count ==2:
                    outfile.write(line)
            if line.startswith("#AVG") or line.startswith("#STD"):
                outfile.write(line)
    '''
    # remove temp file
    for temp_file in [all_data_temp, data_analysis_temp]:
        shutil.move(temp_file, os.path.join(removed_files_folder, os.path.basename(temp_file)))
    '''
    frame = []
    avg = []
    std = []
    with open(all_cycle_data, 'r') as infile:
        for line in infile:
            if line.startswith("#frame"):
                frame = line.strip().split()[1:]
            elif line.startswith("#AVG"):
                avg = line.strip().split()[1:]
            elif line.startswith("#STD"):
                std = line.strip().split()[1:]
    if frame and avg and std:
        if sequence == 0:
            output_lines = ["Results for Configuation"]
        output_lines = [f"Results for Mutant{sequence}"]
        output_lines += [f"{frame[i]}: {avg[i]} +- {std[i]} kcal/mol"
                        for i in range(len(frame))
                       ]
        logging.info("\n".join(output_lines))
    shutil.move("AllData.temp", folders["REMOVED_FILES_FOLDER"])
    shutil.move(data_analysis_temp, folders["REMOVED_FILES_FOLDER"])
    AVG = float(avg[0])
    STD = float(std[0])
    if sequence == 0:
        MUTANT_signal = True
        Stored_AVG = float(avg[0]) # DeltaG(kcal/mol)
        Stored_STD = float(std[0])
        Stored_system_file = protein_infile
        logging.info(f"Finished with Configuration{sequence}")
        
        MonteCarlo_row = pd.DataFrame([{
            "Mutant_Name": f"Wild Type",
            "Binding Energy": AVG,
            "delta_delta_G": " ",
            "Status": " ",
            "Sequence": " "
        }])
        # Append to CSV
        MonteCarlo_row.to_csv(final_csv_file, mode='a', header=not pd.io.common.file_exists(final_csv_file), index=False)
        sequence += 1 
        os.chdir(ROOT_OUTPUT)
        # if FAST == TRue delete removed files folder
        continue


    logging.info("Metropolis algorithm")
    Prob = None
    if not Prob:
        RandNum = random.uniform(0,1)
    else:
        RandNum = float(Prob)
    # Metropolis
    if MP_correction == True:
        MP = math.exp(-(AVG+STD/2-Stored_AVG) / Eff_Metropolis_Temp)
    else:
        MP = math.exp(-(AVG-Stored_AVG) / Eff_Metropolis_Temp)
    # new G < old G
    if MP >= 1:
        MP = 1
    delta_delta_G = AVG - Stored_AVG



    logging.info(f"Random Number: {RandNum}  Metropolis Prob: {MP}  AVG: {AVG}  Stored AVG: {Stored_AVG} delta_delta_G:{delta_delta_G}")
    Metropolis_flag = 1 if RandNum < MP else 0


    if Metropolis_flag == 1:
        logging.info("New Configuration Accepted")
        Stored_AVG = AVG
        Stored_STD = STD
        Stored_system_file = protein_infile
        Consecutive_DISCARD_Count = 0
        Eff_Metropolis_Temp = Metropolis_temp
    else:
        logging.info("New Configuration Declined")
        protein_infile = Stored_system_file
        Consecutive_DISCARD_Count += 1
        '''
        if Consecutive_DISCARD_Count > 5:
            if Eff_Metropolis_Temp < Metropolis_Temp_cap:
                Eff_Metropolis_Temp += 0.5*(Consecutive_DISCARD_Count - 5)
        if Eff_Metropolis_Temp > Metropolis_Temp_cap:
            Eff_Metropolis_Temp = Metropolis_Temp_cap
        '''
    
    # log in MonteCarlo Process
    MonteCarlo_row = pd.DataFrame([{
        "Mutant_Name": f"Mutant{sequence}",
        "Binding Energy": AVG,
        "delta_delta_G": delta_delta_G,
        "Status": "Acc" if Metropolis_flag == 1 else "Dec",
        "Sequence": sequence_log  
    }])
    # Append to CSV
    MonteCarlo_row.to_csv(final_csv_file, mode='a', header=not pd.io.common.file_exists(final_csv_file), index=False)


    MUTANT_signal = True
    logging.info(f"Finished Mutant{sequence}")
    sequence += 1
    os.chdir(ROOT_OUTPUT)



logging.info("ALL DONE.")
    

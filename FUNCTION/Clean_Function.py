import os
import shutil
import logging



#  NO SUCH FILE NEED TO IMPROVE
def clean_up(removed_files_folder):
    logging.info("Cleaning files.")
    try:
        for pattern in ["*#*", "*~*", "*temp*"]:
            for filename in os.listdir('.'):
                if filename.startswith(pattern):
                    file_path = os.path.join('.', filename)
                    if os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                    else:
                        os.remove(file_path)
        
        with open(os.path.join(removed_files_folder, 'rm.out'), 'w') as f:
            f.write("Cleanup completed.")
    except Exception as e:
        logging.error(f"Error during cleanup: {e}")

# move and copy
def move_and_copy_files(cycle_number, repository_folder, cycle_number_MD_folder, root_name):
    logging.info("Moving files.")
    try:
        # move
        for filename in os.listdir('.'):
            if filename.startswith(f'RUN1_cycle{cycle_number}_prot') or filename.endswith('.out') or filename == root_name:
                shutil.move(filename, os.path.join(repository_folder, filename))
        
        # copy
        #shutil.copy('system_Compl_MDstart.gro', os.path.join(repository_folder, f'system_cycle{cycle_number}_MDstart.gro'))
        if os.path.isdir(cycle_number_MD_folder):
            shutil.copytree(cycle_number_MD_folder, os.path.join(repository_folder, f'cycle{cycle_number}_MD'))
        
    except Exception as e:
        logging.error(f"Error during file operations: {e}")

def clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_folder, root_name, removed_files_folder, current_conf_path):
    clean_up(removed_files_folder)
    move_and_copy_files(cycle_number, repository_folder, cycle_number_MD_folder, root_name)



    # move to removed_files_folder
    for filename in os.listdir('.'):
        if os.path.isfile(filename):  
            shutil.move(filename, os.path.join(removed_files_folder, filename))

    # remove
    for filename in os.listdir('.'):
        file_path = os.path.join('.', filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    # delete
    cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")
    if os.path.exists(cycle_number_md_folder):
        shutil.rmtree(cycle_number_md_folder)


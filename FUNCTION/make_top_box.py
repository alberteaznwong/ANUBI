import os
import subprocess
import glob
import logging

def make_top_protein(input_file_path, forcefield, watermodel, protein_outfile, topfile, gmx_path):
    """
    generate GROMACS topology (.top and .itp) from .pdb
    buid box to get .gro

    parameters:
        protein_infile (str): input file PDB name;
        forcefield (str): forcrfild used;
        watermodel (str): water model used;
        protein_outfile (str): output filr GRO name (no extension);
        topfile (str): TOP name (no extension);
    """
   
    # output logging 
    out_file = "MakeTOP_protein.out"


    # use pdb2gmx to generate topology
    pdb2gmx_cmd = f"{gmx_path} pdb2gmx -f {input_file_path} -o system.pdb -p {topfile}.top -ignh -ff {forcefield} -water {watermodel}"
    try:
        subprocess.run(pdb2gmx_cmd, shell=True, check=True, stdout=open(out_file, 'a'), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with pdb2gmx!")
        return


    #print("Topology generation completed successfully.")
    logging.info("Topology generation completed successfully.")

    # input and ouput for buiding box
    input_file = "system.pdb"
    output_file = f"{protein_outfile}.gro"

    # box type
    #editconf_option = "-bt triclinic -d 1.5"  
    editconf_option = "-c -bt cubic -d 1.5"

    try:
        # run editconf to buid box
        editconf_cmd = f"{gmx_path} editconf -f {input_file} {editconf_option} -o {output_file}"
        # log
        with open(out_file, 'a') as log:
            subprocess.run(editconf_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with editconf!")

    #print("Simulation box definition completed successfully.")
    logging.info("Simulation box definition completed successfully.")
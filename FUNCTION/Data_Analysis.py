import csv
import os
import shutil
import logging
import numpy as np
import pandas as pd
import math

def extract_delta_energy_terms(input_filename, output_filename):
    try:
        with open(input_filename, 'r') as infile:
            lines = infile.readlines()  
    except FileNotFoundError:
        print(f"Error: {input_filename} file not found.")
        return

    # search Delta Energy Terms 
    delta_start = False
    delta_lines = []

    for line in lines:
        # search Delta Energy Terms 
        if "Delta Energy Terms" in line:
            delta_start = True
            continue  
        # get data
        if delta_start:
            if line.strip() == "":  # end
                break
            delta_lines.append(line.strip().split(','))  # split by ','

    if not delta_lines:
        print("Error: 'Delta Energy Terms' data not found.")
        return

    # write in new csv 
    try:
        with open(output_filename, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(delta_lines)  # write data
        print(f"Delta Energy Terms extracted and saved to {output_filename}")
    except Exception as e:
        print(f"Error saving the file: {e}")

def append_results(cycle_number, results_folder):
    
    output_file = os.path.join(results_folder, "MoleculesResults.dat")
    
    input_file = os.path.join(results_folder, f"cycle{cycle_number}_results.dat")
    
    
    with open(output_file, "a") as outfile:
        
        outfile.write(f"{cycle_number:<5}\t")
        
        
        if os.path.exists(input_file):
            with open(input_file, "r") as infile:
                for line in infile:
                    if "#AVG" in line:
                        
                        columns = line.strip().split("\t")[1:]
                        outfile.write("\t".join(columns) + "\n")

def Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe = "all"):
    input_filename = "gmx_MMPBSA_plot.csv"  
    output_filename = "delta_energy_terms.csv"  
    extract_delta_energy_terms(input_filename, output_filename)
    # read the csv file
    df = pd.read_csv('delta_energy_terms.csv')

    # create df
    output_data = []
    for index, row in df.iterrows():
        # for each line
    
        # Frame #
        frame = row['Frame #']  
    
        # DeltaG(kcal/mol)
        delta_g = row['TOTAL']  
    
        # Coul(kcal/mol) = EEL + EEL14
        EEL = row['EEL']  
        EEL14 = row['1-4 EEL']  
        coul = EEL + EEL14
    
        # vdW(kcal/mol) = VDW + VDW14
        VDW = row['VDWAALS']  
        VDW14 = row['1-4 VDW']  
        vdW = VDW + VDW14
    
        # PolSol(kcal/mol) = EPB + ENPOLAR (if have)
        # for this infile the NpoSol is the value of [ENPOLAR], so it is included in the PolSol

        if pd.notna(row['EPB']) and pd.notna(row['ENPOLAR']):
            EPB = row['EPB']
            ENPOLAR = row['ENPOLAR']
            pol_sol = EPB + ENPOLAR
        else:
            pol_sol = row['EGB']  # if have
    
        # NpoSol(kcal/mol) = EDISPER or ESURF (if have)
        non_pol_solv = 0
        if pd.notna(row['EDISPER']):
            non_pol_solv = row['EDISPER']
        elif pd.notna(row['ESURF']):
            non_pol_solv = row['ESURF']
    
        # add the values 
        output_data.append({
            '# frame': frame,
            'DeltaG(kcal/mol)': delta_g,
            'Coul(kcal/mol)': coul,
            'vdW(kcal/mol)': vdW,
            'PolSol(kcal/mol)': pol_sol,
            'NpoSol(kcal/mol)': non_pol_solv
        })

    # transfer
    output_df = pd.DataFrame(output_data)

    # write it in
    output_df.to_csv('energy_plot_temp.csv', sep='\t', index=False)

    
    #NUMframe = "all"  

    #logging.info("Cleaning files.")
    # energy_plot_temp 
    with open('./energy_plot_temp.csv', 'r') as f:
        lines = f.readlines()

   
    #header = lines[0]  

    if NUMframe == "all":
        selected_lines = lines[1:]  # if "all"，save all last 150 frames, no the frame at 2ns, so it should be lines[2:], if with boundary, it should be lines[1:]
    else:
        selected_lines = lines[int(NUMframe):]  


    with open('../energy_plot_temp.csv', 'w') as f:
        f.writelines(selected_lines)



    # clean files
    for filename in os.listdir('.'):
        if filename.startswith('_GMXMMPBSA'):
            shutil.move(filename, os.path.join(REMOVED_FILES_FOLDER, filename))
        elif '_temp' in filename:
            shutil.move(filename, os.path.join(REMOVED_FILES_FOLDER, filename))

  
    for filename in os.listdir('.'):
        if filename.endswith(('#', '~', '.ff')):
            os.remove(filename)

 
    with open(os.path.join(REMOVED_FILES_FOLDER, 'rm.out'), 'w') as rm_out_file:
        rm_out_file.write("Removed backup and .ff files.")

   
    try:
        os.chdir(cycle_number_MD_FOLDER)
    except FileNotFoundError:
        raise SystemExit(f"Cannot enter '{cycle_number_MD_FOLDER}' folder")


def Data_Analysis_Cal_child(input_file, output_file, Data_Analysis_Signal = True):
    
    
    #df = pd.read_csv(input_file, delim_whitespace=True, header=None)
    df = pd.read_csv(input_file, sep=r'\s+', header=None)
    if Data_Analysis_Signal == True:
        df.columns = ['frame', 'DeltaG', 'Coul', 'VdW', 'PolSol', 'NpoSol']
    else:
        #df.columns = ['frame', 'DeltaG', 'Coul', 'VdW', 'PolSol', 'NpoSol', 'SF1', 'SF2', 'Canonical_AVG', 'MedianDG', 'DeltaG_2s']
        df.columns = ['frame', 'DeltaG', 'Coul', 'VdW', 'PolSol', 'NpoSol', 'SF1', 'SF2', 'MedianDG', 'DeltaG_2s']

   
    Population = len(df)
    DeltaG = df['DeltaG'].sum()
    Coul = df['Coul'].sum()
    VdW = df['VdW'].sum()
    PolSol = df['PolSol'].sum()
    NpoSol = df['NpoSol'].sum()

    # calculate SF1 , SF2
    df['SF1'] = (df['Coul'] / 10) - (df['PolSol'] / 10) + (df['NpoSol'] * 10)
    df['SF2'] = (3 * df['Coul']) + df['PolSol']

    '''
    Canonical_AVG = 0.0
    Canonical_AVG_w = 0.0

    for _, row in df.iterrows():
        
        #DeltaG_temp = row[1]  # get deltaG_temp
        DeltaG_temp = row.iloc[1]
        # calculate weight 1 KT =2.479 KJ/mol ????
        weight = math.exp(-int(DeltaG_temp / 2.479))
        
        Canonical_AVG += DeltaG_temp * weight
        Canonical_AVG_w += weight

    # 
    if Canonical_AVG_w != 0:
        Canonical_AVG /= Canonical_AVG_w
    '''
    # 
    mean_DeltaG = DeltaG / Population
    mean_Coul = Coul / Population
    mean_VdW = VdW / Population
    mean_PolSol = PolSol / Population
    mean_NpoSol = NpoSol / Population

    mean_SF1 = mean_Coul / 10 - mean_PolSol / 10 + mean_NpoSol * 10
    mean_SF2 = (3 * mean_Coul) + mean_PolSol


    # std
    std_DeltaG = np.std(df['DeltaG'], ddof=1)
    std_Coul = np.std(df['Coul'], ddof=1)
    std_VdW = np.std(df['VdW'], ddof=1)
    std_PolSol = np.std(df['PolSol'], ddof=1)
    std_NpoSol = np.std(df['NpoSol'], ddof=1)
    std_SF1 = std_Coul / 10 - std_PolSol / 10 + std_NpoSol * 10
    std_SF2 = (3 * std_Coul) + std_PolSol

    # calculate  DeltaG within 2sigma
    DeltaG_2s = df[(np.abs(df['DeltaG'] - mean_DeltaG) < 2 * std_DeltaG)]['DeltaG'].mean()

    #
    median_DeltaG = df['DeltaG'].median()
    # 
    std_DeltaG = np.std(df['DeltaG'], ddof=1)  # ddof=1 delta degree of freedom; devide by n-1; sample standard deviation


    
    warnings = []

    # check if it > 2 sigma
    for i, row in df.iterrows():
        var = np.abs(row['DeltaG'] - mean_DeltaG)
        if var >= 2 * std_DeltaG:
            warnings.append(f"# WARNING: frame {row['frame']} is out of 2 sigma!!")

    
    with open(output_file, 'w') as f:
       
        f.write("# SF1=Coulomb/10-PolarSolvation/10+Non-PolarSolvation*10\n")
        f.write("# SF2=3*Coulomb+PolarSolvation\n")
        #f.write("# C_AVG=norm(SUM Gi*e^BGi)\n")
        f.write(f"#frame\tDeltaG(kcal/mol)\tCoul(kcal/mol)\tVdW(kcal/mol)\tPolSol(kcal/mol)\tNpoSol(kcal/mol)\tSF1\tSF2\n")

        
        for _, row in df.iterrows():
            f.write(f"{row['frame']:<10}{row['DeltaG']:>12.3f}{row['Coul']:>13.3f}{row['VdW']:>13.3f}{row['PolSol']:>13.3f}{row['NpoSol']:>13.3f}{row['SF1']:>13.3f}{row['SF2']:>13.3f}\n")
    
         
        for warning in warnings:
            f.write(f"{warning}\n")   
   
        
        f.write("\n# FINAL RESULTS\n")
        #f.write(f"#frame\t{'DeltaG(kcal/mol)':>15}\t{'Coul(kcal/mol)':>15}\t{'VdW(kcal/mol)':>15}\t{'PolSol(kcal/mol)':>15}\t{'NpoSol(kcal/mol)':>15}\t{'SF1':>15}\t{'SF2':>15}\t{'Canonical_AVG':>15}\t{'MedianDeltaG(kcal/mol)':>15}\t{'DeltaG_2s(kcal/mol)':>15}\n")
        f.write(f"#frame\t{'DeltaG(kcal/mol)':>15}\t{'Coul(kcal/mol)':>15}\t{'VdW(kcal/mol)':>15}\t{'PolSol(kcal/mol)':>15}\t{'NpoSol(kcal/mol)':>15}\t{'SF1':>15}\t{'SF2':>15}\t{'MedianDeltaG(kcal/mol)':>15}\t{'DeltaG_2s(kcal/mol)':>15}\n")
        #f.write(f"#AVG\t{mean_DeltaG:>15.1f}\t{mean_Coul:>15.1f}\t{mean_VdW:>15.1f}\t{mean_PolSol:>15.1f}\t{mean_NpoSol:>15.1f}\t{mean_SF1:>15.1f}\t{mean_SF2:>15.1f}\t{Canonical_AVG:>15.1f}\t{median_DeltaG:>15.1f}\t{DeltaG_2s:>15.1f}\n")
        #f.write(f"#STD\t{std_DeltaG:>15.1f}\t{std_Coul:>15.1f}\t{std_VdW:>15.1f}\t{std_PolSol:>15.1f}\t{std_NpoSol:>15.1f}\t{std_SF1:>15.1f}\t{std_SF2:>15.1f}\t{'nan':>15}\t{'nan':>15}\t{std_DeltaG:>15.1f}\n")    
        f.write(f"#AVG\t{mean_DeltaG:>15.1f}\t{mean_Coul:>15.1f}\t{mean_VdW:>15.1f}\t{mean_PolSol:>15.1f}\t{mean_NpoSol:>15.1f}\t{mean_SF1:>15.1f}\t{mean_SF2:>15.1f}\t{median_DeltaG:>15.1f}\t{DeltaG_2s:>15.1f}\n")
        f.write(f"#STD\t{std_DeltaG:>15.1f}\t{std_Coul:>15.1f}\t{std_VdW:>15.1f}\t{std_PolSol:>15.1f}\t{std_NpoSol:>15.1f}\t{std_SF1:>15.1f}\t{std_SF2:>15.1f}\t{'nan':>15}\t{std_DeltaG:>15.1f}\n")    



def Data_Analysis_Cal(cycle_number, results_folder):
    
    logging.info("Data Analysis.")
    
    input_file = 'energy_plot_temp.csv'  
    # NEED TO CHANGE cycle${Cycle_Number}_results.dat
    output_file = f'cycle{cycle_number}_results.dat'
    Data_Analysis_Cal_child(input_file, output_file)
    

    try:
        with open(f'cycle{cycle_number}_results.dat', 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        logging.error(f"File not found: cycle{cycle_number}_results.dat")
        exit()

    
    header = None
    avg_values = None
    std_values = None

    
    for line in lines:
        if line.startswith("#frame"):
            header = line.strip().lstrip("#").split()
        elif line.startswith("#AVG"):
            avg_values = line.strip().lstrip("#").split()
        elif line.startswith("#STD"):
            std_values = line.strip().lstrip("#").split()

   
    if not header or not avg_values or not std_values:
        logging.error("Missing required data (header, AVG, or STD).")
        exit()

   
    log_message = (
        f"Results for cycle{cycle_number}:\n"
        f"Headers: {' | '.join(header)}\n"
        f"AVG Values: {' | '.join(avg_values)}\n"
        f"STD Values: {' | '.join(std_values)}")

   
    logging.info(log_message)
    # move outpufile (cycle1_results.dat) to results folder
    shutil.move(output_file, os.path.join(results_folder, output_file))
    # append cycle data on MoleculeResults.dat
    append_results(cycle_number, results_folder)
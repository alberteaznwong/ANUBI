import string
from itertools import product
from modeller import *
from modeller.automodel import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
import os
from Bio.PDB import PDBParser, PDBIO, Chain, Residue, Atom
from Bio.PDB.Structure import Structure
from Bio import PDB
import json
from FUNCTION.Structure_Build import process_pdb_safe

# function to change chain id (make it A B C D)

def right_order_chain(input_file,output_file,last_chain_id=''):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", input_file)
    model = structure[0]

    # Create the final chain ID (A, B, C...)

    available_ids = list(string.ascii_uppercase)
    if len(list(model.get_chains())) > 26:
        available_ids += [''.join(p) for p in product(string.ascii_uppercase, repeat=2)]

    #  Step 1: Temporarily replace all chain IDs with non-conflicting temporary IDs (to avoid conflicts)
    temp_id_map = {}
    for idx, chain in enumerate(model.get_chains()):
        temp_id = f"TMP{idx}"
        temp_id_map[chain.id] = temp_id
        chain.id = temp_id  

    # Step 2: Replace the temporary ID with the target chain ID (A, B, C...)
    chain_objects = list(model.get_chains())
    final_map = {}
    for chain_obj, new_id in zip(chain_objects, available_ids):
        old_id = chain_obj.id
        final_map[old_id] = new_id
        chain_obj.id = new_id

    
    # Output mapping relationships
    print("Chain ID mapping (final):")
    for k, v in final_map.items():
        print(f"{k} ➜ {v}")

    last_chain_id = chain_objects[-1].id
    # save new pdb
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)
    return last_chain_id



def right_order_chain(input_file, output_file, last_chain_id=''):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("model", input_file)
    model = structure[0]

    # Create the final chain ID (A, B, C...)
    available_ids = list(string.ascii_uppercase)
    if len(list(model.get_chains())) > 26:
        available_ids += [''.join(p) for p in product(string.ascii_uppercase, repeat=2)]

    #  Step 1: Temporarily replace all chain IDs with non-conflicting temporary IDs (to avoid conflicts)
    original_to_tmp = {}
    temp_id_map = {}
    for idx, chain in enumerate(model.get_chains()):
        temp_id = f"TMP{idx}"
        original_to_tmp[chain.id] = temp_id
        chain.id = temp_id  

    # Step 2: Replace the temporary ID with the target chain ID (A, B, C...)
    final_map = {}
    chain_objects = list(model.get_chains())
    for chain_obj, new_id in zip(chain_objects, available_ids):
        old_id = chain_obj.id
        final_map[old_id] = new_id
        chain_obj.id = new_id

    # Output chain ID mappings
    print("Chain ID mappings:")
    print("Original to Temporary:", original_to_tmp)
    print("Temporary to Final:", final_map)

    # Save mappings to JSON
    mappings = {
        "original_to_tmp": original_to_tmp,
        "tmp_to_original": {v: k for k, v in original_to_tmp.items()},
    }

    with open("chain_map.json", "w") as f:
        json.dump(mappings, f, indent=4)

    last_chain_id = chain_objects[-1].id
    # Save the modified PDB
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    return last_chain_id

def optimize(atmsel, sched):
    # conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    # md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


# molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (iterations, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                       (200, 600,
                                        (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                        max_iterations=iterations, equilibrate=equil)
            init_vel = False

def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    sel = Selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(sel, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(sel, restraint_type=typ + '_dihedral', spline_range=4.0,
                 spline_dx=0.3, spline_min_points=5, aln=aln,
                 spline_on_site=True)
        

# function to add ALA by using Modeller
# setting part, create ali 
def add_ALA_setting (modelname,chain,num_ala,p_choice):

    # Initializing the Modeller Environment
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read the original structure
    mdl = Model(env, file=f'{modelname}.pdb')

    # Getting the target chain and inserting the ALA residue
    try:
        target_chain = mdl.chains[chain]
    except KeyError:
        raise ValueError(f"Chain '{chain}' not found in {modelname}.pdb")

    # Get the sequence of the chain, ready to add the ALA
    #original_seq = ''.join(res.code for res in target_chain.residues)
    added_seq = 'A' * num_ala  # Add the number of ALAs

    # Creating Alignment Objects
    aln = Alignment(env)

    # Adding Template Models to Aligned Objects
    aln.append_model(mdl, atom_files=modelname, align_codes=modelname)

    # Write the extracted sequence to a .seq file
    aln.write(file=f'{modelname}_origin.seq')
    with open(f'{modelname}_origin.seq', 'r') as f:
        lines = f.read().splitlines()

    # find 'structureX:' 
    for i, line in enumerate(lines):
        if line.strip().startswith('structureX:'):
            temp_seq = ''.join(lines[i + 1:])  # Merge from the next line
            break
    temp_seq = temp_seq.replace('*','')
    print(f"original sequence :{temp_seq}")
    # Find the sequence after the last '/'
    if '/' in temp_seq:
        prefix, old_tail = temp_seq.rsplit('/', 1)
        if p_choice == 'ADD_FIRST':
            new_tail = added_seq + old_tail
        if p_choice == 'ADD_END':
            new_tail = old_tail + added_seq
        temp_seq = prefix + '/' + new_tail

    print(f"new sequence :{temp_seq}")
    ######## finished add ########
    with open(f'{modelname}_origin.seq', 'r') as f:
        lines = f.read().splitlines()

    for line in lines:
        if line.startswith('structureX:'):
            # Split, take out the target part
            parts = line.split(':')
            original_info = ':'.join(parts[2:])  # get"25:A:+116:B:::-1.00:-1.00"

            # Split in order to modify +116
            info_parts = original_info.split(':')
        
            # Item 3 is '+116', we take out the number part, add one and spell it back together
            num = int(info_parts[2].lstrip('+'))
            num += 1
            info_parts[2] = f'+{num}'

            # Combine to final string
            new_info = ':'.join(info_parts)
            break  # exit
    print(f"check the info:{new_info}")
    # 自动生成对齐文件内容

    alignment_ali = f""">P1;
sequence::{new_info}
{temp_seq}*
    """
    with open(f'{modelname}_origin.seq', 'r') as f:
        origin_content = f.read()  # Remove empty lines at the end
    combined_ali = origin_content + alignment_ali

    # save .ali
    with open('combined.ali', 'w') as f:
        f.write(combined_ali)

    aln.read(file='combined.ali')
    mdl.clear_topology()
    mdl.generate_topology(aln[-1])
    mdl.transfer_xyz(aln)

    # Build the remaining unknown coordinates
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.write(file=f'{modelname}_mutated.pdb')
    
# function to create pdb, with optimazation
def all_ALA_opt(modelname,last_chain_id,p_choice):
    # Initializing the Modeller Environment
    env = Environ()
    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')
   
    aln = Alignment(env)
    aln.read(file='combined.ali')
    # loading structure
    #mdl1 = Model(env, file='reumbered_output.pdb')  
    #f'{modelname}_mutated.pdb'
    mdl1 = Model(env, file=f'{modelname}_mutated.pdb')
    mdl1.clear_topology()
    mdl1.generate_topology(aln[-1])
    print(aln[-1])
    # select ALA
    if p_choice == 'ADD_FIRST':
        #s = Selection(mdl1.chains['B'].residues[-1:])
        s = Selection(mdl1.chains[last_chain_id].residues[0])
    else:
        s = Selection(mdl1.chains[last_chain_id].residues[-1])
    # clearing the previously set restraints
    make_restraints(mdl1, aln)
    # a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms = 1

    sched = autosched.loop.make_for_model(mdl1)
    print(s.energy())
    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

    #feels environment (energy computed on pairs that have at least one member
    #in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)
    print(s.energy())
    mdl1.write(file=f'{modelname}_finished.pdb')


def restore_chain_ids(input_file, output_file, mapping_file="chain_map.json"):
    # read map
    with open(mapping_file, "r") as f:
        mappings = json.load(f)

    tmp_to_original = mappings["tmp_to_original"]           # TMPx ➜ origin ID
    original_to_tmp = mappings["original_to_tmp"]           # origin ID ➜ TMPx
    tmp_to_final = {}                                       # TMPx ➜ now ID
    final_to_tmp = {}                                       # now ID ➜ TMPx（reverse）

    # create tmp_to_final and final_to_tmp
    tmp_indices = list(tmp_to_original.keys())              # TMP0, TMP1...
    final_chain_ids = list(PDB.PDBParser(QUIET=True).get_structure("model", input_file)[0].get_chains())
    for tmp, chain in zip(tmp_indices, final_chain_ids):
        tmp_to_final[tmp] = chain.id
        final_to_tmp[chain.id] = tmp

    #  restore id
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("model", input_file)
    model = structure[0]

    print("Restoring chain IDs:")
    for chain in model.get_chains():
        current_id = chain.id
        if current_id in final_to_tmp:
            tmp_id = final_to_tmp[current_id]
            if tmp_id in tmp_to_original:
                original_id = tmp_to_original[tmp_id]
                print(f"{current_id} ➜ {tmp_id} ➜ {original_id}")
                chain.id = original_id
            else:
                print(f"No original mapping for TMP id: {tmp_id}")
        else:
            print(f"No TMP mapping found for chain ID: {current_id}")

    # save
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    print(f"Chain IDs restored and saved to: {output_file}")

# function to delete residue
def remove_terminal_residue(input_pdb, output_pdb, peptide_chain_id, mode):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb)

    for model in structure:
        for chain in model:
            if chain.id == peptide_chain_id:
                residues = list(chain)
                if not residues:
                    print(f"Chain {peptide_chain_id} has no residues.")
                    return

                if mode == "REMOVE_LAST":
                    chain.detach_child(residues[-1].id)
                    print(f"Deleted last residue: {residues[-1].id}")
                elif mode == "REMOVE_FIRST":
                    chain.detach_child(residues[0].id)
                    print(f"Deleted first residue: {residues[0].id}")
                else:
                    raise ValueError("Please check your mode choice")

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Saved modified structure to {output_pdb}")

def get_last_chain_id(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    last_chain_id = None
    
    # Iterate through all models, residues, and chains
    for model in structure:
        for chain in model:
            last_chain_id = chain.id  # Update last_chain_id with the current chain ID
    
    return last_chain_id

def peptide_mode(input_file,num_ala,p_choice,sequence):
    input_file_NAME = os.path.basename(input_file)
    input_file_basename = os.path.splitext(input_file_NAME)[0]
    if p_choice in  ['ADD_FIRST','ADD_END']:
        #input_file = 'HLA.pdb'
        output_file = f"{input_file_basename}_RightOrder.pdb"
        # make the chain ID (A B C D ) get f"{input_file_basename}_RightOrder.pdb"
        last_chain_id = right_order_chain(input_file,output_file)

        output_file_basename = os.path.splitext(output_file)[0]
        modelname = output_file_basename
        #chain = 'D'  # target chian id 
        chain = last_chain_id
        num_ala = num_ala  # ALA number
        p_choice = p_choice
        # add ALA, get combined.ali and f'{modelname}_mutated.pdb'
        add_ALA_setting (modelname,chain,num_ala,p_choice)
        # optimazition,get '{modelname}_finished.pdb'
        all_ALA_opt(modelname,chain,p_choice)
        # reorder, get f'{modelname}_reorder.pdb'
        process_pdb_safe(f"{input_file_basename}_RightOrder.pdb",f'{modelname}_finished.pdb',f'{modelname}_reorder.pdb')
        # restore id, get f'{modelname}_RESTORE.pdb'
        #restore_chain_ids(f'{modelname}_reorder.pdb',f'{modelname}_RESTORE.pdb')
        restore_chain_ids(f'{modelname}_reorder.pdb',f'Mutant{sequence}.pdb')
    else:
        chain = get_last_chain_id(input_file)
        #remove_terminal_residue(input_file, f'{input_file_basename}_RESTORE.pdb', chain, p_choice)
        remove_terminal_residue(input_file, f'Mutant{sequence}.pdb', chain, p_choice)

# peptide_mode('HLA.pdb',num_ala,choice) choice:ADD_FIRST/ADD_END/REMOVE_LAST/REMOVE_FIRST

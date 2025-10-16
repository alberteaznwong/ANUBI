import argparse
import os
import random
from os import access, R_OK
from os.path import isfile

from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched


#from compute_weights import compute_weights
#from FUNCTION.compute_weights import compute_weights

#
# from: https://salilab.org/modeller/wiki/Mutate_model
#  mutate_model.py
#
#     Usage:   python mutate_model.py modelname respos resname chain > logfile
#
#     Example: python mutate_model.py 1t29 1699 LEU A > 1t29.log
#
#
#  Creates a single in silico point mutation to sidechain type and at residue position
#  input by the user, in the structure whose file is modelname.pdb
#  The conformation of the mutant sidechain is optimized by conjugate gradient and
#  refined using some MD.
#
#  Note: if the model has no chain identifier, specify "" for the chain argument.
#

# six residue types for he 20 amino acids - dictionary of tuples (that are immutable) -
RESIDUE_TYPES = {'HydPhobic': ('LEU', 'VAL', 'ILE', 'MET', 'PHE', 'TYR', 'TRP'),
                 'Neg': ('GLU', 'ASP'),
                 'Pos': ('ARG', 'LYS'),
                 'HydPhilic': ('SER', 'THR', 'ASN', 'GLN', 'HIS'),
                 'Other': ('ALA', 'CYS', 'PRO'),
                 'Gly': ('GLY',)
                 }

# Real random number generator (slower)
RAND = random.SystemRandom()


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


# use homologs and dihedral library for dihedral angle restraints
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


def mutate_residue(modelname, respos, chain, new_restype):
    # first argument
    # modelname, respos, new_restype, chain, = sys.argv[1:]

    # log.verbose()
    log.none()
    # Set a different value for rand_seed to get a different final model
    env = Environ(rand_seed=random.randint(-500000, 500000))

    env.io.hetatm = True
    # soft sphere potential
    env.edat.dynamic_sphere = False
    # lennard-jones potential (more accurate)
    env.edat.dynamic_lennard = True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = Model(env, file=modelname)
    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    # Saving the starting residue type
    start_res_type = mdl1.residues[str(respos) + ':' + chain].pdb_name

    # set up the mutate-residue selection segment
    select = Selection(mdl1.chains[chain].residues[respos])

    # perform the mutate-residue operation
    select.mutate(residue_type=new_restype)
    # get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    # here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    # yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = Model(env, file=modelname)

    # required to do a transfer_res_numb
    # ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    # transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2, ali)

    # It is usually necessary to write the mutated sequence out and read it in
    # before proceeding, because not all sequence related information about MODEL
    # is changed by this command (e.g., internal coordinates, charges, and atom
    # types and radii are not updated).
    mdl1.write(file=modelname + new_restype + str(respos) + '.tmp')
    mdl1.read(file=modelname + new_restype + str(respos) + '.tmp')

    # set up restraints before computing energy
    # we do this a second time because the model has been written out and read in,
    # clearing the previously set restraints
    make_restraints(mdl1, ali)

    # a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms = 1

    sched = autosched.loop.make_for_model(mdl1)

    # only optimize the selected residue (in first pass, just atoms in selected
    # residue, in second pass, include non-bonded neighboring atoms)
    # set up the mutate-residue selection segment
    select = Selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(select)

    select.energy()

    select.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms = 2
    optimize(select, sched)

    # feels environment (energy computed on pairs that have at least one member
    # in the selected)
    mdl1.env.edat.nonbonded_sel_atoms = 1
    optimize(select, sched)

    select.energy()
    # give a proper name
    mdl1.write(file=modelname[0:-4] + '_' + chain + '_' + start_res_type + str(respos) + new_restype + '.pdb')
    # delete the temporary file
    os.remove(modelname + new_restype + str(respos) + '.tmp')
    # return the name of the new model
    return modelname[0:-4] + '_' + chain + '_' + start_res_type + str(respos) + new_restype + '.pdb'


def make_new_mutation(pdb_file,res_pos_list,
                      new_restype_list, output_name, 
                      system=None, verbose=False):
    # Create a Modeller environment
    log.none()
    env = Environ(rand_seed=1037)
    # Read the PDB file into a Modeller model
    mdl = Model(env, file=pdb_file)
    new_model_name = ""
    res_weight_list = None
    sequence_mutant_log = "OUTPUT_MUTANT.log"

    if res_pos_list:
        # If I get the res-pos-list I need to select one random residue from it.
        # Random mutation in a random residue from a list of res-type
        #                     conserve the hydration properties of the protein (Carol method)
        # Randomly select one residue among the residue:chain of the list
        residue_pos_list = list(res_pos_list)

        if verbose:
            print("\nresidue_pos_list: " + str(residue_pos_list))
            #print("res_weight_files: " + str(res_weight_files))
            #print("res_weight_list: " + str(res_weight_list))
        # print("residue_pos_list"+str(residue_pos_list)+" res_weight_list"+str(args.res_weight_list))
        # if not assigned, res_weight_list=None that is the Default value
        random_res_position = RAND.choices(residue_pos_list, weights=res_weight_list)[0]
        
        res_position = random_res_position.split(':')[0]
        chain = random_res_position.split(':')[1]
        res_type = mdl.residues[str(res_position) + ':' + chain].pdb_name
        if verbose:
            print("random_res_position: " + str(random_res_position))
            print("\nres_position: " + str(res_position))
            print("chain: " + str(chain))
            print("res_type: " + str(res_type))

        if new_restype_list:
            res_type = mdl.residues[str(res_position) + ':' + chain].pdb_name
            new_restype_list = list(new_restype_list)
            if res_type in new_restype_list:
                new_restype_list.remove(res_type)
            new_restype = RAND.choice(new_restype_list)

            new_model_name = mutate_residue(pdb_file, res_position, chain, new_restype)
            if verbose:
                # Check the residues that have been mutated and print them
                mut_mdl = Model(env, file=new_model_name)
                res_type_mut = mut_mdl.residues[str(res_position) + ':' + chain].pdb_name
                print("Selecting a residue among: " + str(new_restype_list))
                print("Mutate_rand_residue-> RES=" + res_type + str(res_position) + ":" + chain + " to " + new_restype)

                log_entry = f"Mutate_rand_residue-> RES={res_type}{res_position}:{chain} to {new_restype}"   
                with open(sequence_mutant_log, 'a') as file:  
                    file.write(log_entry)

                print("Check on the mutant->     " + res_type_mut + str(res_position) + ":" + chain)
            else:
                print(chain + "_" + res_type + str(res_position) + new_restype)
        else:
            print("ERROR: you must provide a non-empty new_restype_list")
            exit(1)
    else:
        print("ERROR: you must provide a list of residue from which choose (list with res_pos:chain)")
        exit(1)

    if output_name:
        # change the output name to the one set by the user
        os.rename(new_model_name, output_name + '.pdb')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates a single in silico point mutation to sidechain type and at '
                                                 'residue position input by the user, in the structure whose file is '
                                                 'modelname.pdb The conformation of the mutant sidechain is optimized '
                                                 'by conjugate gradient and refined using some MD.')
    parser.add_argument('pdb_file', type=str, help='PDB file with the starting model')
    parser.add_argument('-rl', '--res_pos_list', type=str, nargs='+',
                        help='list of Residue "position:chain" of the residue you want to mutate')
#    parser.add_argument('-rw', '--res_weight_files', type=str, nargs='+', default=[],
#                       help='Input files name with the per-residue energy decomposition. I will extract the average.'
#                           '(Default: None)')
    parser.add_argument('-nl', '--new_restype_list', type=str, nargs='+',
                        default=['ALA', 'ARG', 'ASN', 'ASP', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MET'],
                        help='list of 3letters name of the new possible residues')
    parser.add_argument('-o', '--output_name', type=str, help='output file name (with no extention)', default="")
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('-s', '--system', default=None,
                        help='System that we are using.. used for adhoc_fixing on compute_weights.py')
    args = parser.parse_args()
    # verbose = args.verbose

    make_new_mutation(pdb_file=args.pdb_file,
                      res_pos_list=args.res_pos_list,
                      new_restype_list=args.new_restype_list,
                      #res_weight_files=args.res_weight_files,
                      output_name=args.output_name,
                      system=args.system,
                      verbose=args.verbose)

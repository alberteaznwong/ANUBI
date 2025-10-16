
from Bio.PDB import PDBParser, PDBIO, Model, Chain, Residue, Atom
from Bio.PDB.Structure import Structure


def get_chain_start_residues(structure):
    chain_start_ids = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    chain_start_ids[chain.id] = residue.id[1]
                    break
    return chain_start_ids

def copy_and_renumber_structure(structure, chain_start_ids):
    new_structure = Structure("renumbered")
    model = Model.Model(0)
    new_structure.add(model)

    for chain in structure[0]:
        chain_id = chain.id
        if chain_id not in chain_start_ids:
            raise ValueError(f"Chain {chain_id} not in original PDB")

        new_chain = Chain.Chain(chain_id)
        start_res_id = chain_start_ids[chain_id]
        i = 0

        for residue in chain:
            if residue.id[0] != ' ':
                continue

            new_res_id = (' ', start_res_id + i, ' ')
            new_residue = Residue.Residue(new_res_id, residue.resname, residue.segid)

            for atom in residue:
                new_atom = Atom.Atom(
                    atom.name, atom.coord, atom.bfactor, atom.occupancy,
                    atom.altloc, atom.fullname, atom.serial_number, atom.element
                )
                new_residue.add(new_atom)

            new_chain.add(new_residue)
            i += 1

        model.add(new_chain)
    return new_structure

def process_pdb_safe(original_pdb, new_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    original = parser.get_structure("original", original_pdb)
    new = parser.get_structure("new", new_pdb)

    chain_starts = get_chain_start_residues(original)
    safe_structure = copy_and_renumber_structure(new, chain_starts)

    io = PDBIO()
    io.set_structure(safe_structure)
    io.save(output_pdb)
    print(f"Renumbering completed: {output_pdb}")



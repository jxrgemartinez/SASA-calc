from Bio.PDB import PDBParser

class Atom:
    def __init__(self, name, atom_type, coord):
        self.name = name
        self.type = atom_type
        self.x = coord[0]
        self.y = coord[1]
        self.z = coord[2]
        self.rad_vdw = ATOM_VdW.get(atom_type)

    def __str__(self):
        return f"atom {self.name}, type {self.type}, Van der Waals radius = {self.rad_vdw:.1f} A, coordinates ({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"

    def calc_distance(self, other):
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2) ** 0.5
    

class Protein:
    def __init__(self, name):
        self.name = name
        self.list_atoms = []
        
    def add_atom(self, atom):
        if isinstance(atom, Atom):
            self.list_atoms.append(atom)
    
    def build_mlc_from_pdb(self, filename):
        parser = PDBParser()
        structure = parser.get_structure(self.name, filename)
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        name = atom.get_name()
                        atom_type = atom.element
                        coord = atom.get_coord()
                        new_atom = Atom(name, atom_type, coord)
                        self.add_atom(new_atom)
from Bio.PDB import PDBParser
import numpy as np

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
                        
class Sphere:
    def __init__(self, n_points=100):
        self.n_points = n_points
        self.unit_sphere_points = self._generate_unit_sphere_points()
    
    def _generate_unit_sphere_points(self): # Sphere of radius 1 to reescale after
        points = []
        inc = np.pi * (3 - np.sqrt(5))  # Calculate the golden angle for point distribution
        offset = 2 / self.n_points  # Determine vertical spacing between points

        for i in range(self.n_points):
            y = i * offset - 1 + (offset / 2)  
            r = np.sqrt(1 - y * y)  
            phi = i * inc  
            x = np.cos(phi) * r  
            z = np.sin(phi) * r  
    
            points.append([x, y, z])  
        
        return np.array(points)

    
    def get_adjusted_points(self, center_atom): # Rescale and relocate sphere
        adjusted_points = []
        radius = center_atom.rad_vdw
        cx, cy, cz = center_atom.x, center_atom.y, center_atom.z
        for x, y, z in self.unit_sphere_points:
            adj_x = x * radius + cx
            adj_y = y * radius + cy
            adj_z = z * radius + cz
            adjusted_points.append([adj_x, adj_y, adj_z])
        return np.array(adjusted_points)
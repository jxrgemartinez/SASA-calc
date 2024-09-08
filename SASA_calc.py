from Bio.PDB import PDBParser
import numpy as np

class Atom:
    def __init__(self, serial, name, atom_type, aa_name, aa_number, chain_id, coord):
        self.serial = serial  # Atom number 
        self.name = name  # Atom name
        self.type = atom_type  # Element type
        self.aa_name = aa_name  # Amino acid name
        self.aa_number = aa_number  # Amino acid sequence number
        self.chain_id = chain_id  # Chain identifier
        self.coord = np.array(coord)  
        self.rad_vdw = ATOM_VdW.get(atom_type, 1.0)
        self.neighbors = []

    def __str__(self):
        return (f"Atom {self.serial}, {self.name}, Type {self.type}, Amino Acid: {self.aa_name} {self.aa_number}, "
                f"Chain: {self.chain_id}, VdW Radius: {self.rad_vdw:.1f}, "
                f"Coordinates: ({self.coord[0]:.3f}, {self.coord[1]:.3f}, {self.coord[2]:.3f})")

    def calc_distance(self, other):
        return np.linalg.norm(self.coord - other.coord)

    def add_neighbor(self, other):
        self.neighbors.append(other)
    
    
class Protein:
    def __init__(self, filepath):
        self.list_atoms = []
        self._build_mlc_from_pdb(filepath)
        self.distance_matrix = self.calculate_distances() 
        
    def add_atom(self, atom):
        if isinstance(atom, Atom):
            self.list_atoms.append(atom)
    
    def _build_mlc_from_pdb(self, filepath):
        parser = PDBParser(QUIET=True)
        structure_id = filepath.split('/')[-1].split('.')[0]
        structure = parser.get_structure(structure_id, filepath)
        for model in structure:
            for chain in model:
                for residue in chain:
                    aa_name = residue.resname  # Get the amino acid name
                    aa_number = residue.get_id()[1]  # Get the sequence number of the amino acid
                    chain_id = chain.id  # Chain identifier
                    for atom in residue.get_atoms():
                        serial = atom.serial_number
                        name = atom.get_name()
                        atom_type = atom.element
                        coord = atom.get_coord()
                        new_atom = Atom(serial, name, atom_type, aa_name, aa_number, chain_id, coord)
                        self.add_atom(new_atom)

    def calculate_distances(self):
        n = len(self.list_atoms)
        distance_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                dist = self.list_atoms[i].calc_distance(self.list_atoms[j])
                distance_matrix[i][j] = dist
                distance_matrix[j][i] = dist
        return distance_matrix
    
    def determine_neighbors(self, atom_index, cutoff=15):
        for j in range(len(self.list_atoms)):
            if j != atom_index and self.distance_matrix[atom_index][j] < cutoff:
                self.list_atoms[atom_index].add_neighbor(self.list_atoms[j])
                        
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

    def rescale_sphere(self, center_atom): # Rescale and relocate sphere
        adjusted_points = []
        
        radius = center_atom.rad_vdw
        cx, cy, cz = center_atom.coord
        
        adjusted_points = (self.unit_sphere_points * radius) + np.array([cx, cy, cz])
        
        return adjusted_points
 
def calculate_distances_from_sphere_to_neighbors(atom, sphere, protein, sonde):
    rescaled_sphere_points = sphere.rescale_sphere(atom)
    accessible_points = []
    for point in rescaled_sphere_points:
        is_accessible = True
        for neighbor in atom.neighbors:
            distance_to_neighbor = np.linalg.norm(point - neighbor.coord)
            if distance_to_neighbor < (neighbor.rad_vdw + sonde): 
                is_accessible = False
                break 

        if is_accessible:
            accessible_points.append(point)

    return accessible_points, {}


def SASA(filepath):
    protein = Protein(filepath)
    sphere = Sphere(100)
    
    for i, atom in enumerate(protein.list_atoms):
        protein.determine_neighbors(i, cutoff=15)
        
        accessible_points, neighbor_distances = calculate_distances_from_sphere_to_neighbors(atom, sphere, protein, sonde=1.4)
        
        print(f"Atom {atom.name}: {len(accessible_points)} accessible points")
        

ATOM_VdW = {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
    
if __name__ == "__main__":
    SASA("/Users/jorge/Library/Mobile Documents/com~apple~CloudDocs/Downloads/1a5d.pdb")
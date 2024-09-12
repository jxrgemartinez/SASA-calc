import os
from multiprocessing import Pool
import numpy as np
from Bio.PDB import PDBParser

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

MAX_ASA ={
    'ALA': 121.0,
    'ARG': 265.0,
    'ASN': 187.0,
    'ASP': 187.0,
    'CYS': 148.0,
    'GLU': 214.0,
    'GLN': 214.0,
    'GLY': 97.0,
    'HIS': 216.0,
    'ILE': 195.0,
    'LEU': 191.0,
    'LYS': 230.0,
    'MET': 203.0,
    'PHE': 228.0,
    'PRO': 154.0,
    'SER': 143.0,
    'THR': 163.0,
    'TRP': 264.0,
    'TYR': 255.0,
    'VAL': 165.0
}

class Atom:
    """
    Represents an atom within a molecular structure, encapsulating properties
    relevant to the atom itself and its context within a molecule.
    """
    
    def __init__(self, serial, name, atom_type, aa_name, aa_number, chain_id, coord):
        """
        Initialize a new Atom instance with detailed properties and default values for derived attributes.

        Parameters:
        - serial (int): Atom number.
        - name (str): Atom name.
        - atom_type (str): Element type.
        - aa_name (str): Amino acid name.
        - aa_number (int): Amino acid sequence number.
        - chain_id (str): Chain identifier.
        - coord (list or tuple): 3D coordinates of the atom as a list [x, y, z].
        """
        self.serial = serial 
        self.name = name
        self.type = atom_type
        self.aa_name = aa_name
        self.aa_number = aa_number
        self.chain_id = chain_id
        self.coord = np.array(coord)
        self.rad_vdw = ATOM_VdW.get(atom_type, 1.0)
        self.neighbors = []
        self.abs_sasa = 0
        self.rel_sasa = 0
    
    def add_neighbor(self, other):
        """
        Adds another atom to the list of this atom's neighbors.

        Parameters:
        - other (Atom): The atom to be added as a neighbor.
        """
        self.neighbors.append(other)
 
 
def distance_pair(coords):
    """
    Calculates the Euclidean distance between two 3D points.

    Parameters:
    - coords (tuple of np.ndarray): A tuple containing two numpy arrays representing the coordinates
      of the two points (coord1, coord2).

    Returns:
    - float: The Euclidean distance between the two points.
    """
    coord1, coord2 = coords
    return np.linalg.norm(coord1 - coord2)

   
class Protein:
    """
    Represents a protein structure by parsing a PDB file and computing distances
    between all protein atoms.
    """
    
    def __init__(self, filepath, model_index=0):
        """
        Initializes the Protein object by building the molecular structure from a PDB file
        and computing the distance matrix.

        Parameters:
        - filepath (str): Path to the PDB file.
        - model_index (int): Index of the model to parse from the PDB file.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File does not exist: {filepath}")
        
        self.list_atoms = []
        self.atom_index_map = {}
        self._build_mlc_from_pdb(filepath, model_index)
        
        self.distance_matrix = self.compute_distance_matrix()

    def add_atom(self, atom):
        """
        Adds an atom to the protein structure.

        Parameters:
        - atom (Atom): The atom to be added.
        """
        if isinstance(atom, Atom):
            self.list_atoms.append(atom)
            self.atom_index_map[atom.serial] = len(self.list_atoms) - 1
    
    def _build_mlc_from_pdb(self, filepath, model_index):
        """
        Private method to parse the PDB file and store the atom coordinates.

        Parameters:
        - filepath (str): Path to the PDB file.
        - model_index (int): Index of the model to parse.
        """
        parser = PDBParser(QUIET=True)
        structure_id = filepath.split('/')[-1].split('.')[0]
        
        structure = parser.get_structure(structure_id, filepath)
        
        # Check if model index is valid
        if model_index >= len(structure):
            raise IndexError(f"Model index {model_index} is out of bounds. Structure only has {len(structure)} models.")
        model = structure[model_index]
        
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'HOH' or residue.get_resname() == 'HETATM' or residue.get_resname() == 'WAT':
                    continue  # Skip water molecules and heteroatoms
                
                aa_name = residue.resname
                aa_number = residue.get_id()[1]
                
                chain_id = chain.id
                
                for atom in residue.get_atoms():
                    serial = atom.serial_number
                    name = atom.get_name()
                    atom_type = atom.element
                    coord = atom.get_coord()
                    
                    new_atom = Atom(serial, name, atom_type, aa_name, aa_number, chain_id, coord)
                    self.add_atom(new_atom)
   
    def compute_distance_matrix(self):
        """
        Computes the distance matrix for all atoms in the protein using parallel processing.

        Returns:
        - np.ndarray: A symmetric matrix of distances between all atoms.
        """
        num_atoms = len(self.list_atoms)
        tasks = [(self.list_atoms[i].coord, self.list_atoms[j].coord) for i in range(num_atoms) for j in range(i, num_atoms)]

        with Pool() as pool:
            results = pool.map(distance_pair, tasks)

        distance_matrix = np.zeros((num_atoms, num_atoms))
        result_index = 0
        for i in range(num_atoms):
            for j in range(i, num_atoms):
                distance = results[result_index]
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance
                result_index += 1

        return distance_matrix

    def determine_neighbors(self, atom_index, cutoff=15):
        """
        Determines close neighbors for a specified atom within a given cutoff distance.

        Parameters:
        - atom_index (int): Index of the atom to find neighbors for.
        - cutoff (float): Distance cutoff for determining neighbors.
        """
        for j in range(len(self.list_atoms)):
            if j != atom_index and self.distance_matrix[atom_index][j] < cutoff:
                self.list_atoms[atom_index].add_neighbor(self.list_atoms[j])


class Sphere:
    """
    Represents a geometric sphere used for calculating solvent accessible surface area (SASA).
    This class generates a unit sphere with the Fibonacci lattice method and can relocate and
    rescale it based on the coordinates and van der Waals radius of a center atom.
    """
    
    def __init__(self, n_points=100):
        """
        Initialize the Sphere with a specified number of points.

        Parameters:
        - n_points (int): The number of points to generate on the unit sphere.
        """
        self.n_points = n_points
        self.unit_sphere_points = self._generate_unit_sphere_points()
    
    def _generate_unit_sphere_points(self):
        """
        Generates points on a unit sphere using the Fibonacci lattice method.

        Returns:
        - np.ndarray: Coordinates of points on the unit sphere.
        """
        points = []
        for k in range(self.n_points):
            theta = np.arccos(1 - 2 * (k + 0.5) / self.n_points)  # polar angle
            phi = np.pi * (1 + np.sqrt(5)) * k  # azimuthal angle (golden angle approximation)
            
            x = np.sin(theta) * np.cos(phi)  # x coordinate
            y = np.sin(theta) * np.sin(phi)  # y coordinate
            z = np.cos(theta)  # z coordinate
            
            points.append([x, y, z])
        
        return np.array(points)

    def rescale_sphere(self, center_atom, probe_radius):
        """
        Rescales and relocates the sphere based on the center atom's coordinates and radius.

        Parameters:
        - center_atom (Atom): The atom at the center of the sphere.
        - probe_radius (float): The radius of the probe used in SASA calculation.

        Returns:
        - np.ndarray: The coordinates of the rescaled sphere points.
        """
        radius = center_atom.rad_vdw + probe_radius
        adjusted_points = (self.unit_sphere_points * radius) + center_atom.coord
        
        return adjusted_points
 
    
class SASACalculator:
    """
    Calculates the Solvent Accessible Surface Area (SASA) of a protein using a sphere-based
    approximation method. This class handles the computation across all atoms in the protein
    and provides detailed SASA reports.
    """
    
    def __init__(self, protein, n_points=100, probe_radius=1.4):
        """
        Initializes a SASA calculator for a given protein.

        Parameters:
        - protein (Protein): The protein for which SASA will be calculated.
        - n_points (int): Number of points on the sphere used in the calculation.
        - probe_radius (float): Radius of the probe sphere used in SASA calculation.
        """
        self.protein = protein
        self.sphere = Sphere(n_points)
        self.probe_radius = probe_radius
        self.protein.total_absolute_sasa = 0
        self.protein.total_rel_sasa = 0
        
    def compute_atom_sasa(self, atom):
        """
        Computes the SASA for a single atom.

        Parameters:
        - atom (Atom): The atom for which SASA is being calculated.

        Returns:
        - float: The solvent accessible surface area of the atom.
        """
        atom_index = self.protein.atom_index_map[atom.serial]
        self.protein.determine_neighbors(atom_index)

        rescaled_sphere_points = self.sphere.rescale_sphere(atom, self.probe_radius)
        accessible_points = 0

        for point in rescaled_sphere_points:
            if all(np.linalg.norm(point - neighbor.coord) >= (neighbor.rad_vdw + self.probe_radius) for neighbor in atom.neighbors):
                accessible_points += 1

        sphere_area = 4 * np.pi * (atom.rad_vdw + self.probe_radius) ** 2
        accessible_area = (accessible_points / len(rescaled_sphere_points)) * sphere_area
        return accessible_area

    def run(self):
        """
        Executes the SASA calculation for all atoms in the protein and updates their SASA values.
        """
        with Pool() as pool:
            results = pool.map(self.compute_atom_sasa, self.protein.list_atoms)

        total_absolute_sasa = sum(results)
        for atom, abs_sasa in zip(self.protein.list_atoms, results):
            atom.abs_sasa = abs_sasa
        
        self.protein.total_absolute_sasa = total_absolute_sasa
    
    def print_sasa_report(self, output_type="complete"):
        """
        Prints a report of the SASA calculations.

        Parameters:
        - output_type (str): Specifies the type of SASA report ('total', 'atomic', 'residue', 'complete').
        """
        possible_arguments = ["total", "atomic", "residue", "complete"]
        
        if output_type not in possible_arguments:
            raise ValueError(f"Invalid output type. Possible types are {possible_arguments}")
        
        print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8} {:>8}".format("TYPE", "ATM_N", "ATM", "RES", "CHN", "RES_N", "ABS_SASA", "REL_SASA"))
        
        chain_sasa = {}
        amino_acid_sasa = {}
        self.protein.total_rel_sasa = 0
        
        for atom in self.protein.list_atoms:
            chain_id = atom.chain_id
            aa_key = (atom.aa_name, atom.aa_number, chain_id)
            
            if aa_key in amino_acid_sasa:
                amino_acid_sasa[aa_key]["abs"] += atom.abs_sasa
            else:
                amino_acid_sasa[aa_key] = {"abs": atom.abs_sasa}
        
        for aa_key, sasa_value in amino_acid_sasa.items():
            max_sasa = MAX_ASA.get(aa_key[0], 0)
            rel_sasa = (sasa_value["abs"] / max_sasa) * 100 if max_sasa > 0 else 0
            sasa_value["rel"] = rel_sasa 
            
            chain_id = aa_key[2]
            if chain_id in chain_sasa:
                chain_sasa[chain_id]["abs"] += sasa_value["abs"]
                chain_sasa[chain_id]["rel"] += rel_sasa
            else:
                chain_sasa[chain_id] = {"abs": sasa_value["abs"], "rel": rel_sasa}
            
            self.protein.total_rel_sasa += rel_sasa
            
        if output_type == "atomic":
            for atom in self.protein.list_atoms:
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8}".format("ATOM", atom.serial, atom.name, atom.aa_name, atom.chain_id, atom.aa_number, atom.abs_sasa, ""))

        if output_type in ["residue", "complete"]:
            for aa_key, sasa_value in amino_acid_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("RES", "", "", aa_key[0], aa_key[2], aa_key[1], sasa_value["abs"], sasa_value["rel"]))
        
        if output_type in ["residue", "complete", "total"]:
            for chain, sasa in chain_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("CHAIN", "", "", "", chain, "", sasa["abs"], sasa["rel"]))
            print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("TOTAL", "", "", "", "", "", self.protein.total_absolute_sasa, self.protein.total_rel_sasa))
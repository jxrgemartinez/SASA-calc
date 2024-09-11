import os
from Bio.PDB import PDBParser
import numpy as np
from multiprocessing import Pool

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

RESIDUE_SASA = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0
}


class Atom:
    """
    Represents an atom within a molecular structure, containing
    properties relevant to both the atom itself and its context
    within a larger molecule
    """
    def __init__(self, serial, name, atom_type, aa_name, aa_number, chain_id, coord):
        """
        Initializes a new instance of the Atom class.

        Args:
        serial (int): Atom number.
        name (str): Atom name.
        atom_type (str): Element type.
        aa_name (str): Amino acid name.
        aa_number (int): Amino acid sequence number.
        chain_id (str): Chain identifier.
        coord (list): 3D coordinates of the atom as a list [x, y, z].
        """
        self.serial = serial  # Atom number 
        self.name = name  # Atom name
        self.type = atom_type  # Element type
        self.aa_name = aa_name  # Amino acid name
        self.aa_number = aa_number  # Amino acid sequence number
        self.chain_id = chain_id  # Chain
        self.coord = np.array(coord)  # Coordinates
        self.rad_vdw = ATOM_VdW.get(atom_type, 1.0) # Van der Waals radius
        self.neighbors = []
        self.abs_sasa = 0 # Absolute solvent accessible surface area
        self.rel_sasa = 0 # Relative solvent accessible surface area
    
    def add_neighbor(self, other):
        """
        Adds another atom to the list of this atom's neighbors.

        Args:
        other (Atom): Another atom to be added as a neighbor.
        """
        self.neighbors.append(other)
 
 
def distance_pair(coords):
    """
    Calculates the Euclidean distance between two points in 3D space.

    Args:
    coords (tuple): A tuple containing two numpy arrays representing the coordinates of the two points.

    Returns:
    float: The Euclidean distance between the two points.
    """
    coord1, coord2 = coords
    return np.linalg.norm(coord1 - coord2)

   
class Protein:
    """
    Represents a protein structure, capable of parsing PDB files,
    storing atoms, and computing distances between them.
    """
    def __init__(self, filepath, model_index=0):
        """
        Initializes a Protein object by parsing a PDB file.

        Args:
        filepath (str): Path to the PDB file.
        model_index (int): Index of the model in the PDB file to be parsed.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File does not exist: {filepath}")
        
        self.list_atoms = []
        self.atom_index_map = {}
        self._build_mlc_from_pdb(filepath, model_index)
        
        self.distance_matrix = self.compute_distance_matrix()

    def add_atom(self, atom):
        """
        Adds an Atom object to the protein structure.

        Args:
        atom (Atom): The Atom object to be added.
        """
        if isinstance(atom, Atom):
            self.list_atoms.append(atom)
            self.atom_index_map[atom.serial] = len(self.list_atoms) - 1
    
    def _build_mlc_from_pdb(self, filepath, model_index):
        """
        Private method to parse a PDB file and build the molecular structure.

        Args:
        filepath (str): Path to the PDB file.
        model_index (int): Index of the model to parse.
        """
        parser = PDBParser(QUIET=True)
        structure_id = filepath.split('/')[-1].split('.')[0]
        
        # Load structure
        structure = parser.get_structure(structure_id, filepath)
        
        # Check if model index is valid
        if model_index >= len(structure):
            raise IndexError(f"Model index {model_index} is out of bounds. Structure only has {len(structure)} models.")
        model = structure[model_index]
        
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'HOH':
                    continue  # Skip water molecules
                
                aa_name = residue.resname  # Get the amino acid name
                aa_number = residue.get_id()[1]  # Get the sequence number of the amino acid
                
                chain_id = chain.id  # Get chain
                
                for atom in residue.get_atoms():
                    serial = atom.serial_number
                    name = atom.get_name()
                    atom_type = atom.element
                    coord = atom.get_coord()
                    
                    new_atom = Atom(serial, name, atom_type, aa_name, aa_number, chain_id, coord)
                    self.add_atom(new_atom)
   
    def compute_distance_matrix(self):
        """
        Computes the distance matrix for all atoms in the protein.

        Returns:
        np.ndarray: A symmetric matrix of distances between all atom pairs.
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
        Determines neighboring atoms for a given atom within a specified cutoff distance.

        Args:
        atom_index (int): Index of the atom in list_atoms.
        cutoff (float): Distance cutoff to consider an atom as a neighbor.
        """
        for j in range(len(self.list_atoms)):
            if j != atom_index and self.distance_matrix[atom_index][j] < cutoff:
                self.list_atoms[atom_index].add_neighbor(self.list_atoms[j])


class Sphere:
    """
    Represents a geometric sphere used for calculating solvent accessible surface area (SASA).
    This class generates points uniformly distributed over a unit sphere and can rescale them
    based on the van der Waals radius of atoms.
    """
    def __init__(self, n_points=100):
        """
        Initializes a Sphere object with a specified number of points.

        Args:
        n_points (int): Number of points to generate on the sphere.
        """
        self.n_points = n_points
        self.unit_sphere_points = self._generate_unit_sphere_points()
    
    def _generate_unit_sphere_points(self): # Sphere of radius 1 to reescale after
        """
        Generates points uniformly distributed over the surface of a unit sphere
        using the golden angle.

        Returns:
        np.ndarray: Array of points on the unit sphere.
        """
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

    def rescale_sphere(self, center_atom, sonde_radius):
        """
        Rescales and relocates the sphere points according to the atom's radius and the sonde radius.

        Args:
        center_atom (Atom): The atom around which the sphere is centered.
        sonde_radius (float): Radius of the probe used in SASA calculation.

        Returns:
        np.ndarray: Array of rescaled and relocated sphere points.
        """
        # Rescale and relocate the sphere
        radius = center_atom.rad_vdw + sonde_radius
        adjusted_points = (self.unit_sphere_points * radius) + center_atom.coord
        
        return adjusted_points
 
    
class SASACalculator:
    """
    Calculates the Solvent Accessible Surface Area (SASA) of a protein using a sphere-based
    approximation method. This class handles the computation across all atoms in the protein
    and provides detailed SASA reports.
    """
    def __init__(self, protein, n_points=100, sonde_radius=1.4):
        """
        Initializes a SASA calculator for a given protein.

        Args:
        protein (Protein): The protein for which SASA will be calculated.
        n_points (int): Number of points on the sphere used in the calculation.
        sonde_radius (float): Radius of the probe sphere used in SASA calculation.
        """
        self.protein = protein
        self.sphere = Sphere(n_points)
        self.sonde_radius = sonde_radius
        self.protein.total_absolute_sasa = 0
        
    def compute_atom_sasa(self, atom):
        """
        Computes the SASA for a single atom.

        Args:
        atom (Atom): The atom for which SASA is being calculated.

        Returns:
        float: The solvent accessible surface area of the atom.
        """
        atom_index = self.protein.atom_index_map[atom.serial]
        self.protein.determine_neighbors(atom_index)

        rescaled_sphere_points = self.sphere.rescale_sphere(atom, self.sonde_radius)
        accessible_points = 0

        for point in rescaled_sphere_points:
            if all(np.linalg.norm(point - neighbor.coord) >= (neighbor.rad_vdw + self.sonde_radius) for neighbor in atom.neighbors):
                accessible_points += 1

        sphere_area = 4 * np.pi * (atom.rad_vdw + self.sonde_radius) ** 2
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
            atom.rel_sasa = (abs_sasa / total_absolute_sasa * 100) if total_absolute_sasa > 0 else 0
        
        self.protein.total_absolute_sasa = total_absolute_sasa
    
    def print_sasa_report(self, output_type="copmlete"):
        """
        Prints a report of the SASA calculations.

        Args:
        output_type (str): Specifies the type of SASA report ('total', 'atomic', 'residue', 'complete').
        """
        possible_arguments = ["total", "atomic", "residue", "complete"]
        if output_type not in possible_arguments:
            raise ValueError(f"Invalid output type. Possible types are {possible_arguments}")
        
        if output_type == "atomic":
            for atom in self.protein.list_atoms:
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("ATOM", atom.serial, atom.name, atom.aa_name, atom.chain_id, atom.aa_number, atom.abs_sasa, atom.rel_sasa))
        
        if output_type == "total":
            print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8}".format("TOTAL", "", "", "", "", "", self.protein.total_absolute_sasa, ""))        
        if output_type == "residue":
            chain_sasa = {}
            amino_acid_sasa = {}
            
            for atom in self.protein.list_atoms:
                chain_id = atom.chain_id
                aa_key = (atom.aa_name, atom.aa_number, atom.chain_id)
                
                if chain_id in chain_sasa:
                    chain_sasa[chain_id]["abs"] += atom.abs_sasa
                    chain_sasa[chain_id]["rel"] += atom.rel_sasa
                else:
                    chain_sasa[chain_id] = {"abs": atom.abs_sasa, "rel": atom.rel_sasa}
            
                if aa_key in amino_acid_sasa:
                    amino_acid_sasa[aa_key]["abs"] += atom.abs_sasa
                    amino_acid_sasa[aa_key]["rel"] += atom.rel_sasa
                else:
                    amino_acid_sasa[aa_key] = {"abs": atom.abs_sasa, "rel": atom.rel_sasa}
            
            for aa_key, sasa_value in amino_acid_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("RES", "", "", aa_key[0], aa_key[2], aa_key[1], sasa_value["abs"], sasa_value["rel"]))
                
            for chain, sasa_value in chain_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("CHAIN", "", "", "", chain, "", sasa_value["abs"], sasa_value["rel"]))

            print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("TOTAL", "", "", "", "", "", self.protein.total_absolute_sasa, ""))
            
        if output_type == "complete":
            for atom in self.protein.list_atoms:
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("ATOM", atom.serial, atom.name, atom.aa_name, atom.chain_id, atom.aa_number, atom.abs_sasa, atom.rel_sasa))
                
            chain_sasa = {}
            amino_acid_sasa = {}
            
            for atom in self.protein.list_atoms:
                chain_id = atom.chain_id
                aa_key = (atom.aa_name, atom.aa_number, atom.chain_id)
                
                if chain_id in chain_sasa:
                    chain_sasa[chain_id]["abs"] += atom.abs_sasa
                    chain_sasa[chain_id]["rel"] += atom.rel_sasa
                else:
                    chain_sasa[chain_id] = {"abs": atom.abs_sasa, "rel": atom.rel_sasa}
            
                if aa_key in amino_acid_sasa:
                    amino_acid_sasa[aa_key]["abs"] += atom.abs_sasa
                    amino_acid_sasa[aa_key]["rel"] += atom.rel_sasa
                else:
                    amino_acid_sasa[aa_key] = {"abs": atom.abs_sasa, "rel": atom.rel_sasa}
            
            for aa_key, sasa_value in amino_acid_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("RES", "", "", aa_key[0], aa_key[2], aa_key[1], sasa_value["abs"], sasa_value["rel"]))
                
            for chain, sasa_value in chain_sasa.items():
                print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8.2f}".format("CHAIN", "", "", "", chain, "", sasa_value["abs"], sasa_value["rel"]))

            print("{:<5} {:>5} {:>4} {:>4} {:>4} {:>4} {:>8.2f} {:>8}".format("TOTAL", "", "", "", "", "", self.protein.total_absolute_sasa, ""))
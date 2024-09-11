import os
from Bio.PDB import PDBParser
import numpy as np
from multiprocessing import Pool
import time

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

class Atom:
    def __init__(self, serial, name, atom_type, aa_name, aa_number, chain_id, coord):
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
        self.neighbors.append(other)
    
class Protein:
    def __init__(self, filepath, model_index=0):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File does not exist: {filepath}")
        
        self.list_atoms = []
        self.atom_index_map = {}
        self._build_mlc_from_pdb(filepath, model_index)
        
        self.distance_matrix = self.compute_distance_matrix()

    def add_atom(self, atom):
        if isinstance(atom, Atom):
            self.list_atoms.append(atom)
            self.atom_index_map[atom.serial] = len(self.list_atoms) - 1
    
    def _build_mlc_from_pdb(self, filepath, model_index):
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
        # Take the row of each atom in the distance matrix and add coundt as neighbors the ones that are closer than the cutoff
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

    def rescale_sphere(self, center_atom, sonde_radius):
        # Rescale and relocate the sphere
        radius = center_atom.rad_vdw + sonde_radius
        adjusted_points = (self.unit_sphere_points * radius) + center_atom.coord
        
        return adjusted_points
 
    
class SASACalculator:
    def __init__(self, protein, n_points=100, sonde_radius=1.4):
        self.protein = protein
        self.sphere = Sphere(n_points)
        self.sonde_radius = sonde_radius
        self.protein.total_absolute_sasa = 0
        
    def compute_atom_sasa(self, atom):
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
        with Pool() as pool:
            results = pool.map(self.compute_atom_sasa, self.protein.list_atoms)

        total_absolute_sasa = sum(results)
        for atom, abs_sasa in zip(self.protein.list_atoms, results):
            atom.abs_sasa = abs_sasa
            atom.rel_sasa = (abs_sasa / total_absolute_sasa * 100) if total_absolute_sasa > 0 else 0
        
        self.protein.total_absolute_sasa = total_absolute_sasa
    
    def print_sasa_report(self, output_type="total"):
        possible_arguments = ["total", "atomic", "residue"]
        if output_type not in possible_arguments:
            raise ValueError(f"Invalid output type. Possible types are {possible_arguments}")
        
        if output_type == "atomic":
            print("{:>4} {:>4} {:>4} {:>4} {:>5} {:>5} {:>8} {:>8}".format("ATOM", "NAME", "TYPE", "RES", "NUM", "CHAIN", "ABS_SASA", "REL_SASA"))
            
            for atom in self.protein.list_atoms:
                print("{:>4} {:>4} {:>4} {:>4} {:>5} {:>5} {:>8.2f} {:>8.1f}".format(atom.serial, atom.name, atom.type, atom.aa_name, atom.aa_number, atom.chain_id, atom.abs_sasa, atom.rel_sasa))
        
        if output_type == "total":
            print("{:>5} {:>3} {:>4} {:>4} {:>5} {:>5} {:>8.2f}".format("TOTAL", "", "", "", "", "", self.protein.total_absolute_sasa))
        
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
            
            print("{:>4} {:>4} {:>5} {:>8} {:>8}".format("RES", "NUM", "CHAIN", "ABS", "REL"))
            for aa_key, sasa_value in amino_acid_sasa.items():
                print("{:>4} {:>4} {:>5} {:>8.2f} {:>8.1f}".format(aa_key[0], aa_key[1], aa_key[2], sasa_value["abs"], sasa_value["rel"]))
            
            print("\nAbsolute sums over single chains surface areas")
            for chain, sasa_value in chain_sasa.items():
                print("{:>5} {:>4} {:>8.2f}".format("CHAIN", chain, sasa_value["abs"]))

            print("\nAbsolute sums over all chains")
            print("{:>5} {:>4} {:>8.2f}".format("TOTAL", "", self.protein.total_absolute_sasa))
  
    
def distance_pair(coords):
    coord1, coord2 = coords
    return np.linalg.norm(coord1 - coord2)
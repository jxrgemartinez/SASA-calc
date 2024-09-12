# Protein Solvent Accessible Surface Area (SASA) Calculation Tool

## Overview
This Python application, developed as part of academic research at [University/Institute Name], calculates and reports the Solvent Accessible Surface Area (SASA) of proteins from PDB files. Designed to support both educational and research activities in computational biology, this tool computes SASA using precise geometric methods to determine the exposure of atoms within protein structures based on their atomic coordinates.

### Objective
The goal of this project is to offer a robust and efficient tool for the bioinformatics community to facilitate the analysis of protein solvation properties, crucial for understanding molecular interactions and biological functions. This tool stands out by providing detailed features for customizing the analysis, comparing results with established methods like NACCESS, and supporting various output formats for diverse research needs.

## Features
* Computes absolute and relative SASA for each atom in the protein
* Offers different output types: atom by atom, residue by residue and chains, total SASA and a complete report
* Supports selecting specific models within multi-model PDB files
* Adjustable resolution for sphere point density and probe size.

## Setup

To install the algorithm and its dependencies, you need to perform the following steps:

### Clone the repository

```bash
git clone https://github.com/jxrgemartinez/SASA-calc.git

cd SASA-calc
```

### Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Create a Conda environment

```bash
conda env create -f environment.yml
```

### Activate the Conda environment

```bash
conda activate sasa_env
```

## Usage
To run the tool, use the following command from the terminal:

```bash
python src/SASA_calc.py <pdb_file> > <output_file>
```

> The output format depends on the --output flag and is printed directly to the console, consider redirecting the output to a file

### Optional command line Arguments
* `--model`: Model index to use from the PDB file. Defaul is 0, representing the first model
* `--points`: Number of points on the calculation sphere. Default is 100
* `--sonde`: Radius of the probe sphere in Å. Default is 1.4 (typical size of a water molecule)
* `--output`: Output format. Options are 'residue', 'atomic', 'total' and 'complete'. Defaul is 'complete'

### Output
The output format depends on the --output flag and is printed directly to the console

#### Atomic
The atomic output lists SASA details for each atom in the protein

**Exapmle of Atomic Output:**
```bash
ATOM   1    N     MET    A    1     0.35     1.13
ATOM   2    CA    MET    A    1     1.47     1.13
ATOM   3    C     MET    A    1     0.76     1.13
...
```
*Columns represent Atom ID, Atom Name, Residue Name, Chain Identifier, Residue Number, Absolute SASA, Relative SASA*

#### Residue
The residue output provides a summary of SASA for each residue, chain and the total protein

**Example of Residue Output:**
```bash
RES               MET    A    1   195.12     1.13
RES               TYR    A    2   195.12     1.13
RES               TYR    A    3   195.12     1.13
...
CHAIN                    A       8617.87    49.96
CHAIN                    B       8632.99    50.04
TOTAL                           17250.86         
```

#### Total
This output summarizes the total SASA for the entire protein

```bash
TOTAL                           17250.86         
```

#### Complete
The complete output combines all the above details into a fully detailed report
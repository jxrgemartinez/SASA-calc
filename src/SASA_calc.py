import argparse
import classes

def parse_arguments():
    description = ('Calculate and report the Solvent Accessible Surface Area '
                   '(SASA) for a given protein structure file.')
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument(
        'pdb_file', 
        type=str, 
        help='Path to the PDB file containing the protein structure.'
    )
    
    parser.add_argument(
        '--model', 
        type=int, 
        default=0,
        help='Index of the desired model of the protein structure. '
             'Default is 0.'
    )
    
    parser.add_argument(
        '--points', 
        type=int, 
        default=100, 
        help='Number of points on the sphere used in the SASA calculation. '
             'Default is 100.'
    )
    
    parser.add_argument(
        '--sonde', 
        type=float, 
        default=1.4, 
        help='Radius of the probe sphere used in SASA calculation. '
             'Default is 1.4 Angstroms.'
    )
    
    parser.add_argument(
        '--output', 
        type=str, 
        default='complete',
        choices=['residue', 'atomic', 'total', 'complete'],
        help='Type of output for SASA report. Options are "residue", "atomic", "total" or "copmlete. '
             'Default is "complete".'
    )
    
    return parser.parse_args()

def main():
    """
    Main function to parse arguments and perform SASA calculation on a protein structure.
    """
    args = parse_arguments()
    
    try:
        protein = classes.Protein(args.pdb_file, model_index=args.model)
        analysis = classes.SASACalculator(protein, n_points=args.points, sonde_radius=args.sonde)
        analysis.run()
        analysis.print_sasa_report(output_type=args.output)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
if __name__ == '__main__':
    main()
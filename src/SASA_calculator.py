import classes
         
def run_calculator(protein):
    calculator = classes.SASACalculator(protein)
    calculator.run()
    return calculator

if __name__ == '__main__':
    protein = classes.Protein("/Users/jorge/Library/Mobile Documents/com~apple~CloudDocs/Downloads/1a5d.pdb")
    total_sasa = run_calculator(protein)
    
    total_sasa.print_sasa_report(output_type="residue")
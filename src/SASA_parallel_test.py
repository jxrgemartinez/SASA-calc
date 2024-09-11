import classes_parallel as classes

if __name__ == '__main__':
    protein = classes.Protein("data/6fqf.pdb")

    analysis = classes.SASACalculator(protein)
    analysis.run()
    
    analysis.print_sasa_report(output_type="residue")
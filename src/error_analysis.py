import os
import re
import numpy as np
import csv

def extract_number_from_file(file_path, pattern="TOTAL"):
    """
    Extracts the first floating-point number from the specified line of a file.

    Parameters:
    - file_path (str): Path to the file.
    - pattern (str): Pattern to identify the relevant line.

    Returns:
    - float: The first floating-point number on the identified line, or None if not found.
    """
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if pattern in line:
                    numbers = re.findall(r'\d+\.\d+', line)
                    return float(numbers[0]) if numbers else None
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None

def calculate_accuracy(rsa_value, txt_value):
    """
    Calculates the percentage accuracy between two values.

    Parameters:
    - rsa_value (float): Reference value.
    - txt_value (float): Test value.

    Returns:
    - float: Percentage accuracy if both values are not None, otherwise None.
    """
    if rsa_value is None or txt_value is None:
        return None
    return (abs(txt_value - rsa_value) / rsa_value) * 100

def save_results_to_csv(results, filename='res/errors.csv'):
    """
    Saves the results to a CSV file.

    Parameters:
    - results (list of tuples): List containing results data.
    - filename (str): Path to the output CSV file.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['PDB ID', 'NACCESS', 'CUSTOM', 'ERROR (%)'])
        writer.writerows(results)

def print_results_table(results):
    """
    Prints a formatted table of results.

    Parameters:
    - results (list of tuples): List containing results data.
    """
    header = f"{'PDB ID':<10} {'NACCESS':<10} {'CUSTOM':<10} {'ERROR (%)':<10}"
    print(header)
    print("=" * len(header))
    for pdb_id, rsa_value, txt_value, error in results:
        print(f"{pdb_id:<10} {rsa_value:<10.2f} {txt_value:<10.2f} {error:<10.2f}")

def main(directory):
    """
    Main function to analyze results from a directory.

    Parameters:
    - directory (str): Directory containing the files to analyze.
    """
    rsa_files = {f.split('.')[0]: os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.rsa')}
    txt_files = {f.split('.')[0]: os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.txt')}
    
    results = []
    accuracies = []

    for pdb_id in rsa_files:
        if pdb_id in txt_files:
            rsa_value = extract_number_from_file(rsa_files[pdb_id])
            txt_value = extract_number_from_file(txt_files[pdb_id])
            if rsa_value is not None and txt_value is not None:
                accuracy = calculate_accuracy(rsa_value, txt_value)
                accuracies.append(accuracy)
                results.append((pdb_id, rsa_value, txt_value, accuracy))

    if results:
        print_results_table(results)

        average_accuracy = sum(accuracies) / len(accuracies)
        max_error = max(accuracies)
        min_error = min(accuracies)
        std_dev = np.std(accuracies)
        
        print(f"\nAverage error: {average_accuracy:.2f}%")
        print(f"Max error: {max_error:.2f}%, Min error: {min_error:.2f}%, Std Dev: {std_dev:.2f}%")
        
        save_results_to_csv(results)

    else:
        print("No valid data to calculate accuracy.")

if __name__ == "__main__":
    directory = 'res'
    main(directory)

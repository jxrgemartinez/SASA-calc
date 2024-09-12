import os
import re
import numpy as np
import csv

def extract_number_from_rsa(file_path):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if "TOTAL" in line:
                    numbers = re.findall(r'\d+\.\d+', line)
                    if numbers:
                        return float(numbers[0])
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None

def extract_number_from_txt(file_path):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if "TOTAL" in line:
                    numbers = re.findall(r'\d+\.\d+', line)
                    if numbers:
                        return float(numbers[0])
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None

def calculate_accuracy(rsa_value, txt_value):
    if rsa_value is None or txt_value is None:
        return None
    return ((txt_value - rsa_value) / rsa_value) * 100

def save_results_to_csv(results, filename='errors.csv'):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['PDB ID', 'NACCESS', 'CUSTOM', 'ERROR (%)'])
        for pdb_id, rsa_value, txt_value, error in results:
            writer.writerow([pdb_id, rsa_value, txt_value, error])

def print_results_table(results):
    print(f"{'PDB ID':<10} {'NACCESS':<10} {'CUSTOM':<10} {'ERROR (%)':<10}")
    print("=" * 50)
    for pdb_id, rsa_value, txt_value, error in results:
        print(f"{pdb_id:<10} {rsa_value:<10.2f} {txt_value:<10.2f} {error:<10.2f}")

def main(directory):
    rsa_files = {f.split('.')[0]: os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.rsa')}
    txt_files = {f.split('.')[0]: os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.txt')}
    
    results = []
    accuracies = []

    for pdb_id in rsa_files:
        if pdb_id in txt_files:
            rsa_value = extract_number_from_rsa(rsa_files[pdb_id])
            txt_value = extract_number_from_txt(txt_files[pdb_id])
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
        
        # Save results to CSV
        save_results_to_csv(results)

    else:
        print("No valid data to calculate accuracy.")

if __name__ == "__main__":
    directory = 'res'
    main(directory)

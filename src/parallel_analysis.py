import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('res/times.csv')

# Calculate speedup
data['Speedup'] = data['No Parallel'] / data['Parallel']

# Prepare results for printing
results = []
for index, row in data.iterrows():
    pdb_id = row['File'].split('.')[0]  # Extract PDB ID without extension
    no_parallel = row['No Parallel']
    parallel = row['Parallel']
    speedup = row['Speedup']
    results.append((pdb_id, no_parallel, parallel, speedup))

# Print the formatted results
print(f"\n{'PDB ID':<10} {'No Parallel':<12} {'Parallel':<12} {'Speedup':<12}")
print("=" * 50)
for pdb_id, no_parallel, parallel, speedup in results:
    print(f"{pdb_id:<10} {no_parallel:<12.2f} {parallel:<12.2f} {speedup:<12.2f}")

# Calculate overall statistics
speedups = [speedup for _, _, _, speedup in results]
average_speedup = sum(speedups) / len(speedups)
max_speedup = max(speedups)
min_speedup = min(speedups)
std_dev = pd.Series(speedups).std()

# Print overall statistics
print(f"\nAverage speedup: {average_speedup:.2f}")
print(f"Max speedup: {max_speedup:.2f}, Min speedup: {min_speedup:.2f}, Std Dev: {std_dev:.2f}")
# Visualize speedup
plt.figure(figsize=(10, 5))
plt.bar(data['File'], data['Speedup'], color='#F19963')
plt.title('Speedup of Parallel Execution')
plt.ylabel('Speedup')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

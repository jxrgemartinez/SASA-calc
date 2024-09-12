import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('res/times.csv')

# Calculate speedup
data['Speedup'] = data['No Parallel'] / data['Parallel']

# Prepare results for printing
results = data.apply(lambda row: (row['File'].split('.')[0], row['No Parallel'], row['Parallel'], row['Speedup']), axis=1).tolist()

# Print the formatted results
header = f"{'PDB ID':<10} {'No Parallel':<12} {'Parallel':<12} {'Speedup':<12}"
print(f"\n{header}")
print("=" * len(header))
for pdb_id, no_parallel, parallel, speedup in results:
    print(f"{pdb_id:<10} {no_parallel:<12.2f} {parallel:<12.2f} {speedup:<12.2f}")

# Calculate overall statistics using pandas Series for speedups
speedups = data['Speedup']
average_speedup = speedups.mean()
max_speedup = speedups.max()
min_speedup = speedups.min()
std_dev = speedups.std()

# Print overall statistics
print(f"\nAverage speedup: {average_speedup:.2f}")
print(f"Max speedup: {max_speedup:.2f}, Min speedup: {min_speedup:.2f}, Std Dev: {std_dev:.2f}")

# Visualize speedup
plt.figure(figsize=(10, 5))
plt.bar(data['File'], data['Speedup'], color='#F19963')
plt.title('Speedup of Parallel Execution')
plt.ylabel('Speedup')
plt.xticks(rotation=45)
plt.tight_layout()  # Adjust layout to make room for rotated x-axis labels
plt.show()

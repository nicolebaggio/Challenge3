import matplotlib.pyplot as plt
import os

# Folder containing the results
data_folder = 'data'

# Collect results
grid_sizes = []
execution_times = []
l2_norm_errors = []

for filename in os.listdir(data_folder):
    if filename.endswith('.txt'):
        with open(os.path.join(data_folder, filename), 'r') as file:
            lines = file.readlines()
            grid_size = int(lines[0].strip().split(':')[1])
            execution_time = float(lines[1].strip().split(':')[1])
            l2_norm_error = float(lines[2].strip().split(':')[1])
            
            grid_sizes.append(grid_size)
            execution_times.append(execution_time)
            l2_norm_errors.append(l2_norm_error)

# Plot the results
plt.figure(figsize=(12, 6))

# Plot Execution Time
plt.subplot(1, 2, 1)
plt.plot(grid_sizes, execution_times, marker='o')
plt.xlabel('Number of Processors')
plt.ylabel('Execution Time (seconds)')
plt.title('Execution Time vs Number of Processors')
plt.grid(True)

# Plot L2 Norm Error
plt.subplot(1, 2, 2)
plt.plot(grid_sizes, l2_norm_errors, marker='o')
plt.xlabel('Number of Processors')
plt.ylabel('L2 Norm Error')
plt.title('L2 Norm Error vs Number of Processors')
plt.grid(True)

# Save the plots as image files
plt.tight_layout()
plt.savefig('execution_time_vs_processors.png')
plt.close()

# Print a message indicating where the plots are saved
print("Plots saved as 'execution_time_vs_processors.png' and 'l2_norm_error_vs_processors.png'.")

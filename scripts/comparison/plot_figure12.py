import numpy as np
import matplotlib.pyplot as plt


openmp = [1, 2, 4, 8, 16, 32, 64]

openmp = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]

# Original Timings 
original_setup = np.array([11.88, 6.8, 5.1, 4.245, 3.394, 2.977, 2.55, 2.334, 2.118, 2.011, 1.916, 1.867, 1.834])
original_total = np.array([37.43, 21.327, 16.877, 13.286, 9.91, 8.64, 6.69, 5.94, 5.324, 5.78, 6.475, 8.778, 16.7696])
original_solve = original_total - original_setup

refactored_solve = np.array([7.63, 3.89, 2.66, 2.03, 1.415, 1.09, 0.789, 0.692, 0.4773, 0.399, 0.33, 0.2949, 0.277])

old_refactored_solve = np.array([7.63, 4.24, 2.858, 2.233, 1.59, 1.3, 1.08, 1.05, 1.468, 2.10, 3.52, 5.089, 13.352])


# Compute speedup for each number of threads
original_speed_up = original_solve[0] / original_solve
refactored_speed_up = refactored_solve[0] / refactored_solve
old_refactored_speed_up = old_refactored_solve[0] / old_refactored_solve

### Plotting the execution times
plt.figure(figsize=(10, 6))
plt.rcParams.update({'font.size': 14})

plt.plot(openmp, original_solve, label='Original (task-based)', marker='s', color='red', linestyle='-')
plt.plot(openmp, old_refactored_solve, label='Refactored (task-based)', marker='s', color='green', linestyle='-')
plt.plot(openmp, refactored_solve, label='Refactored (loop-parallelism)', marker='s', color='blue', linestyle='-')

plt.xscale('log')
plt.yscale('log')

# Set x-ticks explicitly to the OpenMP values
plt.xticks(ticks=openmp, labels=openmp)

plt.xlabel('Number of cores')
plt.ylabel('Execution time [s]')
# plt.title('Parallel Performance - Czarny Geometry (2048x2048)')
plt.legend(loc='lower left')
plt.grid(False)
plt.show()


plt.savefig("speed1.png", dpi=300, bbox_inches="tight", pad_inches=1)

## Plotting the speedup
plt.figure(figsize=(10, 6))
plt.rcParams.update({'font.size': 14})

plt.plot(openmp, openmp, label='Ideal Speedup', color='black', linestyle='--')

plt.plot(openmp, original_speed_up, label='Original (task-based)', marker='s', color='red', linestyle='-')
plt.plot(openmp, old_refactored_speed_up, label='Refactored (task-based)', marker='s', color='green', linestyle='-')
plt.plot(openmp, refactored_speed_up, label='Refactored (loop-parallelism)', marker='s', color='blue', linestyle='-')


plt.xscale('log')
plt.yscale('log')

# Set x-ticks explicitly to the OpenMP values
plt.xticks(ticks=openmp, labels=openmp)
plt.yticks(ticks=openmp, labels=openmp)

plt.xlabel('Number of cores')
plt.ylabel('Speedup')
# plt.title('Parallel Speedup - Czarny Geometry (2048x2048)')
plt.legend()
plt.grid(False)
plt.show()
plt.savefig("speedup.png", dpi=300, bbox_inches="tight", pad_inches=1)

import matplotlib.pyplot as plt
import numpy as np

# Function to read timings from the file
def read_timings(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    timings = {}
    current_level = -1

    for line in lines:
        line = line.strip()
        if line.startswith("Level"):
            current_level = int(line.split()[1])
            timings[current_level] = []
        elif line and current_level != -1:
            timings[current_level].append(int(line) / 1000000.0)

    return timings

# Read timings from the file
timings = read_timings('timings2.txt')
num_threads = list(range(1, len(next(iter(timings.values()))) + 1))

# Plot the timings for each level
plt.figure(figsize=(12, 8))

for level, timing in timings.items():
    plt.plot(num_threads, timing, marker='o', linestyle='-', label=f'Level {level}')

plt.title('ApplyA depend(in, out), (16385x32768), Culham geometry')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (seconds)')
plt.xscale('log')
ticks = np.arange(0, 101, step=5)
ticks[0] = 1
plt.xticks(ticks)
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.show()

import re
import matplotlib.pyplot as plt

### ------------------------ ###
### Input test configuration ###
### ------------------------ ###

folder = "../example_data/data_give_1537x2048"
nr = 1537
ntheta = 2048

### Define the number of sockets and socket size ###
sockets = 4
socket_size = 14
cores = 56

### Version 1: Custom core list ###
# core_list=[1, 2, 4, 8, 16, 32, 56]

### Version 2: Test all core configurations ###
core_list = [i for i in range(1, cores + 1)]

### ---------------- ###
### Used for x-ticks ###
### ---------------- ###

total_threads = sockets * socket_size 

def get_divisors(n):
    """Return the list of divisors of n."""
    divisors = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisors.append(i)
    return divisors

divisors = get_divisors(total_threads)

### ----------------------- ###
### Extract data from table ###
### ----------------------- ###

def parse_data_from_file(filename, num_cores):

    with open(filename, 'r') as file:
        lines = file.readlines()

    headers = []
    current_line = 0
    if num_cores == 1:
        for line in lines:
            if ('Metric' in line and 'Core' in line) or ('Metric' in line and 'HWThread' in line):
                headers.append(current_line)
            current_line += 1
    else:
        for line in lines:
            if 'Metric' in line and 'Sum' in line:
                headers.append(current_line)
            current_line += 1

    matrix = [[''] * 3 for _ in range(10)] 
    matrix[0][0] = 'Metric'
    matrix[1][0] = 'Runtime (RDTSC) [s]'
    matrix[2][0] = 'Runtime unhalted [s]'
    matrix[3][0] = 'CPI'
    matrix[4][0] = 'Energy [J]'
    matrix[5][0] = 'Power [W]'
    matrix[6][0] = 'DP [MFLOP/s]'
    matrix[7][0] = 'Memory bandwidth [MBytes/s]'
    matrix[8][0] = 'Memory data volume [GBytes]'
    matrix[9][0] = 'Operational intensity'
    offset = [0, 2, 3, 5, 6, 7, 10, 18, 19, 20] # Identifies the row in the table


    matrix[0][1] = 'Setup'
    if num_cores==1:
        matrix[1][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[1]]).group(0))
        matrix[2][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[2]]).group(0))
    else:
        matrix[1][1] = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[1]])[2])
        matrix[2][1] = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[2]])[2])

    matrix[3][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[3]]).group(0))
    matrix[4][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[4]]).group(0))
    matrix[5][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[5]]).group(0))
    matrix[6][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[6]]).group(0))
    matrix[7][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[7]]).group(0))
    matrix[8][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[8]]).group(0))
    matrix[9][1] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[0] + offset[9]]).group(0))


    matrix[0][2] = 'Solve'
    if num_cores==1:
        matrix[1][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[1]]).group(0))
        matrix[2][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[2]]).group(0))
    else:
        matrix[1][2] = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[1]])[2])
        matrix[2][2] = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[2]])[2])

    matrix[3][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[3]]).group(0))
    matrix[4][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[4]]).group(0))
    matrix[5][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[5]]).group(0))
    matrix[6][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[6]]).group(0))
    matrix[7][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[7]]).group(0))
    matrix[8][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[8]]).group(0))
    matrix[9][2] = float(re.search(r"[-+]?\d*\.\d+|\d+", lines[headers[1] + offset[9]]).group(0))

    return matrix

### --------- ###
### Plot data ###
### --------- ###

def main():
    ## --------- ##
    ## Load Data ##
    ## --------- ##

    compact_data = []
    for num_cores in core_list:
        filename_compact_mem = folder + '/COMPACT_MEM_DP_%d.txt' % num_cores
        matrix_compact_mem = parse_data_from_file(filename_compact_mem, num_cores)
        compact_data.append(matrix_compact_mem)

        # print(f"Results for {filename_compact_mem}:")
        # for row in matrix_compact_mem:
        #     print(row)
        # print("\n" + "="*50 + "\n") 

    spread_data = []
    for num_cores in core_list:
        filename_spread_mem = folder + '/SPREAD_MEM_DP_%d.txt' % num_cores
        matrix_spread_mem = parse_data_from_file(filename_spread_mem, num_cores)
        spread_data.append(matrix_spread_mem)

        # print(f"Results for {filename_spread_mem}:")
        # for row in matrix_spread_mem:
        #     print(row)
        # print("\n" + "="*50 + "\n") 

    ## ------------- ##
    ## Plot: Runtime ##
    ## ------------- ##

    setup_compact = [sublist[1][1] for sublist in compact_data]
    setup_spread = [sublist[1][1] for sublist in spread_data]
    solve_compact = [sublist[1][2] for sublist in compact_data]
    solve_spread = [sublist[1][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve [s] - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve [s] - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup [s] - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup [s] - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    initial_time = solve_compact[0]  # First value in TotalTime (with the fewest threads)
    optimal_times = [initial_time / t for t in core_list]  # Ideal scaling times
    # Plot the optimal line
    plt.plot(core_list, optimal_times, 'k--', label='_nolegend_')

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('Execution time [s]')
    plt.title(f'Parallel Performance - Czarny Geometry ({nr} x {ntheta})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('1_Runtime.png')

    ## ---------------------------- ##
    ## Plot: Cycles per instruction ##
    ## ---------------------------- ##

    setup_compact = [sublist[3][1] for sublist in compact_data]
    setup_spread = [sublist[3][1] for sublist in spread_data]
    solve_compact = [sublist[3][2] for sublist in compact_data]
    solve_spread = [sublist[3][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve: CPI - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve: CPI - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup: CPI - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup: CPI - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    initial_time = solve_compact[0]  # First value in TotalTime (with the fewest threads)
    optimal_times = [initial_time * t for t in core_list]  # Ideal scaling times
    # Plot the optimal line
    plt.plot(core_list, optimal_times, 'k--', label='_nolegend_')

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    # Labels and title
    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('Cycles per Instruction')
    plt.title(f'Cycles per Instruction - Czarny Geometry ({nr} x {ntheta})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('2_CyclesPerInstruction.png')

    ## ------------------------------------------------- ##
    ## Plot: Double Floating Point Operations Per Second ##
    ## ------------------------------------------------- ##

    setup_compact = [sublist[6][1] for sublist in compact_data]
    setup_spread = [sublist[6][1] for sublist in spread_data]
    solve_compact = [sublist[6][2] for sublist in compact_data]
    solve_spread = [sublist[6][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve: DP [MFLOP/s] - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve: DP [MFLOP/s] - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup: DP [MFLOP/s] - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup: DP [MFLOP/s] - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    initial_time = solve_compact[0]  # First value in TotalTime (with the fewest threads)
    optimal_times = [initial_time * t for t in core_list]  # Ideal scaling times
    # Plot the optimal line
    plt.plot(core_list, optimal_times, 'k--', label='_nolegend_')

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    # Labels and title
    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('DP [MFLOP/s]')
    plt.title(f'DP [MFLOP/s] - Czarny Geometry ({nr} x {ntheta})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('3_FLOPS_DP.png')

    ## ---------------------- ##
    ## Plot: Memory Bandwidth ##
    ## ---------------------- ##

    setup_compact = [sublist[7][1] for sublist in compact_data]
    setup_spread = [sublist[7][1] for sublist in spread_data]
    solve_compact = [sublist[7][2] for sublist in compact_data]
    solve_spread = [sublist[7][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve: Memory bandwidth [MBytes/s] - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve: Memory bandwidth [MBytes/s] - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup: Memory bandwidth [MBytes/s] - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup: Memory bandwidth [MBytes/s] - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    initial_time = solve_compact[0]  # First value in TotalTime (with the fewest threads)
    optimal_times = [initial_time * t for t in core_list]  # Ideal scaling times
    # Plot the optimal line
    plt.plot(core_list, optimal_times, 'k--', label='_nolegend_')

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    # Labels and title
    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('Memory bandwidth [MBytes/s]')
    plt.title(f'Memory bandwidth [MBytes/s] - Czarny Geometry ({nr} x {ntheta})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('4_Memory_Bandwidth.png')

    ## ------------------------ ##
    ## Plot: Memory Data Volume ##
    ## ------------------------ ##

    setup_compact = [sublist[8][1] for sublist in compact_data]
    setup_spread = [sublist[8][1] for sublist in spread_data]
    solve_compact = [sublist[8][2] for sublist in compact_data]
    solve_spread = [sublist[8][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve: Memory bandwidth [MBytes/s] - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve: Memory bandwidth [MBytes/s] - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup: Memory bandwidth [MBytes/s] - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup: Memory bandwidth [MBytes/s] - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    # Labels and title
    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('Memory data volume [GBytes]')
    plt.title(f'Memory data volume [GBytes] - Czarny Geometry ({nr} x {ntheta})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('5_Memory_DataVolume.png')

    ## --------------------------- ##
    ## Plot: Operational Intensity ##
    ## --------------------------- ##

    setup_compact = [sublist[9][1] for sublist in compact_data]
    setup_spread = [sublist[9][1] for sublist in spread_data]
    solve_compact = [sublist[9][2] for sublist in compact_data]
    solve_spread = [sublist[9][2] for sublist in spread_data]

    plt.figure(figsize=(10, 6))

    plt.plot(core_list, solve_compact, marker='o', markersize = 5, label='Solve: Operational intensity - compact pinning', color='red', linewidth=2)
    plt.plot(core_list, solve_spread, marker='o', markersize = 5, label='Solve: Operational intensity - spread pinning', color='blue', linewidth=2)
    plt.plot(core_list, setup_compact, marker='o', markersize = 5, label='Setup: Operational intensity - compact pinning', color='orangered', linestyle='--', linewidth=2)
    plt.plot(core_list, setup_spread, marker='o', markersize = 5, label='Setup: Operational intensity - spread pinning', color='dodgerblue', linestyle='--', linewidth=2)

    plt.xscale('log')
    plt.xticks(divisors, labels=[str(d) for d in divisors])
    # Add vertical lines at multiples of socket_size
    plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
    for i in range(1, sockets + 1):  # 1 to number of sockets
        x_position = i * socket_size
        plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

    plt.ylim(bottom= 0.7 * min(solve_compact + setup_compact), top= 1.3 * max(solve_compact + setup_compact))

    # Labels and title
    plt.yscale('log')
    plt.xlabel('Number of threads')
    plt.ylabel('Operational intensity')
    plt.title(f'Operational intensity - Czarny Geometry ({nr} x {ntheta})')

    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.savefig('6_Operational_Intensity.png')


if __name__ == '__main__':
    main()
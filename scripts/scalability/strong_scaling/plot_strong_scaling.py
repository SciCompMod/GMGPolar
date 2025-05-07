

# import pandas as pd
# import matplotlib.pyplot as plt

# # Load the strong scaling results from CSV
# data = pd.read_csv('strong_scaling_results.csv')

# # Define the number of sockets and socket size
# sockets = 4
# socket_size = 14

# # Calculate the total number of threads
# total_threads = sockets * socket_size  # 2 * 24 = 48

# # Function to calculate divisors
# def get_divisors(n):
#     """Return the list of divisors of n."""
#     divisors = []
#     for i in range(1, n + 1):
#         if n % i == 0:
#             divisors.append(i)
#     return divisors

# # Get all divisors of total_threads
# divisors = get_divisors(total_threads)

# nr = data["nr"][0]
# ntheta = data["ntheta"][0]
# geometry_name = data["geometry"][0]

# # Create the plot
# plt.figure(figsize=(12, 8))

# # Define the colors for each component
# colors = {
#     'Total': '#d62728',            # Total in red
#     'Smoothing': '#ff7f0e',        # Smoothing in orange
#     'Residual': '#1f77b4',         # Residual in blue
#     'Direct Solver': '#2ca02c'     # Direct Solver in green
# }
# # Plot the times with the chosen colors for each component
# plt.plot(data['Threads'], data['t_avg_MGC_total'], marker='o', label='Multigrid Cycle (Total)', color=colors['Total'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_preSmoothing'] + data['t_avg_MGC_postSmoothing'], marker='o', label='Smoothing', color=colors['Smoothing'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_residual'], marker='o', label='Residual', color=colors['Residual'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_directSolver'], marker='o', label='Direct Solver', color=colors['Direct Solver'], linewidth=2)

# # Calculate and plot the optimal scaling line
# initial_time = data['t_avg_MGC_total'][0]  # First value in TotalTime (with the fewest threads)
# optimal_times = [initial_time / t for t in data['Threads']]  # Ideal scaling times

# # Plot the optimal line
# plt.plot(data['Threads'], optimal_times, 'k--', label='__no_legend__')

# # Set x-axis to logarithmic scale
# plt.xscale('log')

# # Set x-axis ticks to the calculated divisors, limiting them to a few key ones
# plt.xticks(divisors, labels=[str(d) for d in divisors if d <= total_threads])

# # Add vertical lines at multiples of socket_size
# plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
# for i in range(1, sockets + 1):  # 1 to number of sockets
#     x_position = i * socket_size
#     plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

# # Labels and title
# plt.yscale('log')
# plt.xlabel('Number of Threads')
# plt.ylabel('Execution Time [s]')
# plt.title(f'Parallel Performance - {geometry_name} Geometry ({nr} x {ntheta})')
# plt.legend(loc='upper right', frameon=False)

# plt.grid(True, which="both", ls="--", linewidth=0.5)
# plt.savefig('image_strong_scaling.png')

# # Show the plot
# plt.show()




















import pandas as pd
import matplotlib.pyplot as plt

# Load the strong scaling results from CSV
data = pd.read_csv('strong_scaling_results.csv')

# Define the number of sockets and socket size
sockets = 2
socket_size = 32

# Calculate the total number of threads
total_threads = sockets * socket_size  # 2 * 24 = 48

# Function to calculate divisors
def get_divisors(n):
    """Return the list of divisors of n."""
    divisors = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisors.append(i)
    return divisors

# Get all divisors of total_threads
divisors = get_divisors(total_threads)

nr = data["nr"][0]
ntheta = data["ntheta"][0]
geometry_name = data["geometry"][0]

### General Plot ###

plt.figure(figsize=(12, 8))

colors = {
    'Total Solve': '#d62728', 
    'Multigrid Iterations': '#ff7f0e',     
    'Check Convergence': '#1f77b4',
    'Total Setup': '#2ca02c',
    'Initial Approximation': 'gray'
}

plt.plot(data['Threads'], data['t_solve_total'], marker='o', label='Solve (Total)', color=colors['Total Solve'], linewidth=2)
plt.plot(data['Threads'], data['t_solve_multigrid_iterations'], marker='o', label='Multigrid Iterations', color=colors['Multigrid Iterations'], linewidth=2)
plt.plot(data['Threads'], data['t_check_convergence'], marker='o', label='Check Convergence', color=colors['Check Convergence'], linewidth=2)
plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors['Total Setup'], linewidth=2)
if data["FMG"][0] > 0:
    plt.plot(data['Threads'], data['t_solve_initial_approximation'], marker='o', label='Initial Approximation', color=colors['Initial Approximation'], linewidth=2)

initial_time = data['t_solve_total'][0]  # First value in TotalTime (with the fewest threads)
optimal_times = [initial_time / t for t in data['Threads']]  # Ideal scaling times

# Plot the optimal line
plt.plot(data['Threads'], optimal_times, 'k--', label='__no_legend__')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')

# Set x-axis ticks to the calculated divisors, limiting them to a few key ones
plt.xticks(divisors, labels=[str(d) for d in divisors if d <= total_threads])

# Add vertical lines at multiples of socket_size
plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
for i in range(1, sockets + 1):  # 1 to number of sockets
    x_position = i * socket_size
    plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

# Labels and title
plt.yscale('log')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Strong Scalability - {geometry_name} Geometry ({nr} x {ntheta}) - Method: {data["stencil_method"][0]}')

plt.legend(loc='upper right', frameon=False)

plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_strong_scaling.png')

### Multigrid Cycle Plot ###

plt.figure(figsize=(12, 8))

colors = {
    'Total': '#d62728',
    'Smoothing': '#ff7f0e',
    'Residual': '#1f77b4',
    'Direct Solver': '#2ca02c'
}

plt.plot(data['Threads'], data['t_avg_MGC_total'], marker='o', label='Multigrid Cycle (Total)', color=colors['Total'], linewidth=2)
plt.plot(data['Threads'], data['t_avg_MGC_preSmoothing'] + data['t_avg_MGC_postSmoothing'], marker='o', label='Smoothing', color=colors['Smoothing'], linewidth=2)
plt.plot(data['Threads'], data['t_avg_MGC_residual'], marker='o', label='Residual', color=colors['Residual'], linewidth=2)
plt.plot(data['Threads'], data['t_avg_MGC_directSolver'], marker='o', label='Direct Solver', color=colors['Direct Solver'], linewidth=2)

initial_time = data['t_avg_MGC_total'][0]  # First value in TotalTime (with the fewest threads)
optimal_times = [initial_time / t for t in data['Threads']]  # Ideal scaling times

# Plot the optimal line
plt.plot(data['Threads'], optimal_times, 'k--', label='__no_legend__')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')

# Set x-axis ticks to the calculated divisors, limiting them to a few key ones
plt.xticks(divisors, labels=[str(d) for d in divisors if d <= total_threads])

# Add vertical lines at multiples of socket_size
plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
for i in range(1, sockets + 1):  # 1 to number of sockets
    x_position = i * socket_size
    plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

# Labels and title
plt.yscale('log')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Strong Scalability: Multigrid Cycle - {geometry_name} Geometry ({nr} x {ntheta}) - Method: {data["stencil_method"][0]}')

plt.legend(loc='upper right', frameon=False)

plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_strong_scaling_MGC.png')

### Setup Plot ###

plt.figure(figsize=(12, 8))

colors = {
    'Total': '#d62728',
    'Create Levels': '#ff7f0e',
    'Smoother': '#1f77b4',
    'Direct Solver': '#2ca02c' 
}

plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors['Total'], linewidth=2)
plt.plot(data['Threads'], data['t_setup_createLevels'], marker='o', label='Create Levels', color=colors['Create Levels'], linewidth=2)
plt.plot(data['Threads'], data['t_setup_smoother'], marker='o', label='Smoother', color=colors['Smoother'], linewidth=2)
plt.plot(data['Threads'], data['t_setup_directSolver'], marker='o', label='Direct Solver', color=colors['Direct Solver'], linewidth=2)

initial_time = data['t_setup_total'][0]  # First value in TotalTime (with the fewest threads)
optimal_times = [initial_time / t for t in data['Threads']]  # Ideal scaling times

# Plot the optimal line
plt.plot(data['Threads'], optimal_times, 'k--', label='__no_legend__')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')

# Set x-axis ticks to the calculated divisors, limiting them to a few key ones
plt.xticks(divisors, labels=[str(d) for d in divisors if d <= total_threads])

# Add vertical lines at multiples of socket_size
plt.axvline(x=1, color='darkgray', linestyle='--', linewidth=1.5)
for i in range(1, sockets + 1):  # 1 to number of sockets
    x_position = i * socket_size
    plt.axvline(x=x_position, color='darkgray', linestyle='--', linewidth=1.5)

# Labels and title
plt.yscale('log')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Strong Scalability: Setup - {geometry_name} Geometry ({nr} x {ntheta}) - Method: {data["stencil_method"][0]}')

plt.legend(loc='upper right', frameon=False)

plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_strong_scaling_Setup.png')




# plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

# plt.xlabel('Number of Threads')
# plt.ylabel('Execution Time [s]')
# plt.title(f'Weak Scaling - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

# plt.legend(loc='lower right', frameon=False)
# plt.grid(True, which="both", ls="--", linewidth=0.5)
# plt.savefig('image_weak_scaling.png')

# ### Multigrid Cycle Plot ###

# plt.figure(figsize=(12, 8))

# colors = {
#     'Total': '#d62728',
#     'Smoothing': '#ff7f0e',
#     'Residual': '#1f77b4',
#     'Direct Solver': '#2ca02c'
# }

# plt.plot(data['Threads'], data['t_avg_MGC_total'], marker='o', label='Multigrid Cycle (Total)', color=colors['Total'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_preSmoothing'] + data['t_avg_MGC_postSmoothing'], marker='o', label='Smoothing', color=colors['Smoothing'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_residual'], marker='o', label='Residual', color=colors['Residual'], linewidth=2)
# plt.plot(data['Threads'], data['t_avg_MGC_directSolver'], marker='o', label='Direct Solver', color=colors['Direct Solver'], linewidth=2)

# plt.xscale('log')
# plt.yscale('log')

# plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

# plt.xlabel('Number of Threads')
# plt.ylabel('Execution Time [s]')
# plt.title(f'Weak Scaling: Multigrid Cycle - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

# plt.legend(loc='lower right', frameon=False)
# plt.grid(True, which="both", ls="--", linewidth=0.5)
# plt.savefig('image_weak_scaling_MGC.png')

# ### Setup Plot ###

# plt.figure(figsize=(12, 8))

# colors = {
#     'Total': '#d62728',
#     'Create Levels': '#ff7f0e',
#     'Smoother': '#1f77b4',
#     'Direct Solver': '#2ca02c' 
# }

# plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors['Total'], linewidth=2)
# plt.plot(data['Threads'], data['t_setup_createLevels'], marker='o', label='Create Levels', color=colors['Create Levels'], linewidth=2)
# plt.plot(data['Threads'], data['t_setup_smoother'], marker='o', label='Smoother', color=colors['Smoother'], linewidth=2)
# plt.plot(data['Threads'], data['t_setup_directSolver'], marker='o', label='Direct Solver', color=colors['Direct Solver'], linewidth=2)

# plt.xscale('log')
# plt.yscale('log')

# plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

# plt.xlabel('Number of Threads')
# plt.ylabel('Execution Time [s]')
# plt.title(f'Weak Scaling: Setup - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

# plt.legend(loc='lower right', frameon=False)
# plt.grid(True, which="both", ls="--", linewidth=0.5)
# plt.savefig('image_weak_scaling_Setup.png')


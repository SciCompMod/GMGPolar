import pandas as pd
import matplotlib.pyplot as plt

# Load the weak scaling results from CSV
data = pd.read_csv('weak_scaling_results.csv')

nr_small = data["nr"][0]
ntheta_small = data["ntheta"][0]
nr_large = data["nr"].iloc[-1]
ntheta_large = data["ntheta"].iloc[-1]
geometry_name = data["geometry"][0]

### General Plot ###

plt.figure(figsize=(12, 8))

colors = {
    'Total Solve': '#d62728', 
    'Multigrid Iterations': '#ff7f0e',     
    'Check Convergence': '#2ca02c', 
    'Total Setup': '#1f77b4',
    'Initial Approximation': 'gray'
}

plt.plot(data['Threads'], data['t_solve_total'], marker='o', label='Solve (Total)', color=colors['Total Solve'], linewidth=2)
plt.plot(data['Threads'], data['t_solve_multigrid_iterations'], marker='o', label='Multigrid Iterations', color=colors['Multigrid Iterations'], linewidth=2)
plt.plot(data['Threads'], data['t_check_convergence'], marker='o', label='Check Convergence', color=colors['Check Convergence'], linewidth=2)
plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors['Total Setup'], linewidth=2)
if data["FMG"][0] > 0:
    plt.plot(data['Threads'], data['t_solve_initial_approximation'], marker='o', label='Initial Approximation', color=colors['Initial Approximation'], linewidth=2)


plt.xscale('log')
plt.yscale('log')

plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Weak Scaling - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

plt.legend(loc='lower right', frameon=False)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_weak_scaling.png')

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

plt.xscale('log')
plt.yscale('log')

plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Weak Scaling: Multigrid Cycle - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

plt.legend(loc='lower right', frameon=False)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_weak_scaling_MGC.png')

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

plt.xscale('log')
plt.yscale('log')

plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])

plt.xlabel('Number of Threads')
plt.ylabel('Execution Time [s]')
plt.title(f'Weak Scaling: Setup - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')

plt.legend(loc='lower right', frameon=False)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.savefig('image_weak_scaling_Setup.png')


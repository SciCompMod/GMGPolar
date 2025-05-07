import pandas as pd
import matplotlib.pyplot as plt

# Load the weak scaling results from CSV
data = pd.read_csv('weak_scaling_results.csv')

nr_small = data["nr"][0]
ntheta_small = data["ntheta"][0]
nr_large = data["nr"].iloc[-1]
ntheta_large = data["ntheta"].iloc[-1]
geometry_name = data["geometry"][0]

# Define common colors
colors_general = {
    'Total Solve': '#d62728', 
    'Multigrid Iterations': '#ff7f0e',     
    'Check Convergence': '#2ca02c', 
    'Total Setup': '#1f77b4',
    'Initial Approximation': 'gray'
}

colors_mgc = {
    'Total': '#d62728',
    'Smoothing': '#ff7f0e',
    'Residual': '#1f77b4',
    'Direct Solver': '#2ca02c'
}

colors_setup = {
    'Total': '#d62728',
    'Create Levels': '#ff7f0e',
    'Smoother': '#1f77b4',
    'Direct Solver': '#2ca02c' 
}

### General Plot - Original and Clean ###
def plot_general(clean=False):
    plt.figure(figsize=(12, 8))
    
    plt.plot(data['Threads'], data['t_solve_total'], marker='o', label='Solve (Total)', color=colors_general['Total Solve'], linewidth=2)
    plt.plot(data['Threads'], data['t_solve_multigrid_iterations'], marker='o', label='Multigrid Iterations', color=colors_general['Multigrid Iterations'], linewidth=2)
    plt.plot(data['Threads'], data['t_check_convergence'], marker='o', label='Check Convergence', color=colors_general['Check Convergence'], linewidth=2)
    plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors_general['Total Setup'], linewidth=2)
    if data["FMG"][0] > 0:
        plt.plot(data['Threads'], data['t_solve_initial_approximation'], marker='o', label='Initial Approximation', color=colors_general['Initial Approximation'], linewidth=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])
    plt.xlabel('Number of Threads')
    plt.ylabel('Execution Time [s]')
    plt.legend(loc='lower right', frameon=False)  # Legend kept for both versions
    
    if not clean:
        plt.title(f'Weak Scaling - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')
        plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.savefig('image_weak_scaling_with_title.png')
    else:
        plt.savefig('clean_image_weak_scaling.png')
    
    plt.close()

### Multigrid Cycle Plot - Original and Clean ###
def plot_mgc(clean=False):
    plt.figure(figsize=(12, 8))
    
    plt.plot(data['Threads'], data['t_avg_MGC_total'], marker='o', label='Multigrid Cycle (Total)', color=colors_mgc['Total'], linewidth=2)
    plt.plot(data['Threads'], data['t_avg_MGC_preSmoothing'] + data['t_avg_MGC_postSmoothing'], marker='o', label='Smoothing', color=colors_mgc['Smoothing'], linewidth=2)
    plt.plot(data['Threads'], data['t_avg_MGC_residual'], marker='o', label='Residual', color=colors_mgc['Residual'], linewidth=2)
    plt.plot(data['Threads'], data['t_avg_MGC_directSolver'], marker='o', label='Direct Solver', color=colors_mgc['Direct Solver'], linewidth=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])
    plt.xlabel('Number of Threads')
    plt.ylabel('Execution Time [s]')
    plt.legend(loc='right', frameon=False)  # Legend kept for both versions
    
    if not clean:
        plt.title(f'Weak Scaling: Multigrid Cycle - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')
        plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.savefig('image_weak_scaling_MGC_with_title.png')
    else:
        plt.savefig('clean_image_weak_scaling_MGC.png')
    
    plt.close()

### Setup Plot - Original and Clean ###
def plot_setup(clean=False):
    plt.figure(figsize=(12, 8))
    
    plt.plot(data['Threads'], data['t_setup_total'], marker='o', label='Setup (Total)', color=colors_setup['Total'], linewidth=2)
    plt.plot(data['Threads'], data['t_setup_createLevels'], marker='o', label='Create Levels', color=colors_setup['Create Levels'], linewidth=2)
    plt.plot(data['Threads'], data['t_setup_smoother'], marker='o', label='Smoother', color=colors_setup['Smoother'], linewidth=2)
    plt.plot(data['Threads'], data['t_setup_directSolver'], marker='o', label='Direct Solver', color=colors_setup['Direct Solver'], linewidth=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(data['Threads'].unique(), labels=[str(int(t)) for t in data['Threads'].unique()])
    plt.xlabel('Number of Threads')
    plt.ylabel('Execution Time [s]')
    plt.legend(loc='right', frameon=False)  # Legend kept for both versions
    
    if not clean:
        plt.title(f'Weak Scaling: Setup - Geometry: {geometry_name} - Method: {data["stencil_method"][0]} - Mesh Sizes: ({nr_small}x{ntheta_small}) → ({nr_large}x{ntheta_large})')
        plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.savefig('image_weak_scaling_Setup_with_title.png')
    else:
        plt.savefig('clean_image_weak_scaling_Setup.png')
    
    plt.close()

# Generate all plots
plot_general(clean=False)  # Original version with title and grid
plot_general(clean=True)   # Clean version without title and grid

plot_mgc(clean=False)      # Original version with title and grid
plot_mgc(clean=True)       # Clean version without title and grid

plot_setup(clean=False)    # Original version with title and grid
plot_setup(clean=True)     # Clean version without title and grid
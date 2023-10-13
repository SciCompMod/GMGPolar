import pandas as pd
import numpy as np
import sys
# Use backend to not plot on UI
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import join, exists, dirname
from os import makedirs
import os.path

### Plots scaling of FLOPS and Caches (saturation) scaling from 0 to n Cores
### as read from data frame

colors  =   [
            [1.00,  0.49,  0.06],
            [0.17,  0.63,  0.18],
            [0.83,  0.15,  0.16],
            [0.13, 0.47, 0.69],
            [0.58,  0.40,  0.74],
            [0.53,  0.35,  0.27],
            [0.92,  0.46,  0.77],
            [0.50,  0.50,  0.50],
            [0.66,  0.85,  0.06],                
            [0.06,  0.85,  0.85],
            [0.85,  0.15,  0.85],
            [0.75,  0.75,  0.75],
            [0.80,  0.67,  0.00],
            [0.30,  0.65,  1.00],
            [0.40,  0.07,  0.00],
            [1.00, 0.80, 0.87],
            [1.00, 0.87, 0.20],
            [1.00, 0.80, 0.70],
            [0.80, 0.87, 1.00]];  

def plot_scaling(path_out, fname, benchname, df, title, ylabel, saturation_limit=0, logplot = True, colors=colors):
    '''
    Plot different timings or other benchmarks against the number of cores used.
    @param path_out Path to output files.
    @param fname First part of filename.
    @param benchname Benchmark name and postfix for filename.
    @param df Data frame with results per core scaling.
    @param title Plot title.
    @param ylabel Y-axis label.
    @param saturation_limit Saturation limit to be plotted for MEM_DP scaling.
    '''
    fontsize = 20

    fig = plt.figure(figsize=(12, 10))
    # fig.subplots_adjust(bottom=0.3) # use for interactive plotting, for savefig, legend+figure is adjusted with bbox_inches='tight'
    ax = fig.add_subplot() 
    start_optimal_line = df.loc[0,[timer for timer in df.iloc[:,1:].columns if 'Total' in timer]].max()
    optimal_line = 1.0/df['Cores'].values * start_optimal_line
    if logplot:
        if benchname == 'Timings':
            plt.loglog(df['Cores'], optimal_line, linewidth=2, linestyle='dashed', color='black', label='Optimal')
        for i in range(len(df.iloc[:,1:].columns)):
            plt.loglog(df['Cores'], df.iloc[:,i+1], linewidth=2, label=df.columns[i+1], color=colors[i]) # Cores is assumed to be in first column
    else:
        if benchname == 'Timings':
            plt.plot(df['Cores'], optimal_line, linewidth=2, linestyle='dashed', color='black', label='Optimal')
        for i in range(len(df.iloc[:,1:].columns)):
            plt.plot(df['Cores'], df.iloc[:,i+1], linewidth=2, label=df.columns[i+1], color=colors[i]) # Cores is assumed to be in first column

    if benchname == 'MEM_DP':
        if saturation_limit > 0:
            plt.plot(df['Cores'], saturation_limit * np.ones(len(df['Cores'])), linestyle='dotted', linewidth=3, color=[0, 0, 0]) 
            ax.text(1, saturation_limit + 3, 'Memory bandwith (AXPY) (' + str(saturation_limit) + ' GBytes/s)', fontsize=14)
            ax.set_ylim(0, saturation_limit + 10)

    ax.legend(bbox_to_anchor=(0, 0, 1, -0.1), mode="expand", ncols=2)
    ax.set_title(title, fontsize=fontsize+6)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xlabel('Number of cores', fontsize=fontsize)
    ax.set_xticks(df['Cores'])
    ax.set_xticklabels(df['Cores'])

    path_out = join(path_out, 'figures')
    if not exists(path_out):
        makedirs(path_out)
    plt.savefig(join(path_out, fname + '_' + benchname.lower()), dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()   


def main(benchmarks=['FLOPS_DP']):

    file_prefix = 'example' # provide correct slurm job id
    if file_prefix != '':
        file_prefix = file_prefix + '-'
    file_postfix = ''
    if file_postfix != '':
        file_postfix = '-' + file_postfix
    problem = 6
    divideBy2 = 7 # steps of division of initial grid
    nr_exp = 4
    mod_pk = 1
    smoother = 3
    extrapolation = 1

    nodes = 1
    ranks = 1
    maxCores = 128 # maxCores simulated in scaling

    maxCoresPlot = 128 # maxCores to plot
    plot_counter = {} # dict on which counter to plot
    plot_counter['Total setup'] = 0
    plot_counter['Building system matrix A and RHS'] = 1
    plot_counter['Factorization of coarse operator Ac'] = 1
    plot_counter['Building intergrid operators (e.g. projections)'] = 1
    plot_counter['Building smoothing operators A_sc'] = 1
    plot_counter['Factorizing smoothing operators A_sc'] = 1
    plot_counter['Total multigrid cycle'] = 1
    plot_counter['Complete smoothing'] = 1
    plot_counter['Computing residual'] = 1
    plot_counter['Applying restriction'] = 1
    plot_counter['Solve coarse system'] = 1
    plot_counter['Applying prolongation (+ coarse grid correction)'] = 1
    plot_counter['Computing residual on finest level'] = 1
    plot_counter['Computing final error'] = 1
    plot_counter['Total application of A'] = 1
    plot_counter['Evaluation of arr, art, and att'] = 1
    plot_counter['Evaluation of alpha and beta'] = 1
    plot_counter['Computing determinant of Jacobian of inverse mapping'] = 1
    plot_counter['Computing exact solution'] = 1
    plot_counter['Total execution time'] = 1

    plot_regions = {} # dict on which likwid region to plot
    plot_regions['Setup'] = 1
    plot_regions['Iteration'] = 1

    bench_to_unit = {}
    bench_to_unit['FLOPS_DP'] = 'MFLOP/s'
    bench_to_unit['MEM_DP'] = 'MBytes/s'   

    ## saturation_limit for MEM_DP scaling (node specific and needs to be adapted).
    saturation_limit = 80

    fname = file_prefix + 'p' + str(problem) + '-r' + str(nr_exp) + '-dbt' + str(divideBy2) + '-mpk' + str(mod_pk) + '-s' + str(
        smoother) + '-e' + str(extrapolation) + '--N' + str(nodes) + '-R' + str(ranks) + '-maxC' + str(maxCores) + file_postfix
    
    path_to_files_rel = os.path.join('..', os.path.join('..', 'scaling_output')) # relative from read_output call
    path_to_files = os.path.join(os.path.dirname(__file__), os.path.join(path_to_files_rel))


    # Problem setting columns
    setting_cols = ['Problem', 'rExp', 'divB2', 'ModPK', 'Extrapolation', 'Nodes', 'Ranks']


    df = []
    for bench in benchmarks:
        df.append(pd.read_csv(
            join(path_to_files, fname + '_' + bench + '.csv'),
            dtype={'Problem': int, 'rExp': int, 'divB2': int, 'ModPK': int,
                    'Extrapolation': int, 'Nodes': int, 'Ranks': int,
                    'Cores': int, 'its': int}))
        
        # Check that the different number of threads/cores where only conducted
        # on one particular setting of the above columns
        if np.max(df[0].loc[:,setting_cols].nunique()) > 1:
            sys.exit('Error in input data, more than one setting found.')   
        # check that only one row per number of cores is available
        if df[0]['Cores'].nunique() != len(df[0]['Cores']):
            sys.exit('Error: Multiple times computed with the same number of cores.')                 

    # check content
    for i in range(len(df)):
        # plot timings
        plot_scaling(path_to_files, fname, 'Timings', df[i].loc[df[i]['Cores']<=maxCoresPlot, ['Cores'] + [k for k, v in plot_counter.items() if v==1]], 'Strong scaling of GMGPolar', 'Time (sec)')
        # plot likwid benchmark
        plot_scaling(path_to_files, fname, benchmarks[i], df[i].loc[df[i]['Cores']<=maxCoresPlot, ['Cores'] + [col for col in df[i].columns if bench_to_unit[benchmarks[i]] in col]], 'Strong scaling of GMGPolar', bench_to_unit[benchmarks[i]])


if __name__ == '__main__':
    main(benchmarks=['FLOPS_DP'])
    
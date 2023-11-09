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


# This file has been written and used to plot the read weak scaling results for
# Leleux, Schwarz, Kühn, Kruse, Rüde - Complexity analysis and scalability of a matrix-free extrapolated geometric multigrid (2023)
# In contrast to the other files with prefix lskkr23_*, it is written quick and dirty and only 
# for weak scaling data collected from different strong scaling experiments.

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

def plot_weak_scaling(path_out, fname, df, ylabel, logplot = True, colors=colors):
    '''
    Plot different timings or other benchmarks against the number of cores used.
    @param path_out Path to output files.
    @param fname First part of filename.
    @param df Data frame with results per core scaling.
    @param ylabel Y-axis label.
    '''
    fontsize = 24

    fig = plt.figure(figsize=(12, 10))
    # fig.subplots_adjust(bottom=0.3) # use for interactive plotting, for savefig, legend+figure is adjusted with bbox_inches='tight'
    ax = fig.add_subplot()
    start_plot = list(df.columns).index('Total multigrid cycle')
    start_optimal_line = df.loc[0,[timer for timer in df.iloc[:,start_plot:].columns if 'Total' in timer]].max()
    optimal_line = [start_optimal_line for i in df['Cores'].values]
    if logplot:
        plt.loglog(df['Cores'], optimal_line, linewidth=2, linestyle='dashed', color='black', label='Optimal')
        for i in range(len(df.iloc[:,start_plot:].columns)):
            col_label = df.columns[start_plot+i]
            col_label = col_label.replace('Computing residual on finest level', 'Computing residual (finest)')
            col_label = col_label.replace('Applying prolongation (+ coarse grid correction)', 'Applying prolongation (+ CGC)')
            plt.loglog(df['Cores'], df.iloc[:,start_plot+i], linewidth=2, label=col_label, color=colors[i]) # Cores is assumed to be in first column
    else:
        plt.plot(df['Cores'], optimal_line, linewidth=2, linestyle='dashed', color='black', label='Optimal')
        for i in range(len(df.iloc[:,start_plot:].columns)):
            col_label = df.columns[start_plot+i]
            col_label = col_label.replace('Computing residual on finest level', 'Computing residual (finest)')
            col_label = col_label.replace('Applying prolongation (+ coarse grid correction)', 'Applying prolongation (+ CGC)')            
            plt.plot(df['Cores'], df.iloc[:,start_plot+i], linewidth=2, label=col_label, color=colors[i]) # Cores is assumed to be in first column


    ax.legend(bbox_to_anchor=(0, 0, 1, -0.1), mode="expand", ncols=2,  fontsize=fontsize-6)
    # ax.set_title(title, fontsize=fontsize+6)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xlabel('Number of cores', fontsize=fontsize)
    ax.set_xticks(df['Cores'])
    ax.set_xticklabels(df['Cores'])

    plt.rcParams['xtick.labelsize']=fontsize-6
    plt.rcParams['ytick.labelsize']=fontsize-6

    path_out = join(path_out, 'figures')
    if not exists(path_out):
        makedirs(path_out)
    plt.savefig(join(path_out, fname), dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()


def main(problem=7, divideBy2=[4, 5, 6, 7], mod_pk=1):

    file_prefix = 'caro-paper-weak-scaling' 
    if file_prefix != '':
        file_prefix = file_prefix + '-'
    file_postfix = 'C1toC64_FLOPS_DP'
    if file_postfix != '':
        file_postfix = '-' + file_postfix
    nr_exp = 4
    smoother = 3
    extrapolation = 1

    nodes = 1
    ranks = 1

    path_to_perf_files_rel = os.path.join('..', os.path.join('..', 'scaling_output')) # relative from read_output call
    path_to_perf_files = os.path.join(os.path.dirname(__file__), os.path.join(path_to_perf_files_rel))

    filename_input = 'p' + str(problem) + '-r' + str(nr_exp) + '-dbt' + str(divideBy2[0]) + '-' + str(divideBy2[-1]) + '-mpk' + str(mod_pk) + '-s' + str(
        smoother) + '-e' + str(extrapolation) + '--N' + str(nodes) + '-R' + str(ranks)
    
    try: # allow deletion of empty error files
        df_weak =  pd.read_csv(os.path.join(path_to_perf_files, file_prefix + filename_input + file_postfix + '.csv'), sep=';')

        plot_weak_scaling(path_to_perf_files, file_prefix + filename_input + file_postfix, df_weak.iloc[:, list(df_weak.columns).index('Cores'):], 'Time (sec)')
    except FileNotFoundError:
        err = 0 # no error file = no error    

if __name__ == "__main__":
    # weak scaling
    main(problem=7, divideBy2=[4,5,6,7], mod_pk=1)
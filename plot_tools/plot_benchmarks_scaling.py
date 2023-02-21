import pandas as pd
import numpy as np
import sys
# Use backend to not plot on UI
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import join, exists
from os import makedirs

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
            [0.75,  0.75,  0.75]];  

def plot_perf_per_core(path_out, fname, df, benchname, saturation_limit, colors=colors):  
    fontsize = 16

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot()
    plt.plot(df['Cores'], df[benchname])

    if benchname == 'CACHES':
        plt.plot(df['Cores'], saturation_limit * np.ones(len(df['Cores'])), linestyle='dotted', linewidth=3, color=[0, 0, 0]) 
        ax.text(1, saturation_limit+3, 'Memory bandwith (AXPY) (' + str(saturation_limit) + ' GBytes/s)', fontsize=14)
        ax.set_ylim(0, 90)
        ax.set_title('Memory bandwidth saturation', fontsize=fontsize+6)
        ax.set_ylabel('Memory bandwidth (GBytes/s)', fontsize=fontsize)

    elif benchname == 'FLOPS_DP':
        ax.set_title('Flop performance in Multi-Threading', fontsize=fontsize+6)
        ax.set_ylabel('Flops (GFlops/s)', fontsize=fontsize)

    ax.set_xlabel('Number of cores used', fontsize=fontsize)


    path_out = join(path_out, 'figures')
    if not exists(path_out):
        makedirs(path_out)
    if benchname == 'CACHES':
        plt.savefig(join(path_out, fname + '_cache_saturation.jpg'), bbox_inches='tight')
    elif benchname == 'FLOPS_DP':
        plt.savefig(join(path_out, fname + '_flop_perf.jpg'), bbox_inches='tight')
    plt.close()   


def main():

    problem = 5
    nr_exp = 4
    mod_pk = 1
    smoother = 3
    extrapolation = 1

    nodes = 1
    ranks = 1
    maxCores = 14

    ## saturation_limit is node specific and needs to be adapted.
    saturation_limit = 80

    fname = 'p' + str(problem) + '-r' + str(nr_exp) + '-mpk' + str(mod_pk) + '-s' + str(
        smoother) + '-e' + str(extrapolation) + '--N' + str(nodes) + '-R' + str(ranks) + '-maxC' + str(maxCores)   
    path_to_files = 'plot_tools/'

    df = pd.read_csv(
        join(path_to_files, fname + '_benchmarks.csv'),
        dtype={'Problem': int, 'rExp': int, 'ModPK': int,
                'Extrapolation': int, 'Nodes': int, 'Ranks': int,
                'Cores': int, 'its': int})

    # benchmarks
    benchmarks = ['FLOPS_DP', 'CACHES']
    setting_cols = ['Problem', 'rExp', 'ModPK', 'Extrapolation', 'Nodes', 'Ranks']

    # check content
    for bench in benchmarks:
        bench_rows = np.where(df[bench].isnull()!=True)[0]
        if len(bench_rows) > 0:
            df_subframe = df.iloc[bench_rows].copy()

            # Check that the different number of threads/cores where only conducted
            # on one particular setting of the above columns
            if np.max(df_subframe.loc[:,setting_cols].nunique()) > 1:
                sys.exit('Error in input data, more than one setting found.')

            # TODO: If the same run was done multiple times, it is not clear which line to
            # take. This should be extended to take the sum of all these lines then
            # and ensure before that we have a minimum number of n identic runs...
            cores_used = df_subframe['Cores'].unique()
            if len(cores_used) != len(df_subframe['Cores']):
                sys.exit('Error: Multiple times computed with the same number of threads.')

            plot_perf_per_core(path_to_files, fname, df_subframe, bench, saturation_limit=saturation_limit)

            # TODO: Plot particular timings from table
            # ...




if __name__ == '__main__':
    main()
    
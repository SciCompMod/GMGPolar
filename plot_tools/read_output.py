import pandas as pd
import numpy as np
import os
import sys


def main():

    file_prefix = 'caro_likwid-perctr-flopsdp-compare-perfctr-withoutlikwid' # provide correct slurm job id
    file_postfix = 'zones copy'
    problem = 7
    divideBy2 = 5 # steps of division of initial grid
    nr_exp = 4
    mod_pk = 1
    smoother = 3
    extrapolation = 1

    nodes = 1
    ranks = 1
    maxCores = 64
    cores_used = [] # will be filled automatically

    path_to_perf_files_rel = os.path.join('..', os.path.join('..', 'output_scripts')) # relative from read_output call
    path_to_perf_files = os.path.join(os.path.dirname(__file__), os.path.join(path_to_perf_files_rel))

    cols_problem = ['Problem', 'rExp', 'divB2', 'ModPK', 'Extrapolation']

    cols_cluster = ['Nodes', 'Ranks', 'Cores']

    cols_convergence = ['its', '2-norm of error', 'inf-norm of error']

    cols_time = []
    cols_time += ['Setup', 'Building A and RHS', 'Factorization of coarse operator Ac',
                  'Building intergrid operators (e.g. projections)', 'Building smoothing operators A_sc', 
                  'Factorizing smoothing operators A_sc']
    cols_time += ['Total multigrid cycle', 'Complete smoothing', 'Computing residual',
                  'Applying restriction', 'Solve coarse system', 'Applying prolongation (+ coarse grid correction)']
    cols_time += ['Computing residual on finest level', 'Computing final error', 'Total application of A']
    cols_time += ['Evaluation of a^\{rr\}, a^\{rt\}, and a^\{tt\}', 'Evaluation of alpha and beta', 'Computing determinant of Jacobian of inverse mapping',
                  'Computing exact solution']
    cols_time += ['Total execution time']

    cols = cols_problem + cols_cluster + cols_convergence + cols_time + ['Benchmark']
    cols_dtypes = {col: 'int' for col in cols_problem + cols_cluster}
    cols_dtypes.update(
        {col: 'double' for col in cols_convergence + cols_time + ['Benchmark']})

    filename_input = 'p' + str(problem) + '-r' + str(nr_exp) + '-dbt' + str(divideBy2) + '-mpk' + str(mod_pk) + '-s' + str(
        smoother) + '-e' + str(extrapolation) + '--N' + str(nodes) + '-R' + str(ranks) + '-maxC' + str(maxCores)

    err = 0
    multi_region_run = True
    regions = []
    search_terms_likwid = {}
    search_terms_likwid['FLOPS_DP'] = 'DP [MFLOP/s]'
    search_terms_likwid['MEM_DP'] = 'Memory bandwidth [MBytes/s]'    
    benchmarks = []

    try: # allow deletion of empty error files
        with open(os.path.join(path_to_perf_files, file_prefix + '-' + filename_input + '-' + file_postfix + '.err')) as f:
            line = f.readline()
            while line and err == 0:
                line = f.readline()
                if ' error' in line:
                    err = 1
                    print('Error in job script.\n')
                    print(line)
    except FileNotFoundError:
        err = 0 # no error file = no error

    if err == 0:
        with open(os.path.join(path_to_perf_files, file_prefix + '-' + filename_input + '-' + file_postfix + '.out')) as f:

            lines = f.readlines()
            i = 0
            while i < len(lines)-1:

                # search for next program execution
                if '--> Iteration' in lines[i]:
                    while '--> Iteration' in lines[i] and '--> Iteration' in lines[i+1]:
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')
                    # iteration scheme ended:
                    # split on empty spaces and take number of iterations in third place
                    its = int(lines[i].split()[2][:-1])
                    while '-norm' not in lines[i]:
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')
                    norm2 = -1
                    if '2-norm' in lines[i]:
                        norm2 = float(lines[i].split()[-1])
                    else:
                        sys.exit('Error. 2-Norm information not available.')
                    norminf = -1
                    if 'inf-norm' in lines[i+1]:
                        norminf = float(lines[i+1].split()[-1])
                    else:
                        sys.exit('Error. Inf-Norm information not available.')

                    while 'Total setup' not in lines[i]:
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding timing output.')
                    # global time print out, split line by line of information and add to dataframe
                    time_dict = dict()
                    while 'Total execution time' not in lines[i-1]:
                        timings = lines[i].split(':') # here we split the timers based on empty spaces
                        # convert timing to floating point
                        if len(timings) > 1:
                            time_dict[timings[0].replace('\t', '')] = float(timings[1])
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')

                    while multi_region_run: # iterate over potentially multiple LIKWID regions

                        x=15 # TODO: check if all regions treated then set false
                    
                        while 'Group 1: ' not in lines[i]: # iterate until LIKWID output
                            i = i+1
                            if i > len(lines)-1:
                                sys.exit(
                                    'End of file reached without finding output.')
                        benchmark_line = lines[i].split()
                        if 'Region' in lines[i]:
                            regions.append(benchmark_line[1].replace(',', ''))
                        else:   
                            multi_region_run = False
                        benchmark = benchmark_line[-1]

                        # create save dataframe for previous benchmark and create new one
                        if benchmark not in benchmarks:
                            x=15 #TODO
                            if len(benchmarks) > 0:
                                perf_results_df.to_csv(
                                    os.path.join(path_to_perf_files, 
                                                 file_prefix + '-' + filename_input + '-' + file_postfix + '_' + benchmarks[-1] + '.csv'),
                                    index=False)

                            benchmarks.append(benchmark)
                            # create empty data frame from column set
                            perf_results_df = pd.DataFrame(columns=cols)
                            # set data types for new data frame as defined above
                            for k, v in cols_dtypes.items():
                                perf_results_df[k] = perf_results_df[k].astype(v)

                        # number of LIKWID outputs gives number of used cores
                        cores = str.count(lines[i+2], 'HWThread')

                        cores_used.append(cores)

                        search_postfix = ''
                        if cores > 1:
                            search_postfix = ' STAT'

                        while all([val + search_postfix not in lines[i] for val in search_terms_likwid.values()]):
                            i = i+1
                            if i > len(lines)-1:
                                sys.exit(
                                    'End of file reached without finding output.')

                        benchmark_result = float(lines[i].split('|')[2].split()[0])

                        # with the number of used cores determined, check if
                        # line for this result is already present
                        row_setting_avail = (
                            perf_results_df[cols_problem[0]] == problem) & (
                            perf_results_df[cols_problem[1]] == nr_exp) & (
                            perf_results_df[cols_problem[2]] == divideBy2) & (
                            perf_results_df[cols_problem[3]] == mod_pk) & (
                            perf_results_df[cols_cluster[0]] == nodes) & (
                            perf_results_df[cols_cluster[1]] == ranks) & (
                            perf_results_df[cols_cluster[2]] == cores)

                        if row_setting_avail.values.sum() <= 1:  # add new line
                            # add to all rows fix identifiers for:
                            #   problem, exponent in r, modified coordinates, extrapolation
                            #   number of nodes and ranks of the job
                            # as well as variable numbers for other columns
                            perf_results_df.loc[len(perf_results_df.index),
                                                cols_problem + cols_cluster +
                                                cols_convergence + 'Benchmark'] = [
                                problem, nr_exp, mod_pk, extrapolation,
                                nodes, ranks, cores, its, norm2, norminf,
                                benchmark_result]
                            # add timings
                            for time_col, time_val in time_dict.items():
                                perf_results_df.loc[len(
                                    perf_results_df.index)-1, time_col] = time_val
                        elif row_setting_avail.values.sum() > 1:  # error
                            sys.exit(
                                'Error. More than two lines corresponds to criterion.'
                                ' We can only have one line for each benchmark: FLOPS_DP and MEM_DP.')

                i += 1
            # end while

        # close file


if __name__ == "__main__":
    main()

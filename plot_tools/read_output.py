import pandas as pd
import numpy as np
import os
import sys


def main():

    jobid = '12425'
    problem = 5
    nr_exp = 4
    mod_pk = 1
    smoother = 3
    extrapolation = 1

    nodes = 1
    ranks = 1
    maxCores = 14
    cores_used = []

    path_to_perf_files = 'plot_tools/'

    cols_problem = ['Problem', 'rExp', 'ModPK', 'Extrapolation']

    cols_cluster = ['Nodes', 'Ranks', 'Cores']

    cols_convergence = ['its', '2-norm of error', 'inf-norm of error']

    cols_time = []
    cols_time += ['t_setup', 't_build', 't_facto_Ac',
                  't_build_P', 't_build_Asc', 't_facto_Asc']
    cols_time += ['t_total(fine)', 't_smoothing', 't_residual',
                  't_restriction', 't_Ac', 't_prolongation', 't_fine_residual',
                  't_error']
    cols_time += ['t_applyA']
    cols_time += ['t_coeff', 't_arr_art_att', 't_sol', 't_detDFinv', 't_trafo']
    cols_time += ['Total_execution_time']

    cols_benchmark = ['FLOPS_DP', 'CACHES']

    cols = cols_problem + cols_cluster + cols_convergence + cols_time + cols_benchmark
    cols_dtypes = {col: 'int' for col in cols_problem + cols_cluster}
    cols_dtypes.update(
        {col: 'double' for col in cols_convergence + cols_time + cols_benchmark})

    filename_input = 'p' + str(problem) + '-r' + str(nr_exp) + '-mpk' + str(mod_pk) + '-s' + str(
        smoother) + '-e' + str(extrapolation) + '--N' + str(nodes) + '-R' + str(ranks) + '-maxC' + str(maxCores)

    err = 0
    with open(os.path.join(path_to_perf_files, 'slurm-' + str(jobid) + '-' + filename_input + '.err')) as f:
        line = f.readline()
        while line and err == 0:
            line = f.readline()
            if ' error' in line:
                err = 1
                print('Error in job script.\n')
                print(line)

    if err == 0:

        with open(os.path.join(path_to_perf_files, 'slurm-' + str(jobid) + '-' + filename_input + '.out')) as f:

            # create empty data frame from column set
            perf_results_df = pd.DataFrame(columns=cols)
            # set data types for new data frame as defined above
            for k, v in cols_dtypes.items():
                perf_results_df[k] = perf_results_df[k].astype(v)

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

                    while 't_setup' not in lines[i]:
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')
                    # global time print out, split line by line of information and add to dataframe
                    time_dict = dict()
                    while 'Total_execution_time' not in lines[i-1]:
                        timings = lines[i].split()
                        for j in range(0, len(timings), 2):
                            # remove last char (':') from string and convert timing to floating point
                            time_dict[timings[j][:-1]
                                      ] = float(timings[j+1].replace(',', ''))
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')

                    # iterate until LIKWID printout
                    while 'Group 1: ' not in lines[i]:
                        i = i+1
                        if i > len(lines)-1:
                            sys.exit(
                                'End of file reached without finding output.')
                    benchmark = lines[i].split()[-1]
                    # number of LIKWID outputs gives number of used cores
                    cores = str.count(lines[i+2], 'HWThread')

                    cores_used.append(cores)

                    search_term_likwid = ''
                    if benchmark == 'FLOPS_DP':
                        search_term_likwid = 'DP [MFLOP/s]'
                    elif benchmark == 'CACHES':
                        search_term_likwid = 'Memory bandwidth [MBytes/s]'
                    else:
                        sys.exit('Error. Benchmark type not supported.')

                    if cores > 1:
                        search_term_likwid = search_term_likwid + ' STAT'

                    while search_term_likwid not in lines[i]:
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
                        perf_results_df[cols_problem[2]] == mod_pk) & (
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
                                            cols_convergence + [benchmark]] = [
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
                            ' We can only have one line for each benchmark: FLOPS_DP and CACHES.')

                i += 1
            # end while

            perf_results_df.to_csv(
                os.path.join(
                    path_to_perf_files, filename_input + '_benchmarks.csv'),
                index=False)
        # close file


if __name__ == "__main__":
    main()

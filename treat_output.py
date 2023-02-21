#read the txt.file with the terminal outputs and process them

file = 'output.txt'
open(file, 'r')

lines = []

with open(file) as f:
    lines = f.readlines()

f.close()

#create new file for latex code
latex = open('latex.txt', 'w')




count = 0
number = 0

for line in lines:
    count += 1

    problem = line.find('Problem size')
    if (problem != -1):
        number += 1
        print(f'{line[:-1]}')       #-1 to leave out the break at the end of the line

        latex.write("\\begin{table} \n \\begin{tabular}[c]{c|ccccc} \n ")
        latex.write("\multicolumn{5}{c}")


        if(number == 1):
            latex.write("{r0=1e-8, Circle geometry, Diriclet boundary condition}")
        elif(number == 2):
            latex.write("{r0=1e-8, Circle geometry, Across-the-origin discretization}")
        elif (number == 3):
            latex.write("{r0=1e-8, Deformed geometry, Diriclet boundary condition}")
        elif (number == 4):
            latex.write("{r0=1e-8, Deformed geometry, Across-the-origin discretization}")
        elif (number == 5):
            latex.write("{r0=1e-5, Circle geometry, Diriclet boundary condition}")
        elif (number == 6):
            latex.write("{r0=1e-5, Circle geometry, Across-the-origin discretization}")
        elif (number == 7):
            latex.write("{r0=1e-5, Deformed geometry, Diriclet boundary condition}")
        elif (number == 8):
            latex.write("{r0=1e-5, Deformed geometry, Across-the-origin discretization}")
        elif (number == 9):
            latex.write("{r0=1e-2, Circle geometry, Diriclet boundary condition}")
        elif (number == 10):
            latex.write("{r0=1e-2, Circle geometry, Across-the-origin discretization}")
        elif (number == 11):
            latex.write("{r0=1e-2, Deformed geometry, Diriclet boundary condition}")
        elif (number == 12):
            latex.write("{r0=1e-2, Deformed geometry, Across-the-origin discretization}")


        latex.write("\\\\ \n \cline{1-6} \n")
        latex.write(" $n_r\\times n_\\theta$  & its & $\widehat{\\rho}$ & $\|err\|_{\ell_2}$ & $\|err\|_{\infty}$ & time \\\\ \n")
        latex.write("\cline{1-6} \n")
        latex.write("  5$\\times$8  & ")


    convergence = line.find('Convergence after')
    if (convergence != -1):
        #print(f'{line}')
        #number_of_conv is found at: convergence+28, until the end
        conv = line[convergence+28:-1]   #take a substring from the line
        #print("conv = ", conv)
        latex.write(conv)
        latex.write(" & ")

    residual = line.find('rho =')
    if (residual != -1):
        #print(f'{line}')
        #mean residual reduction factor is found at: residual+38, until the end
        rho = line[residual+6:-1]
        #print("rho = ", rho)
        latex.write(rho)
        latex.write(" & ")

    norm_2 = line.find('2-norm of error =')
    if (norm_2 != -1):
        #print(f'{line}')
        #error is found at: norm_2+18, until the end
        error_2_norm = line[norm_2+18:-1]
        #print("error_2_norm = ", error_2_norm)
        latex.write(error_2_norm)
        latex.write(" & ")

    norm_inf = line.find('inf-norm of error =')
    if (norm_inf != -1):
        #print(f'{line}')
        #error is found at: norm_inf+20
        error_norm_inf = line[norm_inf+20:-1]
        #print("error_inf_norm = ", error_norm_inf)
        latex.write(error_norm_inf)
        latex.write(" & ")

    time = line.find('Total_execution_time:')
    if (time != -1):
        #print(f'{line}')
        #execution time is found at: time+22, until the end
        t = line[time+22:]
        #print("time = ", t)
        latex.write(t)
        latex.write(" \\\\ \n \cline{1-6} \n")

        latex.write("\end{tabular} \n \caption{\\textbf{Multigrid without extrapolation}. Iteration counts $\\textit{its}$, mean residual reduction factor $\widehat{\\rho}$, and errors of iterative solution to exact solution evaluated at the nodes in $\| \cdot \|_{\ell_2}$ and $\|\cdot\|_{\infty}$ norms.} \n")
        latex.write("\end{table}")


print("done!")











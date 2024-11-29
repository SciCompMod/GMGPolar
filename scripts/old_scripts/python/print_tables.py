import glob
import os
import h5py
import numpy as np

class Result:
    def __init__(self, n1, n2, levels, iterations, rho, resid, l2, linf):
        self.n1 = n1
        self.n2 = n2
        self.levels = levels
        self.iterations = iterations
        self.rho = rho
        self.residual = resid
        self.l2 = l2
        self.linf = linf

    def __lt__(self, other):
        return (self.n1, self.n2) < (other.n1, other.n2)

Shafranov  = 1
Triangular = 2

Polar = 6
Cartesian = 7

Uniform = 0
NonUniform = 1

def get_results(folder):
    files = glob.glob(folder+'/job.out*')
    results = {}
    for geom in (Shafranov, Triangular):
        results[geom] = {}
        for prob in (Polar, Cartesian):
            results[geom][prob] = {}
            for uni in (Uniform, NonUniform):
                results[geom][prob][uni] = {}
                for extrapol in (0, 1):
                    results[geom][prob][uni][extrapol] = []

    for f in files:
        with open(f) as fl:
            lines = fl.readlines()
        larray = [l.strip('\n').split() for l in lines]
        larray = [l for l in larray if l]

        l2 = [l[4] for l in larray if l[0] == "2-norm"]
        linf = [l[4] for l in larray if l[0] == "inf-norm"]
        geom = [l[1] for l in larray if l[0] == "mod_pk:"]
        extrapol = [l[1] for l in larray if l[0] == "extrapolation:"]
        prob = [l[1] for l in larray if l[0] == "prob:"]
        aniso = [l[1] for l in larray if l[0] == "fac_ani:"]
        iterations = [l[3] for l in larray if l[0] == "Convergence"]
        rho = [l[6] for l in larray if len(l) >=5 and l[4] == "rho"]
        resid = [l[6].strip(',') for l in larray if len(l) >=2 and l[1] == "Iteration"]
        prob_info = [l for l in larray if l[0] == "Prob:"]
        assert len(l2)==1
        assert len(linf)==1
        assert len(geom)==1
        assert len(extrapol)==1
        assert len(prob)==1
        assert len(aniso)==1
        assert len(rho)==1
        assert len(iterations)<=1
        assert len(prob_info)==1

        l2 = float(l2[0])
        linf = float(linf[0])
        geom = int(geom[0])
        extrapol = int(extrapol[0])
        prob = int(prob[0])
        aniso = int(aniso[0])
        rho = float(rho[0])
        iterations = int(iterations[0]) if iterations else None
        prob_info = prob_info[0]
        levels = int(prob_info[12])
        resid = float(resid[-1])
        n1 = int(prob_info[10][1:-1])
        n2 = int(prob_info[11][:-2])
        res = Result(n1, n2, levels, iterations, rho, resid, l2, linf)
        results[geom][prob][aniso][extrapol].append(res)

    return results

results = get_results('outputs')

for geom in (Shafranov, Triangular):
    for prob in (Polar, Cartesian):
        for uni in (Uniform, NonUniform):
            for extrapol in (True, False):
                results[geom][prob][uni][extrapol].sort()

configs = [(Triangular, Cartesian, Uniform),
           (Triangular, Polar, Uniform),
           (Triangular, Polar, NonUniform),
           (Shafranov, Cartesian, Uniform),
           (Shafranov, Polar, Uniform),
           (Shafranov, Polar, NonUniform)]

letter = ord('A')
for c in configs:
    print(r"""\begin{table}[!ht]
    \centering
    % \footnotesize
    \setlength{\tabcolsep}{3pt}""")
    print(r"\caption{ Test Case ", chr(letter), "}")
    print(r"""\label{tab:par_results}
    \begin{tabular}{cccccccccc} 
        \toprule
        $\bf n_r \times n_{\theta}$ & $\bf m$ & $\bf levels$ & $\bf iteration$s & $\bf \Hat{\rho}$ & $\bf residual$ & $\bf \norm{err}_{l_2}$ & $\bf ord$ & $\bf \norm{err}_{\infty}$ & $\bf ord$\\""")
    for e in (0,1):
        print("\\toprule")
        if e == 0:
            print(r"\multicolumn{10}{c}{\bf No extrapolation}\\")
        else:
            print(r"\multicolumn{10}{c}{\bf Implicit extrapolation}\\")
        print("\\midrule")
        table_result = results[c[0]][c[1]][c[2]][e]

        n_results = len(table_result)
        
        l2 = 0
        
        for r in table_result:
            new_m = r.n1*r.n2
            if l2 == 0:
                o2 = '-'
                oinf = '-'
            else:
                o2 = np.log(l2/r.l2) / np.log(np.sqrt(new_m/m))
                oinf = np.log(linf/r.linf) / np.log(np.sqrt(new_m/m))
            l2 = r.l2
            linf = r.linf
            m = new_m
            tmp='{:.2e}'.format(l2).split('e')
            l2_str = '${:.2f} \\cdot 10^{{{:d}}}$'.format(float(tmp[0]), int(tmp[1]))
            tmp='{:.2e}'.format(linf).split('e')
            linf_str = '${:.2f} \\cdot 10^{{{:d}}}$'.format(float(tmp[0]), int(tmp[1]))
            tmp='{:.2e}'.format(r.residual).split('e')
            resid_str = '${:.2f} \\cdot 10^{{{:d}}}$'.format(float(tmp[0]), int(tmp[1]))
            iter_str = '{:,}'.format(r.iterations) if r.iterations else 'N/A'
            cells = ['${:,} \\times {:,}$'.format(r.n1, r.n2),
                    '{:,}'.format(new_m),
                    '{:,}'.format(r.levels),
                    iter_str,
                    '{:.2f}'.format(r.rho),
                    resid_str,
                    l2_str,
                    o2 if isinstance(o2,str) else '{:.2f}'.format(o2),
                    linf_str,
                    oinf if isinstance(oinf,str) else '{:.2f}'.format(oinf)]
            cells = [s.replace(',','\,') for s in cells]
            print("&".join("{:<25}".format(s) for s in cells)+'\\\\')
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table} ")
    letter += 1

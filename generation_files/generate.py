import os
import shutil
import subprocess

culham_folder = '/home/emily/Code/culham-metric'
current_folder = os.path.dirname(__file__)

python_folder = os.path.join(culham_folder, 'poisson_code', 'tests', 'python_code')

qn_functions_file = os.path.join(python_folder, 'qn_functions.py')

qn_gen_file = os.path.join(python_folder, 'main_qn_analytical.py')

wrap_to_c = os.path.join(current_folder, 'wrap_c_to_class.py')

out_folder = os.path.join(current_folder, 'out')

src_folder = os.path.join(current_folder, '..', 'src')

include_folder = os.path.join(current_folder, '..', 'include')

PYTHON = shutil.which('python3')

problems = {
        'CartesianR2' : '( 1 - s**2 ) * sp.cos( 2*sp.pi*x ) * sp.sin( 2*sp.pi*y )',
        'PolarR6' : '1e-4 * s**6*(s-1)**6/(0.5**12) * sp.cos( 11 * t )',
        'CartesianR6' : '1e-4 * (s+1)**6 * (s-1)**6/(0.5**12) * sp.cos( 2*sp.pi*x ) * sp.sin( 2*sp.pi*y )',
        }

coeffs = {
        'Sonnendrucker' : '2/(2.6 + 3.14)*(1.3 + sp.atan( (1-1.3*s)/0.09 ))',
        'Zoni' : 'sp.exp( - sp.tanh( ( s - 0.5 ) / 0.1 ) )',
        'ZoniShifted' : 'sp.exp( - sp.tanh( ( s - 0.7 ) / 0.05 ) )'
        }

geometry = {
        'Circular': 'circle',
        'Shafranov': 'shafranov',
        'Triangular': 'triangularity'
        }

for c, coeff_a in coeffs.items():
    for beta in (True, False):
        coeff_b = '1/({})'.format(coeff_a) if beta else '0.0'
        for p, phi in problems.items():
            for g, geom in geometry.items():
                classname = p+('Gyro' if beta else '')+c+g
                print("Generating : ",classname)
                with open(qn_functions_file, "w") as f:
                    print("import sympy as sp", file=f)
                    print("", file=f)
                    print("def phi_exact(s,t,x,y):", file=f)
                    print(f"    return {phi}", file=f)
                    print("", file=f)
                    print("def coeffs1(s,t):", file=f)
                    print(f"    return {coeff_a}", file=f)
                    print("", file=f)
                    print("def coeffs2(s,t):", file=f)
                    print(f"    return {coeff_b}", file=f)
                cmd = [PYTHON, qn_gen_file, 'cpp', '--mapping', geom, '--output', os.path.join(out_folder, classname), '--vectorise', '--Rmax']
                exc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
                out, err = exc.communicate()
                assert(exc.returncode==0)

                shutil.move(os.path.join(out_folder, classname+'.h'), os.path.join(include_folder, classname+'.h'))
                shutil.move(os.path.join(out_folder, classname+'.cpp'), os.path.join(src_folder, classname+'.cpp'))

                shutil.rmtree(out_folder)

c = 'Poisson'
coeff_a = '0.0'
coeff_b = '0.0'
for p, phi in problems.items():
    for g, geom in geometry.items():
        classname = p+c+g
        print("Generating : ",classname)
        with open(qn_functions_file, "w") as f:
            print("import sympy as sp", file=f)
            print("", file=f)
            print("def phi_exact(s,t,x,y):", file=f)
            print(f"    return {phi}", file=f)
            print("", file=f)
            print("def coeffs1(s,t):", file=f)
            print(f"    return {coeff_a}", file=f)
            print("", file=f)
            print("def coeffs2(s,t):", file=f)
            print(f"    return {coeff_b}", file=f)
        cmd = [PYTHON, qn_gen_file, 'cpp', '--mapping', geom, '--output', os.path.join(out_folder, classname), '--vectorise', '--Rmax']
        exc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        out, err = exc.communicate()
        assert(exc.returncode==0)

        shutil.move(os.path.join(out_folder, classname+'.h'), os.path.join(include_folder, classname+'.h'))
        shutil.move(os.path.join(out_folder, classname+'.cpp'), os.path.join(src_folder, classname+'.cpp'))

        shutil.rmtree('out')


from argparse import ArgumentParser
import numpy as np
import sys

parser = ArgumentParser(description='Tool for generating points')
parser.add_argument('Nc', type=int, help='Number of cells (must be a multiple of 2)')
parser.add_argument('outfile', type=str, help='File where points should be printed')
parser.add_argument('--non_uniform', action = 'store_true', help='Cutoff = pi/4')
args = parser.parse_args()

N = args.Nc
has_cutoff = args.non_uniform

if has_cutoff:
    cutoff = np.pi/4
else:
    cutoff = 0

percent_in = cutoff/np.pi

n_in = round(2*percent_in*N)
n_out = round((1-percent_in)*N)

q = [*np.linspace(0, cutoff, n_in//2, endpoint=False),
        *np.linspace(cutoff, 2*np.pi-cutoff, n_out, endpoint=False),
        *np.linspace(2*np.pi-cutoff, 2*np.pi, n_in//2, endpoint=False)]

with open('t.txt','w') as f:
    for p in q:
        print(p, file=f)

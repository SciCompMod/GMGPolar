from argparse import ArgumentParser
import numpy as np

parser = ArgumentParser(description='Tool for generating points')
parser.add_argument('Nc', type=int, help='Number of cells (must be a multiple of 2)')
parser.add_argument('outfile', type=str, help='File where points should be printed')
parser.add_argument('--R0', type=float, default=1e-5, help='Minimum radius')
parser.add_argument('--R', type=float, default=1.0, help='Maximum radius')
args = parser.parse_args()


R0 = args.R0
R = args.R
outfile = args.outfile
N = args.Nc

pts = np.linspace(R0, R, N+1)

with open(outfile,'w') as f:
    for p in pts:
        print(p, file=f)

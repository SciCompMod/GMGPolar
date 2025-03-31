from argparse import ArgumentParser
import numpy as np

parser = ArgumentParser(description='Tool for generating points')
parser.add_argument('Nc', type=int, help='Number of cells (must be a multiple of 2)')
parser.add_argument('outfile', type=str, help='File where points should be printed')
args = parser.parse_args()


outfile = args.outfile
N = args.Nc

pts = np.linspace(0.0, 2*np.pi, N, endpoint=False)

with open(outfile,'w') as f:
    for p in pts:
        print(p, file=f)


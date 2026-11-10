# mesh_study.py
# Usage (in the case directory, next to input_bump.txt and geom_bump.txt):
#   python mesh_study.py generate     -> writes inputs + run_mesh.sh
#   bash run_mesh.sh                  -> runs the solver for every mesh
#   python mesh_study.py plot         -> reads the convergence logs and plots

import sys, shutil
import numpy as np
import matplotlib.pyplot as plt
from routines import read_conv

BASE = 'bump'
NI_LIST = [13, 25, 49, 97, 193, 385, 769]   # 12*2^k + 1  ->  12, 24, ... cells
NSTEPS = 5000                              
# Name of the numeric convergence file the solver writes. CHECK plot_conv.py
# for the real pattern; {name} is replaced by the casename.
CONV_PATTERN = 'conv_{name}.csv'

def name_of(ni):
    return f'{BASE}_ni{ni}'

def generate():
    with open(f'input_{BASE}.txt') as f:
        lines = f.read().splitlines()
    ni0, nj0 = [int(x) for x in lines[4].split()]
    cmds = []
    for ni in NI_LIST:
        # keep cell aspect ratio: scale nj with ni, rounded to 4k+1
        nj = int(round((nj0 - 1) * (ni - 1) / (ni0 - 1) / 4)) * 4 + 1
        name = name_of(ni)
        new = lines.copy()
        new[0] = name                    # casename
        new[3] = str(NSTEPS)             # nsteps
        new[4] = f'{ni} {nj}'            # ni nj
        with open(f'input_{name}.txt', 'w') as f:
            f.write('\n'.join(new) + '\n')
        shutil.copy(f'geom_{BASE}.txt', f'geom_{name}.txt')
        cmds.append(f'../Code/solver.x < input_{name}.txt > screen_{name}.txt &')
    with open('run_mesh.sh', 'w', newline='\n') as f:
        f.write('#!/bin/bash\n' + '\n'.join(cmds) + '\nwait\n')
    print('Written inputs and run_mesh.sh')

def plot():
    cells, err = [], []
    for ni in NI_LIST:
        l = read_conv(CONV_PATTERN.format(name=name_of(ni)))
        cells.append(ni - 1)
        err.append(l['dro_avg'][-1])
        print(f'ni={ni:4d}  steps={int(l["nstep"][-1]):6d}  dro_avg={err[-1]:.3e}')
    cells = np.array(cells, float); err = np.array(err)

    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    x = np.array([cells[0] * 0.8, cells[-1] * 1.4])
    for p, lab in zip([1, 2, 3], ['1st order', '2nd order', '3rd order']):
        ax.loglog(x, err[0] * (x / cells[0]) ** (-p), label=lab)
    ax.loglog(cells, err, 'o-', ms=3, color='C0', label='Bump test case')

    ax.set_xlabel('Number i-direction cells')
    ax.set_ylabel('Average residual density error')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()
    fig.savefig('mesh_convergence.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Observed order between successive refinements (slope on log-log axes)
    order = -np.diff(np.log(err)) / np.diff(np.log(cells))
    print('Observed order between meshes:', np.round(order, 2))

if __name__ == '__main__':
    {'generate': generate, 'plot': plot}[sys.argv[1]]()
# sfac_study.py
# Usage (in the Code directory, next to routines.py, input_<BASE>.txt and
# geom_<BASE>.txt):
#   python3 sfac_study.py generate   -> writes input/geom files for every SFAC
#   python3 sfac_study.py run        -> runs each case ONE AT A TIME,
#                                       writes screen_*.txt
#   python  sfac_study.py plot       -> stagnation pressure ratio vs x

import sys, glob, shutil, subprocess
import numpy as np
import matplotlib.pyplot as plt
import routines as rt

# ---------------- CONFIG ----------------
BASE = 'bump'                           # base case name (input_bump.txt ...)
SFAC_LIST = [0.06, 0.1, 0.2, 0.3, 0.4]
SOLVER = '../Code/solver.x'             # on Windows: '../Code/solver.exe'
NSTEPS = None                           # e.g. 20000, or None to keep the base value
NSTEPS_LINE = 3                         # 0-based line of nsteps in the input file

# Input file line 2 is "cfl sfac d_max" (see read_settings in routines.py)
SFAC_LINE = 2
SFAC_TOKEN = 1

# Solver output read with routines.read_case. {name} = casename.
# If this file name is wrong the script searches for "*final*<name>*".
FINAL_PATTERN = 'out_final_{name}.bin'

# What to plot against x:
#   'mass'  -> mass-averaged stagnation pressure on every i = const line
#   'jline' -> stagnation pressure along the grid line j = J_INDEX (0-based;
#              0 = wall, -1 = opposite boundary, nj//2 handled by 'mid')
PROFILE = 'mass'
J_INDEX = 'mid'
XSCALE = 1000.0                         # x in metres -> mm
# ----------------------------------------

def name_of(sfac):
    return f'{BASE}_sfac{sfac:g}'

def generate():
    with open(f'input_{BASE}.txt') as f:
        lines = f.read().splitlines()
    print(f'Using line {SFAC_LINE}, token {SFAC_TOKEN}:',
          repr(lines[SFAC_LINE]))
    for sfac in SFAC_LIST:
        name = name_of(sfac)
        new = lines.copy()
        new[0] = name                                    # casename
        if NSTEPS is not None:
            new[NSTEPS_LINE] = str(NSTEPS)
        tok = new[SFAC_LINE].split()
        tok[SFAC_TOKEN] = f'{sfac:g}'
        new[SFAC_LINE] = ' '.join(tok)
        with open(f'input_{name}.txt', 'w') as f:
            f.write('\n'.join(new) + '\n')
        shutil.copy(f'geom_{BASE}.txt', f'geom_{name}.txt')
    print('Written input and geom files for SFAC =', SFAC_LIST)

def run():
    for sfac in SFAC_LIST:
        name = name_of(sfac)
        print(f'Running SFAC = {sfac:g} ...', flush=True)
        with open(f'input_{name}.txt') as fin, \
             open(f'screen_{name}.txt', 'w') as fout:
            subprocess.run([SOLVER], stdin=fin, stdout=fout, check=True)

def find_final(name):
    fn = FINAL_PATTERN.format(name=name)
    if glob.glob(fn):
        return fn
    hits = glob.glob(f'*final*{name}*')
    if len(hits) == 1:
        print(f'{fn} not found, using {hits[0]}')
        return hits[0]
    sys.exit(f'Cannot find the final solution file for {name}. '
             f'Tried {fn}; matches for *final*{name}*: {hits}. '
             'Set FINAL_PATTERN at the top of the script.')

def profile(name):
    av = rt.read_settings(f'input_{name}.txt')
    g = rt.read_case(find_final(name))
    b = rt.calc_secondary(av, g)
    ni, nj = b['ni'], b['nj']
    if PROFILE == 'mass':
        x = np.zeros(ni); p0 = np.zeros(ni)
        for i in range(ni):
            c = rt.cut_i(b, i)
            p0[i], _ = rt.mass_av(c, 'pstag')
            x[i] = b['x'][i, 0]
    else:
        j = nj // 2 if J_INDEX == 'mid' else J_INDEX
        x = b['x'][:, j]; p0 = b['pstag'][:, j]
    return x * XSCALE, p0 / av['pstag']          # ratio to inlet stagnation pressure

def plot():
    fig, ax = plt.subplots(figsize=(10.8, 5.0))
    for sfac in SFAC_LIST:
        x, p0 = profile(name_of(sfac))
        ax.plot(x, p0, label=f'SFAC = {sfac:g}')
    ax.set_xlabel('x / mm')
    ax.set_ylabel('Stagnation pressure ratio')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig('sfac_stagnation_pressure.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    {'generate': generate, 'run': run, 'plot': plot}[sys.argv[1]]()
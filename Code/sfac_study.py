# sfac_study.py
# Usage (in the case directory, next to input_<BASE>.txt and geom_<BASE>.txt):
#   python3 sfac_study.py generate   -> writes input/geom files for every SFAC
#   python3 sfac_study.py run        -> runs each case ONE AT A TIME,
#                                       writes screen_*.txt
#   python  sfac_study.py plot       -> stagnation pressure ratio vs x

import sys, shutil, subprocess
import numpy as np
import matplotlib.pyplot as plt

# ---------------- CONFIG: check these against your own files ----------------
BASE = 'bump'                           # base case name (input_bump.txt ...)
SFAC_LIST = [0.06, 0.1, 0.2, 0.3, 0.4]
SOLVER = '../Code/solver.x'             # on Windows: '../Code/solver.exe'
NSTEPS = None                           # e.g. 20000, or None to keep the base value
NSTEPS_LINE = 3                         # 0-based line of nsteps in the input file

# 0-based line of the smoothing factor in input_<BASE>.txt. If the line holds
# several numbers, SFAC_TOKEN picks which one to replace.
# None = auto-detect: first line mentioning "sfac" (case-insensitive).
# If that fails, generate prints the numbered input file so you can set it here.
SFAC_LINE = None
SFAC_TOKEN = 0

# Solver output holding the profile. {name} = casename.
# Columns are read by index; the ratio is normalised by the first (inlet) point.
PROFILE_PATTERN = 'profile_{name}.txt'  # <-- CHECK
SKIPROWS = 0
X_COL, P0_COL = 0, 1                    # x and stagnation pressure columns
XSCALE = 1.0                            # 1000 if x is in metres, 1 if already mm
# -----------------------------------------------------------------------------

def name_of(sfac):
    return f'{BASE}_sfac{sfac:g}'

def generate():
    with open(f'input_{BASE}.txt') as f:
        lines = f.read().splitlines()
    global SFAC_LINE
    if SFAC_LINE is None:
        hits = [i for i, l in enumerate(lines) if 'sfac' in l.lower()]
        if hits:
            SFAC_LINE = hits[0]
    if SFAC_LINE is None or not (0 <= SFAC_LINE < len(lines)):
        print(f'Could not locate SFAC in input_{BASE}.txt. File contents '
              '(0-based line numbers):')
        for i, l in enumerate(lines):
            print(f'{i:3d}: {l}')
        sys.exit('Set SFAC_LINE (and SFAC_TOKEN if the line has several '
                 'numbers) at the top of the script.')
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

def plot():
    fig, ax = plt.subplots(figsize=(10.8, 5.0))
    for sfac in SFAC_LIST:
        d = np.loadtxt(PROFILE_PATTERN.format(name=name_of(sfac)),
                       skiprows=SKIPROWS, ndmin=2)
        x = d[:, X_COL] * XSCALE
        p0 = d[:, P0_COL] / d[0, P0_COL]
        ax.plot(x, p0, label=f'SFAC = {sfac:g}')
    ax.set_xlabel('x / mm')
    ax.set_ylabel('Stagnation pressure ratio')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig('sfac_stagnation_pressure.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    {'generate': generate, 'run': run, 'plot': plot}[sys.argv[1]]()
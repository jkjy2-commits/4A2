# mesh_study.py
# Usage (in the case directory, next to input_bump.txt and geom_bump.txt):
#   python3 mesh_study.py generate   -> writes input/geom files for every mesh
#   nohup python3 mesh_study.py run > run.log 2>&1 &
#                                    -> runs each case ONE AT A TIME (Linux),
#                                       writes screen_*.txt and times.txt
#   python mesh_study.py plot        -> convergence, run time, error vs time

import sys, shutil, subprocess
import numpy as np
import matplotlib.pyplot as plt
from routines import read_conv

BASE = 'bump'
NI_LIST = [13, 25, 49, 97, 193, 385, 769]   # 12*2^k + 1  ->  12, 24, ... cells
NSTEPS = 5000
NSTEPS_OVERRIDE = {}                         # e.g. {385: 10000, 769: 20000}
SOLVER = '../Code/solver.x'
TIMES_FILE = 'times.txt'
# Name of the numeric convergence file the solver writes. CHECK plot_conv.py
# for the real pattern; {name} is replaced by the casename.
CONV_PATTERN = 'conv_{name}.csv'

def name_of(ni):
    return f'{BASE}_ni{ni}'

def generate():
    with open(f'input_{BASE}.txt') as f:
        lines = f.read().splitlines()
    ni0, nj0 = [int(x) for x in lines[4].split()]
    for ni in NI_LIST:
        # keep cell aspect ratio: scale nj with ni, rounded to 4k+1
        nj = int(round((nj0 - 1) * (ni - 1) / (ni0 - 1) / 4)) * 4 + 1
        name = name_of(ni)
        new = lines.copy()
        new[0] = name                                    # casename
        new[3] = str(NSTEPS_OVERRIDE.get(ni, NSTEPS))    # nsteps
        new[4] = f'{ni} {nj}'                            # ni nj
        with open(f'input_{name}.txt', 'w') as f:
            f.write('\n'.join(new) + '\n')
        shutil.copy(f'geom_{BASE}.txt', f'geom_{name}.txt')
    print('Written input and geom files for', NI_LIST)

def run():
    # Linux only: "resource" is not available on Windows
    import resource
    open(TIMES_FILE, 'w').close()                # start a fresh times file
    for ni in NI_LIST:
        name = name_of(ni)
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with open(f'input_{name}.txt') as fin, \
             open(f'screen_{name}.txt', 'w') as fout:
            subprocess.run([SOLVER], stdin=fin, stdout=fout)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = ((after.ru_utime - before.ru_utime)
               + (after.ru_stime - before.ru_stime))
        print(f'ni={ni:4d}  cpu={cpu:.2f} s', flush=True)
        with open(TIMES_FILE, 'a') as f:
            f.write(f'{ni} {cpu:.3f}\n')

def plot():
    cells, err = [], []
    for ni in NI_LIST:
        l = read_conv(CONV_PATTERN.format(name=name_of(ni)))
        d = l['dro_avg']; n = len(d)
        late = d[int(0.8 * n):]                  # last 20 % of the run
        drift = (late.max() - late.min()) / late.mean()
        cells.append(ni - 1)
        err.append(d[-1])
        print(f'ni={ni:4d}  steps={int(l["nstep"][-1]):6d}  '
              f'dro_avg={d[-1]:.3e}  drift(last 20%)={drift:.1e}',
              '' if drift < 0.01 else '  <-- STILL CHANGING, add steps')
    cells = np.array(cells, float); err = np.array(err)

    # --- Mesh convergence ---
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    x = np.array([cells[0] * 0.8, cells[-1] * 1.4])
    for p, lab in zip([1, 2, 3], ['1st order', '2nd order', '3rd order']):
        ax.loglog(x, err[0] * (x / cells[0]) ** (-p), label=lab)
    ax.loglog(cells, err, 'o-', ms=3, color='C0', label='Bump test case')
    ax.set_xlabel('Number i-direction cells')
    ax.set_ylabel('Average residual density error')
    ax.grid(True, which='both', alpha=0.3); ax.legend()
    fig.savefig('mesh_convergence.png', dpi=300, bbox_inches='tight')

    order = -np.diff(np.log(err)) / np.diff(np.log(cells))
    print('Observed order between meshes:', np.round(order, 2))

    # --- Run time (only if times.txt exists) ---
    try:
        data = np.loadtxt(TIMES_FILE, ndmin=2)
    except OSError:
        print(f'{TIMES_FILE} not found, skipping run-time plots')
        plt.show(); return
    t_ni, t = data[:, 0], data[:, 1]
    t_cells = t_ni - 1

    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    ax.loglog(t_cells, t, 'o-', ms=3, color='C0')
    ax.set_xlabel('Number i-direction cells')
    ax.set_ylabel('Run time / seconds')
    ax.grid(True, which='both', alpha=0.3)
    fig.savefig('run_time.png', dpi=300, bbox_inches='tight')
    print('Run-time scaling exponent between meshes:',
          np.round(np.diff(np.log(t)) / np.diff(np.log(t_cells)), 2))

    # --- Error against run time (cost of accuracy) ---
    # Match each time to its error by ni, in case some runs are missing
    err_for_t = [err[list(NI_LIST).index(int(n))] for n in t_ni]
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    ax.loglog(t, err_for_t, 'o-', ms=3, color='C0')
    for n, tt, ee in zip(t_cells, t, err_for_t):
        ax.annotate(f'{int(n)}', (tt, ee), textcoords='offset points',
                    xytext=(5, 5), fontsize=10)
    ax.set_xlabel('Run time / seconds')
    ax.set_ylabel('Average residual density error')
    ax.grid(True, which='both', alpha=0.3)
    fig.savefig('error_vs_time.png', dpi=300, bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    {'generate': generate, 'run': run, 'plot': plot}[sys.argv[1]]()
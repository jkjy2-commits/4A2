# plot_residual.py
# Usage: python plot_residual.py casename
# Needs input_<casename>.txt and out_final_<casename>.bin in the current directory
from routines import *
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
RES_VMIN = 5e-6          # lowest value on the residual colour scale
RES_VMAX = None          # None = data maximum; fix it to compare cases
D2_VMIN = None           # None = 1e-3 x data maximum
D2_VMAX = None
ZOOM_MANUAL = None       # or (x0, x1, y0, y1) to override the automatic zoom box
ZOOM_HALF_WIDTH = 0.08   # automatic box half-width as a fraction of the x-range
FIGSIZE = (10, 4.5)
# ---------------------------------------------------------------------------

def second_deriv(f, x, y, axis):
    # Second derivative of f along an index direction, with the physical
    # (non-uniform) node spacing. Returns interior nodes along that axis.
    f = np.moveaxis(f, axis, 0); x = np.moveaxis(x, axis, 0)
    y = np.moveaxis(y, axis, 0)
    h = np.hypot(np.diff(x, axis=0), np.diff(y, axis=0))   # spacing, (n-1, ...)
    h1, h2 = h[:-1], h[1:]
    d2 = 2.0 * ((f[2:] - f[1:-1]) / h2 - (f[1:-1] - f[:-2]) / h1) / (h1 + h2)
    return np.moveaxis(d2, 0, axis)

def max_second_deriv(f, x, y):
    # Larger magnitude of the i- and j-direction second derivatives at each
    # interior node, padded by the nearest value so the field is full size
    d2i = second_deriv(f, x, y, 0)[:, 1:-1]       # (ni-2, nj-2)
    d2j = second_deriv(f, x, y, 1)[1:-1, :]       # (ni-2, nj-2)
    d2 = np.maximum(np.abs(d2i), np.abs(d2j))
    return np.pad(d2, 1, mode='edge')             # (ni, nj)

def add_inset(ax, x, y, field, norm, g, box):
    x0, x1, y0, y1 = box
    axins = ax.inset_axes([1.03, 0.0, 0.6, 1.0])
    axins.pcolormesh(x, y, field, norm=norm, shading='flat', cmap='viridis')
    plot_wall(axins, g, scale=False)
    axins.set_xlim(x0, x1); axins.set_ylim(y0, y1)
    axins.set_aspect('equal'); axins.axis('off')
    ax.indicate_inset_zoom(axins, edgecolor='black')

def main():
    casename = sys.argv[-1]
    av = read_settings('input_' + casename + '.txt')
    g = read_case('out_final_' + casename + '.bin')
    x, y = g['x'], g['y']                          # nodes, (ni, nj)

    # Reference density (inlet stagnation) and reference length (inlet width)
    rho_ref = av['pstag'] / (av['rgas'] * av['tstag'])
    l_ref = np.sum(np.hypot(g['lx_i'][0, :], g['ly_i'][0, :]))

    # ---------------- Figure 1: density ratio ----------------
    density_ratio = cell_av(g['ro']) / rho_ref     # (ni-1, nj-1)
    fig1, ax1 = plt.subplots(figsize=FIGSIZE)
    m1 = ax1.pcolormesh(x, y, density_ratio, shading='flat', cmap='viridis')
    ax1.set_aspect('equal'); ax1.axis('off'); plot_wall(ax1, g)
    cb = fig1.colorbar(m1, ax=ax1, orientation='horizontal', pad=0.05, fraction=0.05)
    cb.set_label(r'Density ratio $\rho/\rho_{01}$')
    fig1.savefig('density_' + casename + '.png', dpi=300, bbox_inches='tight')

    # ---------------- Figure 2: residual with zoom inset ----------------
    res = np.maximum(np.abs(g['dro']) / rho_ref, RES_VMIN)   # (ni-1, nj-1)
    res_norm = LogNorm(vmin=RES_VMIN,
                       vmax=res.max() if RES_VMAX is None else RES_VMAX)

    xc, yc = cell_av(x), cell_av(y)
    ip, jp = np.unravel_index(res.argmax(), res.shape)
    print('Peak residual %.3e at x = %.4f, y = %.4f'
          % (res[ip, jp], xc[ip, jp], yc[ip, jp]))
    if ZOOM_MANUAL is not None:
        box = ZOOM_MANUAL
    else:
        hw = ZOOM_HALF_WIDTH * (x.max() - x.min())
        box = (xc[ip, jp] - hw, xc[ip, jp] + hw,
               yc[ip, jp] - 0.5 * hw, yc[ip, jp] + 0.5 * hw)

    fig2, ax2 = plt.subplots(figsize=FIGSIZE)
    m2 = ax2.pcolormesh(x, y, res, norm=res_norm, shading='flat', cmap='viridis')
    ax2.set_aspect('equal'); ax2.axis('off'); plot_wall(ax2, g)
    cb = fig2.colorbar(m2, ax=ax2, orientation='horizontal', pad=0.05, fraction=0.05)
    cb.set_label(r'Absolute residual density error $|\Delta\rho|/\rho_{01}$')
    add_inset(ax2, x, y, res, res_norm, g, box)
    fig2.savefig('residual_' + casename + '.png', dpi=300, bbox_inches='tight')

    # ---------------- Figure 3: maximum second spatial derivative ----------------
    d2 = max_second_deriv(g['ro'], x, y) * l_ref**2 / rho_ref   # non-dimensional
    d2_max = d2.max()
    d2_min = 1e-3 * d2_max if D2_VMIN is None else D2_VMIN
    d2_norm = LogNorm(vmin=d2_min, vmax=d2_max if D2_VMAX is None else D2_VMAX)
    i2, j2 = np.unravel_index(d2.argmax(), d2.shape)
    print('Peak d2rho/ds2 (non-dim) %.3e at x = %.4f, y = %.4f'
          % (d2_max, x[i2, j2], y[i2, j2]))

    fig3, ax3 = plt.subplots(figsize=FIGSIZE)
    m3 = ax3.pcolormesh(x, y, np.maximum(d2, d2_min), norm=d2_norm,
                        shading='gouraud', cmap='viridis')
    ax3.set_aspect('equal'); ax3.axis('off'); plot_wall(ax3, g)
    cb = fig3.colorbar(m3, ax=ax3, orientation='horizontal', pad=0.05, fraction=0.05)
    cb.set_label(r'Max. second derivative of density, '
                 r'$\max\,|\partial^2\rho/\partial s^2|\,L^2/\rho_{01}$')
    fig3.savefig('d2rho_' + casename + '.png', dpi=300, bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    main()
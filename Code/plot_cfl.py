import numpy as np
import matplotlib.pyplot as plt
from routines import read_conv

cfls = [0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8]
err = []
for cfl in cfls:
    l = read_conv(f'conv_bump_cfl{cfl}.csv')   # replace with your actual log filename
    err.append(l['dro_avg'][-1])               # final average density residual
err = np.array(err); cfls = np.array(cfls)

fig, ax = plt.subplots(figsize=(9.6, 7.2))
ax.loglog(cfls, err, 'o-', ms=3, color='C0', label='Bump test case')

# Reference lines, anchored to the CFL = 0.1 data point
k = np.argmin(abs(cfls - 0.1)); x = np.array([0.009, 1.1])
ax.loglog(x, err[k] * (x / cfls[k])**0, label='0th order')
ax.loglog(x, err[k] * (x / cfls[k])**1, label='1st order')
ax.loglog(x, err[k] * (x / cfls[k])**2, label='2nd order')

ax.set_xlabel('CFL number'); ax.set_ylabel('Average residual density error')
ax.grid(True, which='both', alpha=0.3); ax.legend()
plt.show()
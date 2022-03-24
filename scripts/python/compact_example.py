#!/usr/bin/env python3
"""
Example script for compact library

Created  on & by: 2022/02/16 - J. Kostelecky (j.kostelecky@posteo.de)
Modified on & by: ...        - ...
"""
import numpy as np 
import matplotlib.pyplot as plt
import compact_lib as cl
#---------------------------------------------------------------------------#
# parameters
alpha_modified = 0.52; n = 10; nex = 100

# create grid - last point is dropped (periodicity)
x          = np.linspace(0, 2*np.pi, num=n,   endpoint=False)
x_exact    = np.linspace(0, 2*np.pi, num=nex, endpoint=False)

# functions - exact
u_exact    = np.sin(x_exact)
dudx_exact = np.cos(x_exact)

# functions - fdm
u          = np.sin(x)
dudx       = cl.compactder(u, x, scheme='6t',                       periodic=True)[:,0]
dudxp      = cl.compactder(u, x, scheme='6p', alpha=alpha_modified, periodic=True)[:,0]

# plot - functions
fig, ax1 = plt.subplots(figsize=(12,8), dpi =210)
ax1.set_title('FDM compact 6th-order (periodic + uniform)')
ax1.set_xlabel('x')
ax1.set_ylabel('functions')
ax1.set_xlim(0, 2*np.pi)
ax1.set_ylim(-1.5, 1.5)
ax1.plot(x_exact, u_exact,    label="u  = sin(x)",       color='g'            )
ax1.plot(x_exact, dudx_exact, label="u' = cos(x)",       color='b', alpha=0.7 )
ax1.scatter(x,    u,          label="u  fdm",            color='k', marker='x')
ax1.scatter(x,    dudxp,      label="u' fdm-pentadiag.", color='m', marker='o')
ax1.scatter(x,    dudx,       label="u' fdm-tridiag.",   color='r', marker='x')
ax1.legend()
ax1.grid(True)
plt.show()

# plot - error
fig, ax1 = plt.subplots(figsize=(12,8), dpi =210)
ax1.set_title('FDM compact 6th-order (periodic + uniform) -- absolute error')
ax1.set_xlabel('x')
ax1.set_ylabel('error')
ax1.set_xlim(0,2*np.pi)
ax1.set_ylim(-4e-5, 4e-5)
ax1.scatter(x, np.cos(x) - dudxp, label="error = cos(x) - u' fdm-pentadiag.", color='m', marker='o')
ax1.scatter(x, np.cos(x) - dudx,  label="error = cos(x) - u' fdm-tridiag.",   color='r', marker='x')
ax1.legend()
ax1.grid(True)
plt.show()
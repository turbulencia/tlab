import numpy as np
import matplotlib.pyplot as plt
from   scipy import interpolate
# import numpy.linalg as la
import os
import sys

# folder for figs
figs      = 'figs_splines_xyz'  # save figs here
if not os.path.exists(figs):
    os.makedirs(figs)

# =============================================================================
# in z
# =============================================================================
# 
xk = np.arange(1,12,1)
xk = np.delete(xk,[4,5,6])

yk = np.random.random(11)
yk = np.delete(yk,[4,5,6])
yk[3] = 0
yk[4] = 0

# spline interpolation
tck = interpolate.splrep(xk, yk, s=0, k=3)

# new
xks = np.arange(-1,11.01,0.1)
yks = interpolate.splev(xks, tck, der=0)

# test if possible for extrapolation
xks = np.arange(0,12.01,0.1)
yks = interpolate.splev(xks, tck, der=0)
#sys.exit()


# plot
plt.close('all')
plt.figure(figsize=(12,8), dpi =210)
plt.rcParams.update({'font.size':11})
plt.xlabel('$xk$')
plt.ylabel("$yk$")
plt.title('Spline interpolation in z')
plt.xlim(0,12)
plt.ylim(-0.5,1.5)
plt.xticks(np.arange(0, 12, step=1))
plt.grid(b=True, which='minor', color='grey', linestyle=':')
plt.grid(b=True, which='major', color='grey', linestyle=':')
#
plt.scatter(xk, yk, marker='s', c='b',label='points')
#
plt.plot(xks, yks, label='spline')
#
plt.legend(loc='best')
plt.show()
# sys.exit()
# output_name = figs+'/spline_test.svg'        
# plt.savefig(output_name)

# =============================================================================
# in y
# =============================================================================
# 
xj = np.arange(1,8,1)
xj = np.delete(xj,[1,2])

yj    = np.random.random(7)
yj    = np.delete(yj,[1,2])
yj[0] = 0
yj[1] = 0

# spline interpolation
tck = interpolate.splrep(xj, yj, s=0, k=3)

# new
xjs = np.arange(1,7.01,0.1)
yjs = interpolate.splev(xjs, tck, der=0)

# plot
plt.figure(figsize=(12,8), dpi =210)
plt.rcParams.update({'font.size':11})
plt.xlabel('$xj$')
plt.ylabel("$yj$")
plt.title('Spline interpolation in y')
plt.xlim(0,8)
plt.ylim(-0.5,1.5)
plt.xticks(np.arange(0, 6, step=1))
plt.grid(b=True, which='minor', color='grey', linestyle=':')
plt.grid(b=True, which='major', color='grey', linestyle=':')
#
plt.scatter(xj, yj, marker='s', c='b',label='points')
#
plt.plot(xjs, yjs, label='spline')
#
plt.legend(loc='best')
plt.show()
# sys.exit()
# output_name = figs+'/spline_test.svg'        
# plt.savefig(output_name)

# =============================================================================
# in x
# =============================================================================

# 
xi     = np.zeros((4))
xi[0]  = 1
xi[1]  = 2
xi[-2] = 9
xi[-1] = 10

yi = np.zeros((4))

# spline interpolation
tck = interpolate.splrep(xi, yi, s=0, k=3)

# new
xis = np.arange(1,10.01,0.1)
yis = interpolate.splev(xis, tck, der=0)

# plot
plt.figure(figsize=(12,8), dpi =210)
plt.rcParams.update({'font.size':11})
plt.xlabel('$xj$')
plt.ylabel("$yj$")
plt.title('Spline interpolation in x')
plt.xlim(0,11)
plt.ylim(-0.5,1.5)
plt.xticks(np.arange(0, 11, step=1))
plt.grid(b=True, which='minor', color='grey', linestyle=':')
plt.grid(b=True, which='major', color='grey', linestyle=':')
#
plt.scatter(xi, yi, marker='s', c='b',label='points')
#
plt.plot(xis, yis, label='spline')
#
plt.legend(loc='best')
plt.show()
# sys.exit()
# output_name = figs+'/spline_test.svg'        
# plt.savefig(output_name)

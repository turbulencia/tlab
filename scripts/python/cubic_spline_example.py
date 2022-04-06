#!/usr/bin/env python3
import numpy as np 
import matplotlib.pyplot as plt
from   scipy.interpolate import CubicSpline 
import cubic_spline_lib as csl
# =============================================================================
# Change settings here
# =============================================================================
# data size
nran  = 11  # orig data
nint  = 101 # int data

# validata with random data
random = True #False

# choose BCs here (also mixed BCs possible, not for periodic BCs!)
# bc = ['clamped'] *2; bc_val=[0,0] # first derivative at endpoints is zero
# bc = ['fixed1']  *2; bc_val=[4,4] # first  derivative at endpoints is on fixed value
bc = ['natural'] *2; bc_val=[0,0] # second derivative at endpoints is zero
# bc = ['fixed2']  *2; bc_val=[10,2] # second derivative at endpoints is on fixed value
# bc = ['periodic']*2; bc_val=[0,0]
# =============================================================================
# Do not change below here
# =============================================================================
# generate data
if random:
    x    = np.arange(nran) - np.random.random(nran)
    xint = np.linspace(x.min(), x.max(), nint)     
    y    = np.random.random(nran)
else:
    x    = np.linspace(0,1,nran)
    xint = np.linspace(x.min(), x.max(), nint)     
    y    = np.sin(2*np.pi/x[-1]*x)

# adjust data for scipy routine
if 'periodic'  in bc: y[-1] = y[0]; bc_sci = 'periodic'
elif 'clamped' in bc: bc_sci = bc_type=((1, 0.       ), (1, 0.       ))
elif 'fixed1'  in bc: bc_sci = bc_type=((1, bc_val[0]), (1, bc_val[1]))
elif 'natural' in bc: bc_sci = bc_type=((2, 0.       ), (2, 0.       ))
elif 'fixed2'  in bc: bc_sci = bc_type=((2, bc_val[0]), (2, bc_val[1]))

# spline interpolation
sp   = CubicSpline(x, y , bc_type=bc_sci)            # scipy impl.
yint = csl.cubic_spline(bc, x, y, xint, bc_val)     # owm   impl.
# =============================================================================
# First and second derivative to check
# =============================================================================
dyint   = np.round(np.gradient(yint,xint),2)
dyints  = np.round(np.gradient(sp(xint),xint),2)
ddyint  = np.round(np.gradient(dyint,xint),2)
ddyints = np.round(np.gradient(dyints,xint),2)
print('------------------------------------------------------------')
print('own: 1st deriv. at start = ', dyint[:4])
print('                at end   = ', dyint[-4:])
print('own: 2nd deriv. at start = ', ddyint[:4])
print('                at end   = ', ddyint[-4:])
print('------------------------------------------------------------')
print('sci: 1st deriv. at start = ', dyints[:4])
print('                at end   = ', dyints[-4:])
print('sci: 2nd deriv. at start = ', ddyints[:4])
print('                at end   = ', ddyints[-4:])
print('------------------------------------------------------------')

# error - own - scipy implementation
if sum((yint   - sp(xint))**2)**0.5 < 10e-10:
    print('Validation successfull!')
print('L2   error =', sum((yint   - sp(xint))**2)**0.5) 
# =============================================================================
# Plot
# =============================================================================
plt.close('all')
plt.figure(figsize=(12,8), dpi =210)
plt.title('Validation of cubic splines')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(x.min(), x.max())
plt.scatter(x,y, label='data points')
plt.plot(xint, sp(xint),label= 'scipy - truth')
plt.plot(xint, yint, label='own impl.')
plt.plot(xint, yint - sp(xint), label='error')
plt.legend()
plt.grid(True)
plt.show()
# =============================================================================
# Validating tlab implementation
# =============================================================================
# read .dat files
path = './'
a_int = csl.read_txt(path,'cspline_int.dat',header=0,footer=0)
a_org = csl.read_txt(path,'cspline_org.dat',header=0,footer=0)

# compute spline with BCs from tlab routine
yint  = csl.cubic_spline(['periodic','periodic'], a_org[:,0], a_org[:,1], a_int[:,0], [0,0]) 

# error
print('L2   error =', sum((yint - a_int[:,2])**2)**0.5) 

# plot
plt.figure(figsize=(12,8), dpi =210)
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(a_org[:,0].min(), a_org[:,0].max())
plt.scatter(a_org[:,0], a_org[:,1], label= 'data points')
plt.plot(a_int[:,0], yint, label= 'python own. impl')
plt.plot(a_int[:,0], a_int[:,2], label= 'tlab')
plt.plot(a_int[:,0], yint - a_int[:,2], label= 'error')
plt.legend()
plt.grid(True)
plt.show()
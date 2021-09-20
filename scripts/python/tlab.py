#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 12:35:54 2019

@author: jpmellado
"""

import netCDF4 as nc
import sys
import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt

if __name__ == "__main__":
    B0 = 0.005       # Surface buoyancy flux in simulation units.
    N0 = 3. **0.5    # Buoyancy frequency in the free atmosphere, in simulation units.
    L0 = (B0/N0**3.)**0.5
    #
    sdata = nc.Dataset(sys.argv[1], 'r')
    sdata.set_auto_mask(False)
    t  = sdata.variables['t'][:]
    y  = sdata.variables['y'][:]
    b  = sdata.variables['rS'][:,:]
    b_bg  = N0**2.0 *y
    z_enc = np.sqrt( np.trapz(b[:,:]-b_bg[None,:], y, axis=1) *2.0 /N0**2.0 )
    #
    print(z_enc/L0)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def find_nearest_bisection(array,value):
    """
    Returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.
    """

    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n
    jl = 0                      # Initialize lower
    ju = n-1                    # and upper limits.
    while (ju-jl > 1):          # If we are not yet done,
        jm=(ju+jl) >> 1         # compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm               # and replace either the lower limit
        else:
            ju=jm               # or the upper limit, as appropriate.
                                # Repeat until the test condition is satisfied.
    if (value == array[0]):     # edge cases at bottom
        return 0
    elif (value == array[n-1]): # and top
        return n-1
    else:
        return jl

def running_average(x,f,dx,m):
    """
    Calculate running average over dx at equally spaced m points.
    Assuming the independent variable x is non-decreasing.
    """

    # Defining the final independent variable
    yl = x[0]  +0.5 *dx
    jl = find_nearest_bisection(x,yl)
    yr = x[-1] -0.5 *dx
    jr = find_nearest_bisection(x,yr)
    y = np.linspace(yl,yr,m)

    if   np.ndim(f) == 1:
        avg = np.empty(m)
        std = np.empty(m)
    elif np.ndim(f) == 2:
        avg = np.empty((m,np.shape(f)[1]))
        std = np.empty((m,np.shape(f)[1]))

    for j in range(m):
        yl = y[j] -0.5 *dx
        jl = find_nearest_bisection(x,yl)
        yr = y[j] +0.5 *dx
        jr = find_nearest_bisection(x,yr)

        if   np.ndim(f) == 1:
            avg[j] = np.trapz(  f[jl:jr]            ,x[jl:jr]) /( x[jr-1]-x[jl] )
            std[j] = np.trapz( (f[jl:jr]-avg[j])**2.,x[jl:jr]) /( x[jr-1]-x[jl] )
        elif np.ndim(f) == 2:
            avg[j,:] = np.trapz(  f[jl:jr,:]              ,x[jl:jr], axis=0) /( x[jr-1]-x[jl] )
            std[j,:] = np.trapz( (f[jl:jr,:]-avg[j,:])**2.,x[jl:jr], axis=0) /( x[jr-1]-x[jl] )

    std = np.sqrt( std )

    return y, avg, std

def gradient(x,f,order):
    """
    Calculate 1st-order derivative using 6th order compact scheme
    """

    n = np.size(f)
    b = np.empty(n)

    # defining vector left-hand side
    A = np.ones((3,n))
    A[0,:] = 1. /3.
    A[2,:] = 1. /3.

    A[0, 1] = 2.
    A[0, 2] = 1. /2.
    A[0,-1] = 1. /6.

    A[2, 0] = 1. /6.
    A[2,-3] = 1. /2.
    A[2,-2] = 2.

    # Jacobian
    b[0] = -5./2. *x[0] +2.0 *x[1] +0.5 *x[2]
    b[1] = -5./9. *x[0] -0.5 *x[1] +     x[2] +1./18. *x[3]
    for i in range(2,n-2):
        b[i] = 7./9. *(x[i+1]-x[i-1]) +1./36. *(x[i+2]-x[i-2])
    b[n-2] =-1./18. *x[n-4] -     x[n-3] +0.5 *x[n-2] +5./9. *x[n-1]
    b[n-1] =                -0.5 *x[n-3] -2.0 *x[n-2] +5./2. *x[n-1]

    xp = solve_banded((1, 1), A, b)

    A = A *xp

    # defining right-hand side
    b[0] = -5./2. *f[0] +2.0 *f[1] +0.5 *f[2]
    b[1] = -5./9. *f[0] -0.5 *f[1] +     f[2] +1./18. *f[3]
    for i in range(2,n-2):
        b[i] = 7./9. *(f[i+1]-f[i-1]) +1./36. *(f[i+2]-f[i-2])
    b[n-2] =-1./18. *f[n-4] -     f[n-3] +0.5 *f[n-2] +5./9. *f[n-1]
    b[n-1] =                -0.5 *f[n-3] -2.0 *f[n-2] +5./2. *f[n-1]

    # solve system
    fp = solve_banded((1, 1), A, b)

    if order == 1:
        return fp
    else:
        """
        Calculate 2nd-order derivative using 6th order compact scheme
        """
        # defining vector left-hand side
        A = np.ones((3,n))
        A[0,:] = 2. /11.
        A[2,:] = 2. /11.

        A[0, 1] = 11.
        A[0, 2] = 1. /10.
        A[0,-1] = 1. /10.

        A[2, 0] = 1. /10.
        A[2,-3] = 1. /10.
        A[2,-2] = 11.

        # Jacobian
        b[0] = 13.  *x[0] -27.   *x[1] +15.   *x[2] -x[3]
        b[1] = 6./5.*x[0] -12./5.*x[1] +6./5. *x[2]
        for i in range(2,n-2):
            b[i] = 12./11. *(x[i+1]+x[i-1]) +3./44. *(x[i+2]+x[i-2]) -51./22. *x[i]
        b[n-2] =          6./5.*x[n-3] -12./5.*x[n-2] +6./5.*x[n-1]
        b[n-1] = -x[n-4] +15.  *x[n-3] -27.   *x[n-2] +13.  *x[n-1]

        xpp = solve_banded((1, 1), A, b)
        Ap= A *xpp

        A = A *xp **2.

        # defining right-hand side
        b[0] = 13.  *f[0] -27.   *f[1] +15.   *f[2] -f[3]
        b[1] = 6./5.*f[0] -12./5.*f[1] +6./5. *f[2]
        for i in range(2,n-2):
            b[i] = 12./11. *(f[i+1]+f[i-1]) +3./44. *(f[i+2]+f[i-2]) -51./22. *f[i]
        b[n-2] =          6./5.*f[n-3] -12./5.*f[n-2] +6./5.*f[n-1]
        b[n-1] = -f[n-4] +15.  *f[n-3] -27.   *f[n-2] +13.  *f[n-1]

        i=0; \
            b[i] -=                     Ap[1,i] *fp[i] +Ap[0,i+1] *fp[i+1]
        for i in range(1,n-1):
            b[i] -= Ap[2,i-1] *fp[i-1] +Ap[1,i] *fp[i] +Ap[0,i+1] *fp[i+1]
        i=n-1; \
            b[i] -= Ap[2,i-1] *fp[i-1] +Ap[1,i] *fp[i]

        fpp = solve_banded((1, 1), A, b)
        return fpp


"""
Testing
"""

# n = 100
# x = np.linspace(0.,1.,n)
# x = x **2. +x
# f = np.sin(x)
# df = gradient(x,f,order=2)
# # plt.plot(x,np.exp(x))
# # plt.plot(x,df,'o')
# plt.plot(x,df+np.sin(x),'o')
# plt.show()

#n = 100
#f = np.random.rand(n)
#x = np.linspace(0.,1.,n)
#
#plt.plot(x,f)
#plt.show()
#
#y, avg, std = running_average(x,f,0.1,20)
#plt.plot(y,avg)
#plt.fill_between(y,avg-std,avg+std,alpha=0.2)
#plt.show()

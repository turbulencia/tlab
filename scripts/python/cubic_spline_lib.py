#!/usr/bin/env python3
import numpy as np 
import bisect
import compact_lib as cl
"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
Library for evaluation of cubic splines.

Created  on & by: 2022/03/29 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: ...        - ...
#---------------------------------------------------------------------------#
Contains:
            - wrapper function for cubic splines with 
              different boundary conditions
#############################################################################
"""
def cubic_spline(bcs,x,y,xint,bcs_val=[0,0]):
    """
    Wrapper function for cubic splines with different boundary conditions.

    Parameters
    ----------
    bcs     : boundary conditions for start and end of spline
    x       : grid nodes of original data
    y       : function values of original data
    xint    : grid nodes of interpolated data
    bcs_val : boundary values for first or second derivative at 
              start and end of spline, default is [0,0].

    Returns
    -------
    yint    : interpolated function on xint.

    """
    # step size
    h = np.diff(x)
    # check input data
    check_input(bcs, x, y, xint, bcs_val)
    # diagonals of tridiagonal matrix 
    aa, bb, cc = tridiag(bcs, h)
    # compute rhs 
    rhs = rhs_mat(bcs, h, y, bcs_val)
    # solve tridiagonal matrix with thomas algorithm (m = 2nd deriv.)
    if 'periodic' in bcs:
        # solve circulant tridiagonal matrix
        n = h.size; dd = np.zeros(n+1); ee = np.zeros(n+1)
        m = cl.tdmap(aa[:-1],bb[:-1],cc[:-1],dd[:-1],ee[:-1],rhs[:-1])
        m = np.insert(m, m.size, m[0,0])
    else:
        # solve regular tridiagonal matrix
        m = cl.tdma(aa,bb,cc,rhs)
    # compute spline coefficients
    a, b, c, d = coeff(h,y,m)
    # compute interpolated function values
    yint = spline(a,b,c,d,x,xint)
    return yint

def check_input(bcs,x,y,xint,bcs_val):
    """
    Function to validate intput data. 
    """
    if x.size != y.size:
        raise ValueError('x and y must have same length.')
    else:
        if x.size < 3:
            raise ValueError('At least three data points needed for cubic spline interpolation.')
    if ((xint.min() < x.min()) or (xint.max() > x.max())): 
        raise ValueError('No extrapolation possible, just interpolation - check boarders of xint.')
    if any(i < 0 for i in np.diff(x)):
        raise ValueError('x values must be strictly increasing.')
    if len(bcs) != 2:
        raise ValueError("'bcs' must contain two elements to specify start and end conditions.")
    if 'periodic' in bcs:
        if not('periodic' == bcs[0] and 'periodic' == bcs[1]): 
            raise ValueError("'periodic' `bc_type` is defined for both "
                             "curve ends and cannot be used with other "
                             "boundary conditions.")
        if not np.allclose(y[0], y[-1], rtol=1e-15, atol=1e-15):
            raise ValueError('Start and end point of y needs to be similar for periodic BCs.')
    bcs_list = ['natural','periodic','clamped','fixed2','fixed1']
    for i in range(2):
        if not bcs[i] in bcs_list:
            raise ValueError('Wrong choice of boundary condition.')
        if bcs[i] in ['fixed1', 'fixed2']:
            if all(i==0 for i in bcs_val):
                print('Fixed derivative values at endpoints to zero, similar to clamped or natural Bcs.')

def tridiag(bcs,h): 
    '''
    Building left-hand side, tridiagonal matrix of the linear system. 
    
    Parameters
    ----------
    bc    : boundary conditions for start and end of spline
    h     : grid spacing

    Returns
    -------
    a,b,c : diagonals
    '''
    n = h.size; a = np.zeros(n+1); b = np.zeros(n+1); c = np.zeros(n+1)
    if bcs[0] in ['natural', 'fixed2']:
        a[ 0] = 0 
        c[ 0] = 0 
    if bcs[0] in ['clamped', 'fixed1']: 
        a[ 0] = 0 
        c[ 0] = 1 
    if bcs[0] == 'periodic': 
        a[ 0] = h[n-1] / (h[0] + h[-1])
        c[ 0] = h[  0] / (h[0] + h[-1])
    b[ 0] = 2
    b[-1] = 2
    for i in range(1,n): 
        a[i] = h[i-1] / (h[i-1] + h[i])
        b[i] = 2
        c[i] = h[i  ] / (h[i-1] + h[i])
    if bcs[1] in ['natural', 'fixed2']:
        a[-1] = 0
        c[-1] = 0
    if bcs[1] in ['clamped', 'fixed1']: 
        a[-1] = 1
        c[-1] = 0
    if bcs[1] == 'periodic': 
        a[-2] = h[-2] / (h[-2] + h[-1])
        c[-2] = h[-1] / (h[-2] + h[-1])
    return a, b, c

def rhs_mat(bcs,h,y,bcs_val):
    '''
    Building right-hand side of the linear system. 
    
    Parameters
    ----------
    bc      : boundary conditions
    h       : grid spacing
    y       : function values of original data
    bcs_val : boundary values for first or second derivative at 
              start and end of spline, default is [0,0].

    Returns
    -------
    rhs     : right-hand side of linear system
    '''
    n   = h.size; rhs = np.zeros(n+1)    
    if bcs[0] == 'natural': 
        rhs[ 0] = 0
    if bcs[0] == 'fixed2': 
        rhs[ 0] = 2*bcs_val[0]
    if bcs[0] == 'clamped': 
        rhs[ 0] = (6 / h[ 0]) *                 (y[ 1] - y[ 0]) / h[ 0]
    if bcs[0] == 'fixed1': 
        rhs[ 0] = (6 / h[ 0]) * (-bcs_val[0] + ((y[ 1] - y[ 0]) / h[ 0]))
    if bcs[0] == 'periodic': 
        rhs[ 0] = (y[1] - y[0]) / h[0] - (y[-1] - y[-2]) / h[-1]
        rhs[ 0] = 6 * rhs[ 0] / (h[0] + h[-1])
    for i in range(1,n):
        rhs[i] = (y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]
        rhs[i] = 6 * rhs[i] / (h[i] + h[i-1])
    if bcs[1] == 'natural': 
        rhs[-1] = 0
    if bcs[1] == 'fixed2': 
        rhs[-1] = 2*bcs_val[1]
    if bcs[1] == 'clamped': 
        rhs[-1] = (6 / h[-1]) *                 (y[-2] - y[-1]) / h[-1]
    if bcs[1] == 'fixed1': 
        rhs[-1] = (6 / h[-1]) * ( bcs_val[1] + ((y[-2] - y[-1]) / h[-1]))
    if bcs[1] == 'periodic': 
        rhs[-2] = (y[-1] - y[-2]) / h[-1] - (y[-2] - y[-3]) / h[-2]
        rhs[-2] = 6 * rhs[-2] / (h[-1] + h[-2])
    return rhs  

def coeff(h,y,m):
    '''
    Compute spline coefficients. 
    
    Parameters
    ----------
    h       : grid spacing
    y       : function values of original data
    m       : second derivative (solution of linear system)

    Returns
    -------
    a,b,c,d : spline coefficients
    '''
    n = y.size; a = np.zeros(n-1); b = np.zeros(n-1); c = np.zeros(n-1); d = np.zeros(n-1)
    for i in range(n-1):
        a[i] = (m[i+1] - m[i]) / (6 * h[i])
        b[i] = m[i] / 2
        c[i] = (y[i+1] - y[i]) / h[i] - (m[i+1] + 2*m[i])*(h[i] / 6)
        d[i] = y[i]
    return a, b, c, d

def spline(a,b,c,d,x,xint):
    '''
    Evaluate cubic spline. 
    
    Parameters
    ----------
    a,b,c,d : spline coefficients
    x       : grid nodes of original data
    xint    : grid nodes of interpolated data

    Returns
    -------
    yint    : interpolated function on xint.
    '''
    n = x.size; nint = xint.size; yint = np.zeros(nint)
    for i in range(nint):
        idx     = min(bisect.bisect(x, xint[i])-1, n-2)
        z       = xint[i] - x[idx]
        yint[i] = a[idx]*z**3 + b[idx]*z**2 + c[idx]*z + d[idx] 
    return yint

def read_txt(path,name,header=0,footer=0):
    '''
    Reads file.txt / file.dat to numpy array.

    Parameters
    ----------
    path    : path to file
    name    : name of file
    header  : skip lines if header is present
    footer  : skip lines if footer is present
    
    Returns
    -------
    data    : data
    '''
    file = path + name
    print('--------------------------------------------------')
    print('read txt file to array:  ', file)
    with open(file, 'r') as f:
        raw = f.readlines()
        raw = np.asanyarray(raw)
        #---------------------------------------------#
        # clean up
        skiprows = header 
        raw = np.delete(raw, np.s_[:skiprows], axis=0)
        if footer > 0:
            skiprows = footer-header-1 
            raw = np.delete(raw, np.s_[skiprows:], axis=0)  
        # split and get array size
        i = 0
        col_out = len(raw[0].split())
        lin_out = len(raw)
        data    = np.zeros((lin_out,col_out))
        for line in raw:
            data_line = line.split()
            for j in range(col_out):
                data[i,j] = float(data_line[j])
            i += 1
    print('data array size       :  ', int(lin_out), 'x', int(col_out))
    return data
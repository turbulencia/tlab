#!/usr/bin/env python3
import numpy as np 
"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
Library for higher-order compact finite difference schemes

(for testing/debugging/devolpment purposes, based on the 
 implementation in tlab/src/fdm/)

Created  on & by: 2022/02/15 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: ...        - ...
#---------------------------------------------------------------------------#
Contains:
            - wrapper function for sixth-order compact derivative    schemes
            - wrapper function for sixth-order compact interpolation schemes
            - initialization - Jacobian
            - tri-/pentadigonal non-/interpolatory sixth-order 
              derivative schemes for non-/periodic boundary conditions
            - tridiagonal interpolation schemes of sixth order
            - tri-/pentadiagonal LU-solvers for periodic/non-periodic
              boundary schemes
            - ...
            - ...
#############################################################################
"""

def compactder(u, x, imin_set_zero=0, imax_set_zero=0, scheme='6t', order=1, \
               alpha=1/3, periodic=True, interpol=False, direction='vp'):
    '''
    Wrapper function for sixth-order compact derivative schemes

    Parameters
    ----------
    u             : function u[i,jk] to be differentiated along i 
    x             : grid nodes
    imin_set_zero : for non-periodic boundary conditions [0,1]
    imax_set_zero : for non-periodic boundary conditions [0,1]
    scheme        : sixth-order tridiagonal '6t' & pentadiagonal '6p' schemes
    order         : only 1st derivative schemes are implemented
    alpha         : compact parameter, needed for modified pentadiagonal schemes
                    (cf. Lele 1992)
    periodic      : boundary condition periodic/non-periodic
    interpol      : flag for interpolatory schemes on half grid nodes
    direction     : interpolation direction 
                    (velocity grid <--> pressure grid) 'vp' or 'pv'

    Returns
    -------
    dudx          : derviative of u-field
    '''
    # check input combination
    if all(np.diff(x)[0]==np.diff(x)):
        uniform = True
    else:
        uniform = False
    if not uniform and periodic:
        raise ValueError('Error - Only periodic+uniform or non-periodic+(non-)unifrom possible!') 
    if order !=1:
        raise ValueError('Error - Only first derivative schemes are implemented!') 
    # initialization - jacobian
    jac = gjac(x, scheme, order, alpha, periodic, uniform)
    if interpol:
        if periodic:
            xpre   = x + 0.5*np.diff(x)[0]      # pressure grid (periodic==uniform,    dx==dxpre)
        else:
            xpre   = x[:-1] + 0.5*np.diff(x)[:] # pressure grid (different grid sizes, dx!=dxpre) 
        jacpre = gjac(xpre, scheme, order, alpha, periodic, uniform)
    # 6th order tri-(6t)/pentadiagonal(6p) schemes
    if not interpol:
        if scheme == '6t':
            if periodic:
                a,b,c         = fdm_c1n6p_lhs(jac[:,0])
                a,b,c,d,e     = tforwardp(a,b,c)
                rhs           = fdm_c1n6p_rhs(u)
                dudx          = tbackwardp(a,b,c,d,e,rhs)
            else:
                a,b,c         = fdm_c1n6_lhs(jac[:,0], imin_set_zero, imax_set_zero)
                a,b,c         = tforward(a,b,c)
                rhs           = fdm_c1n6_rhs(u, imin_set_zero, imax_set_zero)
                dudx          = tbackward(a,b,c,rhs)
        elif scheme == '6p':
            if periodic: 
                a,b,c,d,e     = fdm_c1n6mp_lhs(jac[:,0], alpha)
                a,b,c,d,e,f,g = pforwardp(a,b,c,d,e) 
                rhs           = fdm_c1n6mp_rhs(u, alpha)
                dudx          = pbackwardp(a,b,c,d,e,f,g,rhs)
            else:
                a,b,c,d,e     = fdm_c1n6m_lhs(jac[:,0], imin_set_zero, imax_set_zero, alpha)
                a,b,c,d,e     = pforward(a,b,c,d,e)
                rhs           = fdm_c1n6m_rhs(u, imin_set_zero, imax_set_zero, alpha)
                dudx          = pbackward(a,b,c,d,e,rhs)
    elif interpol: 
        if scheme == '6t':
            if direction=='vp':
                if periodic:
                    a,b,c     = fdm_c1int6p_lhs(jacpre[:,0])
                    a,b,c,d,e = tforwardp(a,b,c)
                    rhs       = fdm_c1intvp6p_rhs(u)
                    dudx      = tbackwardp(a,b,c,d,e,rhs)
                else:
                    a,b,c     = fdm_c1intvp6_lhs(jacpre[:,0])
                    a,b,c     = tforward(a,b,c)
                    rhs       = fdm_c1intvp6_rhs(u)
                    dudx      = tbackward(a,b,c,rhs)
            elif direction=='pv':
                if periodic:
                    a,b,c     = fdm_c1int6p_lhs(jac[:,0])
                    a,b,c,d,e = tforwardp(a,b,c)
                    rhs       = fdm_c1intpv6p_rhs(u)
                    dudx      = tbackwardp(a,b,c,d,e,rhs)
                else:
                    a,b,c     = fdm_c1intpv6_lhs(jac[:,0])
                    a,b,c     = tforward(a,b,c)
                    rhs       = fdm_c1intpv6_rhs(u)
                    dudx      = tbackward(a,b,c,rhs)
        elif scheme == '6p':
            raise ValueError('Error - Not implemented yet!')       
    return dudx

def compactint(u, x, direction='vp', periodic=True):
    '''
    Wrapper function for sixth-order compact interpolation schemes

    Parameters
    ----------
    u             : function u[i,jk] to be differentiated along i 
    x             : grid nodes
    direction     : interpolation direction 
                    (velocity grid <--> pressure grid) 'vp' or 'pv'
    periodic      : boundary condition periodic/non-periodic

    Returns
    -------
    uint          : interpolation of u-field
    '''
    # check input combination
    if all(np.diff(x)[0]==np.diff(x)):
        uniform = True
    else:
        uniform = False
    if not uniform and periodic:
        raise ValueError('Error - Only periodic+uniform or non-periodic+(non-)unifrom possible!')
    # 6th order tridiagonal interpolation scheme
    if periodic:
        if direction=='vp':
            a,b,c     = fdm_c0int6p_lhs(x.size)
            a,b,c,d,e = tforwardp(a,b,c)
            rhs       = fdm_c0intvp6p_rhs(u)
            uint      = tbackwardp(a,b,c,d,e,rhs)
        elif direction=='pv':
            a,b,c     = fdm_c0int6p_lhs(x.size)
            a,b,c,d,e = tforwardp(a,b,c)
            rhs       = fdm_c0intpv6p_rhs(u)
            uint      = tbackwardp(a,b,c,d,e,rhs)
    else:
        if direction=='vp':
            a,b,c     = fdm_c0intvp6_lhs(x.size)
            a,b,c     = tforward(a,b,c)
            rhs       = fdm_c0intvp6_rhs(u)
            uint      = tbackward(a,b,c,rhs)
        elif direction=='pv':
            a,b,c     = fdm_c0intpv6_lhs(x.size)
            a,b,c     = tforward(a,b,c)
            rhs       = fdm_c0intpv6_rhs(u)
            uint      = tbackward(a,b,c,rhs)
    return uint

def gjac(x, scheme, order, alpha, periodic, uniform):
    '''
    Jacobian function, initialization stage for compact derivative schemes.

    Parameters
    ----------
    x             : grid nodes
    scheme        : sixth-order tridiagonal '6t' & pentadiagonal '6p' schemes
    order         : only 1st derivative schemes are implemented
    alpha         : compact parameter, needed for modified pentadiagonal schemes
                    (cf. Lele 1992)
    periodic      : boundary condition periodic/non-periodic
    uniform       : distribution of grid nodes (non-/uniform)
        
    Returns
    -------
    gjac          : Jacobian
    '''
    # grid scale
    if periodic:
        scale = x[-1] + np.diff(x)[0] - x[0]
    else:
        scale = x[-1]
    # first derivative - inner points
    gjac = np.zeros((x.size, 4))
    if uniform:
        for i in range(1, x.size-1):
            gjac[i,0] = (x[i+1] - x[i-1]) * 0.5
        # first derivative - boundary points
        if periodic:
            gjac[-1,0] = (x[0] + scale - x[-2]              ) * 0.5
            gjac[ 0,0] = (x[1] - x[0] + x[0] + scale - x[-1]) * 0.5
        else:
            gjac[-1,0] = gjac[-2,0]
            gjac[ 0,0] = gjac[ 1,0]
        # second derivative is zero
        gjac[:,1] = 0
    else:
        # first derivative 
        gjac[:,0] = 1        
        if scheme == '6t':
            a,b,c     = fdm_c1n6_lhs(dx=gjac[:,0], imin_set_zero=0, imax_set_zero=0)
            a,b,c     = tforward(a,b,c)
            gjac[:,0] = fdm_c1n6_rhs(u=x, imin_set_zero=0, imax_set_zero=0)[:,0]
            gjac[:,0] = tbackward(a,b,c,gjac[:,0])[:,0]
        elif scheme == '6p':
            a,b,c,d,e = fdm_c1n6m_lhs(dx=gjac[:,0], imin_set_zero=0, imax_set_zero=0, alpha=alpha)
            a,b,c,d,e = pforward(a,b,c,d,e)
            gjac[:,0] = fdm_c1n6m_rhs(u=x, imin_set_zero=0, imax_set_zero=0, alpha=alpha)[:,0]
            gjac[:,0] = pbackward(a,b,c,d,e,gjac[:,0])[:,0]
        if order > 1:
            raise ValueError('Error - order > 1: gjacu for non-uniform grids not implmented yet!')  
    # Saving operations for the time-stability constraint
    gjac[:,2] = 1 / gjac[:,0]
    gjac[:,3] = gjac[:,2] * gjac[:,2]
    return gjac

def fdm_c1n6p_lhs(dx):
    '''
    1st derivative fdm with 6th order compact schemes, periodic, tridiagonal
    by JCP Lele 1992. Interior points according to Eq. 2.1.7 (\alpha=1/3). 
    System multiplied by 18/14 to eliminate one multiplication in the RHS.
    
    Left-hand side; tridiagonal matrix of the linear system.

    Parameters
    ----------
    dx    : grid spacing

    Returns
    -------
    a,b,c : diagonals
    '''
    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    for i in range(n):
        a[i] = 3/7
        b[i] = 18/14
        c[i] = 3/7
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i]   = b[i]  *dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c

    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    for i in range(n):
        a[i] = 3/7
        b[i] = 18/14
        c[i] = 3/7
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i  ] = b[i  ]*dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c
def fdm_c1n6p_rhs(u):
    '''
    Right-hand side; forcing term.

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # (be careful with the indexing: python <--> fortran)
    imm1 = imax - 1 
    for i in range(1,imax+1):
        im2 = i-2; im2=im2+imm1; im2=int(im2%imax)#+1
        im1 = i-1; im1=im1+imm1; im1=int(im1%imax)#+1
        ip1 = i+1; ip1=ip1+imm1; ip1=int(ip1%imax)#+1
        ip2 = i+2; ip2=ip2+imm1; ip2=int(ip2%imax)#+1
        for j in range(u.shape[1]):
            d[i-1,j] = u[ip1,j] - u[im1,j] + (1/28)*(u[ip2,j] - u[im2,j])  
    return d

def fdm_c1n6_lhs(dx, imin_set_zero, imax_set_zero):
    '''
    1st derivative fdm with 6th order compact schemes, non-periodic, tridiagonal
    by JCP Lele 1992. Interior points according to Eq. 2.1.7 (\alpha=1/3).
    Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of 
    boundary points on stability. However, the scheme 3-4-6-4-3 implemented
    here seems to work better than the 3-5-6-5-3 proposed there. This
    fourth-order scheme is Eq. 2.1.6 with \alpha=1/4.
    The BCs are imposed with biased third-order Eq. 4.1.3 with \alpha=2.
    
    Left-hand side; tridiagonal matrix of the linear system.

    Parameters
    ----------
    dx            : grid spacing
    imin_set_zero : boundary conditions [0,1]
    imax_set_zero : boundary conditions [0,1]

    Returns
    -------
    a,b,c         : diagonals
    '''
    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    # Bcs
    vmult_imin = 1
    vmult_imax = 1
    if imin_set_zero == 1: vmult_imin = 0
    if imax_set_zero == 1: vmult_imax = 0
    # third-order biased
    a[0]  = 0
    b[0]  = 1
    c[0]  = 2 * vmult_imin
    #
    a[-1] = 2 * vmult_imax
    b[-1] = 1
    c[-1] = 0
    #
    a[1]  = 1/6
    b[1]  = 1
    c[1]  = 1/2
    #
    a[-2] = 1/2
    b[-2] = 1
    c[-2] = 1/6
    for i in range(2,n-2):    
        a[i] = 1/3
        b[i] = 1
        c[i] = 1/3
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i  ] = b[i  ]*dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c

def fdm_c1n6_rhs(u, imin_set_zero, imax_set_zero):
    '''
    Right-hand side; forcing term.

    Parameters
    ----------
    u             : function u[i,jk] to be differentiated along i 
    imin_set_zero : boundary conditions [0,1]
    imax_set_zero : boundary conditions [0,1]

    Returns
    -------
    d             : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # Bcs
    vmult_imin = 1
    vmult_imax = 1
    if imin_set_zero == 1: vmult_imin = 0
    if imax_set_zero == 1: vmult_imax = 0
    
    c12dx1  = 0.5 * vmult_imin
    c52dx1  = 5   * c12dx1
    c2dx1   = 2   * vmult_imin
    
    c12dxmx = 0.5 * vmult_imax
    c52dxmx = 5   * c12dxmx
    c2dxmx  = 2   * vmult_imax
    
    c0118 = 1/18
    c0102 = 1/2
    c0509 = 5/9
    
    c1418 = 14/18
    c0136 = 1 /36

    for j in range(u.shape[1]):
        # third-order
        d[0,j]  = - c52dx1 *u[ 0,j] + c2dx1 *u[ 1,j] + c12dx1 *u[ 2,j]
        d[-1,j] =   c52dxmx*u[-1,j] - c2dxmx*u[-2,j] - c12dxmx*u[-3,j]
        # fifth-order biased
        d[1,j]  =   c0118*u[ 3,j] + u[ 2,j] - c0102*u[ 1,j] - c0509*u[ 0,j]
        d[-2,j] = - c0118*u[-4,j] - u[-3,j] + c0102*u[-2,j] + c0509*u[-1,j]
    for i in range(2,imax-2):
        for j in range(u.shape[1]):
            d[i,j] = c1418*(u[i+1,j] - u[i-1,j]) + c0136*(u[i+2,j] - u[i-2,j])
    return d

def fdm_c1n6mp_lhs(dx, alpha):
    '''
    1st derivative fdm with 6th order compact schemes, periodic, pentadiagonal
    by JCP Lele 1992. Interior points according to Eq. 2.1.10. 
    
    Left-hand side; pentadiagonal matrix of the linear system.

    Parameters
    ----------
    dx        : grid spacing
    alpha     : compact parameter of lhs of Eq. 2.1.10.

    Returns
    -------
    a,b,c,d,e : diagonals
    '''
    n = dx.size;     a = np.zeros(n); b = np.zeros(n)
    c = np.zeros(n); d = np.zeros(n); e = np.zeros(n)
    aa,bb,cc,beta = fdm_c1_coeff(alpha) # compact parameters
    for i in range(n):
        a[i] = beta  
        b[i] = alpha 
        c[i] = 1      
        d[i] = alpha 
        e[i] = beta  
    # Jacobian Multiplication
    e[-2] = e[-2]*dx[0]
    d[-1] = d[-1]*dx[0]
    c[ 0] = c[ 0]*dx[0]
    b[ 1] = b[ 1]*dx[0]
    a[ 2] = a[ 2]*dx[0]
    #
    e[-1] = e[-1]*dx[1]
    d[ 0] = d[ 0]*dx[1]
    c[ 1] = c[ 1]*dx[1]
    b[ 2] = b[ 2]*dx[1]
    a[ 3] = a[ 3]*dx[1]
    #
    for i in range(2,n-2):
        e[i-2] = e[i-2]*dx[i]
        d[i-1] = d[i-1]*dx[i]
        c[i  ] = c[i  ]*dx[i]
        b[i+1] = b[i+1]*dx[i]
        a[i+2] = a[i+2]*dx[i]
    #
    e[-4] = e[-4]*dx[-2]
    d[-3] = d[-3]*dx[-2]
    c[-2] = c[-2]*dx[-2]
    b[-1] = b[-1]*dx[-2]
    a[ 0] = a[ 0]*dx[-2]
    #
    e[-3] = e[-3]*dx[-1]
    d[-2] = d[-2]*dx[-1]
    c[-1] = c[-1]*dx[-1]
    b[ 0] = b[ 0]*dx[-1]
    a[ 1] = a[ 1]*dx[-1]
    return a,b,c,d,e

def fdm_c1n6mp_rhs(u, alpha):
    '''
    Right-hand side; forcing term.

    Parameters
    ----------
    u     : function u[i,jk] to be differentiated along i 
    alpha : compact parameter on lhs of Eq. 2.1.10.

    Returns
    -------
    d     : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # computed compact coefficients 
    aa,bb,cc,beta = fdm_c1_coeff(alpha)
    # indexing (be careful with the indexing: python <--> fortran)
    for i in range(1,imax+1):
        im3 = int((i+imax-4)%imax)
        im2 = int((i+imax-3)%imax)
        im1 = int((i+imax-2)%imax)
        ip1 = int((i       )%imax)
        ip2 = int((i+1     )%imax)
        ip3 = int((i+2     )%imax)
        for j in range(u.shape[1]):
            d[i-1,j] = (aa/2)*(u[ip1,j] - u[im1,j]) + (bb/4)*(u[ip2,j] - u[im2,j]) + (cc/6)*(u[ip3,j] - u[im3,j])  
    return d

def fdm_c1_coeff(alpha):
    '''
    Parameters for compact pentadiagonal schemes, periodic by JCP Lele 1992 
    according to Eq. 2.1.10 with similar truncation error like tridiagonal
    scheme Eq. 2.1.7 with (\alpha=1/3). 

    Parameters
    ----------
    alpha : compact parameter of lhs of Eq. 2.1.10.

    Returns
    -------
    beta  : compact parameter of lhs of Eq. 2.1.10.
    a,b,c : compact parameter of rhs of Eq. 2.1.10.
    '''
    beta = (2/5)*(alpha - 1/3)
    # two-parameter family of 6th order and pentadiag lhs 
    a = (1/ 6)*( 9 +    alpha - 20*beta)
    b = (1/15)*(-9 + 32*alpha + 62*beta)
    c = (1/10)*( 1 -  3*alpha + 12*beta)
    # e = (12/math.factorial(7))*(3 - 8*alpha + 20*beta) # truncation error
    return a,b,c,beta#,e

def fdm_c1n6m_lhs(dx, imin_set_zero, imax_set_zero, alpha):
    '''
    1st derivative fdm with 6th order compact schemes, non-periodic, pentadiagonal
    by JCP Lele 1992. Interior points according to Eq. 2.1.7 (\alpha=1/3).
    Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of 
    boundary points on stability.Scheme 3-5-c6--penta6m--c6-5-3 is implmented.
    Interior points according to Eq. 2.1.10.
    
    Left-hand side; pentadiagonal matrix of the linear system.

    Parameters
    ----------
    dx            : grid spacing
    imin_set_zero : boundary conditions [0,1]
    imax_set_zero : boundary conditions [0,1]
    alpha         : compact parameter of lhs of Eq. 2.1.10.

    Returns
    -------
    a,b,c,d,e     : diagonals
    '''
    n = dx.size;     a = np.zeros(n); b = np.zeros(n)
    c = np.zeros(n); d = np.zeros(n); e = np.zeros(n)
    aa,bb,cc,beta = fdm_c1_coeff(alpha)
    # Bcs
    vmult_imin = 1
    vmult_imax = 1
    if imin_set_zero == 1: vmult_imin = 0
    if imax_set_zero == 1: vmult_imax = 0
    # third-order biased
    a[0]  = 0
    b[0]  = 0
    c[0]  = 1
    d[0]  = 2 * vmult_imin
    e[0]  = 0
    #
    a[-1] = 0
    b[-1] = 2 * vmult_imax
    c[-1] = 1
    d[-1] = 0
    e[-1] = 0
    
    # fifth-order biased
    a[1]  = 0
    b[1]  = 1/6
    c[1]  = 1
    d[1]  = 1/2
    e[1]  = 0
    #
    a[-2] = 0
    b[-2] = 1/2
    c[-2] = 1
    d[-2] = 1/6
    e[-2] = 0
    
    # sixth-order centered (alpha=1/3)
    a[2]  = 0
    b[2]  = 1/3
    c[2]  = 1
    d[2]  = 1/3
    e[2]  = 0
    #
    a[-3] = 0
    b[-3] = 1/3
    c[-3] = 1
    d[-3] = 1/3
    e[-3] = 0
    
    # sixth-order modified centered
    for i in range(3,n-3):    
        a[i] = beta
        b[i] = alpha
        c[i] = 1
        d[i] = alpha
        e[i] = beta
    
    # Jacobian Multiplication
    e[-2] = e[-2]*dx[0]
    d[-1] = d[-1]*dx[0]
    c[ 0] = c[ 0]*dx[0]
    b[ 1] = b[ 1]*dx[0]
    a[ 2] = a[ 2]*dx[0]
    #
    e[-1] = e[-1]*dx[1]
    d[ 0] = d[ 0]*dx[1]
    c[ 1] = c[ 1]*dx[1]
    b[ 2] = b[ 2]*dx[1]
    a[ 3] = a[ 3]*dx[1]
    #
    for i in range(2,n-2):
        e[i-2] = e[i-2]*dx[i]
        d[i-1] = d[i-1]*dx[i]
        c[i  ] = c[i  ]*dx[i]
        b[i+1] = b[i+1]*dx[i]
        a[i+2] = a[i+2]*dx[i]
    #
    e[-4] = e[-4]*dx[-2]
    d[-3] = d[-3]*dx[-2]
    c[-2] = c[-2]*dx[-2]
    b[-1] = b[-1]*dx[-2]
    a[ 0] = a[ 0]*dx[-2]
    #
    e[-3] = e[-3]*dx[-1]
    d[-2] = d[-2]*dx[-1]
    c[-1] = c[-1]*dx[-1]
    b[ 0] = b[ 0]*dx[-1]
    a[ 1] = a[ 1]*dx[-1]
    return a,b,c,d,e

def fdm_c1n6m_rhs(u, imin_set_zero, imax_set_zero, alpha):   
    '''
    Right-hand side; forcing term.

    Parameters
    ----------
    u             : function u[i,jk] to be differentiated along i 
    imin_set_zero : boundary conditions [0,1]
    imax_set_zero : boundary conditions [0,1]
    alpha         : compact parameter of lhs of Eq. 2.1.10.

    Returns
    -------
    d             : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # Bcs
    vmult_imin = 1
    vmult_imax = 1
    if imin_set_zero == 1: vmult_imin = 0
    if imax_set_zero == 1: vmult_imax = 0
    
    c12dx1  = 0.5 * vmult_imin
    c52dx1  = 5   * c12dx1
    c2dx1   = 2   * vmult_imin
    
    c12dxmx = 0.5 * vmult_imax
    c52dxmx = 5   * c12dxmx
    c2dxmx  = 2   * vmult_imax
    
    c0118 = 1/18
    c0102 = 1/2
    c0509 = 5/9
    
    c1418 = 14/18
    c0136 = 1 /36
    
    # computed coefficients 
    aa,bb,cc,beta = fdm_c1_coeff(alpha)

    for j in range(u.shape[1]):
        # third-order
        d[0,j]  = - c52dx1 *u[ 0,j] + c2dx1 *u[ 1,j] + c12dx1 *u[ 2,j]
        d[-1,j] =   c52dxmx*u[-1,j] - c2dxmx*u[-2,j] - c12dxmx*u[-3,j]
        # fifth-order biased
        d[1,j]  =   c0118*u[ 3,j] + u[ 2,j] - c0102*u[ 1,j] - c0509*u[ 0,j]
        d[-2,j] = - c0118*u[-4,j] - u[-3,j] + c0102*u[-2,j] + c0509*u[-1,j]
        # 6th-order centered with alpha=(1/3)
        d[2,j]  =   c1418*(u[ 3,j] - u[ 1,j]) + c0136*(u[ 4,j] - u[ 0,j])
        d[-3,j] =   c1418*(u[-2,j] - u[-4,j]) + c0136*(u[-1,j] - u[-5,j])
    # 6th-order modified (7-point stencil)   
    for i in range(3,imax-3):
        for j in range(u.shape[1]):
            d[i,j] =   (aa/2)*(u[i+1,j] - u[i-1,j]) \
                     + (bb/4)*(u[i+2,j] - u[i-2,j]) \
                     + (cc/6)*(u[i+3,j] - u[i-3,j])  
    return d

def fdm_c0int6p_lhs(imax):
    '''
    6th order compact interpolation schemes, periodic, tridiagonal by JCP Lele 1992. 
    Interior points according to Eq. C.1.4. (\alpha=3/10, \beta=0, c=0).
    System multiplied by 4/3 to eliminate one multiplication in the RHS.
    
    Left-hand side; tridiagonal matrix of the linear system.

    Parameters
    ----------
    imax  : number of nodes
    
    Returns
    -------
    a,b,c : diagonals
    '''
    n = imax; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    for i in range(n):
        a[i] = 2/5 # 3/10
        b[i] = 4/3 # 1
        c[i] = 2/5 # 3/10
    return a,b,c

def fdm_c0intvp6p_rhs(u):
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the right (velocity --> pressure grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # indexing (be careful with the indexing: python <--> fortran)
    for i in range(1,imax+1):
        im1 = int((i+imax-2)%imax)
        ip1 = int((i       )%imax)
        ip2 = int((i+1     )%imax)
        for j in range(u.shape[1]):
            d[i-1,j] = (u[ip1,j] + u[i-1,j]) + (1/15)*(u[ip2,j] + u[im1,j])  
            # d[i-1,j] = (3/4)*(u[ip1,j] + u[i-1,j]) + (1/20)*(u[ip2,j] + u[im1,j])  
    return d

def fdm_c0intpv6p_rhs(u):
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the left (pressure --> velocity grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # indexing (be careful with the indexing: python <--> fortran)
    for i in range(1,imax+1):
        im1 = int((i+imax-2)%imax)
        im2 = int((i+imax-3)%imax)
        ip1 = int((i       )%imax)
        for j in range(u.shape[1]):
            d[i-1,j] = (u[i-1,j] + u[im1,j]) + (1/15)*(u[ip1,j] + u[im2,j])  
            # d[i-1,j] = (3/4)*(u[i-1,j] + u[im1,j]) + (1/20)*(u[ip1,j] + u[im2,j])  
    return d

def fdm_c0intvp6_lhs(imax): 
    '''
    6th order compact interpolation schemes, non-periodic,tridiagonal by JCP Lele 1992. 
    Interior points according to Eq. C.1.4. (\alpha=3/10, \beta=0, c=0).
    Different boundary closures can be found in Albin 2010 
    https://doi.org/10.1002/fld.2520 table 6 (typo in Sd2_3i) and 
    Boersma 2005 https://doi.org/10.1016/j.jcp.2005.03.004 .
    
    Left-hand side; tridiagonal matrix of the linear system.
    
    ===> half gridpoint to the right (velocity --> pressure grid)

    Parameters
    ----------
    imax  : number of nodes
    
    Returns
    -------
    a,b,c : diagonals
    '''
    n = imax-1; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    # forth-order biased  (implicit)
    a[0]  = 0
    b[0]  = 1
    c[0]  = 1
    #
    a[-1] = 1
    b[-1] = 1
    c[-1] = 0
    # sixth-order (implicit)
    for i in range(1,n-1):    
        a[i] = 3/10
        b[i] = 1
        c[i] = 3/10  
    # # forth-order (implicit)
    # a[1]  = 1/6
    # b[1]  = 1
    # c[1]  = 1/6
    # #
    # a[-2] = 1/6
    # b[-2] = 1
    # c[-2] = 1/6
    return a,b,c

def fdm_c0intpv6_lhs(imax): 
    '''
    Left-hand side; tridiagonal matrix of the linear system.
    
    ===> half gridpoint to the left (pressure --> velocity grid)

    Parameters
    ----------
    imax  : number of nodes
    
    Returns
    -------
    a,b,c : diagonals
    '''
    n = imax+1; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    # # third-order biased (explicit)
    # a[0]  = 0
    # b[0]  = 1
    # c[0]  = 0
    # #
    # a[-1] = 0
    # b[-1] = 1
    # c[-1] = 0  
    # third-order biased (implicit)
    a[0]  = 0
    b[0]  = 1
    c[0]  = 3
    #
    a[-1] = 3
    b[-1] = 1
    c[-1] = 0   
    # forth-order (implicit)
    a[1]  = 1/6
    b[1]  = 1
    c[1]  = 1/6
    #
    a[-2] = 1/6
    b[-2] = 1
    c[-2] = 1/6
    # sixth-order (implicit)
    for i in range(2,n-2):    
        a[i] = 3/10
        b[i] = 1
        c[i] = 3/10
    return a,b,c

def fdm_c0intvp6_rhs(u): 
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the right (velocity --> pressure grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''    
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    imax = u.shape[0]-1; d = np.zeros((imax,u.shape[1]))
    #
    c0302 = 3/2
    c0104 = 1/4
    c0304 = 3/4
    c0120 = 1/20
    #
    for j in range(u.shape[1]):
        # forth-order biased (implicit)
        d[ 0,j] = c0302*u[ 1,j] + c0104*(u[ 0,j] + u[ 2,j])
        d[-1,j] = c0302*u[-2,j] + c0104*(u[-1,j] + u[-3,j])
    for i in range(1,imax-1):
        for j in range(u.shape[1]):
            d[i,j] = c0304*(u[i+1,j] + u[i,j]) + c0120*(u[i+2,j] + u[i-1,j])
    #     # forth-order biased (implicit)
    #     d[ 1,j] = (2/3)*(u[ 2,j] + u[ 1,j])
    #     d[-2,j] = (2/3)*(u[-3,j] + u[-2,j])
    # for i in range(2,imax-2):
    #     for j in range(u.shape[1]):
    #         d[i,j] = c0304*(u[i+1,j] + u[i,j]) + c0120*(u[i+2,j] + u[i-1,j]) 
    return d

def fdm_c0intpv6_rhs(u): 
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the left (pressure --> velocity grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''    
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    imax = u.shape[0]+1; d = np.zeros((imax,u.shape[1]))
    #
    c0304 = 3/4
    c0120 = 1/20
    #
    for j in range(u.shape[1]):
        # # third-order biased (explicit)
        # d[ 0,j] = 15/8*u[ 0,j] - 5/4*u[ 1,j] + 3/8*u[ 2,j]
        # d[-1,j] = 15/8*u[-1,j] - 5/4*u[-2,j] + 3/8*u[-3,j]
        # third-order biased (implicit) 
        d[ 0,j] = 3*u[ 0,j] + u[ 1,j]
        d[-1,j] = 3*u[-1,j] + u[-2,j]
        # forth-order (implicit)
        d[ 1,j] = 2/3*(u[ 0,j] + u[ 1,j])
        d[-2,j] = 2/3*(u[-1,j] + u[-2,j])
    for i in range(2,imax-2):
        for j in range(u.shape[1]):
            d[i,j] = c0304*(u[i,j] + u[i-1,j]) + c0120*(u[i+1,j] + u[i-2,j]) 
    return d
 
def fdm_c1int6p_lhs(dx):
    '''
    1st interpolatory derivative fdm with 6th order compact schemes, periodic, 
    tridiagonal by JCP Lele 1992. Interior points according to Eq. B.1.1. 
    (\alpha=9/62, \beta=0, c=0). System multiplied by 62/63 to eliminate 
    one multiplication in the RHS
    
    Left-hand side; tridiagonal matrix of the linear system.

    Parameters
    ----------
    dx    : grid spacing

    Returns
    -------
    a,b,c : diagonals
    '''
    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    for i in range(n):
        a[i] =  9/63 # 9/62
        b[i] = 62/63 # 1
        c[i] =  9/63 # 9/62
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i]   = b[i]  *dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c

def fdm_c1intvp6p_rhs(u):
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the right (velocity --> pressure grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # indexing (be careful with the indexing: python <--> fortran)
    for i in range(1,imax+1):
        im1 = int((i+imax-2)%imax)
        ip1 = int((i       )%imax)
        ip2 = int((i+1     )%imax)
        for j in range(u.shape[1]):
            d[i-1,j] = (u[ip1,j] - u[i-1,j]) + (17/189)*(u[ip2,j] - u[im1,j])  
            # d[i-1,j] = (63/62)*(u[ip1,j] - u[i-1,j]) + (1/3)*(17/62)*(u[ip2,j] - u[im1,j])  
    return d

def fdm_c1intpv6p_rhs(u): 
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the left (pressure --> velocity grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    d = np.zeros(u.shape); imax = u.shape[0]
    # indexing (be careful with the indexing: python <--> fortran)
    for i in range(1,imax+1):
        im1 = int((i+imax-2)%imax)
        im2 = int((i+imax-3)%imax)
        ip1 = int((i       )%imax)
        for j in range(u.shape[1]):
            d[i-1,j] = (u[i-1,j] - u[im1,j]) + (17/189)*(u[ip1,j] - u[im2,j])  
            # d[i-1,j] = (63/62)*(u[i-1,j] - u[im1,j]) + (1/3)*(17/62)*(u[ip1,j] - u[im2,j]) 
    return d

def fdm_c1intvp6_lhs(dx):
    '''
    1st interpolatory derivative fdm with 6th order compact schemes, non-periodic, 
    tridiagonal by JCP Lele 1992. Interior points according to Eq. B.1.1. 
    (\alpha=9/62, \beta=0, c=0).System for this scheme is multiplied by 62/63 
    to eliminate one multiplication in the RHS. Different boundary closures 
    can be found in Albin 2010 https://doi.org/10.1002/fld.2520 table 6
    (typo in Sd2_3i) and Boersma 2005 https://doi.org/10.1016/j.jcp.2005.03.004 .
        
    Left-hand side; tridiagonal matrix of the linear system.

    ===> half gridpoint to the right (velocity --> pressure grid)
    
    Parameters
    ----------
    dx    : grid spacing

    Returns
    -------
    a,b,c : diagonals
    '''
    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    # # third-order biased (explicit)
    # a[0]  =  0 
    # b[0]  =  1
    # c[0]  =  0
    # #
    # a[-1] =  0
    # b[-1] =  1
    # c[-1] =  0
    # third-order biased (implicit)
    a[0]  =  0 
    b[0]  =  1
    c[0]  = -1
    #
    a[-1] = -1
    b[-1] =  1
    c[-1] =  0
    # sixth-order (implicit)
    for i in range(1,n-1):
        a[i] =  9/63 # 9/62
        b[i] = 62/63 # 1
        c[i] =  9/63 # 9/62
    # # forth-order (implicit)
    # a[1]  =  1/22 
    # b[1]  =  1
    # c[1]  =  1/22
    # #
    # a[-2] =  1/22
    # b[-2] =  1
    # c[-2] =  1/22
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i]   = b[i]  *dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c

def fdm_c1intpv6_lhs(dx): 
    '''
    Left-hand side; tridiagonal matrix of the linear system.

    ===> half gridpoint to the left (pressure --> velocity grid)
    
    Parameters
    ----------
    dx    : grid spacing

    Returns
    -------
    a,b,c : diagonals
    '''
    n = dx.size; a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
    # # third-order biased (explicit)
    # a[0]  =  0 
    # b[0]  =  1
    # c[0]  =  0
    # #
    # a[-1] =  0
    # b[-1] =  1
    # c[-1] =  0
    # third-order biased (implicit)
    a[0]  =  0 
    b[0]  =  1
    c[0]  =  23
    #
    a[-1] =  23
    b[-1] =  1
    c[-1] =  0
    # forth-order (implicit)
    a[1]  =  1/22
    b[1]  =  1
    c[1]  =  1/22
    #
    a[-2] =  1/22
    b[-2] =  1
    c[-2] =  1/22
    for i in range(2,n-2):
        a[i] =  9/63 # 9/62
        b[i] = 62/63 # 1
        c[i] =  9/63 # 9/62
    # Jacobian Multiplication
    c[-1] = c[-1]*dx[0]
    b[ 0] = b[ 0]*dx[0]
    a[ 1] = a[ 1]*dx[0]
    for i in range(1,n-1):
        c[i-1] = c[i-1]*dx[i]
        b[i]   = b[i]  *dx[i]
        a[i+1] = a[i+1]*dx[i]
    c[-2] = c[-2]*dx[-1]
    b[-1] = b[-1]*dx[-1]
    a[ 0] = a[ 0]*dx[-1]
    return a,b,c

def fdm_c1intvp6_rhs(u):
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the right (velocity --> pressure grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    imax = u.shape[0]-1; d = np.zeros((imax,u.shape[1]))
    #
    for j in range(u.shape[1]):
        # third-order biased (explicit)
        # d[ 0,j] = -(23/24)*u[ 0,j] + (7/8)*u[ 1,j] + (1/8)*u[ 2,j] - (1/24)*u[ 3,j]
        # d[-1,j] =  (23/24)*u[-1,j] - (7/8)*u[-2,j] - (1/8)*u[-3,j] + (1/24)*u[-4,j]
        # third-order biased (implicit)
        d[ 0,j] = -u[ 0,j] + 2*u[ 1,j] - u[ 2,j]
        d[-1,j] =  u[-1,j] - 2*u[-2,j] + u[-3,j]
    for i in range(1,imax-1):
        for j in range(u.shape[1]):
            d[i,j] = (u[i+1,j] - u[i,j]) + (17/189)*(u[i+2,j] - u[i-1,j])  
    #     # forth-order (implicit)
    #     d[ 1,j] =  12/11*(u[ 2,j] - u[ 1,j])
    #     d[-2,j] = -12/11*(u[-3,j] - u[-2,j])
    # for i in range(2,imax-2):
    #     for j in range(u.shape[1]):
    #         d[i,j] = (u[i+1,j] - u[i,j]) + (17/189)*(u[i+2,j] - u[i-1,j])  
    return d

def fdm_c1intpv6_rhs(u):
    '''
    Right-hand side; forcing term.

    ===> half gridpoint to the left (pressure --> velocity grid)

    Parameters
    ----------
    u : function u[i,jk] to be differentiated along i 

    Returns
    -------
    d : right-hand side vector of the linear system
    '''
    if len(u.shape)==1: u=np.expand_dims(u, axis=1) # add dimension 
    imax = u.shape[0]+1; d = np.zeros((imax,u.shape[1]))
    #
    for j in range(u.shape[1]):     
        # third-order biased (explicit)
        # d[ 0,j] = -(23/24)*u[ 0,j] + (7/8)*u[ 1,j] + (1/8)*u[ 2,j] - (1/24)*u[ 3,j]
        # d[-1,j] =  (23/24)*u[-1,j] - (7/8)*u[-2,j] - (1/8)*u[-3,j] + (1/24)*u[-4,j]    
        # third-order biased (implicit)
        d[ 0,j] = -25*u[ 0,j] + 26*u[ 1,j] - u[ 2,j]
        d[-1,j] =  25*u[-1,j] - 26*u[-2,j] + u[-3,j]
        # forth-order (implicit)
        d[ 1,j] =  12/11*(u[ 1,j] - u[ 0,j])
        d[-2,j] = -12/11*(u[-2,j] - u[-1,j])
    for i in range(2,imax-2):
        for j in range(u.shape[1]):
            d[i,j] = (u[i,j] - u[i-1,j]) + (17/189)*(u[i+1,j] - u[i-2,j]) 
    return d

def tdma(a,b,c,rhs):
    '''
    Wrapper function for tridiagonal non-periodic solver based 
    on Thomas algorithm.

    Parameters
    ----------
    a,b,c : diagonals
    rhs   : rhs
    
    Returns
    -------
    x     : solution
    '''
    a,b,c= tforward(a,b,c) 
    x    = tbackward(a,b,c,rhs)
    return x

def tforward(a,b,c):
    '''
    LU factorization stage; L with diagonal unity.

    Parameters
    ----------
    a,b,c : diagonals
    
    Returns
    -------
    a,b,c : factored LHS
    '''
    imax = a.size
    for i in range(1,imax):
        a[i] = a[i] / b[i-1]
        b[i] = b[i] - a[i]*c[i-1]
    # Final operations
    a[:] = - a[:]
    b[:] = 1 / b[:]
    c[:] = - c[:]
    return a,b,c 
 
def tbackward(a,b,c,r): 
    '''
    Backward substitution step in the Thomas algorithm

    Parameters
    ----------
    a,b,c : factored LHS
    r     : rhs 
    
    Returns
    -------
    r     : solution 
    '''
    if len(r.shape)==1: r=np.expand_dims(r, axis=1) # add dimension 
    # Forward sweep
    for i in range(1,a.size):
        dummy1 = a[i]
        for j in range(r.shape[1]):
            r[i,j] = r[i,j] + dummy1*r[i-1,j]                  
    # Backward sweep
    for j in range(r.shape[1]):
        r[-1,j] = r[-1,j]*b[-1]
    for i in range(a.size-2,-1,-1):
        for j in range(r.shape[1]):
            r[i,j] = (r[i,j] + c[i]*r[i+1,j])*b[i]
    return r

def tdmap(a,b,c,d,e,rhs):
    '''
    Wrapper function for tridiagonal periodic/circulant solver based 
    on Thomas algorithm.

    Parameters
    ----------
    a,b,c : diagonals
    rhs   : rhs
    
    Returns
    -------
    x     : solution
    '''
    a,b,c,d,e = tforwardp(a,b,c) 
    x         = tbackwardp(a,b,c,d,e,rhs)
    return x

def tforwardp(a,b,c): 
    '''
    LU factorization stage; L with diagonal unity.

    Parameters
    ----------
    a,b,c     : diagonals
    
    Returns
    -------
    a,b,c,d,e : factored LHS
    '''
    imax = a.size; d = np.zeros(imax); e = np.zeros(imax)
    # Generate first elements of LU
    c[0] = c[0]/b[0]
    e[0] = a[0]/b[0]
    d[0] = c[-1]
    # Generate n=2 to n=n-2 elements of LU
    for i in range(1,imax-2):
        b[i] = b[i] - a[i]*c[i-1]
        c[i] = c[i]/b[i]               
        e[i] =-a[i]*e[i-1]/b[i]
        d[i] =-d[i-1]*c[i-1]
    # Generate n-1 elements
    b[-2] =  b[-2] - a[-2]*c[-3]
    e[-2] = (c[-2] - a[-2]*e[-3])/b[-2]
    d[-2] =  a[-1] - d[-3]*c[-3]
    # Generate the n-th element
    summ = 0
    for i in range(imax-1):
        summ = summ + d[i]*e[i]
    b[-1] = b[-1] - summ
    # Final operations
    for i in range(imax):
        b[i] = 1/b[i]
        a[i] =-a[i]*b[i]
        c[i] =-c[i]
        e[i] =-e[i]
    return a,b,c,d,e

def tbackwardp(a,b,c,d,e,r): 
    '''
    Backward substitution step in the Thomas algorithm

    Parameters
    ----------
    a,b,c,d,e : factored LHS
    r         : rhs 
    
    Returns
    -------
    r         : solution 
    '''
    imax = a.size
    if len(r.shape)==1: r=np.expand_dims(r, axis=1) # add dimension 
    # Forward sweep
    dummy1 = b[0]
    for j in range(r.shape[1]):
        r[0,j] = r[0,j]*dummy1
    for i in range(1,imax-1):    
        dummy1 = a[i]
        dummy2 = b[i]
        for j in range(r.shape[1]):
            r[i,j] = r[i,j]*dummy2 + dummy1*r[i-1,j]
    wrk = np.zeros((r.shape[1]))
    for i in range(imax-1):
        dummy1 = d[i]
        for j in range(r.shape[1]):
            wrk[j] = wrk[j] + dummy1*r[i,j]
    dummy1 = b[-1]
    for j in range(r.shape[1]):
        r[-1,j] = (r[-1,j] - wrk[j])*dummy1
    # Backward sweep
    dummy1 = e[-2]
    for j in range(r.shape[1]):
        r[-2,j] = dummy1*r[-1,j] + r[-2,j]
    for i in range(imax-3,-1,-1):
        dummy1 = c[i]
        dummy2 = e[i]        
        for j in range(r.shape[1]):
            r[i,j] =  r[i,j] + dummy1* r[i+1,j] + dummy2* r[-1,j]
    return r

def pdma(a,b,c,d,e,rhs):
    '''
    Wrapper function for pentadiagonal non-periodic solver based 
    on Thomas algorithm (reverse ordering).

    Parameters
    ----------
    a,b,c,d,e : diagonals
    rhs       : rhs
    
    Returns
    -------
    x         : solution
    '''
    a,b,c,d,e = pforward( a,b,c,d,e) 
    x         = pbackward(a,b,c,d,e,rhs)
    return x

def pforward(a,b,c,d,e): 
    '''
    LU factorization stage.

    Parameters
    ----------
    a,b,c,d,e : diagonals
    
    Returns
    -------
    a,b,c,d,e : factored LHS
    '''
    imax = a.size
    i = imax -1
    e[i] = 1 # padding
    d[i] = 1 # padding
    c[i] = c[i] 
    b[i] = b[i] 
    i = imax-2
    e[i] = 1 # padding
    d[i] =(d[i]              )/c[i+1]
    c[i] = c[i] - d[i]*b[i+1] 
    b[i] = b[i] - d[i]*a[i+1]
    for i in range(imax-3,1,-1):
       e[i] = e[i]                             /c[i+2]
       d[i] =(d[i]               - e[i]*b[i+2])/c[i+1]
       c[i] = c[i] - d[i]*b[i+1] - e[i]*a[i+2]
       b[i] = b[i] - d[i]*a[i+1]
    i = 1
    e[i] = e[i]                             /c[i+2]
    d[i] =(d[i]               - e[i]*b[i+2])/c[i+1]
    c[i] = c[i] - d[i]*b[i+1] - e[i]*a[i+2]
    b[i] = b[i] - d[i]*a[i+1]
    a[i] = 1 # padding
    i = 0
    e[i] = e[i]                             /c[i+2]
    d[i] =(d[i]               - e[i]*b[i+2])/c[i+1]
    c[i] = c[i] - d[i]*b[i+1] - e[i]*a[i+2]
    b[i] = 1 # padding
    a[i] = 1 # padding  
    return a,b,c,d,e

def pbackward(a,b,c,d,e,r): 
    '''
    Backward substitution step in the Thomas algorithm.

    Parameters
    ----------
    a,b,c,d,e : factored LHS
    r         : rhs 
    
    Returns
    -------
    r         : solution 
    '''
    imax = a.size
    if len(r.shape)==1: r=np.expand_dims(r, axis=1) # add dimension 
    # Solve Ly=f, forward
    i = imax-2
    r[i,:] = r[i,:] - r[i+1,:]*d[i]
    for i in range(imax-3,-1,-1):
        r[i,:] = r[i,:] -r[i+1,:]*d[i] -r[i+2,:]*e[i]
    # Solve Ux=y, backward
    i = 0
    r[i,:] = r[i,:] /c[i]
    i = 1
    r[i,:] =(r[i,:] -r[i-1,:]*b[i])/c[i]
    for i in range(2,imax):
        r[i,:] =(r[i,:] -r[i-1,:]*b[i] -r[i-2,:]*a[i])/c[i]
    return r

def pdmap(a,b,c,d,e,rhs):
    '''
    Wrapper function for pentadiagonal periodic/cyclic Toeplitz solver based 
    on Sherman–Morrison–Woodbury Formula  "A new algorithm for solving nearly 
    penta-diagonal Toeplitz linear systems". Algorithm 2.2 is implemented here
    https://doi.org/10.1016/j.camwa.2011.12.044

    Parameters
    ----------
    a,b,c,d,e : diagonals
    rhs       : rhs
    
    Returns
    -------
    x         : solution
    '''
    a,b,c,d,e,f,g = pforwardp( a,b,c,d,e) 
    x             = pbackwardp(a,b,c,d,e,f,g,rhs)
    return x

def pforwardp(a,b,c,d,e): 
    '''
    LU factorization stage.

    Parameters
    ----------
    a,b,c,d,e     : diagonals
    
    Returns
    -------
    a,b,c,d,e,f,g : factored LHS
    '''    
    # build regular modified pentadiagonal matrix A1  
    # save off-digonal entries
    a0 = a[ 0] # upper-right corner
    b0 = b[ 0]  
    en = e[-1] # lower-left corner
    dn = d[-1]   
    # diagonals of A1
    b[1] = b[1] - d[-1]
    c[0] = c[0] - e[-1]
    c[1] = c[1] - e[-1]
    #
    c[-2] = c[-2] - a[ 0]
    c[-1] = c[-1] - a[ 0]
    d[-2] = d[-2] - b[ 0]
    # set offdigonal entries to zero
    a[:2]  = 0
    b[:1]  = 0
    d[-1:] = 0
    e[-2:] = 0
    # regular forward step for A1
    a,b,c,d,e = pforward( a,b,c,d,e)    
    # save offdigonal entries again
    a[ 0] = a0 # upper-right corner
    b[ 0] = b0  
    e[-1] = en # lower-left corner
    d[-1] = dn  
    # define matrix u
    u = np.zeros((a.size,2)) 
    u[ 0,0] = 1
    u[-2,0] = 1
    u[ 1,1] = 1
    u[-1,1] = 1
    # regular backward step for u
    r = pbackward(a,b,c,d,e,u)
    # compute entries of matrix m
    m1  = e[-1] * r[0,-2] + a[ 0] * r[-2,-2] + b[ 0] * r[-1,-2] + 1
    m2  = e[-1] * r[0,-1] + a[ 0] * r[-2,-1] + b[ 0] * r[-1,-1]
    m3  = d[-1] * r[0,-2] + e[-1] * r[ 1,-2] + a[ 0] * r[-1,-2]
    m4  = d[-1] * r[0,-1] + e[-1] * r[ 1,-1] + a[ 0] * r[-1,-1] + 1
    # check if m is invertible
    if((m1*m4 - m2*m3) < 1e-8):
        raise ValueError('Error - Matrix m is not invertible.')    
    # store u solutions in f,g additional diagonals
    f = r[:,0]
    g = r[:,1]
    return a,b,c,d,e,f,g

def pbackwardp(a,b,c,d,e,f,g,r): 
    '''
    Backward substitution step.

    Parameters
    ----------
    a,b,c,d,e,f,g : factored LHS
    r             : rhs 
    
    Returns
    -------
    r             : solution 
    '''
    if len(r.shape)==1: r=np.expand_dims(r, axis=1) # add dimension 
    imax = a.size
    # regular backward step
    r = pbackward(a,b,c,d,e,r)
    # store solution of u in r
    r = np.append(r,np.expand_dims(f, axis=1),axis=1)
    r = np.append(r,np.expand_dims(g, axis=1),axis=1)
    # compute entries of matrix m
    m1  = e[-1] * r[0,-2] + a[ 0] * r[-2,-2] + b[ 0] * r[-1,-2] + 1
    m2  = e[-1] * r[0,-1] + a[ 0] * r[-2,-1] + b[ 0] * r[-1,-1]
    m3  = d[-1] * r[0,-2] + e[-1] * r[ 1,-2] + a[ 0] * r[-1,-2]
    m4  = d[-1] * r[0,-1] + e[-1] * r[ 1,-1] + a[ 0] * r[-1,-1] + 1
    # compute d coefficients
    di  =  1 / (m1*m4    - m2*m3)
    d11 = di * (m4*e[-1] - m2*d[-1])
    d12 = di * (m4*b[ 0] - m2*a[ 0])
    d13 = di *  m4*a[ 0]
    d14 = di *  m2*e[-1]
    d21 = di * (m1*d[-1] - m3*e[-1])
    d22 = di * (m1*a[ 0] - m3*b[ 0])
    d23 = di *  m3*a[ 0]
    d24 = di *  m1*e[-1]
    # solve
    for j in range(r.shape[1]-2):
        dummy1 = d11*r[0,j] + d12*r[-1,j] + d13*r[-2,j] - d14*r[1,j]
        dummy2 = d21*r[0,j] + d22*r[-1,j] - d23*r[-2,j] + d24*r[1,j]  
        for i in range(imax):
            r[i,j] = r[i,j] - dummy1*r[i,-2] - dummy2*r[i,-1]
    
    return r[:,:-2] # exclude u
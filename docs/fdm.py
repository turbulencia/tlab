import numpy as np

def fdm1_c6_A(n):
# Creates the system matrix for the finite-difference approximation
# to the first-order derivative using 6th-order compact schemes over
# a uniform grid
# Input arguments:
#    - n: number of grid points
# Output arguments:
#    - l: vectors with the diagonals of the LHS system matrix
    l = np.empty( (3,n) )
    
    l[0,:] = 1. /3.
    l[1,:] = 1.
    l[2,:] = 1. /3.

    # 2nd grid point from the boundary
    l[0,1] = 1. /6.
    l[2,1] = 1. /2.

    l[0,-2] = 1. /2.
    l[2,-2] = 1. /6.

    # 1st grid point from the boundary
    l[2,0] = 2.

    l[0,-1] = 2.

    return l
    
def fdm1_c6_B(n):
# Creates the system matrix for the finite-difference approximation
# to the first-order derivative using 6th-order compact schemes over
# a uniform grid
# Input arguments:
#    - n: number of grid points
# Output arguments:
#    - a,b,c: vectors with the diagonals of the RHS system matrix
    r = np.empty( (5,n) )
        
    r[0,:] =-1. /36.
    r[1,:] =-7. /9.
    r[2,:] = 0.
    r[3,:] = 7. /9.
    r[4,:] = 1. /36.

    # 2nd grid point from the boundary
    r[4][1] = 1./18.
    r[3][1] = 1.
    r[2][1] =-1./2.
    r[1][1] =-5./9.

    r[0][-2] =-1./18.
    r[1][-2] =-1.
    r[2][-2] = 1./2.
    r[3][-2] = 5./9.

    # 1st grid point from the boundary
    r[4][0] = 1./2.
    r[3][0] = 2.
    r[2][0] =-5./2.

    r[0][-1] =-1./2.
    r[1][-1] =-2.
    r[2][-1] = 5./2.

    return r

def fdm2_c6_A(n, kc):
# Creates the system matrix for the finite-difference approximation
# to the second-order derivative using 6th-order compact schemes over
# a uniform grid
# Input arguments:
#    - n:   number of grid points
#    - kc:  Type of scheme: =48./7. clasical with c=0, =np.pi**2. Lamballais2011
# Output arguments:
#    - l:   vectors with the diagonals of the LHS system matrix
    l = np.empty( (3,n) )
    
    alpha = (272.  -45. *kc ) /( 416.  -90.  *kc )

    l[0,:] = alpha
    l[1,:] = 1.
    l[2,:] = alpha
  
    # 3rd grid point from the boundary
    l[0,2] = 2. /11.
    l[2,2] = 2. /11.

    l[0,-3] = 2. /11.
    l[2,-3] = 2. /11.
    # l[0,2] = 1. /10   # Biased 5th order
    # l[2,2] =-7. /20.

    # l[0,-3] =-7. /20.
    # l[2,-3] = 1. /10.
  
    # 2nd grid point from the boundary
    l[0,1] = 1. /10     # Centered 4th order
    l[2,1] = 1. /10.

    l[0,-2] = 1. /10.
    l[2,-2] = 1. /10.
    # l[0,1] = 1. /10   # Biased 5th order
    # l[2,1] =-7. /20.

    # l[0,-2] =-7. /20.
    # l[2,-2] = 1. /10.
    
    # 1st grid point from the boundary
    l[2,0] = 11.

    l[0,-1] = 11.
    
    return l
    
def fdm2_c6_B(n, kc):
# Creates the system matrix for the finite-difference approximation
# to the second-order derivative using 6th-order compact schemes over
# a uniform grid
# Input arguments:
#    - n: number of grid points
#    - kc:  Type of scheme: =48./7. clasical with c=0, =np.pi**2. Lamballais2011
# Output arguments:
#    - r: vectors with the diagonals of the RHS system matrix
    r = np.empty( (7,n) )
    
    a =     (48.  -135. *kc ) /( 1664. -360. *kc )
    b =     (528.  -81. *kc ) /( 208.  -45.  *kc ) /4.
    c =    -(432.  -63. *kc ) /( 1664. -360. *kc ) /9.
    
    r[0,:] = c
    r[1,:] = b
    r[2,:] = a
    r[3,:] =-2. *(a+b+c)
    r[4,:] = a
    r[5,:] = b
    r[6,:] = c
    
    # 3rd grid point from the boundary
    r[1,2] = 3. /44.
    r[2,2] = 12./11.
    r[3,2] =-2. *( 3./44. +12./11 )
    r[4,2] = 12./11.
    r[5,2] = 3. /44.
    r[6,2] = 0.
 
    r[0,-3] = 0.
    r[1,-3] = 3. /44.
    r[2,-3] = 12./11.
    r[3,-3] =-2. *( 3./44. +12./11 )
    r[4,-3] = 12./11.
    r[5,-3] = 3. /44.

    # r[1,2] = 0.
    # r[2,2] = 99. /80.     # Biased 5th order
    # r[3,2] =-3.
    # r[4,2] = 186. /80.
    # r[5,2] =-3. /5.
    # r[6,2] = 3. /80.
 
    # r[0,-3] = 3. /80.
    # r[1,-3] =-3. /5.
    # r[2,-3] = 186. /80.
    # r[3,-3] =-3.
    # r[4,-3] = 99. /80.
    # r[5,-3] = 0.

    # 2nd grid point from the boundary
    r[2,1] = 6. /5.     # Centered 4th order
    r[3,1] =-2. *( 6. / 5. )
    r[4,1] = 6. /5.
    r[5,1] = 0.
    r[6,1] = 0.
 
    r[0,-2] = 0.
    r[1,-2] = 0.
    r[2,-2] = 6. /5.
    r[3,-2] =-2. *( 6. / 5. )
    r[4,-2] = 6. /5.

    # r[2,1] = 99. /80.     # Biased 5th order
    # r[3,1] =-3.
    # r[4,1] = 186. /80.
    # r[5,1] =-3. /5.
    # r[6,1] = 3. /80.
 
    # r[0,-2] = 3. /80.
    # r[1,-2] =-3. /5.
    # r[2,-2] = 186. /80.
    # r[3,-2] =-3.
    # r[4,-2] = 99. /80.
    
    # 1st grid point from the boundary
    r[3,0] = 13.
    r[4,0] =-27.
    r[5,0] = 15.
    r[6,0] =-1.
 
    r[0,-1] =-1
    r[1,-1] = 15.
    r[2,-1] =-27.
    r[3,-1] = 13.
    
    return r

def fdm2_c6_wavenumber( w, kc ):
    # kc = 48./7.     # Classical Lele scheme with c = 0
    # kc = np.pi **2. # Lamballais2011
    alpha = (272.  -45. *kc ) /( 416.  -90.  *kc )
    a =     (48.  -135. *kc ) /( 1664. -360. *kc )
    b =     (528.  -81. *kc ) /( 208.  -45.  *kc )
    c =    -(432.  -63. *kc ) /( 1664. -360. *kc )
    
    num = 2. *a *( 1. -np.cos( w ) ) + 0.5 *b *( 1. -np.cos( 2. *w ) ) + 2. *c /9.  *( 1. -np.cos( 3. *w ) )
    den = 1. + 2. *alpha *np.cos( w )
    
    return num /den

def fdm1_c6_wavenumber( w ):
    alpha = 1. /3.
    a =     2. /3. *( 2. +alpha     )
    b =     1. /3. *( 4. *alpha -1. )
    c =     0.
    
    num = a *np.sin( w ) + 0.5 *b *np.sin( 2. *w ) + c /3.  *np.sin( 3. *w )
    den = 1. + 2. *alpha *np.cos( w )
    
    return num /den


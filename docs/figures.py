import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import fdm

from matplotlib import rc
rc('text',  usetex=True)
rc('font',  family='serif', size=11)
# rc('axes',  titlesize='medium')
rc('legend',fontsize='small')
rc('axes',  grid=True)
rc('grid',  color='0.9')
rc('lines', solid_capstyle='round')
# rc('figure',     dpi=100)
rc('savefig',    dpi=100)

FiguresToPlot=['spectra','wavenumber','stability']
FiguresToPlot=['convergence']

###############################################################################
tag='spectra'
if tag in FiguresToPlot:
    fig_id = 0

    n = 100
    f = np.linspace(0, 5., num=n)

    fig_id = fig_id +1
    plt.figure( figsize = (4,3) )
    plt.plot( f, np.ones(n),                        label=r'uniform')
    plt.plot( f, f**2. *np.exp(-2. *(f-1.) ),       label=r'quadratic')
    plt.plot( f, f**4. *np.exp(-2. *(f **2. -1) ),  label=r'quartic')
    plt.plot( f, np.exp(-0.5 *(f-1.) **2. *6 **2.), label=r'Gaussian')
    plt.xlabel(r'normalized frequency $f/f_0$')
    plt.ylabel(r'normalized spectra $E/E_0$')
    plt.legend( )
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(tag+str(fig_id)))
    plt.show()

###############################################################################
tag='wavenumber'
if tag in FiguresToPlot:
    fig_id = 0

    # Modified wavenumber
    w = np.linspace( 0., np.pi, num = 100)

    # Modified wavenumber of first order derivative
    fig_id = fig_id +1
    plt.figure( figsize = (4,3) )

    plt.plot( w, w,         label=r'exact',         c='k', lw=1.0 )
    plt.plot( w, np.sin(w), label=r'$\delta_x$ C2', c='C0', alpha=0.25 )
    wm = fdm.fdm1_c6_wavenumber(w)
    plt.plot( w, wm,        label=r'$\delta_x$ C6', c='C0' )
    j = wm.argmax()
    plt.text( w[j]-1., wm[j], r'${:3.2f}\,@\,{:3.2f}\pi$'.format(wm[j],w[j]/np.pi), va='bottom' )
    plt.plot( (w[j],w[j]-1), (wm[j],wm[j]), 'k', ls='-', lw=0.5, marker='o', mfc='w', markevery=2, ms=4 )

    plt.xticks( np.linspace(0., np.pi, num = 7), [r'$0.0$', r'', r'$\pi/3$', r'', r'$2\pi/3$', r'', r'$\pi$'])
    plt.yticks( np.linspace(0., np.pi, num = 7), [r'$0.0$', r'', r'$\pi/3$', r'', r'$2\pi/3$', r'', r'$\pi$'], rotation='vertical')
    plt.xlabel(r'$\omega=\kappa h = 2\pi/\mathrm{PPW}$')
    plt.ylabel(r'Im($\lambda_1)$')
    plt.legend( )
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(tag+str(fig_id)))

    # Relative
    fig_id = fig_id +1
    plt.figure( figsize = (4,3) )

    plt.plot( w, np.sin(w) /w,                  label=r'$\delta_x$ C2', c='C0', alpha=0.25 )
    plt.plot( w, fdm.fdm1_c6_wavenumber(w) /w,  label=r'$\delta_x$ C6', c='C0' )

    plt.xticks( np.linspace(0., np.pi, num = 7), [r'$0.0$', r'', r'$\pi/3$', r'', r'$2\pi/3$', r'', r'$\pi$'])
    plt.xlabel(r'$\omega=\kappa h = 2\pi/\mathrm{PPW}$')
    plt.ylabel(r'$\mathrm{Im}(\lambda_1) /[\mathrm{Im}(\lambda_1)]_\mathrm{e}$')
    plt.legend( )
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(tag+str(fig_id)))

    # Modified wavenumber of second order derivative
    fig_id = fig_id +1
    plt.figure( figsize = (4,3) )

    plt.plot( w, w **2.,                                label=r'exact',                 c='k', lw=1.0 )
    plt.plot( w, fdm.fdm1_c6_wavenumber(w) **2.,        label=r'$(\delta_x$ C6)$^2$',   c='C0' )
    plt.plot( w, 2. *( 1 -np.cos(w) ),                  label=r'$\delta_{xx}$ C2',      c='C1', alpha=0.25 )
    wm = fdm.fdm2_c6_wavenumber(w, 48./7.)
    plt.plot( w, wm,                                    label=r'$\delta_{xx}$ C6',      c='C1', alpha=0.5 )
    j = -1
    plt.text( w[j]-1., wm[j], r'${:3.2f}\,@\,\pi$'.format(wm[j]), va='bottom' )
    plt.plot( (w[j],w[j]-1), (wm[j],wm[j]), 'k', ls='-', lw=0.5, marker='o', mfc='w', markevery=2, ms=4 )
    plt.plot( w, fdm.fdm2_c6_wavenumber(w, np.pi **2.),   label=r'$\delta_{xx}$ C6b',             c='C1' )

    plt.xticks( np.linspace(0., np.pi, num = 7), [r'$0.0$', r'', r'$\pi/3$', r'', r'$2\pi/3$', r'', r'$\pi$'])
    plt.yticks( np.linspace(0., np.pi **2., num = 7), [r'$0.0$', r'', r'$\pi^2/3$', r'', r'$2\pi^2/3$', r'', r'$\pi^2$'], rotation='vertical')
    plt.xlabel(r'$\omega=\kappa h = 2\pi/\mathrm{PPW}$')
    plt.ylabel(r'-Re($\lambda_2)$')
    plt.legend( )
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(tag+str(fig_id)))

    # Relative
    fig_id = fig_id +1
    plt.figure( figsize = (4,3) )

    plt.plot( w, fdm.fdm1_c6_wavenumber(w) **2. /w **2.,        label=r'$(\delta_x$ C6)$^2$',   c='C0' )
    plt.plot( w, 2. *( 1 -np.cos(w) ) /w **2.,                  label=r'$\delta_{xx}$ C2',      c='C1', alpha=0.25 )
    plt.plot( w, fdm.fdm2_c6_wavenumber(w, 48./7.) /w **2.,     label=r'$\delta_{xx}$ C6',      c='C1', alpha=0.5 )
    plt.plot( w, fdm.fdm2_c6_wavenumber(w, np.pi **2.) /w **2., label=r'$\delta_{xx}$ C6b',     c='C1' )

    plt.xticks( np.linspace(0., np.pi, num = 7), [r'$0.0$', r'', r'$\pi/3$', r'', r'$2\pi/3$', r'', r'$\pi$'])
    # plt.yticks( np.linspace(0., np.pi **2., num = 7), [r'$0.0$', r'', r'$\pi^2/3$', r'', r'$2\pi^2/3$', r'', r'$\pi^2$'], rotation='vertical')
    plt.xlabel(r'$\omega=\kappa h = 2\pi/\mathrm{PPW}$')
    plt.ylabel(r'$\mathrm{Re}(\lambda_2)/[\mathrm{Re}(\lambda_2)]_\mathrm{e}$')
    plt.legend( )
    plt.tight_layout(pad=0.1)
    plt.savefig("{}.pdf".format(tag+str(fig_id)))

    plt.show()

###############################################################################
tag='stability'
if tag in FiguresToPlot:
    fig_id = 0

    def PlotBackground(w,r,s,tag,wi,wr,cfl_a,cfl_d):
        colors = [ '#507dbc', '#86A7D3', '#bbd1ea', '#f9b5ac', '#F49690', '#ee7674' ]

        # plt.contourf(np.real(w),np.imag(w),r,[0., 1.],colors=['#aedcc0'],alpha=0.5)
        plt.contourf(np.real(w),np.imag(w),s,[-0.1,-0.01,0.,0.01,0.1],colors=colors,alpha=0.75,extend='both')
        # plt.contour( np.real(w),np.imag(w),abs(s_masked),[1.],linewidths=[1.0],colors='w')
        plt.colorbar(orientation='horizontal',shrink=0.6, label=tag, pad=0.05)
        plt.contour( np.real(w),np.imag(w),r,    [1.],colors=['k'],linewidths=[1.0])
        plt.xlabel(r'Re($\lambda\tau)$',loc='right',labelpad=-2)
        plt.ylabel(r'Im($\lambda\tau)$',loc='bottom',labelpad=-22)
        plt.gca().set_aspect('equal')#,'box')
        plt.gca().spines['left'].set_position(('data', 0))
        plt.gca().spines['bottom'].set_position(('data', 0))
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.grid()
        plt.gca().xaxis.set_ticklabels([])
        plt.gca().yaxis.set_ticklabels([])

        # wi, wr = 3.34, -4.65
        plt.text( -1.5, wi, r'${:3.2f}$'.format(wi), va='bottom' )
        plt.plot( (0.,-1.5), (wi,wi), 'k', ls='-', lw=0.5, marker='o', mfc='w', markevery=2, ms=4 )
        plt.text( wr, -1.5, r'${:3.2f}$'.format(wr), ha='left' )
        plt.plot( (wr,wr), (0,-1.5), 'k', ls='-', lw=0.5, marker='o', mfc='w', markevery=2, ms=4 )

        wi, wr = cfl_a *1.99,-cfl_d *np.pi **2. #6.86
        plt.plot( wr *np.array([0., 1., 1., 0., 0.]), wi *np.array([1., 1., -1., -1., 1.]), 'k', lw='0.5')

        plt.text( wr-1.5, wi, r'${}={:3.2f}$'.format('\mathrm{CFL}_\mathrm{a}',cfl_a),     va='bottom' )
        plt.plot( (wr-0.1,wr-1.5), (wi,wi), 'k', ls='-', lw=0.5 )
        plt.text( wr, -wi-1., r'${}={:3.2f}$'.format('\mathrm{CFL}_\mathrm{d}',cfl_d),     ha='right' )
        plt.plot( (wr,wr), (-wi-0.1,-wi-1.), 'k', ls='-', lw=0.5 )

        wi, wr = cfl_a *np.pi /2., -cfl_d *(np.pi /2.) **2.
        plt.plot( wr *np.array([0., 1., 1., 0., 0.]), wi *np.array([1., 1., -1., -1., 1.]), 'k', lw='0.5')
        plt.text( wr-1.5, -wi, r'$\mathrm{PPW}=4$',     va='bottom' )
        plt.plot( (wr-0.1,wr-1.5), (-wi,-wi), 'k', ls='-', lw=0.5 )

        return

    # Stability Regions
    # x = np.linspace(-3.0,0.6,400)
    # y = np.linspace(-2.75,2.75,400)
    # x, y = np.meshgrid(x,y)
    # w = x + 1j *y
    # r = 1. + w + w **2. /2. + w **3. /6.                                # Runge-Kitta 3
    # wi, wr = 1.73, -2.52
    # cfl_a = 0.6
    # cfl_d = 0.15

    x = np.linspace(-5.0,1.0,400)
    y = np.linspace(-4.0,4.0,400)
    x, y = np.meshgrid(x,y)
    w = x + 1j *y
    r = 1. + w + w **2. /2. + w **3. /6. + w **4. /24. + w **5. /200.   # Runge-Kitta 4-5
    wi, wr = 3.34, -4.65
    cfl_a = 1.2
    cfl_d = 0.3

    s = r/np.exp(w)
    r = np.abs(r)
    s_masked = np.ma.masked_where( r > 1., s )

    # Eigenvalues
    n = 128

    # Periodic case
    w1 = np.linspace( -np.pi, np.pi, num = n, endpoint=False)
    lambdasS = w1 *cfl_a *1j -w1 **2. *cfl_d
    lambdasC = fdm.fdm1_c6_wavenumber(w1) *cfl_a *1j -fdm.fdm2_c6_wavenumber(w1, np.pi **2.) *cfl_d

    # Nonperiodic case, grid step from 1 to 1+h2 as sigmoid function. Set h2 = 0 for uniform grid
    h2 = 0.0 #1.5
    ndelta = 4.
    s = np.linspace( 1., float(n), num=n)
    dummy  = np.exp( -(s -float(n)/2.) /ndelta )
    xp1 = 1. + h2 /( 1. + dummy )                               # First order derivative of x=x(s)
    xp2 = ( h2 /ndelta ) *( 1. /( 1. + dummy ) ) **2. *dummy    # Second-order derivative
    D2 = np.diagflat( xp2 /xp1 **2 )                            # Correction diagonal matrix

    # Nonperiodic case, 1st order derivative
    l1 = fdm.fdm1_c6_A(n)               *xp1
    A1 = np.diagflat(l1[0,1:],-1) +np.diagflat(l1[1,:],0) +np.diagflat(l1[2,:-1],1)

    r1 = fdm.fdm1_c6_B(n)
    B1 = np.diagflat(r1[0][2:],-2) +np.diagflat(r1[1][1:],-1) +np.diagflat(r1[2],0) \
        +np.diagflat(r1[3][:-1],1) +np.diagflat(r1[4][:-2],2)

    # Nonperiodic case, 2nd order derivative
    l2 = fdm.fdm2_c6_A(n, np.pi **2.)   *( xp1 **2. )
    # l2 = fdm.fdm2_c6_A(n, 48./7.) *( xp1 **2. )
    A2 = np.diagflat(l2[0,1:],-1) +np.diagflat(l2[1,:],0) +np.diagflat(l2[2,:-1],1)

    r2 = fdm.fdm2_c6_B(n, np.pi **2.)
    # r2 = fdm.fdm2_c6_B(n, 48./7.)
    B2 = np.diagflat(r2[0][3:],-3) +np.diagflat(r2[1][2:],-2) +np.diagflat(r2[2][1:],-1) +np.diagflat(r2[3],0) \
        +np.diagflat(r2[4][:-1],1) +np.diagflat(r2[5][:-2],2) +np.diagflat(r2[6][:-3],3)

    # Reduced arrays for Dirichlet boundary conditions at j=0
    A1 = A1[1:,1:]
    A1[0,0] -= l1[0,1] *l1[2,0] /l1[1,0]
    B1 = B1[1:,1:]
    B1[0,0] -= l1[0,1] *r1[3,0] /l1[1,0]
    B1[0,1] -= l1[0,1] *r1[4,0] /l1[1,0]

    A2 = A2[1:,1:]
    A2[0,0] -= l2[0,1] *l2[2,0] /l2[1,0]
    B2 = B2[1:,1:]
    B2[0,0] -= l2[0,1] *r2[4,0] /l2[1,0]
    B2[0,1] -= l2[0,1] *r2[5,0] /l2[1,0]
    B2[0,2] -= l2[0,1] *r2[6,0] /l2[1,0]

    D2 = D2[1:,1:]

    # Reduced arrays for Dirichlet boundary conditions at j=n
    A1 = A1[:-1,:-1]
    A1[-1,-1] -= l1[2,-2] *l1[0,-1] /l1[1,-1]
    B1 = B1[:-1,:-1]
    B1[-1,-1] -= l1[2,-2] *r1[1,-1] /l1[1,-1]
    B1[-1,-2] -= l1[2,-2] *r1[0,-1] /l1[1,-1]

    A2 = A2[:-1,:-1]
    A2[-1,-1] -= l2[2,-2] *l2[0,-1] /l2[1,-1]
    B2 = B2[:-1,:-1]
    B2[-1,-1] -= l2[2,-2] *r2[2,-1] /l2[1,-1]
    B2[-1,-2] -= l2[2,-2] *r2[1,-1] /l2[1,-1]
    B2[-1,-3] -= l2[2,-2] *r2[0,-1] /l2[1,-1]

    D2 = D2[:-1,:-1]

    L = -( cfl_a +cfl_d *D2 )*scipy.linalg.solve(A1,B1) + cfl_d *scipy.linalg.solve(A2,B2)
    lambdas = scipy.linalg.eigvals( L )

    # Plot
    fig_id = fig_id +1
    fig, ((f1, f2)) = plt.subplots(nrows=1, ncols=2, figsize=(8,6))

    plt.subplot(f1)
    PlotBackground( w, r, abs(s_masked)-1.,             r'amplitude error $\rho-1$',wi,wr,cfl_a,cfl_d )
    plt.plot(np.real(lambdasS),np.imag(lambdasS),marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='0.5')
    plt.plot(np.real(lambdasC),np.imag(lambdasC),marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='#6a2202')
    plt.plot(np.real(lambdas), np.imag(lambdas), marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='#bc7201')

    plt.subplot(f2)
    PlotBackground( w, r, np.angle(s_masked) /np.pi,    r'phase error $\theta/\pi$',wi,wr,cfl_a,cfl_d )
    plt.plot(np.real(lambdasS),np.imag(lambdasS),marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='0.5')
    plt.plot(np.real(lambdasC),np.imag(lambdasC),marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='#6a2202')
    plt.plot(np.real(lambdas), np.imag(lambdas), marker='o',markersize=5.,markeredgewidth=0.,lw=0.,color='#bc7201')

    plt.tight_layout(pad=0.0)
    plt.savefig("{}.pdf".format(tag+str(fig_id)),bbox_inches='tight')

    plt.show()

###############################################################################
tag='convergence'
if tag in FiguresToPlot:
    fig_id = 0

    # Define funtions
    def exp_(x):
        c = 1.
        return (                        # Parenthesis to write comments at end of line
                'Exp',                  # Name
                np.exp(c *x),           # Function
                np.exp(c *x) *c,        # First-order derivative
                np.exp(c *x) *c **2.    # Second-order derivative
                )

    def sin_(x):
        c = 1.
        return 'Sin', \
               np.sin(c *np.pi *x), \
               np.cos(c *np.pi *x) *c *np.pi, \
              -np.sin(c *np.pi *x) *(c *np.pi) **2.

    def cos_(x):
        c = 1.
        return 'Cos', \
               np.cos(c *np.pi *x), \
              -np.sin(c *np.pi *x) *c *np.pi, \
              -np.cos(c *np.pi *x) *(c *np.pi) **2.

    def gauss_(x):
        c = 40.
        tmp = np.exp(-c *x *x)
        return 'Gaussian', \
               tmp, \
              -tmp *2. *c *x, \
              -tmp *2. *c *( 1. -2. *c *x *x )

    # Error analysis of the second-order derivative
    n = 10
    x = np.linspace(-1.,1.,n)
    h = (x[n-1]-x[0]) /(n-1)

    # Define list of functions to be processed
    fs = [ exp_(x) ]#, sin_(x), cos_(x), gauss_(x)]

    # Calculate FD approximation to the second-order derivative
    l2 = fdm.fdm2_c6_A(n, np.pi **2.)
    # l2 = fdm.fdm2_c6_A(n, 48./7.)
    A2 = np.diagflat(l2[0,1:],-1) +np.diagflat(l2[1,:],0) +np.diagflat(l2[2,:-1],1)

    r2 = fdm.fdm2_c6_B(n, np.pi **2.)
    # r2 = fdm.fdm2_c6_B(n, 48./7.)
    B2 = np.diagflat(r2[0][3:],-3) +np.diagflat(r2[1][2:],-2) +np.diagflat(r2[2][1:],-1) +np.diagflat(r2[3],0) \
        +np.diagflat(r2[4][:-1],1) +np.diagflat(r2[5][:-2],2) +np.diagflat(r2[6][:-3],3)

    fdm2s = [ scipy.linalg.solve(A2,B2@f[1]) /h**2. for f in fs ]

    # e = [ A2@(fdm2s[i]-fs[i][3]) for i in range(len(fdm2s)) ]
    e = [ fdm2s[i]-fs[i][3] for i in range(len(fdm2s)) ]

    # Plot result
    plt.figure( figsize = (4,3))

    for i in range(len(fdm2s)):
        # plt.plot(x, f[3], label=f[0])
        plt.plot(x, e[i], label=fs[i][0])

    plt.title("Second-order derivative")
    plt.xlabel("$x$")
    plt.ylabel("$d^2f/dx^2$")
    plt.legend(loc="best")
    plt.show()

    # # Convergence study: we increment the number of grid points n by factors of 2
    # # between 2**imin and 2**imax

    # h   = []
    # e1s = []
    # e2s = []

    # for n in [ 2**i for i in range(4,11)]:
    #     x = np.linspace(-1.,1.,n)
    #     h.append( (x[n-1]-x[0]) /(n-1) )

    #     # Define list of functions to be processed
    #     fs = [ exp_(x) ]#, sin_(x), cos_(x), gauss_(x)]

    #     # Calculate FD approximation to the first-order derivative and error
    #     l1 = fdm.fdm1_c6_A(n)
    #     A1 = np.diagflat(l1[0,1:],-1) +np.diagflat(l1[1,:],0) +np.diagflat(l1[2,:-1],1)

    #     r1 = fdm.fdm1_c6_B(n)
    #     B1 = np.diagflat(r1[0][2:],-2) +np.diagflat(r1[1][1:],-1) +np.diagflat(r1[2],0) \
    #         +np.diagflat(r1[3][:-1],1) +np.diagflat(r1[4][:-2],2)

    #     fdm1s = [ scipy.linalg.solve(A1,B1@f[1]) /h[-1] for f in fs ]

    #     # e1s.append( [ scipy.linalg.norm(fdm1s[i]-fs[i][2]) /np.sqrt(float(n)) for i in range(len(fdm1s)) ] )
    #     e1s.append( [ np.amax(np.abs(fdm1s[i]-fs[i][2])) for i in range(len(fdm1s)) ] )

    #     # Calculate FD approximation to the second-order derivative
    #     l2 = fdm.fdm2_c6_A(n, np.pi **2.)
    #     # l2 = fdm.fdm2_c6_A(n, 48./7.)
    #     A2 = np.diagflat(l2[0,1:],-1) +np.diagflat(l2[1,:],0) +np.diagflat(l2[2,:-1],1)

    #     r2 = fdm.fdm2_c6_B(n, np.pi **2.)
    #     # r2 = fdm.fdm2_c6_B(n, 48./7.)
    #     B2 = np.diagflat(r2[0][3:],-3) +np.diagflat(r2[1][2:],-2) +np.diagflat(r2[2][1:],-1) +np.diagflat(r2[3],0) \
    #         +np.diagflat(r2[4][:-1],1) +np.diagflat(r2[5][:-2],2) +np.diagflat(r2[6][:-3],3)

    #     fdm2s = [ scipy.linalg.solve(A2,B2@f[1]) /h[-1]**2. for f in fs ]

    #     # e2s.append( [ scipy.linalg.norm(fdm2s[i]-fs[i][3]) /np.sqrt(float(n)) for i in range(len(fdm2s)) ] )
    #     e2s.append( [ np.amax(np.abs(fdm2s[i]-fs[i][3])) for i in range(len(fdm2s)) ] )

    #     legends = [ f[0] for f in fs ]

    # h_ = np.array( h )          # We need arrays to plot

    # # e1s_ = np.array( e1s )      # We need arrays to plot
    # # plt.figure( figsize = (4,3))
    # # for i in range(np.shape(e1s_)[1]):
    # #     plt.plot( h_/h_[0], e1s_[:,i] )#/e1s_[0,i] )
    # # plt.legend(legends,loc='best')
    # # plt.xscale("log")
    # # plt.yscale("log")
    # # plt.title("First-order derivative")
    # # plt.xlabel("Grid spacing $h/h_0$")
    # # plt.ylabel("Global error $e_2/e_{2,0}$")
    # # #plt.ylabel("Global error $e_\infty/e_{\infty,0}$")
    # # plt.show()

    # e2s_ = np.array( e2s )      # We need arrays to plot
    # plt.figure( figsize = (4,3))
    # for i in range(np.shape(e2s_)[1]):
    #     plt.plot( h_/h_[0], e2s_[:,i] )#/e2s_[0,i] )
    # plt.legend(legends,loc='best')
    # plt.xscale("log")
    # plt.yscale("log")
    # # plt.axis([None,None,1e-10,1e0])
    # plt.title(r"Second-order derivative")
    # plt.xlabel(r"Grid spacing $h/h_0$")
    # plt.ylabel(r"Global error $e_2/e_{2,0}$")
    # #plt.ylabel(r"Global error $e_\infty/e_{\infty,0}$")
    # plt.show()

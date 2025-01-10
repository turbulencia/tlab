
#ifndef DNS_CONST_H_INCLUDED
#define DNS_CONST_H_INCLUDED

! write format of real numbers
#define G_FORMAT_R  E13.5E3

! File formats
#define IO_MPIIO                1
#define IO_NETCDF               2
#define IO_NOFILE               3

#define IO_TYPE_SINGLE          1
#define IO_TYPE_DOUBLE          2

! Flow Types
#define DNS_FLOW_SHEAR          1
#define DNS_FLOW_JET            2
#define DNS_FLOW_VORTEX         3
#define DNS_FLOW_ISOTROPIC      4

! Flow Mode
#define DNS_MODE_TEMPORAL       1
#define DNS_MODE_SPATIAL        2

! Equations mode
#define DNS_EQNS_TOTAL              0
#define DNS_EQNS_INTERNAL           1
#define DNS_EQNS_INCOMPRESSIBLE     2
#define DNS_EQNS_ANELASTIC          3

! Equation terms
#define EQNS_NONE                   0

#define EQNS_DIVERGENCE             1
#define EQNS_SKEWSYMMETRIC          2
#define EQNS_CONVECTIVE             3
#define EQNS_EXPLICIT               4

#define EQNS_BOD_HOMOGENEOUS        5
#define EQNS_BOD_LINEAR             6
#define EQNS_BOD_BILINEAR           7
#define EQNS_BOD_QUADRATIC          8
#define EQNS_BOD_NORMALIZEDMEAN     9
#define EQNS_BOD_SUBTRACTMEAN      10

#define EQNS_COR_NORMALIZED        12

#define EQNS_RHS_SPLIT             18
#define EQNS_RHS_COMBINED          19
#define EQNS_RHS_NONBLOCKING       20

#define EQNS_TRANS_POWERLAW             21
#define EQNS_TRANS_SUTHERLAND           22

#define EQNS_RAND_MULTIPLY              31 

! Finite-differences method
#define FDM_COM4_JACOBIAN           4
#define FDM_COM6_JACOBIAN_PENTA     5
#define FDM_COM6_JACOBIAN           6
#define FDM_COM6_JACOBIAN_HYPER     7
#define FDM_COM8_JACOBIAN           8

#define FDM_COM6_DIRECT      16
#define FDM_COM4_DIRECT      17

! Operators
#define OPR_P1                1
#define OPR_P2                2
#define OPR_P2_P1             3
#define OPR_P1_BCS            4
#define OPR_P1_INT_VP         5
#define OPR_P1_INT_PV         6
#define OPR_P0_INT_VP         7
#define OPR_P0_INT_PV         8
#define OPR_P0_IBM            9

! Runge-Kutta method
#define RKM_EXP3              3
#define RKM_EXP4              4
#define RKM_IMP3_DIFFUSION    5
#define RKM_IMP3_SOURCE       6
#define RKM_IMP3_DIFFSOURCE   7

! Boundary conditions
#define DNS_BCS_PERIODIC    0

#define DNS_BCS_NONE        0
#define DNS_BCS_NR          1
#define DNS_BCS_INFLOW      2
#define DNS_BCS_DIRICHLET   3
#define DNS_BCS_NEUMANN     4

! Surface Models
#define DNS_SFC_STATIC      0
#define DNS_SFC_LINEAR      1

! Buffer
#define DNS_BUFFER_NONE     0
#define DNS_BUFFER_RELAX    1
#define DNS_BUFFER_FILTER   2
#define DNS_BUFFER_BOTH     3

! Filters
#define DNS_FILTER_NONE      0
#define DNS_FILTER_COMPACT   1
#define DNS_FILTER_6E        2
#define DNS_FILTER_4E        3
#define DNS_FILTER_ADM       4
#define DNS_FILTER_HELMHOLTZ 5
#define DNS_FILTER_BAND      6
#define DNS_FILTER_ERF       7
#define DNS_FILTER_TOPHAT    8
#define DNS_FILTER_COMPACT_CUTOFF    9

#define DNS_FILTER_BCS_PERIODIC  0
#define DNS_FILTER_BCS_BIASED    1
#define DNS_FILTER_BCS_FREE      2
#define DNS_FILTER_BCS_SOLID     3
#define DNS_FILTER_BCS_DIRICHLET 4
#define DNS_FILTER_BCS_NEUMANN   5
#define DNS_FILTER_BCS_ZERO      6

! Mixture Type
#define MIXT_TYPE_NONE             0
#define MIXT_TYPE_BS               1
#define MIXT_TYPE_PETERS1991       2
#define MIXT_TYPE_PETERS1988       3
#define MIXT_TYPE_UNIDECOMP        4
#define MIXT_TYPE_BSZELDOVICH      5
#define MIXT_TYPE_ONESTEP          6
#define MIXT_TYPE_BILGER1997       7
#define MIXT_TYPE_QUASIBS          8
#define MIXT_TYPE_AIR              9
#define MIXT_TYPE_AIRVAPOR        10
#define MIXT_TYPE_AIRWATER        11
#define MIXT_TYPE_AIRWATER_LINEAR 12
#define MIXT_TYPE_CHEMKIN         13

! Cubic Splines
#define CS_BCS_PERIODIC 0
#define CS_BCS_CLAMPED  1
#define CS_BCS_FIXED_1  2
#define CS_BCS_NATURAL  3
#define CS_BCS_FIXED_2  4

! Obs log-file type
#define OBS_TYPE_NONE  0
#define OBS_TYPE_EKMAN 1

! Pressure Decomposition
#define DCMP_TOTAL      0
#define DCMP_RESOLVED   1
#define DCMP_ADVECTION  2
#define DCMP_ADVDIFF    3
#define DCMP_DIFFUSION  4
#define DCMP_CORIOLIS   5
#define DCMP_BUOYANCY   6

#endif
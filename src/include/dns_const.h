
#ifndef DNS_CONST_H_INCLUDED
#define DNS_CONST_H_INCLUDED

! File formats
#define DNS_FILE_RAWARRAY   1
#define DNS_FILE_RAWSPLIT   2
#define DNS_FILE_NETCDF     3
#define DNS_NOFILE          4

! Flow Types
#define DNS_FLOW_SHEAR       1
#define DNS_FLOW_JET         2
#define DNS_FLOW_VORTEX      3
#define DNS_FLOW_ISOTROPIC   4

! Flow Mode
#define DNS_MODE_TEMPORAL   1
#define DNS_MODE_SPATIAL    2

! Equations mode
#define DNS_EQNS_TOTAL            0
#define DNS_EQNS_INTERNAL         1
#define DNS_EQNS_INCOMPRESSIBLE   2
#define DNS_EQNS_ANELASTIC        3

! Equation terms
#define EQNS_NONE             0

#define EQNS_DIVERGENCE       1
#define EQNS_SKEWSYMMETRIC    2
#define EQNS_CONVECTIVE       3
#define EQNS_EXPLICIT         4

#define EQNS_BOD_HOMOGENEOUS  5
#define EQNS_BOD_LINEAR       6
#define EQNS_BOD_BILINEAR     7
#define EQNS_BOD_QUADRATIC    8

#define EQNS_COR_NORMALIZED        12

#define EQNS_RAD_BULK1D_LOCAL      13
#define EQNS_RAD_BULK1D_GLOBAL     14

#define EQNS_RHS_SPLIT             18
#define EQNS_RHS_COMBINED          19
#define EQNS_RHS_NONBLOCKING       20

#define EQNS_TRANS_POWERLAW        21
#define EQNS_TRANS_SUTHERLAND      22
#define EQNS_TRANS_AIRWATER            23
#define EQNS_TRANS_AIRWATERSIMPLIFIED  24

#define EQNS_CHEM_QUADRATIC          25
#define EQNS_CHEM_QUADRATIC3         26
#define EQNS_CHEM_LAYEREDRELAXATION  27
#define EQNS_CHEM_OZONE              28

#define EQNS_SUB_CONSTANT_LOCAL      29
#define EQNS_SUB_CONSTANT_GLOBAL     30

! Finite-differences method
#define FDM_COM4_JACOBIAN     4
#define FDM_COM6_JACPENTA     5
#define FDM_COM6_JACOBIAN     6
#define FDM_COM8_JACOBIAN     8

#define FDM_COM6_DIRECT      16

! Operators
#define OPR_P1                1
#define OPR_P2                2
#define OPR_P2_P1             3
#define OPR_P1_BCS            4

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

#define DNS_FILTER_BCS_PERIODIC  0
#define DNS_FILTER_BCS_BIASED    1
#define DNS_FILTER_BCS_FREE      2
#define DNS_FILTER_BCS_SOLID     3
#define DNS_FILTER_BCS_DIRICHLET 4
#define DNS_FILTER_BCS_NEUMANN   5
#define DNS_FILTER_BCS_ZERO      6

! Mean profiles
#define PROFILE_NONE         0
#define PROFILE_LINEAR       1
#define PROFILE_TANH         2
#define PROFILE_ERF          3
#define PROFILE_BICKLEY      4
#define PROFILE_GAUSSIAN     5
#define PROFILE_LINEAR_ERF   6
#define PROFILE_EKMAN_U      7
#define PROFILE_EKMAN_V      8
#define PROFILE_EKMAN_U_P    9
#define PROFILE_PARABOLIC   10
#define PROFILE_LINEAR_CROP 11
#define PROFILE_MIXEDLAYER  12
#define PROFILE_ERF_ANTISYM 13
#define PROFILE_ERF_SURFACE 14
#define PROFILE_LINEAR_ERF_SURFACE 15
#define PROFILE_PARABOLIC_SURFACE  16
#define PROFILE_GAUSSIAN_SURFACE   17
#define PROFILE_GAUSSIAN_ANTISYM   18
#define PROFILE_GAUSSIAN_SYM       19
#define PROFILE_TANH_ANTISYM       20
#define PROFILE_TANH_SYM           21

! Chemistry Constants
#define CHEM_NONE           0
#define CHEM_INFINITE       1
#define CHEM_FINITE         2

#define CHEM_NONPREMIXED    0
#define CHEM_PREMIXED       1

! Chemistry Type
#define CHEM_TYPE_BS          1
#define CHEM_TYPE_PETERS1991  2
#define CHEM_TYPE_PETERS1988  3
#define CHEM_TYPE_UNIDECOMP   4
#define CHEM_TYPE_BSZELDOVICH 5
#define CHEM_TYPE_ONESTEP     6
#define CHEM_TYPE_BILGER1997  7
#define CHEM_TYPE_QUASIBS     8

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

! Information for IO subarrays
#define IO_SUBARRAY_VISUALS_XOY  1
#define IO_SUBARRAY_VISUALS_ZOY  2
#define IO_SUBARRAY_VISUALS_XOZ  3
#define IO_SUBARRAY_PLANES_XOY   4
#define IO_SUBARRAY_PLANES_ZOY   5
#define IO_SUBARRAY_PLANES_XOZ   6
#define IO_SUBARRAY_BUFFER_ZOY   7
#define IO_SUBARRAY_BUFFER_XOZ   8
#define IO_SUBARRAY_SPECTRA_X    9
#define IO_SUBARRAY_SPECTRA_Z   10
#define IO_SUBARRAY_SPECTRA_XZ  11
#define IO_SUBARRAY_ENVELOPES   12

#define IO_SUBARRAY_SIZE        12

! Lagrangian Type
#define LAG_TYPE_NONE          0
#define LAG_TYPE_TRACER        1
#define LAG_TYPE_SIMPLE_SETT   2
#define LAG_TYPE_BIL_CLOUD_3   3
#define LAG_TYPE_BIL_CLOUD_4   4

! Lagrangian Trajectories
#define LAG_TRAJECTORY_NONE      0
#define LAG_TRAJECTORY_FIRST     1
#define LAG_TRAJECTORY_LARGEST   2
#define LAG_TRAJECTORY_VORTICITY 3

#endif

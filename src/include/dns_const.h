
#ifndef DNS_CONST_H_INCLUDED
#define DNS_CONST_H_INCLUDED

! write format of real numbers
#define G_FORMAT_R  E13.5E3

! Flow Types
#define DNS_FLOW_SHEAR          1
#define DNS_FLOW_JET            2
#define DNS_FLOW_VORTEX         3
#define DNS_FLOW_ISOTROPIC      4

! Flow Mode
#define DNS_MODE_TEMPORAL       1
#define DNS_MODE_SPATIAL        2

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
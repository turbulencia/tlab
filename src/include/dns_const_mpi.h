
#ifndef DNS_CONST_MPI_H_INCLUDED
#define DNS_CONST_MPI_H_INCLUDED

! MPI structures id
#define DNS_MPI_K_PARTIAL     1
#define DNS_MPI_K_SHEAR       2
#define DNS_MPI_K_POISSON     3
#define DNS_MPI_K_INFLOW      4
#define DNS_MPI_K_OUTBCS      5
#define DNS_MPI_K_TOPBCS      6
#define DNS_MPI_K_NRBCX       7
#define DNS_MPI_K_NRBCY       8
#define DNS_MPI_K_AUX1        9
#define DNS_MPI_K_AUX2       10

#define DNS_MPI_K_IBM_NOB    11    
#define DNS_MPI_K_IBM_NOB_BE 12    

#define DNS_MPI_K_MAXTYPES   12

#define DNS_MPI_I_PARTIAL     1
#define DNS_MPI_I_SHEAR       2
#define DNS_MPI_I_POISSON1    3
#define DNS_MPI_I_POISSON2    4
#define DNS_MPI_I_AUX1        5
#define DNS_MPI_I_AUX2        6

#define DNS_MPI_I_IBM_NOB     7    
#define DNS_MPI_I_IBM_NOB_BE  8    

#define DNS_MPI_I_MAXTYPES    8

#define DNS_MPI_J_PARTIAL     1 
#define DNS_MPI_J_MAXTYPES    1

! Control of MPI Transpositions
#define TLAB_MPI_TRP_NONE         0
#define TLAB_MPI_TRP_ASYNCHRONOUS 1
#define TLAB_MPI_TRP_SENDRECV     2

#endif

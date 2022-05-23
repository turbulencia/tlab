
#ifndef DNS_CONST_MPI_H_INCLUDED
#define DNS_CONST_MPI_H_INCLUDED

! MPI ERROR HANDLING
#define DO_MPI_ERROR_CHECK  IF ( ims_err /= MPI_SUCCESS ) CALL TLAB_MPI_PANIC(__FILE__,ims_err)

! MPI structures id
#define TLAB_MPI_K_PARTIAL   1
#define TLAB_MPI_K_SHEAR     2
#define TLAB_MPI_K_POISSON   3
#define TLAB_MPI_K_INFLOW    4
#define TLAB_MPI_K_OUTBCS    5
#define TLAB_MPI_K_TOPBCS    6
#define TLAB_MPI_K_NRBCX     7
#define TLAB_MPI_K_NRBCY     8
#define TLAB_MPI_K_AUX1      9
#define TLAB_MPI_K_AUX2     10

#define TLAB_MPI_K_MAXTYPES 10

#define TLAB_MPI_I_PARTIAL   1
#define TLAB_MPI_I_SHEAR     2
#define TLAB_MPI_I_POISSON1  3
#define TLAB_MPI_I_POISSON2  4
#define TLAB_MPI_I_AUX1      5
#define TLAB_MPI_I_AUX2      6

#define TLAB_MPI_I_MAXTYPES  6

! Control of MPI Transpositions
#define TLAB_MPI_TRP_NONE         0
#define TLAB_MPI_TRP_ASYNCHRONOUS 1
#define TLAB_MPI_TRP_SENDRECV     2

#endif

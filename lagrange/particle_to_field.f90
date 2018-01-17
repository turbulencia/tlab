#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# DESCRIPTION
!#
!# Interpolate particle information to a field
!#
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD(l_q, particle_property, field_out, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY: imax,jmax,kmax, isize_particle
#ifdef USE_MPI
   USE DNS_MPI, ONLY: ims_err
#endif

   IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(isize_particle,3)    :: l_q
  TREAL, DIMENSION(isize_particle)      :: particle_property
  TREAL, DIMENSION(imax,jmax,kmax)      :: field_out
  TREAL, DIMENSION(*)                   :: wrk2d
  TREAL, DIMENSION(imax+1,jmax,kmax+1)  :: wrk3d

!#######################################################################
!Write data to field
!Field is a extended field with gird, halo1, halo2 and halo3
!#######################################################################
  wrk3d = C_0_R
  CALL PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, wrk3d)

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!#######################################################################
!SEND to the East and RECV from the West
!Sum the most left colum of field
!SEND to the North and RECV from the South
!Sum the lowest row of field
!#######################################################################
 CALL PARTICLE_TO_FIELD_SEND_RECV_EAST(wrk2d(1),wrk2d((jmax*(kmax+1))+1),wrk3d)
 CALL PARTICLE_TO_FIELD_SEND_RECV_NORTH(wrk2d(1),wrk2d(((imax+1)*jmax)+1),wrk3d)

#endif
 
!#######################################################################
!Put wrk3d into continuous memory field_out
!#######################################################################
  field_out(1:imax,1:jmax,1:kmax)=wrk3d(1:imax,1:jmax,1:kmax)

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD

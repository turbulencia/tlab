#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# Interpolate particle information to a field
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE  PARTICLE_TO_FIELD &
    (l_q,particle_property, wrk1d,wrk2d,wrk3d, field_out)

  USE DNS_GLOBAL, ONLY: imax,jmax,kmax, isize_particle
  USE DNS_GLOBAL, ONLY: g
  USE LAGRANGE_GLOBAL, ONLY: jmin_part
#ifdef USE_MPI
   USE DNS_MPI, ONLY: ims_npro, ims_pro, ims_err
#endif

   IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(imax,jmax,kmax)      ::field_out ! will be txc
  TREAL, DIMENSION(isize_particle,3)    ::l_q
  TREAL, DIMENSION(imax+1,jmax,kmax+1)  :: wrk3d
  TREAL, DIMENSION(isize_particle)      :: particle_property ! will be particle_txc
  TREAL, DIMENSION(*)                   ::wrk1d, wrk2d

#ifdef USE_MPI
  
!#######################################################################
!Write data to field
!Field is a extended field with gird, halo1, halo2 and halo3
!#######################################################################
  wrk3d=C_0_R
  wrk1d(1)= g(1)%scale/g(1)%size ! wrk1d 1-3 intervalls
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale/g(3)%size

  CALL RHS_PARTICLE_TO_FIELD(l_q, particle_property, wrk1d, wrk3d)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!#######################################################################
!SEND to the East and RECV from the West
!Sum the most left colum of field
!SEND to the North and RECV from the South
!Sum the lowest row of field
!#######################################################################
 

 CALL PARTICLE_TO_FIELD_SEND_RECV_EAST(wrk2d(1),wrk2d((jmax*(kmax+1))+1),wrk3d)
 CALL PARTICLE_TO_FIELD_SEND_RECV_NORTH(wrk2d(1),wrk2d(((imax+1)*jmax)+1),wrk3d)

!#######################################################################
!Put wrk3d into continious memory field_out
!#######################################################################
  field_out=C_0_R

  field_out(1:imax,1:jmax,1:kmax)=wrk3d(1:imax,1:jmax,1:kmax)

#else


!#######################################################################
!Write data to field
!Field is a extended field with gird, halo1, halo2 and halo3
!#######################################################################
  wrk3d=C_0_R
  wrk1d(1)= g(1)%scale/g(1)%size ! wrk1d 1-3 intervalls
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale/g(3)%size
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation

  CALL RHS_PARTICLE_TO_FIELD(l_q, particle_property, wrk1d, wrk3d)

!#######################################################################
!Put wrk3d into continious memory field_out
!#######################################################################
  field_out=C_0_R

  field_out(1:imax,1:jmax,1:kmax)=wrk3d(1:imax,1:jmax,1:kmax)

#endif

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
SUBROUTINE  FIELD_TO_PARTICLE &
    (nvar, data_in, npar, data_out, l_q,l_hq,l_tags,l_comm, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax, isize_particle
  USE DNS_GLOBAL,     ONLY : g
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI,        ONLY:  ims_err
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nvar, npar
  TYPE(pointers3d_dt), DIMENSION(nvar)               :: data_in
  TYPE(pointers_dt), DIMENSION(nvar)                 :: data_out
  TREAL,             DIMENSION(isize_particle,*)     :: l_q, l_hq
  INTEGER(8),        DIMENSION(isize_particle)       :: l_tags
  TREAL,             DIMENSION(isize_l_comm), TARGET :: l_comm
  TREAL,             DIMENSION(*)                    :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal 
  TINTEGER npar_start
  TINTEGER ip1,ip2,ip3, np1,np2,np3, iv

  TYPE(pointers3d_dt), DIMENSION(nvar) :: data_halo1, data_halo2, data_halo3

!#######################################################################
  IF ( nvar .GT. inb_lag_total_interp ) THEN
     CALL IO_WRITE_ASCII(efile,'FIELD_TO_PARTICLE. Not enough memory.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  np1 = 2*jmax*kmax; ip1 = 1
  np2 = imax*jmax*2; ip2 = 1 +isize_hf_1
  np3 = 2   *jmax*2; ip3 = 1 +isize_hf_1 +isize_hf_2
  DO iv = 1,nvar
     data_halo1(iv)%field(1:2,1:jmax,1:kmax) => l_comm(ip1:ip1+np1-1); ip1 = ip1 +np1
     data_halo2(iv)%field(1:imax,1:jmax,1:2) => l_comm(ip2:ip2+np2-1); ip2 = ip2 +np2
     data_halo3(iv)%field(1:2,1:jmax,1:2)    => l_comm(ip3:ip3+np3-1); ip3 = ip3 +np3
  ENDDO

! #######################################################################
! Setting fields in halo regions
! ######################################################################
#ifdef USE_MPI
  CALL PARTICLE_HALO_K(nvar, data_in, data_halo2(1)%field, wrk3d(1), wrk3d(imax*jmax*nvar+1), &
       wrk2d(1), wrk2d(2*(jmax*nvar+1)))

  CALL PARTICLE_HALO_I(nvar, data_in, data_halo1(1)%field, data_halo3(1)%field, wrk3d(1), wrk3d(jmax*(kmax+1)*nvar+1),&
       wrk2d(1), wrk2d(jmax*nvar+1), wrk2d(2*(jmax*nvar+1)))

#else
  CALL PARTICLE_HALO_SERIAL(nvar, data_in, data_halo1(1)%field, data_halo2(1)%field, data_halo3(1)%field)

#endif
  
!#######################################################################
! Sorting and counting particles for each zone
!#######################################################################
  CALL PARTICLE_SORT_HALO(grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal,&
       l_hq, l_tags, l_q)

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
#endif
  
!#######################################################################
! Interpolating
!#######################################################################
  wrk1d(1)= g(1)%scale / g(1)%size ! wrk1d 1-3 intervalls
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale / g(3)%size

  npar_start = 1
  npar       = grid_zone
  CALL PARTICLE_INTERPOLATION(i0, imax,jmax,kmax,nvar, data_in, data_out, l_q, g(2)%nodes, wrk1d, npar_start, npar)

  IF ( halo_zone_x .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_x
     CALL PARTICLE_INTERPOLATION(i1, i2,jmax,kmax,nvar, data_halo1, data_out, l_q, g(2)%nodes, wrk1d, npar_start,npar)
  END IF
  
  IF ( halo_zone_z .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_z
     CALL PARTICLE_INTERPOLATION(i2, imax,jmax,i2,nvar, data_halo2, data_out, l_q, g(2)%nodes, wrk1d, npar_start,npar)
 END IF

  IF ( halo_zone_diagonal .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_diagonal
     CALL PARTICLE_INTERPOLATION(i3, i2,jmax,i2,nvar, data_halo3, data_out, l_q, g(2)%nodes, wrk1d, npar_start,npar)
  END IF

  DO iv = 1,nvar
     NULLIFY(data_halo1(iv)%field,data_halo2(iv)%field,data_halo3(iv)%field)
  ENDDO
  
  RETURN
END SUBROUTINE FIELD_TO_PARTICLE

!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# Interpolate information from a field to particles
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE  FIELD_TO_PARTICLE_OLD &
    (field, wrk1d,wrk2d,wrk3d, particle_property, l_tags, l_hq, l_q)

  USE DNS_GLOBAL, ONLY: imax,jmax,kmax
  USE DNS_GLOBAL, ONLY: g
  USE DNS_GLOBAL, ONLY: isize_particle, inb_particle
  USE LAGRANGE_GLOBAL, ONLY: jmin_part
#ifdef USE_MPI
  USE DNS_MPI, ONLY: ims_npro, ims_err, ims_pro
#endif

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
  TREAL, DIMENSION(imax,jmax,kmax)                :: field 
  TREAL, DIMENSION(isize_particle,inb_particle)   :: l_q, l_hq
  TREAL, DIMENSION(*)                             :: wrk3d
  TREAL, DIMENSION(isize_particle)                :: particle_property ! will be particle_txc
  TREAL, DIMENSION(*)                             :: wrk1d, wrk2d
  TREAL, DIMENSION(:,:,:), ALLOCATABLE            :: halo_field_1
  TREAL, DIMENSION(:,:,:), ALLOCATABLE            :: halo_field_2
  TREAL, DIMENSION(:,:,:), ALLOCATABLE            :: halo_field_3

  INTEGER(8), DIMENSION(isize_particle)           :: l_tags
  TINTEGER halo_zone_x, halo_zone_z, halo_zone_diagonal 
  TINTEGER  halo_start, halo_end,grid_start ,grid_field_counter
  ALLOCATE(halo_field_1(2,jmax,kmax))
  ALLOCATE(halo_field_2(imax,jmax,2))
  ALLOCATE(halo_field_3(2,jmax,2))  

  particle_property = C_0_R

#ifdef USE_MPI
!#######################################################################
!#######################################################################
  CALL HALO_PLANE_SHIFTING_k_1d(field,halo_field_2,wrk3d(1),wrk3d(imax*jmax+1),wrk2d(1),wrk2d(2*(jmax+1)))
  CALL HALO_PLANE_SHIFTING_i_1d(field,halo_field_1,halo_field_3,wrk3d(1),wrk3d((jmax)*(kmax+1)+1),&
       wrk2d(1),wrk2d(jmax+1),wrk2d(2*(jmax+1)))

!#######################################################################
!#######################################################################
  CALL PARTICLE_SORT_HALO(grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
       l_hq, l_tags, l_q)

!#######################################################################
  wrk1d(1)= g(1)%scale / g(1)%size ! wrk1d 1-3 intervalls
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale / g(3)%size

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_1D(field,l_q,particle_property,g(2)%nodes,wrk1d,grid_start,grid_field_counter)
! RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
     halo_start=grid_field_counter+1
     halo_end=grid_field_counter+halo_zone_x
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1_1D(halo_field_1,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
! RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2_1D(halo_field_2,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
! RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3_1D(halo_field_3,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF

#else

  halo_zone_x=0
  CALL PARTICLE_SORT_HALO(grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
       l_hq, l_tags, l_q)

  wrk1d(1)= g(1)%scale / g(1)%size ! wrk1d 1-3 intervalls
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale / g(3)%size

!#######################################################################
! RHS for particles in normal field
!#######################################################################
  halo_field_1(1,1:jmax,1:kmax)=field(g(1)%size,1:jmax,1:kmax)  
  halo_field_1(2,1:jmax,1:kmax)=field(1,1:jmax,1:kmax)  

  halo_field_2(1:imax,1:jmax,1)=field(1:imax,1:jmax,g(3)%size)
  halo_field_2(1:imax,1:jmax,2)=field(1:imax,1:jmax,1)

  halo_field_3(1,1:jmax,1)=field(g(1)%size,1:jmax,g(3)%size)
  halo_field_3(1,1:jmax,2)=field(g(1)%size,1:jmax,1)
  halo_field_3(2,1:jmax,1)=field(1,1:jmax,g(3)%size)
  halo_field_3(2,1:jmax,2)=field(1,1:jmax,1)

  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_1D(field,l_q,particle_property,g(2)%nodes,wrk1d,grid_start,grid_field_counter)
! RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
     halo_start=grid_field_counter+1
     halo_end=grid_field_counter+halo_zone_x
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1_1D(halo_field_1,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
! RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2_1D(halo_field_2,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
! RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3_1D(halo_field_3,l_q,particle_property,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
#endif

  DEALLOCATE(halo_field_1)
  DEALLOCATE(halo_field_2)
  DEALLOCATE(halo_field_3)

  RETURN
END SUBROUTINE FIELD_TO_PARTICLE_OLD

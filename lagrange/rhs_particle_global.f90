#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
SUBROUTINE RHS_PARTICLE_GLOBAL(q,s, txc, l_q,l_hq,l_txc,l_tags,l_comm, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES,  ONLY : pointers_dt
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_particle
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : visc, radiation
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : thermo_param
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TREAL,      DIMENSION(isize_field,*),    TARGET :: q, s, txc
  TREAL,      DIMENSION(isize_particle,*), TARGET :: l_q, l_hq, l_txc
  INTEGER(8), DIMENSION(isize_particle)           :: l_tags
  TREAL,      DIMENSION(isize_l_comm),     TARGET :: l_comm
  TREAL,      DIMENSION(*)                        :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TREAL dummy, dummy2
  TINTEGER particle_number_local, bcs(2,2), nvar
  TINTEGER i, npar
  TREAL delta_inv0, delta_inv2, delta_inv4
  
  TYPE(pointers_dt), DIMENSION(inb_lag_total_interp) :: data, data_out

! #####################################################################
  bcs = 0
  
! Setting pointers to velocity fields
  nvar = 0
  nvar = nvar+1; data(nvar)%field => q(:,1); data_out(nvar)%field => l_hq(:,1)
  nvar = nvar+1; data(nvar)%field => q(:,2); data_out(nvar)%field => l_hq(:,2)
  nvar = nvar+1; data(nvar)%field => q(:,3); data_out(nvar)%field => l_hq(:,3)
  
! -------------------------------------------------------------------
! Additional terms depending on type of particle evolution equations
! -------------------------------------------------------------------
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN

     dummy2 = -thermo_param(2)
     dummy  = -thermo_param(1)

!LAPLACE Xi
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,1), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,1), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,1), txc(1,4), txc(1,3), wrk2d,wrk3d)

     txc(:,1) = visc *dummy *( txc(:,4) +txc(:,5) +txc(:,6) )

!LAPLACE delta
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,2), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,2), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,2), txc(1,4), txc(1,3), wrk2d,wrk3d)

     txc(:,1) = txc(:,1) +( visc *dummy2* (txc(:,4)+txc(:,5)+txc(:,6)) ) !first eq. without ds/dxi
     txc(:,2) = C_1_R - dummy*s(:,1) - dummy2*s(:,2) !xi field in txc(1,2)

     CALL FI_GRADIENT(imax,jmax,kmax, txc(1,2), txc(1,3),txc(1,4), wrk2d,wrk3d) ! square of chi gradient in txc(1,3)
     txc(:,3) = visc *txc(:,3)

     CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(1)), txc(1,4), wrk1d,wrk3d)
! Radiation *** ATTENTION RADIATION IS MINUS
     txc(:,1) = txc(:,1) + dummy2 *txc(:,4)

! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
     txc(:,4) = dummy2 *txc(:,4)

! Setting pointers
     nvar = nvar+1; data(nvar)%field => txc(:,1); data_out(nvar)%field => l_txc(:,1)
     nvar = nvar+1; data(nvar)%field => txc(:,2); data_out(nvar)%field => l_txc(:,2)
     nvar = nvar+1; data(nvar)%field => txc(:,3); data_out(nvar)%field => l_txc(:,3)
     nvar = nvar+1; data(nvar)%field => txc(:,4); data_out(nvar)%field => l_txc(:,4)
     l_txc(:,1:4) = C_0_R
     
  ENDIF

! -------------------------------------------------------------------
! Interpolating field data into particles
! The interpolated data is added to the existing data, which
!  consitutes already the evolution equation for particle position
! -------------------------------------------------------------------
  CALL FIELD_TO_PARTICLE(nvar, data, npar, data_out, l_q,l_hq,l_tags,l_comm, wrk1d,wrk2d,wrk3d)
  
! -------------------------------------------------------------------
! Completing evolution equations
! -------------------------------------------------------------------
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN
! l_txc(1) = equation without ds/dxi 
! l_txc(2) = xi
! l_txc(3) = evaporation/condensation term without d2s/dxi2 
! l_txc(4) = radiation term without ds/dxi
     
     delta_inv0 =  C_1_R  /thermo_param(1)/thermo_param(3)
     delta_inv2 = -C_05_R /thermo_param(1)/thermo_param(3)
     delta_inv4 = -C_025_R/thermo_param(1)/thermo_param(3)
     
     DO i = 1,npar
        l_hq(i,4) = l_hq(i,4) - l_txc(i,1)/(C_1_R + EXP(l_txc(i,2)*delta_inv0))
        
        l_hq(i,5) = l_hq(i,5) - l_txc(i,4)/(C_1_R + EXP(l_txc(i,2)*delta_inv0)) &
                              - l_txc(i,3)*delta_inv4/(COSH(l_txc(i,2)*delta_inv2)**2) 
     ENDDO
     
  ELSE IF (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN
     
#ifdef USE_MPI
     particle_number_local = ims_size_p(ims_pro+1)
#else
     particle_number_local = INT(particle_number)
#endif

     l_hq(1:particle_number_local,2) = l_hq(1:particle_number_local,2) - lagrange_param(1)

  ENDIF

  RETURN
END SUBROUTINE RHS_PARTICLE_GLOBAL

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2013/11 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!#  All particle interactions inside of the grid (each processor)
!#  -Creating of halo field
!#  -Particle sorting
!#  -Particle RHS-grid, RHS-halo
!#
!########################################################################
SUBROUTINE RHS_PARTICLE_GLOBAL_OLD(q,s, wrk1d,wrk2d,wrk3d,txc,l_q,l_hq, l_tags,l_comm)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : isize_particle
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_GLOBAL, ONLY : radiation
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : thermo_param
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax, *)           :: q
  TREAL, DIMENSION(isize_field, *)              :: s
  TREAL, DIMENSION(isize_particle,*)            :: l_q,l_hq
  TREAL, DIMENSION(isize_l_comm)                :: l_comm
  TREAL, DIMENSION(*)                           :: wrk1d
  TREAL, DIMENSION(*)                           :: wrk2d
  TREAL, DIMENSION(*)                           :: wrk3d
!  TREAL, DIMENSION(isize_txc_field,*)           :: txc
  TREAL, DIMENSION(isize_field,*)               :: txc
  
#ifdef USE_MPI
  TREAL, DIMENSION(:), POINTER                  :: halo_field_1
  TREAL, DIMENSION(:), POINTER                  :: halo_field_2
  TREAL, DIMENSION(:), POINTER                  :: halo_field_3
#else
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE    :: halo_field_1
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE    :: halo_field_2
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE    :: halo_field_3
  TINTEGER alloc_err
#endif

  TARGET :: l_comm

  TREAL dummy, dummy2
  TINTEGER  halo_start, halo_end,grid_start ,grid_field_counter,local_np,ij, bcs(2,2)
  INTEGER(8), DIMENSION(isize_particle)  :: l_tags
  TINTEGER halo_zone_x, halo_zone_z, halo_zone_diagonal 

#ifdef USE_MPI
#else
  ALLOCATE(halo_field_1(2,jmax,kmax,inb_lag_total_interp))
  ALLOCATE(halo_field_2(imax,jmax,2,inb_lag_total_interp))
  ALLOCATE(halo_field_3(2,jmax,2,inb_lag_total_interp))
#endif

#ifdef USE_MPI
    
  halo_field_1(1:isize_hf_1) => l_comm(1:isize_hf_1)
  halo_field_2(1:isize_hf_2) => l_comm(isize_hf_1+1:isize_hf_1+isize_hf_2)
  halo_field_3(1:isize_hf_3) => l_comm(isize_hf_1+isize_hf_2+1:isize_hf_1+isize_hf_2+isize_hf_3)

  bcs = 0
  
! #####################################################################
! Put source terms into txc variables
! #####################################################################
  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !using combination of both versions of equation

     dummy2 = -thermo_param(2)
     dummy = -thermo_param(1)

!LAPLACE Xi
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,1), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,1), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,1), txc(1,4), txc(1,3), wrk2d,wrk3d)

     DO ij = 1,isize_field
        txc(ij,1) =visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
     ENDDO

!LAPLACE delta
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,2), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,2), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,2), txc(1,4), txc(1,3), wrk2d,wrk3d)

     DO ij = 1,isize_field
        txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6))) !first eq. without ds/dxi
     ENDDO

     txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)

     CALL FI_GRADIENT(imax,jmax,kmax, txc(1,2), txc(1,3),txc(1,4), wrk2d,wrk3d) ! square of chi gradient in txc(1,3)

     DO ij = 1,isize_field
        txc(ij,3) = visc*txc(ij,3)
     ENDDO

     CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(1)), txc(1,4), wrk1d,wrk3d)
! Radiation *** ATTENTION RADIATION IS MINUS
     DO ij = 1,isize_field
        txc(ij,1) =txc(ij,1) + dummy2*txc(ij,4)
     ENDDO

! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
     DO ij = 1,isize_field
        txc(ij,4) = dummy2*txc(ij,4)
     ENDDO

  ENDIF

! ######################################################################
! Swapping grid information to halo_fields
! ######################################################################
  CALL HALO_PLANE_SHIFTING_k(q,txc,halo_field_2,wrk3d(1),wrk3d(imax*jmax*inb_lag_total_interp+1),&
       wrk2d(1),wrk2d(2*(jmax*inb_lag_total_interp+1)))

  CALL HALO_PLANE_SHIFTING_i(q,txc,halo_field_1,halo_field_3,wrk3d(1),wrk3d(jmax*(kmax+1)*inb_lag_total_interp+1),&
       wrk2d(1),wrk2d(jmax*inb_lag_total_interp+1),wrk2d(2*(jmax*inb_lag_total_interp+1)))

! #######################################################################
! Sorting of particles -> NO HALO / HALO ZONES
! #######################################################################
! Sorting algorithm:
! Counting of particles for each zone
! Zone: 0=grid 1=x_side 2=z_side 3=diagonal 
  halo_zone_x=0

  CALL PARTICLE_SORT_HALO(g(1)%nodes,g(3)%nodes, grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
       l_hq, l_tags, l_q)

  wrk1d(1)= g(1)%scale/g(1)%size ! wrk1d 1-3 intervalls
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale/g(3)%size

!#######################################################################
! RHS for particles in normal field
!#######################################################################

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION(q,l_q,l_hq,g(2)%nodes,wrk1d,txc,grid_start,grid_field_counter)
!RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
     halo_start=grid_field_counter+1
     halo_end=grid_field_counter+halo_zone_x
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1(halo_field_1,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
!RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2(halo_field_2,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
!RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3(halo_field_3,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF

#else
! #####################################################################
! Put source terms into txc variables
! #####################################################################
  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !using second version of equation

     dummy2 = -thermo_param(2)
     dummy = -thermo_param(1)

!LAPLACE Xi
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,1), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,1), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,1), txc(1,4), txc(1,3), wrk2d,wrk3d)

     DO ij = 1,isize_field
!       txc(ij,1) =txc(ij,1) + visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
        txc(ij,1) = visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
     ENDDO

!LAPLACE delta
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,2), txc(1,6), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,2), txc(1,5), txc(1,3), wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,2), txc(1,4), txc(1,3), wrk2d,wrk3d)

     DO ij = 1,isize_field
        txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6))) !first eq. without ds/dxi
     ENDDO

     txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)

     CALL FI_GRADIENT(imax,jmax,kmax, txc(1,2), txc(1,3),txc(1,4), wrk2d,wrk3d) ! square of chi gradient in txc(1,3)

     DO ij = 1,isize_field
        txc(ij,3) = visc*txc(ij,3)
     ENDDO

     CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(1)), txc(1,4), wrk1d,wrk3d)
! Radiation *** ATTENTION RADIATION IS MINUS
     DO ij = 1,isize_field
        txc(ij,1) =txc(ij,1) + dummy2*txc(ij,4)
     ENDDO

! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
     DO ij = 1,isize_field
        txc(ij,4) = dummy2*txc(ij,4)
     ENDDO

  ENDIF

  halo_zone_x=0
  CALL PARTICLE_SORT_HALO(g(1)%nodes,g(3)%nodes, grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
       l_hq, l_tags, l_q)

  wrk1d(1)= g(1)%scale/g(1)%size ! wrk1d 1-3 intervalls
!  wrk1d(2)= g(2)%scale/g(2)%size ! needed for interpolation
  wrk1d(2)= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  wrk1d(3)= g(3)%scale/g(3)%size

!#######################################################################
! RHS for particles in normal field
!#######################################################################
  CALL HALO_PLANE_SHIFTING_SERIAL(q,txc,halo_field_1,halo_field_2,halo_field_3)

  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION(q,l_q,l_hq,g(2)%nodes,wrk1d,txc,grid_start,grid_field_counter)
!RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
     halo_start=grid_field_counter+1
     halo_end=grid_field_counter+halo_zone_x
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1(halo_field_1,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
!RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2(halo_field_2,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF
!RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
     halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
     halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
     CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3(halo_field_3,l_q,l_hq,g(2)%nodes,wrk1d,halo_start,halo_end)
  END IF

  DEALLOCATE(halo_field_1, STAT=alloc_err)
  DEALLOCATE(halo_field_2, STAT=alloc_err)
  DEALLOCATE(halo_field_3, STAT=alloc_err)

#endif

  IF (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN

#ifdef USE_MPI
     local_np = ims_size_p(ims_pro+1)
#else
     local_np = particle_number
#endif

     DO ij=1,local_np
        l_hq(ij,2) = l_hq(ij,2) - lagrange_param(1)
     END DO
  ENDIF

  RETURN
END SUBROUTINE RHS_PARTICLE_GLOBAL_OLD

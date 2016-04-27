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
!# ARGUMENTS 
!#
!#
!#
!########################################################################
SUBROUTINE RHS_PARTICLE_GLOBAL( &
     x,y,z,dx,dy,dz,q,s, wrk1d,wrk2d,wrk3d,txc,l_q,l_hq,& 
     l_tags,l_comm)


#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field, imax_total,jmax_total, kmax_total
  USE DNS_GLOBAL, ONLY : isize_particle, scalex, scaley, scalez, inb_particle, inb_scal_array
  USE DNS_GLOBAL, ONLY : body_param, imode_fdm, i1bc, j1bc, k1bc, visc, isize_wrk1d
  USE DNS_GLOBAL, ONLY : iunifx,iunify,iunifz,iradiation,rad_param, inb_txc
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TREAL, DIMENSION(*)                           :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(imax,jmax,kmax, *)           :: q
  TREAL, DIMENSION(isize_field, *)              :: s
  TREAL, DIMENSION(isize_particle,*)            :: l_q,l_hq
  TREAL, DIMENSION(isize_l_comm)                :: l_comm
  TREAL, DIMENSION(*)                           :: wrk1d
  TREAL, DIMENSION(*)                           :: wrk2d
  TREAL, DIMENSION(*)                           :: wrk3d
!  TREAL, DIMENSION(isize_txc_field,*)           :: txc
  TREAL, DIMENSION(isize_field,*)           :: txc
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
  TINTEGER  halo_start, halo_end,grid_start ,grid_field_counter,local_np,ij
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

! #####################################################################
! Put source terms into txc variables
! #####################################################################
  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD) THEN
    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s

    ! CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,inb_scal_array), rad_param,&
    !     wrk1d(1:isize_wrk1d), wrk1d(isize_wrk1d+1:2*isize_wrk1d), wrk1d(2*isize_wrk1d+1:3*isize_wrk1d), wrk3d) ! Put radiation in wrk3d
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,1), wrk1d,wrk3d)
    ! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
    DO ij = 1,isize_field
       txc(ij,1) = dummy2*txc(ij,1)
    ENDDO 
    
    

    !LAPLACE Xi
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,1), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,1), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,1), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 
 
    !LAPLACE delta
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,2), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,2), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,2), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6)))
    ENDDO 

    !LAPLACE s
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
         dz, s(1,inb_scal_array), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
         dy, s(1,inb_scal_array), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
         dx, s(1,inb_scal_array), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)

    DO ij = 1,isize_field
       txc(ij,3) = visc*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 

    txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)
  
  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2) THEN !using second version of equation

    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s


    
    
    !LAPLACE s
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
         dz, s(1,inb_scal_array), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
         dy, s(1,inb_scal_array), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
         dx, s(1,inb_scal_array), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)

    DO ij = 1,isize_field
       txc(ij,3) = visc*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 


    txc(1:isize_field,4) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,4)

    CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, txc(1,4), txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d) ! square of chi gradient in txc(1,1)

    DO ij = 1,isize_field
       txc(ij,1) = visc*txc(ij,1)
    ENDDO 
  
    ! CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,inb_scal_array), rad_param,&
    !     wrk1d(1:isize_wrk1d), wrk1d(isize_wrk1d+1:2*isize_wrk1d), wrk1d(2*isize_wrk1d+1:3*isize_wrk1d), wrk3d) ! Put radiation in wrk3d
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,2), wrk1d,wrk3d)
    ! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
    DO ij = 1,isize_field
       txc(ij,2) = dummy2*txc(ij,2)
    ENDDO 

  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !using combination of both versions of equation

    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s

    !LAPLACE Xi
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,1), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,1), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,1), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 
 
    !LAPLACE delta
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,2), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,2), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,2), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6))) !first eq. without ds/dxi
    ENDDO 


   
    txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)



    CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, txc(1,2), txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d) ! square of chi gradient in txc(1,3)

    DO ij = 1,isize_field
       txc(ij,3) = visc*txc(ij,3)
    ENDDO 
  
 

    ! CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,inb_scal_array), rad_param,&
    !     wrk1d(1:isize_wrk1d), wrk1d(isize_wrk1d+1:2*isize_wrk1d), wrk1d(2*isize_wrk1d+1:3*isize_wrk1d), wrk3d) ! Put radiation in wrk3d
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,4), wrk1d,wrk3d)
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

  CALL PARTICLE_SORT_HALO(x,z, grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
                           l_hq, l_tags, l_q)

  wrk1d(1)= scalex/imax_total ! wrk1d 1-3 intervalls
!  wrk1d(2)= scaley/jmax_total ! needed for interpolation
  wrk1d(2)= y(jmin_part+1)-y(jmin_part)
  wrk1d(3)= scalez/kmax_total

!#######################################################################
! RHS for particles in normal field
!#######################################################################

CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION(q,l_q,l_hq,y,wrk1d,txc,grid_start,grid_field_counter)
  !RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
    halo_start=grid_field_counter+1
    halo_end=grid_field_counter+halo_zone_x
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1(halo_field_1,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF
  !RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
    halo_start=grid_field_counter+halo_zone_x+1
    halo_end=grid_field_counter+halo_zone_x+halo_zone_z
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2(halo_field_2,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF
  !RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
    halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
    halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3(halo_field_3,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF

#else
  ! #####################################################################
  ! Put source terms into txc variables
  ! #####################################################################

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD) THEN
    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s


    ! CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,inb_scal_array), rad_param,&
    !     wrk1d(1:isize_wrk1d), wrk1d(isize_wrk1d+1:2*isize_wrk1d), wrk1d(2*isize_wrk1d+1:3*isize_wrk1d), wrk3d) ! Put radiation in wrk3d
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,1), wrk1d,wrk3d)
    ! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
    DO ij = 1,isize_field
       txc(ij,1) = dummy2*txc(ij,1)
    ENDDO 
    
    

    !LAPLACE Xi
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,1), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,1), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,1), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 
 
    !LAPLACE delta
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,2), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,2), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,2), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6)))
    ENDDO 

    !LAPLACE s
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
         dz, s(1,inb_scal_array), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
         dy, s(1,inb_scal_array), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
         dx, s(1,inb_scal_array), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)

    DO ij = 1,isize_field
       txc(ij,3) = visc*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 

    txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)

  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2) THEN !using second version of equation
    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s


    
    
    !LAPLACE s
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
         dz, s(1,inb_scal_array), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
         dy, s(1,inb_scal_array), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
         dx, s(1,inb_scal_array), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)

    DO ij = 1,isize_field
       txc(ij,3) = visc*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 


    txc(1:isize_field,4) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,4)

    CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, txc(1,4), txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d) ! square of chi gradient in txc(1,1)

    DO ij = 1,isize_field
       txc(ij,1) = visc*txc(ij,1)
    ENDDO 
  
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,2), wrk1d,wrk3d)
    ! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
    DO ij = 1,isize_field
       txc(ij,2) = dummy2*txc(ij,2)
    ENDDO 


  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !using second version of equation

    dummy2 = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(6)/body_param(5) ! delta_s
    dummy2 = 1/dummy2 !1/delta_s
    dummy = 1/body_param(3)  ! 1/chi_s

    !LAPLACE Xi
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,1), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,1), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,1), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
!       txc(ij,1) =txc(ij,1) + visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
       txc(ij,1) = visc*dummy*(txc(ij,4)+txc(ij,5)+txc(ij,6))
    ENDDO 
 
    !LAPLACE delta
    CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, s(1,2), txc(1,6), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,2), txc(1,5), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, s(1,2), txc(1,4), i0,i0, i0,i0, txc(1,3), wrk1d,wrk2d,wrk3d)
    
    DO ij = 1,isize_field
       txc(ij,1) =txc(ij,1) + (visc*dummy2*(txc(ij,4)+txc(ij,5)+txc(ij,6))) !first eq. without ds/dxi
    ENDDO 


   
    txc(1:isize_field,2) = C_1_R - dummy*s(1:isize_field,1) - dummy2*s(1:isize_field,2) !xi field in txc(1,2)



    CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, txc(1,2), txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d) ! square of chi gradient in txc(1,3)

    DO ij = 1,isize_field
       txc(ij,3) = visc*txc(ij,3)
    ENDDO 
  
 
    ! CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,inb_scal_array), rad_param,&
    !     wrk1d(1:isize_wrk1d), wrk1d(isize_wrk1d+1:2*isize_wrk1d), wrk1d(2*isize_wrk1d+1:3*isize_wrk1d), wrk3d) ! Put radiation in wrk3d
    CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,4), wrk1d,wrk3d)
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
  CALL PARTICLE_SORT_HALO(x,z, grid_field_counter, halo_zone_x, halo_zone_z,halo_zone_diagonal,&
                           l_hq, l_tags, l_q)
  
  wrk1d(1)= scalex/imax_total ! wrk1d 1-3 intervalls
!  wrk1d(2)= scaley/jmax_total ! needed for interpolation
  wrk1d(2)= y(jmin_part+1)-y(jmin_part)
  wrk1d(3)= scalez/kmax_total

!#######################################################################
! RHS for particles in normal field
!#######################################################################

  CALL HALO_PLANE_SHIFTING_SERIAL(q,txc,halo_field_1,halo_field_2,halo_field_3)

 
  grid_start=1
  CALL RHS_PARTICLE_GLOBAL_INTERPOLATION(q,l_q,l_hq,y,wrk1d,txc,grid_start,grid_field_counter)
  !RHS for particles in halo_zone_x
  IF (halo_zone_x .NE. 0) THEN
    halo_start=grid_field_counter+1
    halo_end=grid_field_counter+halo_zone_x
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_1(halo_field_1,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF
  !RHS for particles in halo_zone_z
  IF (halo_zone_z .NE. 0) THEN
    halo_start=grid_field_counter+halo_zone_x+1
    halo_end=grid_field_counter+halo_zone_x+halo_zone_z
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_2(halo_field_2,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF
  !RHS for particles in halo_zone_diagonal
  IF (halo_zone_diagonal .NE. 0) THEN
    halo_start=grid_field_counter+halo_zone_x+halo_zone_z+1
    halo_end=grid_field_counter+halo_zone_x+halo_zone_z+halo_zone_diagonal
    CALL RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3(halo_field_3,l_q,l_hq,y,wrk1d,halo_start,halo_end)
  END IF

DEALLOCATE(halo_field_1, STAT=alloc_err)
DEALLOCATE(halo_field_2, STAT=alloc_err)
DEALLOCATE(halo_field_3, STAT=alloc_err)

#endif



  IF (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN

#ifdef USE_MPI
    local_np = particle_vector(ims_pro+1)
#else
    local_np = particle_number
#endif

    DO ij=1,local_np
      l_hq(ij,2) = l_hq(ij,2) - lagrange_param(1)
    END DO
  ENDIF


        
  RETURN
END SUBROUTINE RHS_PARTICLE_GLOBAL

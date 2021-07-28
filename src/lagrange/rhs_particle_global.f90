#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
SUBROUTINE RHS_PARTICLE_GLOBAL(q,s, txc, l_q,l_hq,l_txc,l_comm, wrk1d,wrk2d,wrk3d)

  USE TLAB_TYPES,  ONLY : pointers_dt, pointers3d_dt
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, isize_particle
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : visc, radiation
  USE LAGRANGE_VARS, ONLY : l_g, ilagrange, lagrange_param
  USE THERMO_VARS, ONLY : thermo_param

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TREAL, DIMENSION(isize_field,*),    TARGET :: q, s, txc
  TREAL, DIMENSION(isize_particle,*), TARGET :: l_q, l_hq, l_txc
  TREAL, DIMENSION(*)                        :: l_comm
  TREAL, DIMENSION(*)                        :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TREAL dummy, dummy2
  TINTEGER bcs(2,2), nvar
  TINTEGER i
  TREAL delta_inv0, delta_inv2, delta_inv4
  
  TYPE(pointers3d_dt), DIMENSION(7) :: data
  TYPE(pointers_dt),   DIMENSION(7) :: data_out

! #####################################################################
  bcs = 0
  
! Setting pointers to velocity fields
  nvar = 0
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => q(:,1); data_out(nvar)%field => l_hq(:,1)
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => q(:,2); data_out(nvar)%field => l_hq(:,2)
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => q(:,3); data_out(nvar)%field => l_hq(:,3)
  
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
     nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,1); data_out(nvar)%field => l_txc(:,1)
     nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,2); data_out(nvar)%field => l_txc(:,2)
     nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,3); data_out(nvar)%field => l_txc(:,3)
     nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,4); data_out(nvar)%field => l_txc(:,4)
     l_txc(:,1:4) = C_0_R
     
  ENDIF

! -------------------------------------------------------------------
! Interpolating field data into particles
! The interpolated data is added to the existing data, which
!  consitutes already the evolution equation for particle position
! -------------------------------------------------------------------
  CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_g,l_q,l_comm, wrk3d)
  
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
     
     DO i = 1,l_g%np
        l_hq(i,4) = l_hq(i,4) - l_txc(i,1)/(C_1_R + EXP(l_txc(i,2)*delta_inv0))
        
        l_hq(i,5) = l_hq(i,5) - l_txc(i,4)/(C_1_R + EXP(l_txc(i,2)*delta_inv0)) &
                              - l_txc(i,3)*delta_inv4/(COSH(l_txc(i,2)*delta_inv2)**2) 
     ENDDO

  ELSE IF (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN

     l_hq(1:l_g%np,2) = l_hq(1:l_g%np,2) - lagrange_param(1)

  ENDIF

  DO i = 1,nvar
     NULLIFY(data(i)%field,data_out(i)%field)
  ENDDO
  
  RETURN
END SUBROUTINE RHS_PARTICLE_GLOBAL

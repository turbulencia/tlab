#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2011/11/01 - J.P. Mellado
!#              Created
!# 2013/01/20 - J.P. Mellado
!#              Introducing direct formulation of non-uniform grid
!#
!########################################################################
!# DESCRIPTION
!#
!# Apply the non-linear operator N(u) = visc* d^2/dx^2 s - u d/dx s
!# Derived from OPR_PARTIAL_XX to avoid one MPI transposition
!#
!########################################################################
!# ARGUMENTS 
!#
!# i1bc       In   grid structure: 0  periodic
!#                                 1  non-periodic
!# bcs2_imin  In   BC derivative at imin: 0  biased, non-zero
!#                                        1  forced to zero
!# bcs2_imax  In   BC derivative at imax: 0  biased, non-zero
!#                                        1  forced to zero
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# result     Out  Result N(u)
!# tmp       Out  Transpose velocity
!#
!########################################################################
SUBROUTINE OPR_BURGERS(iunif,imode_fdm,is,nxy,nz, i1bc, dx, &
  s,u,bcs1_imin,bcs1_imax,bcs2_imin,bcs2_imax,tmp,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : inb_grid_1, inb_grid_2, inb_grid_3
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER,                 INTENT(IN)    :: is
  TREAL, DIMENSION(nz,*),   INTENT(IN)    :: dx
  TREAL, DIMENSION(nxy*nz), INTENT(IN)    :: s,u
  TREAL, DIMENSION(nxy*nz), INTENT(OUT)   :: tmp
  TREAL, DIMENSION(nxy*nz), INTENT(INOUT) :: wrk3d
  TREAL, DIMENSION(*),      INTENT(INOUT) :: wrk2d

  TINTEGER iunif, imode_fdm, nxy,nz 
  TINTEGER i1bc, bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax,ip,ij

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( i1bc .EQ. 0 ) THEN 
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; 
        CALL FDM_C1N4P_RHS(nz,nxy, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
               imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; 
        CALL FDM_C1N6P_RHS(nz,nxy, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; 
        CALL FDM_C1N8P_RHS(nz,nxy, s, wrk3d)
     ENDIF
     ip = inb_grid_1 - 1
     CALL TRIDPSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),&
          dx(1,ip+4),dx(1,ip+5), wrk3d,wrk2d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; 
        CALL FDM_C2N4P_RHS(nz,nxy, s, tmp)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. & 
               imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; 
        CALL FDM_C2N6P_RHS(nz,nxy, s, tmp)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; 
        CALL FDM_C2N6P_RHS(nz,nxy, s, tmp) !8th not yet developed
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),&
          dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), tmp,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; 
        CALL FDM_C1N4_RHS(nz,nxy, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; 
        CALL FDM_C1N6_RHS(nz,nxy, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; 
        CALL FDM_C1N8_RHS(nz,nxy, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; 
        CALL FDM_C1N6_RHS(nz,nxy, bcs1_imin,bcs1_imax, s, wrk3d) ! not yet implemented
     ENDIF
     ip = inb_grid_1 + (bcs1_imin + bcs1_imax*2)*3 - 1
     CALL TRIDSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), wrk3d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; 
        CALL FDM_C2N4_RHS(iunif, nz,nxy, bcs2_imin,bcs2_imax, &
             dx, s, wrk3d, tmp)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; 
        CALL FDM_C2N6_RHS(iunif, nz,nxy, bcs2_imin,bcs2_imax, &
             dx, s, wrk3d, tmp)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; 
        CALL FDM_C2N6_RHS(iunif, nz,nxy, bcs2_imin,bcs2_imax, &
             dx, s, wrk3d, tmp)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; 
        CALL FDM_C2N6N_RHS(nz,nxy, dx(1,inb_grid_2+3), &
             s, tmp)
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), tmp)

  ENDIF

! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nxy*nz
     tmp(ij) = tmp(ij) - u(ij)*wrk3d(ij) ! diffusivity included in 2.-order derivative
  ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE OPR_BURGERS

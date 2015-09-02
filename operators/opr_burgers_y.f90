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
!# Apply the non-linear operator N(u) = visc* d^2/dy^2 s - u d/dy s
!# Derived from OPR_PARTIAL_YY to avoid one MPI transposition
!#
!########################################################################
!# ARGUMENTS 
!#
!# j1bc       In  grid structure: 0  periodic
!#                                1  non-periodic
!# bcs2_jmin  In  BC derivative at imin: 0  biased, non-zero
!#                                       1  forced to zero
!# bcs2_jmax  In  BC derivative at imax: 0  biased, non-zero
!#                                       1  forced to zero
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# result     Out  Result N(u)
!# tmp1       Out  Transpose velocity
!#
!########################################################################
SUBROUTINE OPR_BURGERS_Y(ivel, is, iunif,imode_fdm, nx,ny,nz, j1bc, dy, s,u1,u2, result, &
     bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax, tmp1, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : jmax_total, inb_grid_1, inb_grid_2, inb_grid_3
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER ivel, is, iunif
  TINTEGER imode_fdm, nx,ny,nz, j1bc, bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax
  TREAL, DIMENSION(jmax_total,*), INTENT(IN)    :: dy
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER ip, nxy, nxz, ij
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst1, p_dst2, p_vel

! ###################################################################
  IF ( bcs2_jmin + bcs2_jmax .GT. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_BURGERS_Y. Only developed for biased BCs.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  IF ( jmax_total .EQ. 1 ) THEN ! Set to zero in 2D case
  result = C_0_R
     
  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Make y direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org  => s
     p_dst1 => tmp1
     p_dst2 => result
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (s, nxy, nxy, nz, tmp1, nz)
#else
     CALL DNS_TRANSPOSE(s, nxy, nz, nxy, tmp1, nz)
#endif
     p_org  => tmp1
     p_dst1 => result
     p_dst2 => wrk3d
  ENDIF

! pointer to velocity
  IF ( ivel .EQ. 0 ) THEN; p_vel => p_org
  ELSE;
     IF ( nz .EQ. 1 ) THEN; p_vel => u1         ! I do not need the transposed
     ELSE;                  p_vel => u2; ENDIF  ! I do     need the transposed
  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( j1bc  .EQ. 0 ) THEN
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(jmax_total,nxz, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(jmax_total,nxz, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(jmax_total,nxz, p_org, p_dst1)
     ENDIF
     ip = inb_grid_1 - 1
     CALL TRIDPSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3),dy(1,ip+4),dy(1,ip+5), p_dst1,wrk2d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C2N4P_RHS(jmax_total,nxz, p_org, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C2N6P_RHS(jmax_total,nxz, p_org, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C2N6P_RHS(jmax_total,nxz, p_org, p_dst2) ! 8th not yet developed
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3),dy(1,ip+4),dy(1,ip+5), p_dst2,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)  ! not yet implemented
     ENDIF
     ip = inb_grid_1 + (bcs1_jmin + bcs1_jmax*2)*3 - 1
     CALL TRIDSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3), p_dst1)
     
! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C2N4_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C2N6N_RHS(jmax_total,nxz, dy(1,inb_grid_2+3), p_org, p_dst2)
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3), p_dst2)

  ENDIF
     
! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nx*ny*nz
     p_dst2(ij) = p_dst2(ij) - p_vel(ij)*p_dst1(ij) ! diffusivity included in 2.-order derivative
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst2, nz, nz, nxy, result, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst2, nz, nxy, nz, result, nxy)
#endif
  ENDIF

  NULLIFY(p_org,p_dst1,p_dst2, p_vel)

  ENDIF

  RETURN
END SUBROUTINE OPR_BURGERS_Y

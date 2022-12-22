#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2010/12/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!########################################################################
! Sinusoidal forcing
!########################################################################
SUBROUTINE FI_FORCING_0(imax,jmax,kmax, time,visc, u,v, h_u,h_v)
    use OPR_PARTIAL
  IMPLICIT NONE

#include "integers.h"
  
  TINTEGER imax,jmax,kmax
  TREAL time,visc
  TREAL, DIMENSION(imax,jmax,kmax) :: u,v, h_u,h_v

! -----------------------------------------------------------------------
  TREAL omega, sigma, amplitude
!  TINTEGER i,j

!########################################################################
  omega = C_2_R*ACOS(-C_1_R)
  sigma = C_2_R*omega*omega*visc

!  amplitude =-( C_1_R + (sigma/omega)**2 )*omega
!  amplitude = amplitude*sin(omega*time)
!  DO j = 1,jmax; DO i = 1,imax
!     u(i,j,1) = SIN(x(i)*omega)*COS(g(2)%nodes(j)*omega)
!     v(i,j,1) =-COS(x(i)*omega)*SIN(g(2)%nodes(j)*omega)
!  ENDDO; ENDDO

  amplitude = sin(omega*time)

  h_u = h_u + amplitude*u
  h_v = h_v + amplitude*v
  
 RETURN
END SUBROUTINE FI_FORCING_0

!########################################################################
! Velocity field with no-slip
!########################################################################
SUBROUTINE FI_FORCING_1(imax,jmax,kmax, time,visc, h1,h2, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imax,jmax,kmax
  TREAL time,visc

  TREAL, DIMENSION(imax*jmax*kmax) :: h1,h2
  TREAL, DIMENSION(imax*jmax*kmax) :: tmp1,tmp2,tmp3,tmp4
  TREAL                            :: wrk2d(*), wrk3d(*)

! -----------------------------------------------------------------------
  TINTEGER ij, i, j, bcs(2,2)
  TREAL pi_loc

  bcs = 0

  pi_loc     = ACOS(-C_1_R)

! #######################################################################
  DO j = 1,jmax; DO i = 1,imax; ij = imax*(j-1) + i
!     tmp1(ij) = sin(C_2_R*pi_loc*g(1)%nodes(i))*       sin(C_4_R*pi_loc*g(2)%nodes(j))
!     tmp2(ij) =-cos(C_2_R*pi_loc*g(1)%nodes(i))*(C_1_R-cos(C_4_R*pi_loc*g(2)%nodes(j)))*C_05_R
     tmp1(ij) = sin(pi_loc*g(1)%nodes(i))*sin(pi_loc*g(1)%nodes(i))*sin(C_2_R*pi_loc*g(2)%nodes(j))
     tmp2(ij) =-sin(C_2_R*pi_loc*g(1)%nodes(i))*sin(pi_loc*g(2)%nodes(j))*sin(pi_loc*g(2)%nodes(j))
  ENDDO; ENDDO

! Time terms
  DO j = 1,jmax; DO i = 1,imax; ij = imax*(j-1) + i
     h1(ij) = h1(ij) - tmp1(ij)*C_2_R*pi_loc*sin(C_2_R*pi_loc*time)
     h2(ij) = h2(ij) - tmp2(ij)*C_2_R*pi_loc*sin(C_2_R*pi_loc*time)
  ENDDO; ENDDO

! velocities
  DO j = 1,jmax; DO i = 1,imax; ij = imax*(j-1) + i
     tmp1(ij) = tmp1(ij)*cos(C_2_R*pi_loc*time)
     tmp2(ij) = tmp2(ij)*cos(C_2_R*pi_loc*time)
  ENDDO; ENDDO

! Diffusion and convection terms in Ox momentum eqn
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), tmp1, tmp4,  tmp3, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) - visc*( tmp4(ij) ) + ( tmp3(ij)*tmp2(ij) )
  ENDDO

  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4,  tmp3, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) - visc*( tmp4(ij) ) + ( tmp3(ij)*tmp1(ij) )
  ENDDO

! Diffusion and convection terms in Oy momentum eqn
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp4,  tmp3, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h2(ij) = h2(ij) - visc*( tmp4(ij) ) + ( tmp3(ij)*tmp2(ij) )
  ENDDO

  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), tmp2, tmp4,  tmp3, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h2(ij) = h2(ij) - visc*( tmp4(ij) ) + ( tmp3(ij)*tmp1(ij) )
  ENDDO

! #######################################################################
  DO j = 1,jmax; DO i = 1,imax; ij = imax*(j-1) + i
!     tmp1(ij) = cos(C_4_R*pi_loc*g(1)%nodes(i))*(C_2_R-cos(C_4_R*pi_loc*g(2)%nodes(j)))/C_8_R &
!          - C_05_R*(sin(C_2_R*pi_loc*g(2)%nodes(j)))**4
     tmp1(ij) = sin(C_2_R*pi_loc*g(1)%nodes(i))*sin(C_2_R*pi_loc*g(2)%nodes(j))
     tmp1(ij) = tmp1(ij)*(cos(C_2_R*pi_loc*time))**2
  ENDDO; ENDDO

! Pressure gradient
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp1, tmp3, wrk3d, wrk2d,wrk3d)

  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) + tmp2(ij)
     h2(ij) = h2(ij) + tmp3(ij)
  ENDDO

  RETURN
END SUBROUTINE FI_FORCING_1


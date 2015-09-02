!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

SUBROUTINE BOUNDARY_BUFFER_FILTER(x, rho,u,v,w,e,z1, txc1,txc2,txc3,txc4,txc5, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax, kmax_total, i1bc, j1bc, k1bc
  USE DNS_GLOBAL, ONLY : icalc_scal, inb_scal
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_LOCAL,  ONLY : buff_nps_jmin, buff_nps_jmax, buff_imax, buff_nps_imax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)                       :: x
  TREAL, DIMENSION(imax,jmax,kmax  )        :: rho, u, v, w, e
  TREAL, DIMENSION(imax,jmax,kmax,*)        :: z1
  TREAL, DIMENSION(buff_nps_imax,jmax,kmax) :: txc1, txc2, txc3, txc4, txc5 
  TREAL, DIMENSION(*)                       :: wrk1d,wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER id, iloc, i, j, k, is, ibc_x(4), ibc_y(4), ibc_z(4)
  TREAL eta, delta, amp, ampr, rho_ratio

! ###################################################################
! BCs for the filters (see routine FILTER)
  ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
  ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

! ###################################################################
! Bottom boundary
! ###################################################################
  IF ( buff_nps_jmin .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter at the bottom not yet implemented.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! ###################################################################
! Top boundary
! ###################################################################
  IF ( buff_nps_jmax .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter at the top not yet implemented.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! ###################################################################
! Outflow boundary
! ###################################################################
  IF ( buff_nps_imax .GT. 1 ) THEN
     id = DNS_MPI_K_OUTBCS

! -------------------------------------------------------------------
! Flow
! -------------------------------------------------------------------
     DO k = 1,kmax
        DO j = 1,jmax
           DO i = buff_imax,imax
              iloc = i-buff_imax+1
              txc1(iloc,j,k) = rho(i,j,k)
              txc2(iloc,j,k) = rho(i,j,k)*u(i,j,k)
              txc3(iloc,j,k) = rho(i,j,k)*v(i,j,k)
              txc4(iloc,j,k) = rho(i,j,k)*w(i,j,k)
              txc5(iloc,j,k) = rho(i,j,k)*e(i,j,k)
!              txc2(iloc,j,k) = u(i,j,k)
!              txc3(iloc,j,k) = v(i,j,k)
!              txc4(iloc,j,k) = w(i,j,k)
!              txc5(iloc,j,k) = e(i,j,k)
           ENDDO
        ENDDO
     ENDDO

     CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
          txc1, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
          txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
          txc3, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
          txc4, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
          txc5, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
     delta = (x(imax)-x(buff_imax))/C_16_R
     DO i = buff_imax,imax
        iloc = i-buff_imax+1

        eta = x(i) - C_05_R*( x(imax)+x(buff_imax) )
        amp = C_05_R*(C_1_R+TANH(C_05_R*eta/delta))
        ampr= C_1_R-amp

        DO k = 1,kmax
           DO j = 1,jmax
              rho_ratio = rho(i,j,k)
              rho(i,j,k) = ampr*rho(i,j,k) + amp*txc1(iloc,j,k)
              rho_ratio = rho_ratio/rho(i,j,k)

              u(i,j,k) = ampr*rho_ratio*u(i,j,k) + amp*txc2(iloc,j,k)/rho(i,j,k)
              v(i,j,k) = ampr*rho_ratio*v(i,j,k) + amp*txc3(iloc,j,k)/rho(i,j,k)
              w(i,j,k) = ampr*rho_ratio*w(i,j,k) + amp*txc4(iloc,j,k)/rho(i,j,k)
              e(i,j,k) = ampr*rho_ratio*e(i,j,k) + amp*txc5(iloc,j,k)/rho(i,j,k)

! store rho_ratio to be used in scalar section
              txc1(iloc,j,k) = rho_ratio

!              rho(i,j,k) = ampr*rho(i,j,k) + amp*txc1(iloc,j,k)
!              u(i,j,k)   = ampr*u(i,j,k)   + amp*txc2(iloc,j,k)
!              v(i,j,k)   = ampr*v(i,j,k)   + amp*txc3(iloc,j,k)
!              w(i,j,k)   = ampr*w(i,j,k)   + amp*txc4(iloc,j,k)
!              e(i,j,k)   = ampr*e(i,j,k)   + amp*txc5(iloc,j,k)

           ENDDO
        ENDDO
     ENDDO

! -------------------------------------------------------------------
! Scalar
! -------------------------------------------------------------------
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal

           DO k = 1,kmax
              DO j = 1,jmax
                 DO i = buff_imax,imax
                    iloc = i-buff_imax+1
                    txc2(iloc,j,k) = rho(i,j,k)*z1(i,j,k,is)
!                    txc2(iloc,j,k) = z1(i,j,k,is)
                 ENDDO
              ENDDO
           ENDDO

           CALL OPR_FILTER(i2, buff_nps_imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
                txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
           delta = (x(imax)-x(buff_imax))/C_16_R
           DO i = buff_imax,imax
              iloc = i-buff_imax+1

              eta = x(i) - C_05_R*( x(imax)+x(buff_imax) )
              amp = C_05_R*(C_1_R+TANH(C_05_R*eta/delta))
              ampr= C_1_R-amp

              DO k = 1,kmax
                 DO j = 1,jmax
                    z1(i,j,k,is) = ampr*txc1(iloc,j,k)*z1(i,j,k,is) + &
                         amp*txc2(iloc,j,k)/rho(i,j,k)
!                    z1(i,j,k,is) = ampr*z1(i,j,k,is) + amp*txc2(iloc,j,k)
                 ENDDO
              ENDDO

           ENDDO

        ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BUFFER_FILTER

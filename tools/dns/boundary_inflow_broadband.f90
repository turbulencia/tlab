#include "types.h"
  
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
SUBROUTINE BOUNDARY_INFLOW_BROADBAND&
     (etime, inf_rhs, x_inf, y_inf, z_inf, q_inf,z1_inf, txc, wrk1d, wrk2d, wrk3d)
  
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : ifrc_mode, ifrc_ifield, frc_adapt
  USE DNS_LOCAL, ONLY : imax_inf,jmax_inf,kmax_inf
  USE DNS_LOCAL, ONLY : scalex_inf

  IMPLICIT NONE
  
  TREAL etime
  TREAL inf_rhs(jmax,kmax,*) 
  TREAL x_inf(imax_inf), y_inf(jmax_inf), z_inf(kmax_total)
  TREAL q_inf(imax_inf,jmax_inf,kmax_inf,*)
  TREAL z1_inf(imax_inf,jmax_inf,kmax_inf,*)

  TREAL txc(*), wrk1d(*), wrk2d(*), wrk3d(*)

  TARGET :: q_inf

! -------------------------------------------------------------------
  TREAL xaux, dx_loc, vmult
  TINTEGER joffset, jglobal, ileft, iright, j, k, n, is
  TREAL BSPLINES3P, BSPLINES3

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:),   POINTER :: u_inf, v_inf, w_inf, p_inf, rho_inf

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_BROADBAND' )
#endif

! Define pointers
  u_inf   => q_inf(:,:,:,1)
  v_inf   => q_inf(:,:,:,2)
  w_inf   => q_inf(:,:,:,3)
  p_inf   => q_inf(:,:,:,4)
  rho_inf => q_inf(:,:,:,5)

! Transient factor
  IF ( frc_adapt .GT. C_0_R .AND. etime .LE. frc_adapt ) THEN
     vmult = etime / frc_adapt
  ELSE
     vmult = C_1_R
  ENDIF

! check if we need to read again inflow data
  IF ( ifrc_mode .EQ. 3 .AND. INT(mean_u*etime/scalex_inf)+1 .NE. ifrc_ifield ) THEN
     CALL BOUNDARY_INFLOW_INIT&
          (etime, x_inf,y_inf,z_inf, q_inf,z1_inf, txc, wrk1d,wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Getting the position
! ###################################################################
  joffset = (jmax-jmax_inf)/2

  xaux = mean_u*etime
! Remove integral length scales of box
  xaux = xaux - INT(xaux/scalex_inf)*scalex_inf
! Set distance from box initial length
  xaux = scalex_inf-xaux

  dx_loc = x_inf(2) - x_inf(1)
! Get left index
  ileft = INT(xaux/dx_loc) +1
! Check bounds
  IF ( ileft .GT. imax_inf ) THEN
     ileft = 1
  ENDIF
! Set right index
  iright = ileft + 1
! Check bounds
  IF ( iright .GT. imax_inf ) THEN
     iright = 1
  ENDIF
! Get relative distance from left point
  xaux = (xaux-(x_inf(ileft)-x_inf(1)))/dx_loc

! ###################################################################
! Sampling the information
! ###################################################################
! -------------------------------------------------------------------
! Periodic
! -------------------------------------------------------------------
  IF ( ifrc_mode .EQ. 2 ) THEN
     DO k = 1,kmax_inf
        DO j = 1,jmax_inf
           jglobal = joffset + j
           inf_rhs(jglobal,k,1) = inf_rhs(jglobal,k,1) + vmult*&
                BSPLINES3P(rho_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,2) = inf_rhs(jglobal,k,2) + vmult*&
                BSPLINES3P(u_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,3) = inf_rhs(jglobal,k,3) + vmult*&
                BSPLINES3P(v_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,4) = inf_rhs(jglobal,k,4) + vmult*&
                BSPLINES3P(w_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,5) = inf_rhs(jglobal,k,5) + vmult*&
                BSPLINES3P(p_inf(1,j,k), imax_inf, ileft, xaux)
           
           IF ( icalc_scal .EQ. 1 ) THEN
              DO is = 1,inb_scal_array
                 inf_rhs(jglobal,k,5+is) = inf_rhs(jglobal,k,5+is) + vmult*&
                   BSPLINES3P(z1_inf(1,j,k,is), imax_inf, ileft, xaux)
              ENDDO
           ENDIF
           
        ENDDO
     ENDDO
    
! -------------------------------------------------------------------
! Sequential
! -------------------------------------------------------------------
  ELSE
     DO k = 1,kmax_inf
        DO j = 1,jmax_inf
           jglobal = joffset + j
           inf_rhs(jglobal,k,1) = inf_rhs(jglobal,k,1) + vmult*&
                BSPLINES3(rho_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,2) = inf_rhs(jglobal,k,2) + vmult*&
                BSPLINES3(u_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,3) = inf_rhs(jglobal,k,3) + vmult*&
                BSPLINES3(v_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,4) = inf_rhs(jglobal,k,4) + vmult*&
                BSPLINES3(w_inf(1,j,k), imax_inf, ileft, xaux)
           inf_rhs(jglobal,k,5) = inf_rhs(jglobal,k,5) + vmult*&
                BSPLINES3(p_inf(1,j,k), imax_inf, ileft, xaux)
           
           IF ( icalc_scal .EQ. 1 ) THEN
              DO is = 1,inb_scal_array
                 inf_rhs(jglobal,k,5+is) = inf_rhs(jglobal,k,5+is) + vmult*&
                   BSPLINES3(z1_inf(1,j,k,is), imax_inf, ileft, xaux)
              ENDDO
           ENDIF

        ENDDO
     ENDDO
           
  ENDIF

! ###################################################################
! Filling the rest
! ###################################################################
  DO n = 1,5+inb_scal_array
     DO k = 1,kmax_inf
        DO j = 1,joffset
           inf_rhs(j,k,n) = C_0_R
        ENDDO
        DO j = jmax-joffset+1,jmax
           inf_rhs(j,k,n) = C_0_R
        ENDDO
     ENDDO
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_BROADBAND' )
#endif
  RETURN
END SUBROUTINE BOUNDARY_INFLOW_BROADBAND

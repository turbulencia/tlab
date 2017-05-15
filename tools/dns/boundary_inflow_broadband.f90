#include "types.h"
  
!########################################################################
!# DESCRIPTION
!#
!# Calculating RHS forcings at the inflow plane in spatially evolving cases
!#
!########################################################################
SUBROUTINE BOUNDARY_INFLOW_BROADBAND(etime, inf_rhs, q_inf,s_inf, txc, wrk2d,wrk3d)
  
  USE DNS_GLOBAL, ONLY : jmax,kmax, inb_flow, inb_scal
  USE DNS_GLOBAL, ONLY : icalc_scal
  USE DNS_GLOBAL, ONLY : qbg
  USE DNS_LOCAL, ONLY  : ifrc_mode, ifrc_ifield, frc_adapt
  USE DNS_LOCAL, ONLY  : g_inf

  IMPLICIT NONE
  
  TREAL etime
  TREAL, DIMENSION(jmax,kmax,inb_flow+inb_scal), INTENT(OUT)   :: inf_rhs
  TREAL, DIMENSION(g_inf(1)%size,&
                   g_inf(2)%size,&
                   g_inf(3)%size,inb_flow),      INTENT(IN)    :: q_inf
  TREAL, DIMENSION(g_inf(1)%size,&
                   g_inf(2)%size,&
                   g_inf(3)%size,inb_scal),      INTENT(IN)    :: s_inf
  TREAL, DIMENSION(g_inf(1)%size*&
                   g_inf(2)%size*&
                   g_inf(3)%size),               INTENT(INOUT) :: txc, wrk3d
  TREAL, DIMENSION(*),                           INTENT(INOUT) :: wrk2d

  TARGET :: q_inf

! -------------------------------------------------------------------
  TREAL xaux, dx_loc, vmult
  TINTEGER joffset, jglobal, ileft, iright, j, k, is, ip
  TREAL BSPLINES3P, BSPLINES3

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_BROADBAND' )
#endif

! Transient factor
  IF ( frc_adapt .GT. C_0_R .AND. etime .LE. frc_adapt ) THEN
     vmult = etime / frc_adapt
  ELSE
     vmult = C_1_R
  ENDIF

! check if we need to read again inflow data
  IF ( ifrc_mode .EQ. 3 .AND. INT(qbg(1)%mean*etime/g_inf(1)%scale)+1 .NE. ifrc_ifield ) THEN
     CALL BOUNDARY_INFLOW_INIT(etime, q_inf,s_inf, txc, wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Getting the position
! ###################################################################
  joffset = (jmax-g_inf(2)%size)/2

  xaux = qbg(1)%mean*etime
! Remove integral length scales of box
  xaux = xaux - INT(xaux/g_inf(1)%scale)*g_inf(1)%scale
! Set distance from box initial length
  xaux = g_inf(1)%scale-xaux

  dx_loc = g_inf(1)%nodes(2) - g_inf(1)%nodes(1)
! Get left index
  ileft = INT(xaux/dx_loc) +1
! Check bounds
  IF ( ileft .GT. g_inf(1)%size ) THEN
     ileft = 1
  ENDIF
! Set right index
  iright = ileft + 1
! Check bounds
  IF ( iright .GT. g_inf(1)%size ) THEN
     iright = 1
  ENDIF
! Get relative distance from left point
  xaux = (xaux-(g_inf(1)%nodes(ileft)-g_inf(1)%nodes(1)))/dx_loc

! ###################################################################
! Sampling the information
! ###################################################################
! -------------------------------------------------------------------
! Periodic
! -------------------------------------------------------------------
  IF ( ifrc_mode .EQ. 2 ) THEN
     DO k = 1,kmax
        DO j = 1,g_inf(2)%size
           jglobal = joffset + j
           DO is = 1,inb_scal
              inf_rhs(jglobal,k,is) = inf_rhs(jglobal,k,is) + vmult *BSPLINES3P(q_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
           ENDDO

           IF ( icalc_scal .EQ. 1 ) THEN
              DO is = 1,inb_scal
                 ip = inb_flow +is
                 inf_rhs(jglobal,k,ip) = inf_rhs(jglobal,k,ip) + vmult *BSPLINES3P(s_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
              ENDDO
           ENDIF
           
        ENDDO
     ENDDO
    
! -------------------------------------------------------------------
! Sequential
! -------------------------------------------------------------------
  ELSE
     DO k = 1,kmax
        DO j = 1,g_inf(2)%size
           jglobal = joffset + j
           DO is = 1,inb_flow
              inf_rhs(jglobal,k,is) = inf_rhs(jglobal,k,is) + vmult *BSPLINES3(q_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
           ENDDO
           
           IF ( icalc_scal .EQ. 1 ) THEN
              DO is = 1,inb_scal
                 ip = inb_flow +is
                 inf_rhs(jglobal,k,ip) = inf_rhs(jglobal,k,ip) + vmult *BSPLINES3(s_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
              ENDDO
           ENDIF

        ENDDO
     ENDDO
           
  ENDIF

! ###################################################################
! Filling the rest
! ###################################################################
  DO j = 1,joffset
     inf_rhs(j,:,:) = inf_rhs(j,:,:) + C_0_R
  ENDDO
  DO j = jmax-joffset+1,jmax
     inf_rhs(j,:,:) = inf_rhs(j,:,:) + C_0_R
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_BROADBAND' )
#endif
  RETURN
END SUBROUTINE BOUNDARY_INFLOW_BROADBAND

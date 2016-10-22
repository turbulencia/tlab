#include "types.h"
#include "dns_const.h"

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
SUBROUTINE BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs)
  
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : frc_length, frc_adapt, frc_delta, ifrcdsc_mode
  USE DNS_LOCAL, ONLY : nx2d, nx3d, nz3d, A2D, A3D, Phix2d, Phix3d, Phiz3d

#ifdef USE_MPI
  USE DNS_MPI
#endif
 
  IMPLICIT NONE

  TREAL etime
  TREAL inf_rhs(jmax,kmax,*)

! -------------------------------------------------------------------
  TINTEGER j, k, jsim, idsp
  TREAL ycenter, fy, fyp, wx, wz, wxloc, wzloc, xaux
  TREAL u2d, v2d, u3d, v3d, w3d, vmult
  TINTEGER inx2d, inx3d, inz3d

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_DISCRETE' )
#endif

#ifdef USE_MPI
  idsp = ims_offset_k 
#else 
  idsp = 0
#endif

  wx = C_2_R * C_PI_R / frc_length
  wz = C_2_R * C_PI_R / scalez
  xaux =-mean_u*etime

! Transient factor
  IF ( frc_adapt .GT. C_0_R .AND. etime .LE. frc_adapt ) THEN
     vmult = etime / frc_adapt
  ELSE
     vmult = C_1_R
  ENDIF

! ###################################################################
! Forcing for shear 
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

     DO j = 1, jmax

        ycenter = g(2)%nodes(j) - g(2)%scale *ycoor_u - g(2)%nodes(1)
        fy  = EXP(-(ycenter/(C_2_R*frc_delta))**2)*frc_delta
        fyp =-ycenter*fy/(C_2_R * frc_delta**2)

        ! 2D perturbation

        DO inx2d = 1,nx2d

           wxloc = M_REAL(inx2d)*wx

           DO k = 1,kmax
              u2d = A2d(inx2d) * wxloc*        COS(wxloc*xaux+Phix2d(inx2d)) * fyp 
              v2d = A2d(inx2d) * wxloc*wxloc * SIN(wxloc*xaux+Phix2d(inx2d)) * fy
              inf_rhs(j,k,2) = inf_rhs(j,k,2) - vmult*mean_u*u2d
              inf_rhs(j,k,3) = inf_rhs(j,k,3) - vmult*mean_u*v2d
           ENDDO

        ENDDO

        ! 3D perturbation

        IF (kmax .GT. 1) THEN

           DO inx3d = 1, nx3d
              DO inz3d = 1, nz3d

                 wxloc = M_REAL(inx3d)*wx
                 wzloc = M_REAL(inz3d)*wz

                 DO k=1, kmax
                    u3d = A3d(inx3d)*wxloc*SIN(wxloc*xaux+Phix3d(inx3d)) * &
                         SIN(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fyp
                    v3d = A3d(inx3d)*wxloc*COS(wxloc*xaux+Phix3d(inx3d)) * &
                         SIN(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fy * &
                         (wxloc+wzloc)
                    w3d =-A3d(inx3d)*wxloc*SIN(wxloc*xaux+Phix3d(inx3d)) * &
                         COS(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fyp
                    inf_rhs(j,k,2) = inf_rhs(j,k,2) - vmult*mean_u*u3d
                    inf_rhs(j,k,3) = inf_rhs(j,k,3) - vmult*mean_u*v3d
                    inf_rhs(j,k,4) = inf_rhs(j,k,4) - vmult*mean_u*w3d
                 ENDDO

              ENDDO
           ENDDO

        ENDIF

     ENDDO

! ###################################################################
! Forcing for jet
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN

     DO j = 1,jmax/2

        jsim = jmax - j + 1

        ycenter = g(2)%nodes(j) - g(2)%scale *ycoor_u + diam_u/C_2_R - g(2)%nodes(1)
        fy  = EXP(-(ycenter/(C_2_R*frc_delta))**2)*frc_delta
        fyp =-ycenter*fy/(C_2_R * frc_delta**2)

        ! 2D perturbation

        DO inx2d = 1,nx2d

           wxloc = M_REAL(inx2d)*wx

           DO k = 1,kmax
              u2d = A2d(inx2d) * wxloc *       COS(wxloc*xaux+Phix2d(inx2d)) *fyp 
              v2d = A2d(inx2d) * wxloc*wxloc * SIN(wxloc*xaux+Phix2d(inx2d)) *fy
              inf_rhs(j,k,2) = inf_rhs(j,k,2) - vmult*mean_u*u2d
              inf_rhs(j,k,3) = inf_rhs(j,k,3) - vmult*mean_u*v2d
              !          varicose
              IF (ifrcdsc_mode .EQ. 1) THEN
                 inf_rhs(jsim,k,2) = inf_rhs(jsim,k,2) - vmult*mean_u*u2d
                 inf_rhs(jsim,k,3) = inf_rhs(jsim,k,3) + vmult*mean_u*v2d
                 !          sinuous
              ELSE
                 inf_rhs(jsim,k,2) = inf_rhs(jsim,k,2) + vmult*mean_u*u2d
                 inf_rhs(jsim,k,3) = inf_rhs(jsim,k,3) - vmult*mean_u*v2d
              ENDIF

           ENDDO

        ENDDO

        ! 3D perturbation

        IF (kmax .GT. 1) THEN

           DO inx3d = 1, nx3d
              DO inz3d = 1, nz3d

                 wxloc = M_REAL(inx3d)*wx
                 wzloc = M_REAL(inz3d)*wz

                 DO k=1, kmax
                    u3d = A3d(inx3d)*wxloc*SIN(wxloc*xaux+Phix3d(inx3d)) * &
                         SIN(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fyp
                    v3d = A3d(inx3d)*wxloc*COS(wxloc*xaux+Phix3d(inx3d)) * &
                         SIN(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fy * &
                         (wxloc+wzloc)
                    w3d =-A3d(inx3d)*wxloc*SIN(wxloc*xaux+Phix3d(inx3d)) * &
                         COS(wzloc*g(3)%nodes(k)+Phiz3d(inz3d)) * fyp
                    inf_rhs(j,k,2) = inf_rhs(j,k,2) - vmult*mean_u*u3d
                    inf_rhs(j,k,3) = inf_rhs(j,k,3) - vmult*mean_u*v3d
                    inf_rhs(j,k,4) = inf_rhs(j,k,4) - vmult*mean_u*w3d
                    !             varicose
                    IF (ifrcdsc_mode .EQ. 1) THEN
                       inf_rhs(jsim,k,2) = inf_rhs(jsim,k,2) - vmult*mean_u*u3d
                       inf_rhs(jsim,k,3) = inf_rhs(jsim,k,3) + vmult*mean_u*v3d
                       inf_rhs(jsim,k,4) = inf_rhs(jsim,k,4) - vmult*mean_u*w3d
                       !             sinuous
                    ELSE
                       inf_rhs(jsim,k,2) = inf_rhs(jsim,k,2) + vmult*mean_u*u3d
                       inf_rhs(jsim,k,3) = inf_rhs(jsim,k,3) - vmult*mean_u*v3d
                       inf_rhs(jsim,k,4) = inf_rhs(jsim,k,4) + vmult*mean_u*w3d
                    ENDIF

                 ENDDO

              ENDDO
           ENDDO

        ENDIF

     ENDDO

  ELSE 

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_DISCRETE' )
#endif

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_DISCRETE

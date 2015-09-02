!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/07/11 - J.P. Mellado
!#              Cleaned. 
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE PDF2V3D(ilim, imax, jmax, kmax, umin_ext, umax_ext, vmin_ext, vmax_ext, &
     u, v, nxp, nyp, x, y, pdf, wrk)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ilim, imax, jmax, kmax
  TINTEGER nxp, nyp
  TREAL u(*)
  TREAL v(*)
  TREAL x(nxp,nyp)
  TREAL y(nxp,nyp)
  TREAL pdf(nxp,nyp)
  TREAL wrk(*)
  TREAL umin_ext, umax_ext, vmin_ext, vmax_ext

! -------------------------------------------------------------------
  TREAL umin, umax, vmin, vmax
  TREAL xstep, ystep, pnorm
  TINTEGER ip, jp, np, i
#ifdef USE_MPI
  INTEGER ims_err, itmp
  TINTEGER np_total
#endif

! ###################################################################
! Calculate Minimum and Maximum

  np = 0

  IF ( ilim .EQ. 0 ) THEN
     CALL MINMAX(imax,jmax,kmax, u, umin,umax)
     CALL MINMAX(imax,jmax,kmax, v, vmin,vmax)

  ELSE
     umin = umin_ext
     umax = umax_ext
     vmin = vmin_ext
     vmax = vmax_ext
  ENDIF

! Calculate Step in Histogram

  DO jp=1, nyp
     DO ip=1, nxp
        pdf(ip,jp) = C_0_R
     ENDDO
  ENDDO

  xstep = (umax-umin)/M_REAL(nxp)
  ystep = (vmax-vmin)/M_REAL(nyp)

  IF ( ABS(xstep) .LT. C_1EM6_R .OR. ABS(ystep) .LT. C_1EM6_R ) THEN

     DO jp = 1,nyp
        DO ip = 1,nxp
           x(ip,jp) = C_0_R
        ENDDO
     ENDDO

     DO jp = 1,nyp
        DO ip = 1,nxp
           y(ip,jp) = C_0_R
        ENDDO
     ENDDO

  ELSE

! Calculate x coordinate of histogram
! set y coordinate to 0.0

     DO jp = 1,nyp
        DO ip = 1,nxp
           x(ip,jp) = umin + xstep*(2*ip-1)/C_2_R
        ENDDO
     ENDDO

     DO jp = 1,nyp
        DO ip = 1,nxp
           y(ip,jp) = vmin + ystep*(2*jp-1)/C_2_R
        ENDDO
     ENDDO

! Calculate Histogram

     DO i=1, imax*jmax*kmax
        ip = INT((u(i)-umin)/xstep) + 1
        IF ( ip .LE. nxp .AND. ip .GE. 1 ) THEN
           jp = INT((v(i)-vmin)/ystep) + 1
           IF ( jp .LE. nyp .AND. jp .GE. 1 ) THEN
              pdf(ip,jp) = pdf(ip,jp) + C_1_R
              np = np + 1
           ENDIF
        ENDIF
     ENDDO

#ifdef USE_MPI
     CALL MPI_ALLREDUCE(np, np_total, 1, MPI_INTEGER4, MPI_SUM, &
          MPI_COMM_WORLD, ims_err)

     itmp = nxp*nyp
     CALL MPI_ALLREDUCE(pdf, wrk, itmp, MPI_REAL8, MPI_SUM, &
          MPI_COMM_WORLD, ims_err)

     DO ip = 1, nyp*nxp
        pdf(ip,1) = wrk(ip)
     ENDDO
     pnorm = C_1_R/(M_REAL(np_total)*xstep*ystep)
#else 
     pnorm = C_1_R/(M_REAL(np)*xstep*ystep)
#endif
     DO ip = 1, nyp*nxp
        pdf(ip,1) = pdf(ip,1)*pnorm
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE PDF2V3D


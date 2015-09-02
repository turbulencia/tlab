      SUBROUTINE PDF2V2D(imax, jmax, kmax, kmax_total, j, u, v, &
           nxp, nyp, x, y, pdf, wrk, ierr)

      IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

      TINTEGER imax, jmax, kmax, kmax_total
      TINTEGER i,j,k, ierr
      TINTEGER nxp, nyp
      TREAL u(imax,jmax,kmax)
      TREAL v(imax,jmax,kmax)
      TREAL x(nxp)
      TREAL y(nyp)
      TREAL pdf(nxp, nyp)
      TREAL wrk(nxp, nyp)

      TREAL umin, umax, vmin, vmax
      TREAL xstep, ystep, pnorm
      TINTEGER ip, jp
#ifdef USE_MPI
      INTEGER ims_err, itmp
      TREAL umin_p, umax_p, vmin_p, vmax_p
#endif

      ierr = 0

! Calculate Minimum and Maximum

      umin = u(1,j,1)
      umax = u(1,j,1)
      vmin = v(1,j,1)
      vmax = v(1,j,1)

      DO k=1, kmax
         DO i=1, imax
            umin = MIN(umin, u(i,j,k))
            umax = MAX(umax, u(i,j,k))
            vmin = MIN(vmin, v(i,j,k))
            vmax = MAX(vmax, v(i,j,k))
         ENDDO
      ENDDO

#ifdef USE_MPI
      CALL MPI_ALLREDUCE(vmin, vmin_p, 1, MPI_REAL8, MPI_MIN, &
           MPI_COMM_WORLD, ims_err)
      vmin = vmin_p
      CALL MPI_ALLREDUCE(vmax, vmax_p, 1, MPI_REAL8, MPI_MAX, &
           MPI_COMM_WORLD, ims_err)
      vmax = vmax_p
      CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, &
           MPI_COMM_WORLD, ims_err)
      umin = umin_p
      CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, &
           MPI_COMM_WORLD, ims_err)
      umax = umax_p
#endif

! Calculate Step in Histogram

      DO jp=1, nyp
         DO ip=1, nxp
            pdf(ip,jp) = C_0_R
         ENDDO
      ENDDO

      xstep = (umax-umin)/M_REAL(nxp)
      ystep = (vmax-vmin)/M_REAL(nyp)

      IF ( ABS(xstep) .LT. C_1EM6_R .OR. ABS(ystep) .LT. C_1EM6_R ) THEN

         DO ip=1, nxp
            x(ip) = C_0_R
         ENDDO

         DO jp=1, nyp
            y(jp) = C_0_R
         ENDDO

         ierr = 1

      ELSE
         
! Calculate x coordinate of histogram
! set y coordinate to 0.0

         DO ip=1, nxp
            x(ip) = umin + xstep*(2*ip-1)/C_2_R
         ENDDO

         DO jp=1, nyp
            y(jp) = vmin + ystep*(2*jp-1)/C_2_R
         ENDDO

! Calculate Histogram

         DO k=1, kmax
            DO i=1, imax
               ip = INT((u(i,j,k)-umin)/xstep) + 1
               ip = MAX(1,MIN(ip,nxp))
               jp = INT((v(i,j,k)-vmin)/ystep) + 1
               jp = MAX(1,MIN(jp,nyp))
               pdf(ip,jp) = pdf(ip,jp) + C_1_R
            ENDDO
         ENDDO

#ifdef USE_MPI
         itmp = nxp*nyp
         CALL MPI_ALLREDUCE(pdf, wrk, itmp, MPI_REAL8, MPI_SUM, &
              MPI_COMM_WORLD, ims_err)
         
         DO jp=1, nyp
            DO ip=1, nxp
               pdf(ip,jp) = wrk(ip,jp)
            ENDDO
         ENDDO
         
         pnorm = C_1_R/(M_REAL(imax*kmax_total)*xstep*ystep)
#else 
         pnorm = C_1_R/(M_REAL(imax*kmax)*xstep*ystep)
#endif

         DO jp=1, nyp
            DO ip=1, nxp
               pdf(ip,jp) = pdf(ip,jp)*pnorm
            ENDDO
         ENDDO

      ENDIF

      RETURN
      END


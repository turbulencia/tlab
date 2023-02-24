#include "types.h"

!#########################################################
!# DESCRIPTION
!#
!# Double top-hat filter coefficients initialization
!# Midpoint(trapezoidal) rule
!# The array cf is transposed wrt to the notes for later
!# better memory location
!# 
!#########################################################

SUBROUTINE FILTH2INI(iunif, ibc, imax, nx0, nx1, scalex, x, cf, wrk2d, wrk1d)

  IMPLICIT NONE

  TINTEGER iunif, ibc, imax
  TINTEGER nx0, nx1
  TREAL scalex
  TREAL x(imax)
  TREAL cf(*)
  TREAL wrk2d(imax,imax,3)
  TREAL wrk1d(imax,*)

! -----------------------------------------------------------------------
  TINTEGER i, j, ii, ip, jm, jj

! #######################################################################
  wrk2d = C_0_R

  IF ( iunif .EQ. 0 ) THEN

     ! ################
     ! # Uniform mesh #
     ! ################

     ! define matrix level 0
     CALL FILTHINI(iunif, ibc, imax, nx0, scalex, x,&
          wrk2d(1,1,3), wrk1d)
     DO i = 1,nx0/2
        ! boundary points
        ip = (i-1)*2 + 1
        wrk2d(i,1,1) = wrk2d(ip,1,3) 
        wrk2d(i,2,1) = wrk2d(ip+1,1,3)
        ! inner points
        DO ii = 3,i+nx0/2-1
           wrk2d(i,ii,1) = C_1_R
        ENDDO
        IF ( nx0 .GT. 2 ) THEN
           ii = i+nx0/2
           wrk2d(i,ii,1) = C_05_R
        ENDIF
     ENDDO

     DO i = nx0/2+1,imax-nx0/2
        wrk2d(i,i-nx0/2,1) = C_05_R
        DO ii = i-nx0/2+1,i+nx0/2-1
           wrk2d(i,ii,1) = C_1_R
        ENDDO
        wrk2d(i,i+nx0/2,1) = C_05_R
     ENDDO

     DO i = imax-nx0/2+1,imax
        ! inner points
        IF ( nx0 .GT. 2 ) THEN
           ii = i-nx0/2
           wrk2d(i,ii,1) = C_05_R
        ENDIF
        DO ii = i-nx0/2+1,imax-2
           wrk2d(i,ii,1) = C_1_R
        ENDDO
        ! boundary points
        wrk2d(i,imax,1)   = wrk2d(imax-i+i1,1,1)
        wrk2d(i,imax-1,1) = wrk2d(imax-i+i1,2,1)
     ENDDO

     ! define matrix level1
     CALL FILTHINI(iunif, ibc, imax, nx1, scalex, x,&
          wrk2d(1,1,3), wrk1d)
     DO i = 1,nx1/2
        ! boundary points
        ip = (i-1)*2 + 1
        wrk2d(i,1,2) = wrk2d(ip,1,3) 
        wrk2d(i,2,2) = wrk2d(ip+1,1,3)
        ! inner points
        DO ii = 3,i+nx1/2-1
           wrk2d(i,ii,2) = C_1_R
        ENDDO
        IF ( nx1 .GT. 2 ) THEN
           ii = i+nx1/2
           wrk2d(i,ii,2) = C_05_R
        ENDIF
     ENDDO

     DO i = nx1/2+1,imax-nx1/2
        wrk2d(i,i-nx1/2,2) = C_05_R
        DO ii = i-nx1/2+1,i+nx1/2-1
           wrk2d(i,ii,2) = C_1_R
        ENDDO
        wrk2d(i,i+nx1/2,2) = C_05_R
     ENDDO

     DO i = imax-nx1/2+1,imax
        ! inner points
        IF ( nx1 .GT. 2 ) THEN
           ii = i-nx1/2
           wrk2d(i,ii,2) = C_05_R
        ENDIF
        DO ii = i-nx1/2+1,imax-2
           wrk2d(i,ii,2) = C_1_R
        ENDDO
        ! boundary points
        wrk2d(i,imax,2)   = wrk2d(imax-i+i1,1,2)
        wrk2d(i,imax-1,2) = wrk2d(imax-i+i1,2,2)
     ENDDO

     ! multiply matrices
     IF ( ibc .EQ. 0 ) THEN
        DO ip = 1,nx0+nx1+1
           cf(ip) = C_0_R
           DO ii = 1,imax
              cf(ip) = cf(ip)+wrk2d(nx0/2+nx1/2+1,ii,2)*&
                   wrk2d(ii,ip,1)
           ENDDO
           cf(ip) = cf(ip)/M_REAL(nx0*nx1)            
        ENDDO
     ELSE
        DO i = 1,nx0/2+nx1/2+1
           DO j = 1,nx0+nx1+1
              ip = (i-1)*(nx0+nx1+1) + j
              cf(ip) = C_0_R
              DO ii = 1,imax
                 cf(ip) = cf(ip)+wrk2d(i,ii,2)*wrk2d(ii,j,1)
              ENDDO
              cf(ip) = cf(ip)/M_REAL(nx0*nx1)
           ENDDO
        ENDDO
     ENDIF

  ELSE

     ! ###################
     ! # Nonuniform mesh #
     ! ###################

     ! level0
     CALL FILTHINI(iunif, ibc, imax, nx0, scalex, x,&
          wrk2d(1,1,3), wrk1d)
     DO i = 1,imax
        DO j = 1,nx0+1
           jj = i-nx0/2-1+j
           jm = jj+imax-1
           jm = MOD(jm,imax)+1
           ip = (i-1)*(nx0+1) + j
           wrk2d(i,jm,1) = wrk2d(ip,1,3)
        ENDDO
     ENDDO
     ! level1
     CALL FILTHINI(iunif, ibc, imax, nx1, scalex, x,&
          wrk2d(1,1,3), wrk1d)
     DO i = 1,imax
        DO j = 1,nx1+1
           jj = i-nx1/2-1+j
           jm = jj+imax-1
           jm = MOD(jm,imax)+1
           ip = (i-1)*(nx1+1) + j
           wrk2d(i,jm,2) = wrk2d(ip,1,3)
        ENDDO
     ENDDO
     ! multiply matrices
     DO i = 1,imax
        DO j = 1,imax
           wrk2d(i,j,3) = C_0_R
           DO ii = 1,imax
              wrk2d(i,j,3) = wrk2d(i,j,3) + &
                   wrk2d(i,ii,2)*wrk2d(ii,j,1)
           ENDDO
        ENDDO
     ENDDO
     DO i = 1,imax
        DO j = 1,nx0+nx1+1
           jj = i-(nx0+nx1)/2-1+j
           jm = jj+imax-1
           jm = MOD(jm,imax)+1
           !           ip = (j-1)*imax + i
           ip = (i-1)*(nx0+nx1+1) + j
           cf(ip) = wrk2d(i,jm,3)
        ENDDO
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE FILTH2INI



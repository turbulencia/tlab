#include "types.h"

!#########################################################
!# DESCRIPTION
!#
!# Top-hat filter coefficients initialization
!# Midpoint(trapezoidal) rule
!# The array cf is transposed wrt to the notes for later
!# better memory location
!# 
!#########################################################

SUBROUTINE FILTHINI(iunif, ibc, imax, nx, scalex, x, cf, wrk1d)

  IMPLICIT NONE

  TINTEGER iunif, ibc, imax
  TINTEGER nx
  TREAL scalex
  TREAL x(imax)
  TREAL cf(*)
  TREAL wrk1d(imax,3)

! -----------------------------------------------------------------------
  TREAL dum
  TINTEGER i, ii, im, ic, ip

#include "integers.h"

! #######################################################################

  wrk1d = C_0_R

  IF ( iunif .EQ. 0 ) THEN

     ! ################
     ! # Uniform mesh #
     ! ################

     IF ( ibc .NE. 0 ) THEN
        DO i = 1,nx/2
           im = nx/2-i+1
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ip = (i-1)*2 + 1
           cf(ip) = dum+C_05_R*M_REAL(im+1)
           !           if ( cf(ip) .lt. C_0_R ) then
           !              print *, iunif, ibc, imax, cf(ip)
           !              stop
           !           endif
           IF ( nx .EQ. 2 ) THEN
              cf(ip+1) = C_05_R-dum+C_05_R*M_REAL(im)
           ELSE
              cf(ip+1) = C_1_R -dum+C_05_R*M_REAL(im)
           ENDIF
           !           if ( cf(ip+1) .lt. C_0_R ) then
           !              print *, iunif, ibc, imax, cf(ip+1)
           !              stop
           !           endif
        ENDDO
     ENDIF
  ELSE

     ! ###################
     ! # Nonuniform mesh #
     ! ###################

     ! delta(i)
     DO i = 1,imax-1
        wrk1d(i,1) = x(i+1)-x(i)
     ENDDO
     IF ( ibc .EQ. 0 ) THEN
        wrk1d(imax,1) = scalex-(x(imax)-x(1))
     ELSE
        wrk1d(imax,1) = wrk1d(imax-1,1)
     ENDIF
     ! deltasum(i) = delta(i-1)+delta(i)
     DO i = 2,imax
        wrk1d(i,2) = wrk1d(i-1,1)+wrk1d(i,1)
     ENDDO
     IF ( ibc .EQ. 0 ) THEN
        wrk1d(1,2) = wrk1d(imax,1)+wrk1d(1,1)
     ELSE
        wrk1d(1,2) = wrk1d(1,1)*C_2_R
     ENDIF
     ! deltaf(i)
     IF ( ibc .EQ. 0 ) THEN
        DO i = 1,imax
           wrk1d(i,3) = C_0_R
           DO ii = i-nx/2,i+nx/2-1
              im = ii+imax-1
              im = MOD(im,imax)+1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(im,1)
           ENDDO
        ENDDO
     ELSE
        DO i = 1,nx/2
           wrk1d(i,3) = wrk1d(1,1)*M_REAL(nx/2-i+1)
           DO ii = 1,i+nx/2-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
        DO i = 1+nx/2,imax-nx/2
           DO ii = i-nx/2,i+nx/2-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
        DO i = imax-nx/2+1,imax
           wrk1d(i,3) = wrk1d(imax,1)*M_REAL(nx/2-(imax-i))
           DO ii = i-nx/2,imax-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
     ENDIF

     !    construct the coefficients array
     DO i = 1,imax
        ii = i-nx/2
        im = ii+imax-1
        im = MOD(im,imax)+1
        ic = ii-i+nx/2+1
        ip = (i-1)*(nx+1) + ic
        cf(ip) = C_05_R*wrk1d(im,1)/wrk1d(i,3)
        if ( cf(ip) .lt. C_0_R ) then
           print *, iunif, ibc, imax, cf(ip), '(1)'
           !           stop
        endif
        DO ii = i-nx/2+1,i+nx/2-1
           im = ii+imax-1
           im = MOD(im,imax)+1
           ic = ii-i+nx/2+1
           ip = (i-1)*(nx+1) + ic
           cf(ip) = C_05_R*wrk1d(im,2)/wrk1d(i,3)
           if ( cf(ip) .lt. C_0_R ) then
              print *, iunif, ibc, imax, cf(ip), '(2)'
              !              stop
           endif
        ENDDO
        ii = i+nx/2
        !        im = ii+imax-1
        im = (ii-1)+imax-1
        im = MOD(im,imax)+1
        !       if ( im .eq. 1 ) im = 2
        ic = ii-i+nx/2+1
        ip = (i-1)*(nx+1) + ic
        !        cf(ip) = C_05_R*wrk1d(im-1,1)/wrk1d(i,3)
        cf(ip) = C_05_R*wrk1d(im,1)/wrk1d(i,3)
        !        if ( cf(ip) .lt. C_0_R ) then
        !           print *, iunif, ibc, imax, cf(ip), '(3)'
        !           stop
        !        endif
     ENDDO

     !    modification in case of free boundary
     IF ( ibc .NE. 0 ) THEN
        DO i = 1,nx/2
           im = nx/2-i+1
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ic = nx/2-i+2
           ip = (i-1)*(nx+1) + ic
           cf(ip) = cf(ip) + &
                C_05_R*(wrk1d(1,2)*(dum-C_1_R) + &
                wrk1d(1,1)*M_REAL(im+1))/wrk1d(i,3)
           !           if ( cf(ip) .lt. C_0_R ) then
           !             print *, iunif, ibc, imax, cf(ip), '(4)'
           !              stop
           !           endif
           ic = nx/2-i+3
           ip = (i-1)*(nx+1) + ic
           cf(ip) = cf(ip) - &
                C_05_R*(wrk1d(1,2)*(dum-M_REAL(im)) +&
                wrk1d(1,1)*M_REAL(im))/wrk1d(i,3)
           !           if ( cf(ip) .lt. C_0_R ) then
           !              print *, iunif, ibc, imax, cf(ip), '(5)'
           !              stop
           !           endif
           ! pad with ceros
           DO ic=1,nx/2-i+1
              ip = (i-1)*(nx+1) + ic
              cf(ip) = C_0_R
           ENDDO
        ENDDO
        DO i = imax-nx/2+1,imax
           im = i-imax+nx/2
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ic = nx+1-im-1
           ip = (i-1)*(nx+1) + ic
           cf(ip) = cf(ip) - &
                C_05_R*(wrk1d(imax,2)*(dum-M_REAL(im)) +&
                wrk1d(imax,1)*M_REAL(im))/wrk1d(i,3)
           !           if ( cf(ip) .lt. C_0_R ) then
           !              print *, iunif, ibc, imax, cf(ip), '(6)'
           !              stop
           !           endif
           ic = nx+1-im
           ip = (i-1)*(nx+1) + ic
           cf(ip) = cf(ip) + &
                C_05_R*(wrk1d(imax,2)*(dum-C_1_R) + &
                wrk1d(imax,1)*M_REAL(im+1))/wrk1d(i,3)
           !           if ( cf(ip) .lt. C_0_R ) then
           !              print *, iunif, ibc, imax, cf(ip), '(7)'
           !              stop
           !           endif
           ! pad with ceros
           DO ic=nx+1-im+1,nx+1
              ip = (i-1)*(nx+1) + ic
              cf(ip) = C_0_R
           ENDDO
        ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FILTHINI

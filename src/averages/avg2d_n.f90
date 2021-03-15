#include "types.h"
#include "dns_error.h"

#define NMOMS_MAX 10
! the next macro is for the format (NVARS x NMOMS)
#define L_FORMAT_MAX 40

!########################################################################
!# HISTORY
!#
!# 2008/04/10 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the first nmom moments of the nvar fields defined by the
!# pointers array vars
!#
!########################################################################
SUBROUTINE AVG2D_N(fname, rtime, imax,jmax,kmax, &
     nvar, nmom, vars, igate, gate, y, avg_loc)

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TREAL rtime
  TINTEGER,                 INTENT(IN) :: imax,jmax,kmax
  TINTEGER,                 INTENT(IN) :: nvar  ! Number of variables to process
  TINTEGER,                 INTENT(IN) :: nmom  ! Number of moments to consider in the analysis
  INTEGER(1),               INTENT(IN) :: igate ! Gate level in gate array to be used. If 0, no intermittency considered
  INTEGER(1), DIMENSION(*), INTENT(IN) :: gate  ! Array with the mask field corresponding to the different gate levels
  TYPE(pointers_dt), DIMENSION(nvar)   :: vars ! Array of pointer to the fields to be processed
  TREAL y(jmax)
  TREAL avg_loc(nmom,nvar,2)

! -------------------------------------------------------------------
  TINTEGER k, j, i, iv, im
  TREAL AVG1V2D, AVG1V2D1G, AVG1V3D, AVG1V3D1G
  TREAL moment(NMOMS_MAX), num, den

  CHARACTER*32 str
  CHARACTER*1024 line1
  CHARACTER*2048 line2

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

! ###################################################################
  IF ( nvar*nmom .GT. L_FORMAT_MAX ) THEN
     CALL IO_WRITE_ASCII(efile,'AVG2D_N. Format length too short.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! -------------------------------------------------------------------
! TkStat file; header
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     OPEN(unit=21,file=fname)

     WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
     WRITE(21, '(A7,I8)') 'IMAX = ', i1
     WRITE(21, '(A7,I8)') 'JMAX = ', jmax

     line2 = 'I J Y'

     line1 = ' '
     DO iv = 1,nvar
        DO im = 1,nmom
           WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(vars(iv)%tag))//'Mom'//TRIM(ADJUSTL(str))
        ENDDO
     ENDDO
     WRITE(21,'(A)') 'GROUP = Plane '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = ' '
     DO iv = 1,nvar
        DO im = 1,nmom
           WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(vars(iv)%tag))//'Mom'//TRIM(ADJUSTL(str))//'Vol'
        ENDDO
     ENDDO
     WRITE(21,'(A)') 'GROUP = Volume '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     WRITE(21,'(A)') TRIM(ADJUSTL(line2))

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Volume averages; storage in array avg_vol
! ###################################################################
  DO iv = 1,nvar
     DO im = 1,nmom
        IF ( igate .GT. 0 ) THEN
           moment(im) = AVG1V3D1G(imax,jmax,kmax, igate, im, vars(iv)%field, gate)
        ELSE
           moment(im) = AVG1V3D(imax,jmax,kmax, im, vars(iv)%field)
        ENDIF

! -------------------------------------------------------------------
! calculate central moments using binomial theorem
! -------------------------------------------------------------------
        IF ( im .GT. 1 ) THEN
           avg_loc(im,iv,2) = C_0_R
           DO k = 0,im-1
              num = C_1_R
              DO i = im,im-k+1,-1
                 num = num*M_REAL(i)
              ENDDO
              den = C_1_R
              DO i = k,1,-1
                 den = den*M_REAL(i)
              ENDDO
              avg_loc(im,iv,2) = avg_loc(im,iv,2) + num/den*moment(im-k)*((-moment(1))**k)
           ENDDO
           avg_loc(im,iv,2) = avg_loc(im,iv,2) + (-moment(1))**im

        ELSE
           avg_loc(im,iv,2) = moment(im)

        ENDIF

     ENDDO
  ENDDO

! ###################################################################
! Average calculation of 1 variable in along planes
! ###################################################################
  DO j = 1,jmax

     DO iv = 1,nvar
        DO im = 1,nmom
           IF ( igate .GT. 0 ) THEN
              moment(im) = AVG1V2D1G(imax,jmax,kmax, j, igate, im, vars(iv)%field, gate)
           ELSE
              moment(im) = AVG1V2D(imax,jmax,kmax, j, im, vars(iv)%field)
           ENDIF

! -------------------------------------------------------------------
! calculate central moments using binomial theorem
! -------------------------------------------------------------------
           IF ( im .GT. 1 ) THEN
              avg_loc(im,iv,1) = C_0_R
              DO k = 0,im-1
                 num = C_1_R
                 DO i = im,im-k+1,-1
                    num = num*M_REAL(i)
                 ENDDO
                 den = C_1_R
                 DO i = k,1,-1
                    den = den*M_REAL(i)
                 ENDDO
                 avg_loc(im,iv,1) = avg_loc(im,iv,1) + num/den*moment(im-k)*((-moment(1))**k)
              ENDDO
              avg_loc(im,iv,1) = avg_loc(im,iv,1) + (-moment(1))**im

           ELSE
              avg_loc(im,iv,1) = moment(im)

           ENDIF

        ENDDO
     ENDDO

! -------------------------------------------------------------------
! TkStat output
! -------------------------------------------------------------------
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif

        IF ( j .EQ. jmax/2 ) THEN
           WRITE(21,1020) 1, j, y(j), (avg_loc(im,1,1),im=1,nmom*nvar), &
                (avg_loc(im,1,2),im=1,nmom*nvar)
        ELSE
           WRITE(21,1010) 1, j, y(j), (avg_loc(im,1,1),im=1,nmom*nvar)
        ENDIF

#ifdef USE_MPI
     ENDIF
#endif

  ENDDO

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CLOSE(21)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(I5,(1X,I5),L_FORMAT_MAX(1X,G_FORMAT_R))
1020 FORMAT(I5,(1X,I5),L_FORMAT_MAX(1X,G_FORMAT_R),L_FORMAT_MAX(1X,G_FORMAT_R))

END SUBROUTINE AVG2D_N

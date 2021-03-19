#include "types.h"
#include "dns_error.h"

#define NMOMS_MAX 10
! the next macro is for the format (NVARS x NMOMS)
#define L_FORMAT_MAX 40

!########################################################################
!#
!# Calculating the first nm moments of the nv fields defined by the pointers array vars
!#
!########################################################################
SUBROUTINE AVG2D_N(fname, rtime, nx,ny,nz, nv, nm, vars, igate,gate, y, avg)

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile, lfile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TREAL rtime
  TINTEGER,          INTENT(IN   ) :: nx,ny,nz, nv, nm ! Number of moments to consider in the analysis
  TYPE(pointers_dt), INTENT(IN   ) :: vars(nv)         ! Array of pointer to the fields to be processed
  INTEGER(1),        INTENT(IN   ) :: gate(*), igate   ! discrete conditioning criteria
  TREAL,             INTENT(IN   ) :: y(ny)            ! heights of each plane
  TREAL,             INTENT(  OUT) :: avg(nm,nv,ny)

  ! -------------------------------------------------------------------
  TINTEGER j,iv,im
  TREAL AVG1V2D, AVG1V2D1G

  CHARACTER*32 str
  CHARACTER*1024 line1
  CHARACTER*2048 line2

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  IF ( nv*nm > L_FORMAT_MAX ) THEN
    CALL IO_WRITE_ASCII(efile,'AVG2D_N. Format length too short.')
    CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  END IF

  DO j = 1,ny                   ! calculation in planes
    DO iv = 1,nv
      DO im = 1,nm
        IF ( igate > 0 ) THEN
          avg(im,iv,j) = AVG1V2D1G(nx,ny,nz, j, igate, im, vars(iv)%field, gate)
        ELSE
          avg(im,iv,j) = AVG1V2D(nx,ny,nz, j, im, vars(iv)%field)
        END IF
      END DO
      CALL RAW_TO_CENTRAL( nm, avg(1,iv,j) )
    END DO
  END DO

  IF ( ny > 1 ) THEN            ! calculation in whole volume, saved as plane j=ny+1
    DO iv = 1,nv
      DO im = 1,nm
        IF ( igate > 0 ) THEN
          ! avg(im,iv,j) = AVG1V3D1G(nx,ny,nz, igate, im, vars(iv)%field, gate)
          avg(im,iv,j) = AVG1V2D1G(nx*ny,1,nz, 1, igate, im, vars(iv)%field, gate)
        ELSE
          ! avg(im,iv,j) = AVG1V3D(nx,ny,nz, im, vars(iv)%field)
          avg(im,iv,j) = AVG1V2D(nx*ny,1,nz, 1, im, vars(iv)%field)
        END IF
      END DO
      CALL RAW_TO_CENTRAL( nm, avg(1,iv,j) )
    END DO
  END IF

  ! ###################################################################
  ! -------------------------------------------------------------------
  ! TkStat file
  ! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif
    OPEN(unit=21,file=fname)

    WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
    WRITE(21, '(A7,I8)') 'IMAX = ', i1
    WRITE(21, '(A7,I8)') 'JMAX = ', ny

    line2 = 'I J Y'

    line1 = ' '
    DO iv = 1,nv
      DO im = 1,nm
        WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(vars(iv)%tag))//'Mom'//TRIM(ADJUSTL(str))
      END DO
    END DO
    WRITE(21,'(A)') 'GROUP = Plane '//TRIM(ADJUSTL(line1))
    line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

    line1 = ' '
    DO iv = 1,nv
      DO im = 1,nm
        WRITE(str,*) im; line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(vars(iv)%tag))//'Mom'//TRIM(ADJUSTL(str))//'Vol'
      END DO
    END DO
    WRITE(21,'(A)') 'GROUP = Volume '//TRIM(ADJUSTL(line1))
    line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

    WRITE(21,'(A)') TRIM(ADJUSTL(line2))

    DO j = 1,ny
      IF ( j == ny/2 ) THEN
        WRITE(21,1020) 1, j, y(j), (avg(im,1,j),im=1,nm*nv), &
            (avg(im,1,ny+1),im=1,nm*nv)
      ELSE
        WRITE(21,1010) 1, j, y(j), (avg(im,1,j),im=1,nm*nv)
      END IF
    END DO

    CLOSE(21)
#ifdef USE_MPI
  END IF
#endif

  RETURN

1010 FORMAT(I5,(1X,I5),L_FORMAT_MAX(1X,G_FORMAT_R))
1020 FORMAT(I5,(1X,I5),L_FORMAT_MAX(1X,G_FORMAT_R),L_FORMAT_MAX(1X,G_FORMAT_R))

END SUBROUTINE AVG2D_N

! ###################################################################
SUBROUTINE RAW_TO_CENTRAL( nm, moments )
  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: nm
  TREAL,    INTENT(INOUT) :: moments(nm)

  TINTEGER im, i,k
  TREAL aux(NMOMS_MAX), num, den

  ! -------------------------------------------------------------------
  aux(1:nm) = moments(1:nm)

  DO im = 2,nm
    moments(im) = C_0_R
    DO k = 0,im-1
      num = C_1_R
      DO i = im,im-k+1,-1
        num = num*M_REAL(i)
      END DO
      den = C_1_R
      DO i = k,1,-1
        den = den*M_REAL(i)
      END DO
      moments(im) = moments(im) + num/den*aux(im-k)*((-aux(1))**k)
    END DO
    moments(im) = moments(im) + (-aux(1))**im
  END DO

  RETURN
END SUBROUTINE RAW_TO_CENTRAL

!########################################################################
!#
!# Intermittency factors, i.e., area fraction occupied by each gate level
!#
!########################################################################
SUBROUTINE INTER2D_N(fname, parname, rtime, nx,ny,nz, npar, y, gate, inter)

  USE DNS_CONSTANTS, ONLY : efile, lfile

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname, parname(npar)
  TREAL rtime
  TINTEGER,   INTENT(IN   ) :: nx,ny,nz, npar ! npar is the number of partitions in gate field
  TREAL,      INTENT(IN   ) :: y(ny)          ! heights of each plane
  INTEGER(1), INTENT(IN   ) :: gate(*)        ! field with partitions
  TREAL,      INTENT(  OUT) :: inter(npar,ny) ! intermittency factor

  ! -------------------------------------------------------------------
  TINTEGER ip, j
  TREAL INTER1V2D
  INTEGER(1) gate_level

  CHARACTER*512 line1

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  IF ( npar > L_FORMAT_MAX ) THEN
    CALL IO_WRITE_ASCII(efile,'INTER2D_N. Format length too short.')
    CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  END IF

  DO j = 1,ny
    DO ip = 1,npar
      gate_level = INT(ip,KIND=1)
      inter(ip,j) = INTER1V2D(nx,ny,nz, j, gate_level, gate)
    ENDDO
  END DO

  ! ###################################################################
  ! -------------------------------------------------------------------
  ! TkStat file
  ! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
    OPEN(unit=21,file=fname)

    WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
    WRITE(21, '(A7,I8)') 'IMAX = ', 1
    WRITE(21, '(A7,I8)') 'JMAX = ', ny

    line1 = 'I J Y'
    DO ip = 1,npar
      line1 = TRIM(ADJUSTL(line1))//' '//TRIM(ADJUSTL(parname(ip)))
    ENDDO
    WRITE(21,'(A)') TRIM(ADJUSTL(line1))

    DO j = 1,ny
      WRITE(21,1030) 1, j, y(j), (inter(ip,j),ip=1,npar)
    ENDDO

    CLOSE(21)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1030 FORMAT(I5,(1X,I5),(1X,G_FORMAT_R),L_FORMAT_MAX(1X,G_FORMAT_R))

END SUBROUTINE INTER2D_N

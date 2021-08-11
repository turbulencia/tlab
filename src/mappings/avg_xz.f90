#include "types.h"
#include "dns_error.h"

#define NMOMS_MAX 10

!########################################################################
!#
!# Calculating the first nm moments of the nv fields defined by the pointers array vars
!#
!########################################################################
SUBROUTINE AVG_N_XZ(fname, itime,rtime, nx,ny,nz, nv,nm, vars, igate,gate, y, avg)

  USE TLAB_TYPES,     ONLY : pointers_dt
  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_PROCS

  IMPLICIT NONE

  CHARACTER*(*) fname
  TINTEGER itime
  TREAL rtime
  TINTEGER,          INTENT(IN   ) :: nx,ny,nz, nv, nm ! Number of moments to consider in the analysis
  TYPE(pointers_dt), INTENT(IN   ) :: vars(nv)         ! Array of pointer to the fields to be processed
  INTEGER(1),        INTENT(IN   ) :: gate(*), igate   ! discrete conditioning criteria
  TREAL,             INTENT(IN   ) :: y(ny)            ! heights of each plane
  TREAL,             INTENT(  OUT) :: avg(ny,nm,nv)

  ! -------------------------------------------------------------------
  TINTEGER j,iv,im
  TREAL AVG1V2D, AVG1V2D1G, moments(nm)

  CHARACTER*32 str
  CHARACTER*32 groupname(1) ! 3 groups for consistency with old TkStat format
  CHARACTER*250 varname(1)  ! to be reduce to just 1 group in the future

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO j = 1,ny
    DO iv = 1,nv
      DO im = 1,nm
        IF ( igate > 0 ) THEN
          moments(im) = AVG1V2D1G( nx,ny,nz, j, igate, im, vars(iv)%field, gate )
        ELSE
          moments(im) = AVG1V2D( nx,ny,nz, j, im, vars(iv)%field )
        END IF
      END DO
      IF ( nm > 1 ) CALL RAW_TO_CENTRAL( nm, moments )
      avg(j,1:nm,iv) = moments(1:nm)
    END DO
  END DO

  ! ###################################################################
  ! Output
  ! ###################################################################
  varname(:) = ''
  groupname(:) = ''
  DO iv = 1,nv
    DO im = 1,nm
      varname(1) = TRIM(ADJUSTL(varname(1)))//' '//TRIM(ADJUSTL(vars(iv)%tag))
      IF ( im > 1 ) THEN
        WRITE(str,*) im; varname(1) = TRIM(ADJUSTL(varname(1)))//'.'//TRIM(ADJUSTL(str))
      ENDIF
    END DO
  END DO
  ! for consistency with old format I pass 2 groups; no longer (Cedrick Jul 3rd 2021)
  CALL IO_WRITE_AVERAGES( fname, itime,rtime, ny,nv*nm,1, y, varname, groupname, avg )

  RETURN
END SUBROUTINE AVG_N_XZ

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
SUBROUTINE INTER_N_XZ(fname, itime,rtime, nx,ny,nz, np, parname, gate, y, inter)

  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_PROCS

  IMPLICIT NONE

  CHARACTER*(*) fname, parname(np)
  TINTEGER itime
  TREAL rtime
  TINTEGER,   INTENT(IN   ) :: nx,ny,nz, np   ! npar is the number of partitions in gate field
  TREAL,      INTENT(IN   ) :: y(ny)          ! heights of each plane
  INTEGER(1), INTENT(IN   ) :: gate(*)        ! field with partitions
  TREAL,      INTENT(  OUT) :: inter(ny,np)   ! intermittency factor

  ! -------------------------------------------------------------------
  TINTEGER ip, j
  TREAL INTER1V2D
  INTEGER(1) gate_level

  CHARACTER*32 groupname(1) ! Two groups for consistency with old TkStat format
  CHARACTER*250 varname(1)  ! to be reduce to just 1 group in the future

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO j = 1,ny
    DO ip = 1,np
      gate_level = INT(ip,KIND=1)
      inter(j,ip) = INTER1V2D(nx,ny,nz, j, gate_level, gate)
    END DO
  END DO

  ! ###################################################################
  ! Output
  ! ###################################################################
  varname(:) = ''
  groupname(:) = ''
  DO ip = 1,np
    varname(1) = TRIM(ADJUSTL(varname(1)))//' '//TRIM(ADJUSTL(parname(ip)))
  ENDDO
  ! for consistency with old format I pass 2 groups
  CALL IO_WRITE_AVERAGES( fname, itime,rtime, ny,np,1, y, varname, groupname, inter )

  RETURN
END SUBROUTINE INTER_N_XZ

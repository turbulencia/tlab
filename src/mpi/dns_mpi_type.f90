#include "types.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2010/04/12 - J. P. Mellado
!#              Modified to set homogeneous decomposition
!#
!########################################################################
!# DESCRIPTION
!#
!# Pointers and types for transposition across ims_npro processors
!#
!########################################################################
SUBROUTINE DNS_MPI_TYPE_I(ims_npro_i, imax, npage, nd, md, n1, n2, &
     nsize, sdisp, rdisp, stype, rtype)

  USE DNS_CONSTANTS, ONLY : efile
  USE TLAB_CORE

  IMPLICIT NONE

#include "mpif.h"

  INTEGER ims_npro_i
  TINTEGER npage, imax, nsize
  TINTEGER nd, md, n1, n2
  TINTEGER sdisp(*), rdisp(*)
  INTEGER, INTENT(OUT) ::  stype, rtype

! -----------------------------------------------------------------------
  TINTEGER i
  INTEGER ims_ss, ims_rs, ims_err
  INTEGER ims_tmp1, ims_tmp2, ims_tmp3
  CHARACTER*64 str, line

! #######################################################################
  IF ( MOD(npage,ims_npro_i) .EQ. 0 ) THEN
     nsize = npage/ims_npro_i
  ELSE
     CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_TYPE_I. Ratio npage/ims_npro_i not an integer.')
     CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
  ENDIF

! Calculate Displacements in Forward Send/Receive
  sdisp(1) = 0
  rdisp(1) = 0
  DO i = 2,ims_npro_i
     sdisp(i) = sdisp(i-1) + imax *nd *nsize
     rdisp(i) = rdisp(i-1) + imax *md
  ENDDO

! #######################################################################
  !DO i = 1,1 !ims_npro_i

     ims_tmp1 = nsize *n1 ! count
     ims_tmp2 = imax  *n2 ! block
     ims_tmp3 = ims_tmp2  ! stride = block because things are together
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
     CALL MPI_TYPE_COMMIT(stype, ims_err)

     ims_tmp1 = nsize           *n1 ! count
     ims_tmp2 = imax            *n2 ! block
     ims_tmp3 = imax*ims_npro_i *n2 ! stride is a multiple of imax_total=imax*ims_npro_i
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
     CALL MPI_TYPE_COMMIT(rtype, ims_err)

     CALL MPI_TYPE_SIZE(stype, ims_ss, ims_err)
     CALL MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

     IF ( ims_ss .NE. ims_rs ) THEN
        WRITE(str, *) ims_ss; WRITE(line,*) ims_rs
        line='Send size '//TRIM(ADJUSTL(str))//'differs from recv size '//TRIM(ADJUSTL(line))
        WRITE(str, *) 1  ! i
        line=TRIM(ADJUSTL(line))//' in message '//TRIM(ADJUSTL(str))
        CALL TLAB_WRITE_ASCII(efile, line)
        CALL TLAB_STOP(DNS_ERROR_MPITYPECHECK)
     ENDIF

  !ENDDO

  RETURN
END SUBROUTINE DNS_MPI_TYPE_I

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TYPE_K(ims_npro, nmax, npage, nd, md, n1, n2, &
     nsize, sdisp, rdisp, stype, rtype)

  USE DNS_CONSTANTS, ONLY : efile
  USE TLAB_CORE

  IMPLICIT NONE

#include "mpif.h"

  INTEGER ims_npro
  TINTEGER npage, nmax, nsize
  TINTEGER nd, md, n1, n2
  TINTEGER sdisp(*), rdisp(*)
  INTEGER, INTENT(OUT) ::  stype, rtype

! -----------------------------------------------------------------------
  TINTEGER i
  INTEGER ims_ss, ims_rs, ims_err
  INTEGER ims_tmp1, ims_tmp2, ims_tmp3

! #######################################################################
  IF ( MOD(npage,ims_npro) .EQ. 0 ) THEN
     nsize = npage/ims_npro
  ELSE
     CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_TYPE_K. Ratio npage/ims_npro not an integer.')
     CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
  ENDIF

! Calculate Displacements in Forward Send/Receive
  sdisp(1) = 0
  rdisp(1) = 0
  DO i = 2,ims_npro
     sdisp(i) = sdisp(i-1) + nsize *nd
     rdisp(i) = rdisp(i-1) + nsize *md *nmax
  ENDDO

! #######################################################################
  !DO i = 1, 1 !ims_npro

     ims_tmp1 = nmax  *n1 ! count
     ims_tmp2 = nsize *n2 ! block
     ims_tmp3 = npage *n2 ! stride
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
     CALL MPI_TYPE_COMMIT(stype, ims_err)

     ims_tmp1 = nmax  *n1 ! count
     ims_tmp2 = nsize *n2 ! block
     ims_tmp3 = ims_tmp2  ! stride = block to put things together
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
     CALL MPI_TYPE_COMMIT(rtype, ims_err)

     CALL MPI_TYPE_SIZE(stype, ims_ss, ims_err)
     CALL MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

     IF ( ims_ss .NE. ims_rs ) THEN
        PRINT *, 'Message   : ', 1, ' size is wrong' ! i
        PRINT *, 'Send size : ', ims_ss
        PRINT *, 'Recv size : ', ims_rs
        CALL TLAB_STOP(DNS_ERROR_MPITYPECHECK)
     ENDIF

  ! ENDDO

  RETURN
END SUBROUTINE DNS_MPI_TYPE_K

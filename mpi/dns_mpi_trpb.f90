#include "types.h"

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created 
!# 2009/09/01 - J.P. Mellado
!#              Debugged
!# 2015/03/14 - J.P. Mellado
!#              Using local communicator
!#
!########################################################################
!# DESCRIPTION
!#
!# MPI_ALLTOALLW has been tested but did not work yet
!#
!########################################################################
SUBROUTINE DNS_MPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)

  USE DNS_MPI, ONLY : ims_time_trans
  USE DNS_MPI, ONLY : ims_npro_k, ims_pro_k
  USE DNS_MPI, ONLY : ims_comm_z
  USE DNS_MPI, ONLY : ims_tag, ims_err

  IMPLICIT NONE
  
#include "mpif.h"
  
  TREAL,    DIMENSION(*),          INTENT(IN)  :: b
  TREAL,    DIMENSION(*),          INTENT(OUT) :: a
  TINTEGER, DIMENSION(ims_npro_k), INTENT(IN)  :: dsend, drecv
  INTEGER,  DIMENSION(ims_npro_k), INTENT(IN)  :: tsend, trecv
  
! -----------------------------------------------------------------------
  TINTEGER n, l
  INTEGER status(MPI_STATUS_SIZE,2*ims_npro_k)
  INTEGER mpireq(                2*ims_npro_k)
  INTEGER ip

#ifdef PROFILE_ON
  TREAL time_loc_1, time_loc_2
#endif  

! #######################################################################
#ifdef PROFILE_ON
  time_loc_1 = MPI_WTIME()
#endif

! #######################################################################
! Same processor
! #######################################################################
  ip = ims_pro_k; n = ip + 1
  CALL MPI_ISEND(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(1), ims_err)
  CALL MPI_IRECV(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 2
  DO n = 1,ims_npro_k
     ip = n - 1
     IF ( ip .NE. ims_pro_k ) THEN
        l = l + 1
        CALL MPI_ISEND(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        l = l + 1
        CALL MPI_IRECV(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro_k*2-2, mpireq(3), status(1,3), ims_err)

  CALL DNS_MPI_TAGUPDT

#ifdef PROFILE_ON
  time_loc_2 = MPI_WTIME()
  ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif
  
  RETURN
END SUBROUTINE DNS_MPI_TRPB_K

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TRPB_I(b, a, dsend, drecv, tsend, trecv)

  USE DNS_MPI, ONLY : ims_npro_i, ims_pro_i
  USE DNS_MPI, ONLY : ims_comm_x
  USE DNS_MPI, ONLY : ims_tag, ims_err

  IMPLICIT NONE
  
#include "mpif.h"
  
  TREAL,    DIMENSION(*),          INTENT(IN)  :: b
  TREAL,    DIMENSION(*),          INTENT(OUT) :: a
  TINTEGER, DIMENSION(ims_npro_i), INTENT(IN)  :: dsend, drecv ! displacements
  INTEGER,  DIMENSION(ims_npro_i), INTENT(IN)  :: tsend, trecv ! types
  
! -----------------------------------------------------------------------
  TINTEGER n, l
  INTEGER status(MPI_STATUS_SIZE,2*ims_npro_i)
  INTEGER mpireq(                2*ims_npro_i)
  INTEGER ip

! #######################################################################
! Same processor
! #######################################################################
  ip = ims_pro_i; n = ip + 1
  CALL MPI_ISEND(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(1), ims_err)
  CALL MPI_IRECV(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 2
  DO n = 1,ims_npro_i
     ip = n-1
     IF ( ip .NE. ims_pro_i ) THEN
        l = l + 1
        CALL MPI_ISEND(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        l = l + 1
        CALL MPI_IRECV(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro_i*2-2, mpireq(3), status(1,3), ims_err)

  CALL DNS_MPI_TAGUPDT

  RETURN
END SUBROUTINE DNS_MPI_TRPB_I

#include "types.h"
#include "dns_const_mpi.h" 

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created 
!# 2009/09/01 - J.P. Mellado
!#              Debugged
!# 2015/03/14 - J.P. Mellado
!#              Using local communicator
!# 2019/05/23 - C. Ansorge
!#              Ring-shape transposes 
!#              Control to switch between synchronous and sendrecv implementation
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
  USE DNS_MPI, ONLY : ims_plan_trps_k, ims_plan_trpr_k 
  USE DNS_MPI, ONLY : ims_trp_mode_k 

  IMPLICIT NONE
  
#include "mpif.h"
  
  TREAL,    DIMENSION(*),          INTENT(IN)  :: b
  TREAL,    DIMENSION(*),          INTENT(OUT) :: a
  TINTEGER, DIMENSION(ims_npro_k), INTENT(IN)  :: dsend, drecv
  INTEGER,  DIMENSION(ims_npro_k), INTENT(IN)  :: tsend, trecv
  
! -----------------------------------------------------------------------
  TINTEGER l,m,ns,nr,ips,ipr
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

  l = 0
  DO m=1,ims_npro_k 
     ns=ims_plan_trps_k(m)+1; ips=ns-1
     nr=ims_plan_trpr_k(m)+1; ipr=nr-1  
     IF ( ims_trp_mode_k .EQ. DNS_MPI_TRP_ASYNCHRONOUS ) THEN 
        l = l + 1
        CALL MPI_ISEND(b(drecv(nr)+1), 1, trecv(nr), ipr, ims_tag, ims_comm_z, mpireq(l), ims_err)
        l = l + 1
        CALL MPI_IRECV(a(dsend(ns)+1), 1, tsend(ns), ips, ims_tag, ims_comm_z, mpireq(l), ims_err) 
     ELSEIF ( ims_trp_mode_k .EQ. DNS_MPI_TRP_SENDRECV ) THEN 
        CALL MPI_SENDRECV(& 
             b(drecv(nr)+1), 1, trecv(nr), ipr, ims_tag,  & 
             a(dsend(ns)+1), 1, tsend(ns), ips, ims_tag, ims_comm_z, status(1,m), ims_err)   
     ELSE; CONTINUE   ! No transpose
     ENDIF
  ENDDO

  IF ( ims_trp_mode_k .EQ. DNS_MPI_TRP_ASYNCHRONOUS ) & 
       CALL MPI_WAITALL(ims_npro_k*2, mpireq(1:), status(1,1), ims_err)

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
  USE DNS_MPI, ONLY : ims_trp_mode_i
  USE DNS_MPI, ONLY : ims_plan_trps_i, ims_plan_trpr_i
  IMPLICIT NONE
  
#include "mpif.h"
  
  TREAL,    DIMENSION(*),          INTENT(IN)  :: b
  TREAL,    DIMENSION(*),          INTENT(OUT) :: a
  TINTEGER, DIMENSION(ims_npro_i), INTENT(IN)  :: dsend, drecv ! displacements
  INTEGER,  DIMENSION(ims_npro_i), INTENT(IN)  :: tsend, trecv ! types
  
! -----------------------------------------------------------------------
  TINTEGER l,m,ns,nr,ips,ipr
  INTEGER status(MPI_STATUS_SIZE,2*ims_npro_i)
  INTEGER mpireq(                2*ims_npro_i)
  INTEGER ip

  l = 0
  DO m = 1,ims_npro_i
     ns=ims_plan_trps_i(m)+1; ips=ns-1
     nr=ims_plan_trpr_i(m)+1; ipr=nr-1  
     IF ( ims_trp_mode_i .EQ. DNS_MPI_TRP_ASYNCHRONOUS ) THEN 
        l = l + 1
        CALL MPI_ISEND(b(drecv(nr)+1), 1, trecv(nr), ipr, ims_tag, ims_comm_x, mpireq(l), ims_err)
        l = l + 1
        CALL MPI_IRECV(a(dsend(ns)+1), 1, tsend(ns), ips, ims_tag, ims_comm_x, mpireq(l), ims_err)
     ELSEIF ( ims_trp_mode_i .EQ. DNS_MPI_TRP_SENDRECV ) THEN 
        CALL MPI_SENDRECV(& 
             b(drecv(nr)+1), 1, trecv(nr), ipr, ims_tag,&
             a(dsend(ns)+1), 1, tsend(ns), ips, ims_tag, ims_comm_x, status(1,m), ims_err)  
     ELSE; CONTINUE    ! No transpose
     ENDIF
  ENDDO

  IF ( ims_trp_mode_i .EQ. DNS_MPI_TRP_ASYNCHRONOUS ) & 
       CALL MPI_WAITALL(ims_npro_i*2, mpireq(1:), status(1,1), ims_err) 

  CALL DNS_MPI_TAGUPDT

  RETURN
END SUBROUTINE DNS_MPI_TRPB_I

! from dns_const.h
#define TREAL      REAL(8)    ! user-defined types
#define TINTEGER   INTEGER(4)

! from dns_const_mpi.h
#define DNS_MPI_K_PARTIAL   1 ! tags and sizes for MPI data
#define DNS_MPI_I_PARTIAL   1

#define DNS_MPI_K_MAXTYPES 10
#define DNS_MPI_I_MAXTYPES  6

MODULE DNS_MPI
  IMPLICIT NONE
  SAVE
  
  TINTEGER, PARAMETER :: imax = 64 ! number of grid points per task
  TINTEGER, PARAMETER :: jmax = 96
  TINTEGER, PARAMETER :: kmax = 64

  INTEGER, PARAMETER :: ims_npro_i = 2 ! number of tasks in Ox and Oz (no decomposition along Oy)
  INTEGER, PARAMETER :: ims_npro_k = 2 
  
  TINTEGER, PARAMETER :: nmax = 10     ! number of repetitions of operations

! Data below should not be changed
  INTEGER  :: ims_pro, ims_pro_i, ims_pro_j, ims_pro_k ! task positioning

  INTEGER  :: ims_comm_xz,     ims_comm_x,     ims_comm_z      ! communicators
  INTEGER  :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

  INTEGER  :: ims_err, ims_tag

  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_i
  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_i
  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_i, ims_dr_i
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_i, ims_tr_i

  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_k
  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_k
  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_k, ims_dr_k
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_k, ims_tr_k

END MODULE DNS_MPI

!########################################################################
! Main program to test forwards and backwards transposition
!########################################################################
PROGRAM VMPI

  USE DNS_MPI
  
  IMPLICIT NONE
  
#include "mpif.h"

  TREAL, DIMENSION(:,:), ALLOCATABLE :: a
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk3d
  
! -------------------------------------------------------------------
  TREAL residual                                      ! Control
  TINTEGER t_srt,t_end,t_dif, PROC_CYCLES, MAX_CYCLES ! Time
  TINTEGER n
  CHARACTER*64 str
  CHARACTER*256 line

  INTEGER ims_npro
  TREAL dummy
  TINTEGER idummy, id

! ###################################################################
  call MPI_INIT(ims_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ims_npro,ims_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro, ims_err)
  
  IF ( ims_npro_i*ims_npro_k .NE. ims_npro ) THEN ! check
     IF ( ims_pro .EQ. 0 ) THEN
        WRITE(*,'(a)') 'Inconsistency in total number of PEs'
     ENDIF     
     CALL MPI_FINALIZE(ims_err)
     STOP
  ENDIF
  
  CALL DNS_MPI_INITIALIZE

  ALLOCATE(a    (imax*jmax*kmax,18)) ! Number of 3d arrays commonly used in the code
  ALLOCATE(wrk3d(imax*jmax*kmax   ))

! ###################################################################
! ###################################################################
! Create random array
  CALL RANDOM_NUMBER(a(1:imax*jmax*kmax,1))

  DO n = 1,nmax
     
! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
     IF ( ims_npro_i .GT. 1 ) THEN
        id = DNS_MPI_I_PARTIAL
        
        CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)

        CALL DNS_MPI_TRPF_I(a(1,1), wrk3d, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
        CALL DNS_MPI_TRPB_I(wrk3d, a(1,2), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

        CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)
        
        idummy = t_end-t_srt
        CALL MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        WRITE(str,'(E13.5E3)') REAL(t_dif)/PROC_CYCLES
        
        dummy = MAXVAL(ABS(a(1:imax*jmax*kmax,1)-a(1:imax*jmax*kmax,2)))
        CALL MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        WRITE(line,'(E13.5E3)') residual

        line = 'Checking MPI transposition for Ox derivatives: Residual '&
             //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        IF ( ims_pro .EQ. 0 ) THEN
           WRITE(*,'(a)') TRIM(ADJUSTL(line))
        ENDIF
        
     ENDIF
     
! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
     IF ( ims_npro_k .GT. 1 ) THEN
        id = DNS_MPI_K_PARTIAL
        
        CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)
        
        CALL DNS_MPI_TRPF_K(a(1,1), wrk3d, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
        CALL DNS_MPI_TRPB_K(wrk3d, a(1,2), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

        CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)
        
        idummy = t_end-t_srt
        CALL MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err) 
        WRITE(str,'(E13.5E3)') REAL(t_dif)/PROC_CYCLES 
        
        dummy = MAXVAL(ABS(a(1:imax*jmax*kmax,1)-a(1:imax*jmax*kmax,2)))
        CALL MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        WRITE(line,'(E13.5E3)') residual

        line = 'Checking MPI transposition for Oz derivatives: Residual '&
             //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        IF ( ims_pro .EQ. 0 ) THEN
           WRITE(*,'(a)') TRIM(ADJUSTL(line))
        ENDIF
        
     ENDIF

  ENDDO
  
  CALL MPI_FINALIZE(ims_err)

END PROGRAM VMPI

! #######################################################################
! Rest of routines
! #######################################################################
SUBROUTINE DNS_MPI_INITIALIZE

  USE DNS_MPI

  IMPLICIT NONE
  
#include "integers.h"
#include "mpif.h"

! -----------------------------------------------------------------------
  TINTEGER id, ip, npage
  TINTEGER dims(2)
  LOGICAL period(2), remain_dims(2), reorder

! #######################################################################
  ALLOCATE(ims_map_i(ims_npro_i))
  ALLOCATE(ims_size_i(DNS_MPI_I_MAXTYPES))
  ALLOCATE(ims_ds_i(ims_npro_i,DNS_MPI_I_MAXTYPES))
  ALLOCATE(ims_dr_i(ims_npro_i,DNS_MPI_I_MAXTYPES))
  ALLOCATE(ims_ts_i(ims_npro_i,DNS_MPI_I_MAXTYPES))
  ALLOCATE(ims_tr_i(ims_npro_i,DNS_MPI_I_MAXTYPES))

  ALLOCATE(ims_map_k(ims_npro_k))
  ALLOCATE(ims_size_k(DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_ds_k(ims_npro_k,DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_dr_k(ims_npro_k,DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_ts_k(ims_npro_k,DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_tr_k(ims_npro_k,DNS_MPI_K_MAXTYPES))

! #######################################################################
  ims_pro_i = MOD(ims_pro,ims_npro_i) ! Starting at 0
  ims_pro_k =     ims_pro/ims_npro_i  ! Starting at 0
  
  ims_map_i(1) = ims_pro_k*ims_npro_i
  DO ip = 2,ims_npro_i
     ims_map_i(ip) = ims_map_i(ip-1) + 1
  ENDDO
  
  ims_map_k(1) = ims_pro_i
  DO ip = 2,ims_npro_k
     ims_map_k(ip) = ims_map_k(ip-1) + ims_npro_i
  ENDDO

! #######################################################################
! Communicators
! #######################################################################
! the first index in the grid corresponds to k, the second to i
  dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .true.; reorder = .false.
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

!  CALL MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
!  coord(1) is ims_pro_k, and coord(2) is ims_pro_i

  remain_dims(1) = .false.; remain_dims(2) = .true.
  CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

  remain_dims(1) = .true.;  remain_dims(2) = .false.
  CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

! #######################################################################
! Derived MPI types to deal with the strides when tranposing data
! #######################################################################
  IF ( ims_npro_i .GT. 1 ) THEN
!  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox derivatives.')
     id = DNS_MPI_I_PARTIAL
     npage = kmax*jmax
     CALL DNS_MPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
  
  IF ( ims_npro_k .GT. 1 ) THEN
!  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Oz derivatives.')
     id = DNS_MPI_K_PARTIAL
     npage = imax*jmax
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
  
  CALL DNS_MPI_TAGRESET

  RETURN
END SUBROUTINE DNS_MPI_INITIALIZE

! ###################################################################
! ###################################################################
SUBROUTINE DNS_MPI_TYPE_I(ims_npro, imax, npage, nd, md, n1, n2, &
     nsize, sdisp, rdisp, stype, rtype)

  USE DNS_MPI, ONLY : ims_pro
  
  IMPLICIT NONE

#include "mpif.h"

  INTEGER ims_npro
  TINTEGER npage, imax, nsize
  TINTEGER nd, md, n1, n2
  TINTEGER sdisp(*), rdisp(*)
  INTEGER stype(*), rtype(*)

! -----------------------------------------------------------------------
  TINTEGER i
  INTEGER ims_ss, ims_rs, ims_err
  INTEGER ims_tmp1, ims_tmp2, ims_tmp3
!  CHARACTER*64 str, line

! #######################################################################
  IF ( MOD(npage,ims_npro) .EQ. 0 ) THEN
     nsize = npage/ims_npro
  ELSE
     IF ( ims_pro .EQ. 0 ) THEN
        WRITE(*,'(a)') 'Ratio npage/ims_npro_i not an integer'
     ENDIF     
     CALL MPI_FINALIZE(ims_err)
     STOP
  ENDIF

! Calculate Displacements in Forward Send/Receive
  sdisp(1) = 0
  rdisp(1) = 0
  DO i = 2,ims_npro
     sdisp(i) = sdisp(i-1) + imax *nd *nsize
     rdisp(i) = rdisp(i-1) + imax *md
  ENDDO

! #######################################################################
  DO i = 1,ims_npro

     ims_tmp1 = nsize *n1 ! count
     ims_tmp2 = imax  *n2 ! block
     ims_tmp3 = ims_tmp2  ! stride = block because things are together
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype(i), ims_err)
     CALL MPI_TYPE_COMMIT(stype(i), ims_err)

     ims_tmp1 = nsize           *n1 ! count
     ims_tmp2 = imax            *n2 ! block
     ims_tmp3 = imax*ims_npro   *n2 ! stride is a multiple of imax_total=imax*ims_npro_i
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype(i), ims_err)
     CALL MPI_TYPE_COMMIT(rtype(i), ims_err)

     CALL MPI_TYPE_SIZE(stype(i), ims_ss, ims_err)
     CALL MPI_TYPE_SIZE(rtype(i), ims_rs, ims_err)

  ENDDO

  RETURN
END SUBROUTINE DNS_MPI_TYPE_I

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TYPE_K(ims_npro, nmax, npage, nd, md, n1, n2, &
     nsize, sdisp, rdisp, stype, rtype)

  USE DNS_MPI, ONLY : ims_pro

  IMPLICIT NONE

#include "mpif.h"

  INTEGER ims_npro
  TINTEGER npage, nmax, nsize
  TINTEGER nd, md, n1, n2
  TINTEGER sdisp(*), rdisp(*)
  INTEGER stype(*), rtype(*)

! -----------------------------------------------------------------------
  TINTEGER i
  INTEGER ims_ss, ims_rs, ims_err
  INTEGER ims_tmp1, ims_tmp2, ims_tmp3

! #######################################################################
  IF ( MOD(npage,ims_npro) .EQ. 0 ) THEN
     nsize = npage/ims_npro
  ELSE
     IF ( ims_pro .EQ. 0 ) THEN
        WRITE(*,'(a)') 'Ratio npage/ims_npro_k not an integer'
     ENDIF     
     CALL MPI_FINALIZE(ims_err)
     STOP
  ENDIF

! Calculate Displacements in Forward Send/Receive
  sdisp(1) = 0
  rdisp(1) = 0
  DO i = 2,ims_npro
     sdisp(i) = sdisp(i-1) + nsize *nd
     rdisp(i) = rdisp(i-1) + nsize *md *nmax
  ENDDO

! #######################################################################
  DO i = 1,ims_npro

     ims_tmp1 = nmax  *n1 ! count
     ims_tmp2 = nsize *n2 ! block
     ims_tmp3 = npage *n2 ! stride
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype(i), ims_err)
     CALL MPI_TYPE_COMMIT(stype(i), ims_err)

     ims_tmp1 = nmax  *n1 ! count
     ims_tmp2 = nsize *n2 ! block
     ims_tmp3 = ims_tmp2  ! stride = block to put things together
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype(i), ims_err)
     CALL MPI_TYPE_COMMIT(rtype(i), ims_err)

     CALL MPI_TYPE_SIZE(stype(i), ims_ss, ims_err)
     CALL MPI_TYPE_SIZE(rtype(i), ims_rs, ims_err)

  ENDDO

  RETURN
END SUBROUTINE DNS_MPI_TYPE_K

! ###################################################################
! ###################################################################
SUBROUTINE DNS_MPI_TRPF_K(a, b, dsend, drecv, tsend, trecv)
  
  USE DNS_MPI, ONLY : ims_npro_k, ims_pro_k
  USE DNS_MPI, ONLY : ims_comm_z
  USE DNS_MPI, ONLY : ims_tag, ims_err

  IMPLICIT NONE
  
#include "mpif.h"

  TREAL,    DIMENSION(*),          INTENT(IN)  :: a
  TREAL,    DIMENSION(*),          INTENT(OUT) :: b
  TINTEGER, DIMENSION(ims_npro_k), INTENT(IN)  :: dsend, drecv ! displacements
  INTEGER,  DIMENSION(ims_npro_k), INTENT(IN)  :: tsend, trecv ! types
  
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
  CALL MPI_ISEND(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(1), ims_err)  
  CALL MPI_IRECV(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 2
  DO n = 1,ims_npro_k
     ip = n - 1
     IF ( ip .NE. ims_pro_k ) THEN
        l = l + 1      
        CALL MPI_ISEND(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        l = l + 1
        CALL MPI_IRECV(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)        
     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro_k*2-2, mpireq(3), status(1,3), ims_err)

  CALL DNS_MPI_TAGUPDT

  RETURN
END SUBROUTINE DNS_MPI_TRPF_K

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TRPF_I(a, b, dsend, drecv, tsend, trecv)
  
  USE DNS_MPI, ONLY : ims_npro_i, ims_pro_i
  USE DNS_MPI, ONLY : ims_comm_x
  USE DNS_MPI, ONLY : ims_tag, ims_err

  IMPLICIT NONE
  
#include "mpif.h"
  
  TREAL,    DIMENSION(*),          INTENT(IN)  :: a
  TREAL,    DIMENSION(*),          INTENT(OUT) :: b
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
  CALL MPI_ISEND(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(1), ims_err)  
  CALL MPI_IRECV(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 2
  DO n = 1,ims_npro_i
     ip = n-1 
     IF ( ip .NE. ims_pro_i ) THEN
        l = l + 1
        CALL MPI_ISEND(a(dsend(n)+1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        l = l + 1 
        CALL MPI_IRECV(b(drecv(n)+1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro_i*2-2, mpireq(3), status(1,3), ims_err)

  CALL DNS_MPI_TAGUPDT

  RETURN
END SUBROUTINE DNS_MPI_TRPF_I

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)

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

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TAGUPDT
  
  USE DNS_MPI, ONLY : ims_tag

  IMPLICIT NONE
  
  ims_tag = ims_tag+1
  
  IF ( ims_tag .GT. 32000 ) THEN
     CALL DNS_MPI_TAGRESET
  ENDIF
  
  RETURN
END SUBROUTINE DNS_MPI_TAGUPDT

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TAGRESET
  
  USE DNS_MPI, ONLY : ims_tag

  IMPLICIT NONE
  
  ims_tag = 0
  
  RETURN
END SUBROUTINE DNS_MPI_TAGRESET
    

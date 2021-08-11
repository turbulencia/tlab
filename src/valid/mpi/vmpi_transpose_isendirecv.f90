!mpif90 -fpp  -nbs -save-temps -xHost -simd -vec-threshold50 -unroll-aggressive    -axcommon-avx512,SSE4.2  -qopt-prefetch -O3 vmpi_transpose.f90 

MODULE DNS_MPI
  IMPLICIT NONE
  SAVE
  
  INTEGER(4), PARAMETER :: imax = 84   ! number of grid points per task
  INTEGER(4), PARAMETER :: jmax = 480
  INTEGER(4), PARAMETER :: kmax = 56

  INTEGER, PARAMETER :: ims_npro_i = 64 ! number of tasks in Ox and Oz (no decomposition along Oy)
  INTEGER, PARAMETER :: ims_npro_k = 96
  
  INTEGER(4), PARAMETER :: nmax = 20000  ! number of repetitions of operations

! Data below should not be changed
  INTEGER  :: ims_pro, ims_pro_i, ims_pro_k        ! task positioning
  INTEGER  :: ims_comm_xz, ims_comm_x, ims_comm_z  ! communicators

! transposition along X
  INTEGER,  DIMENSION(:), ALLOCATABLE :: ims_map_i
  INTEGER(4), DIMENSION(:), ALLOCATABLE :: ims_size_i
  INTEGER(4), DIMENSION(:), ALLOCATABLE :: ims_ds_i, ims_dr_i

! transposition along Z
  INTEGER,  DIMENSION(:), ALLOCATABLE :: ims_map_k
  INTEGER(4), DIMENSION(:), ALLOCATABLE :: ims_size_k
  INTEGER(4), DIMENSION(:), ALLOCATABLE :: ims_ds_k, ims_dr_k

! Information
  INTEGER  :: ims_err, ims_tag
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: status
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: mpireq

END MODULE DNS_MPI

!########################################################################
! Main program to test forwards and backwards transposition
!########################################################################
PROGRAM VMPI

  USE TLAB_MPI_VARS
  
  IMPLICIT NONE
  
#include "mpif.h"

  REAL(8), DIMENSION(:,:), ALLOCATABLE :: a
  REAL(8), DIMENSION(:),   ALLOCATABLE :: b
  
! -------------------------------------------------------------------
  INTEGER(4) it, ip, l, n
  CHARACTER*64 str

  INTEGER ims_npro

  INTEGER(4) dims(2)
  LOGICAL period(2), remain_dims(2), reorder

! ###################################################################
! Allocate arrays
  ALLOCATE(a(imax*jmax*kmax,18)) ! Number of 3d arrays commonly used in the code
  ALLOCATE(b(imax*jmax*kmax   ))

! ###################################################################
! Initialize parallel part
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
  
  ALLOCATE(ims_map_i(ims_npro_i))
  ALLOCATE(ims_size_i(ims_npro_i))
  ALLOCATE(ims_ds_i(ims_npro_i))
  ALLOCATE(ims_dr_i(ims_npro_i))

  ALLOCATE(ims_map_k(ims_npro_k))
  ALLOCATE(ims_size_k(ims_npro_k))
  ALLOCATE(ims_ds_k(ims_npro_k))
  ALLOCATE(ims_dr_k(ims_npro_k))

  ALLOCATE(status(MPI_STATUS_SIZE,2*MAX(ims_npro_k,ims_npro_i)))
  ALLOCATE(mpireq(                2*MAX(ims_npro_k,ims_npro_i)))

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
! Displacements and sizes
! #######################################################################
  IF ( ims_npro_i .GT. 1 ) THEN
     ! Calculate size
     ims_size_i(:) = imax*jmax*kmax /ims_npro_i
     ! Calculate Displacements in Forward Send/Receive
     ims_ds_i(1) = 0
     ims_dr_i(1) = 0
     DO ip = 2,ims_npro_i
        ims_ds_i(ip) = ims_ds_i(ip-1) + ims_size_i(ip-1)
        ims_dr_i(ip) = ims_dr_i(ip-1) + ims_size_i(ip-1)
     ENDDO
  ENDIF
  
  IF ( ims_npro_k .GT. 1 ) THEN
     ! Calculate size
     ims_size_k(:) = imax*jmax*kmax /ims_npro_k
     ! Calculate Displacements in Forward Send/Receive
     ims_ds_k(1) = 0
     ims_dr_k(1) = 0
     DO ip = 2,ims_npro_k
        ims_ds_k(ip) = ims_ds_k(ip-1) + ims_size_k(ip-1)
        ims_dr_k(ip) = ims_dr_k(ip-1) + ims_size_k(ip-1)
     ENDDO
  ENDIF
  
! ###################################################################
! ###################################################################
! Create random array
  CALL RANDOM_NUMBER(a(1:imax*jmax*kmax,1))
  
  DO it = 1,nmax
     WRITE(str,'(I)') it
     
! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
     IF ( ims_npro_i .GT. 1 ) THEN
        ! CALL MPI_ALLTOALLV(a(1,1), ims_size_i, ims_ds_i, MPI_REAL8, &
        !                    b,      ims_size_i, ims_dr_i, MPI_REAL8, ims_comm_x, ims_err)

        ! Same processor
        ip = ims_pro_i; n = ip + 1
        CALL MPI_ISEND(a(ims_ds_i(n)+1,1), ims_size_i(n), MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(1), ims_err)  
        CALL MPI_IRECV(b(ims_dr_i(n)+1),   ims_size_i(n), MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(2), ims_err)
        CALL MPI_WAITALL(2, mpireq, status, ims_err)

        ! Different processors
        l = 2
        DO n = 1,ims_npro_i
           ip = n-1 
           IF ( ip .NE. ims_pro_i ) THEN
              l = l + 1
              CALL MPI_ISEND(a(ims_ds_i(n)+1,1), ims_size_i(n), MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
              l = l + 1 
              CALL MPI_IRECV(b(ims_dr_i(n)+1),   ims_size_i(n), MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
           ENDIF
        ENDDO
        CALL MPI_WAITALL(ims_npro_i*2-2, mpireq(3:), status(1,3), ims_err)

        ! Output
        IF ( ims_pro .EQ. 0 ) THEN
           WRITE(*,'(a)') 'It '//TRIM(ADJUSTL(str))//'. Tranpose along I.'
        ENDIF
        
     ENDIF
     
! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
     IF ( ims_npro_k .GT. 1 ) THEN
        ! CALL MPI_ALLTOALLV(a(1,1), ims_size_k, ims_ds_k, MPI_REAL8, &
        !                    b,      ims_size_k, ims_dr_k, MPI_REAL8, ims_comm_z, ims_err)

        ! Same processor
        ip = ims_pro_k; n = ip + 1
        CALL MPI_ISEND(a(ims_ds_k(n)+1,1), ims_size_k(n), MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(1), ims_err)  
        CALL MPI_IRECV(b(ims_dr_k(n)+1),   ims_size_k(n), MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(2), ims_err)
        CALL MPI_WAITALL(2, mpireq, status, ims_err)

        ! Different processors
        l = 2
        DO n = 1,ims_npro_k
           ip = n - 1
           IF ( ip .NE. ims_pro_k ) THEN
              l = l + 1      
              CALL MPI_ISEND(a(ims_ds_k(n)+1,1), ims_size_k(n), MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
              l = l + 1
              CALL MPI_IRECV(b(ims_dr_k(n)+1),   ims_size_k(n), MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)        
           ENDIF
        ENDDO
        CALL MPI_WAITALL(ims_npro_k*2-2, mpireq(3:), status(1,3), ims_err)

        ! Output
        IF ( ims_pro .EQ. 0 ) THEN
           WRITE(*,'(a)') 'It '//TRIM(ADJUSTL(str))//'. Tranpose along K.'
        ENDIF

     ENDIF

  ENDDO
  
  CALL MPI_FINALIZE(ims_err)

END PROGRAM VMPI

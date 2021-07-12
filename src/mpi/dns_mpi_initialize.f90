#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

SUBROUTINE DNS_MPI_INITIALIZE

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax,g
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_dimx
  USE DNS_GLOBAL, ONLY : imode_sim, ifourier
  USE DNS_GLOBAL, ONLY : imode_ibm
  USE DNS_IBM,    ONLY : xbars_geo
  USE DNS_CONSTANTS, ONLY : lfile
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
  ALLOCATE(ims_ts_i(1,DNS_MPI_I_MAXTYPES))   ! ims_npro_i, DMS_MPI_I_MAXTYPES
  ALLOCATE(ims_tr_i(1,DNS_MPI_I_MAXTYPES))   ! ims_npro_i, DNS_MPI_I_MAXTYPES
  ALLOCATE(ims_plan_trps_i(ims_npro_i))
  ALLOCATE(ims_plan_trpr_i(ims_npro_i)) 

  ALLOCATE(ims_map_k(ims_npro_k))
  ALLOCATE(ims_size_k(DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_ds_k(ims_npro_k,DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_dr_k(ims_npro_k,DNS_MPI_K_MAXTYPES))
  ALLOCATE(ims_ts_k(1,DNS_MPI_K_MAXTYPES))  ! ims_npro_k,DNS_MPI_K_MAXTYPES
  ALLOCATE(ims_tr_k(1,DNS_MPI_K_MAXTYPES))  ! ims_npro_k,DNS_MPI_K_MAXTYPES
  ALLOCATE(ims_plan_trps_k(ims_npro_k))
  ALLOCATE(ims_plan_trpr_k(ims_npro_k)) 

  ALLOCATE(ims_size_j(DNS_MPI_J_MAXTYPES)) ! IBM

  ALLOCATE(ims_size_p(ims_npro)) ! Particle information

! #######################################################################
  ims_pro_i = MOD(ims_pro,ims_npro_i) ! Starting at 0
  ims_pro_k =     ims_pro/ims_npro_i  ! Starting at 0
  
  ims_offset_i = ims_pro_i *imax
  ims_offset_j = 0
  ims_offset_k = ims_pro_k *kmax
  
  ims_map_i(1) = ims_pro_k*ims_npro_i
  DO ip = 2,ims_npro_i
     ims_map_i(ip) = ims_map_i(ip-1) + 1
  ENDDO
  
  ims_map_k(1) = ims_pro_i
  DO ip = 2,ims_npro_k
     ims_map_k(ip) = ims_map_k(ip-1) + ims_npro_i
  ENDDO

! #######################################################################
  CALL IO_WRITE_ASCII(lfile,'Initializing MPI communicators.')

! the first index in the grid corresponds to k, the second to i
  dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .true.; reorder = .false.
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

!  CALL MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
!  coord(1) is ims_pro_k, and coord(2) is ims_pro_i

  remain_dims(1) = .false.; remain_dims(2) = .true.
  CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

  remain_dims(1) = .true.;  remain_dims(2) = .false.
  CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

  ! ip = ims_pro
  ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_x, ims_err)
  ! print*, 'P:', ims_pro, 'Sum along X', id
  ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_z, ims_err)
  ! print*, 'P:', ims_pro, 'Sum along Z', id

! #######################################################################
! Main
! #######################################################################
  IF ( ims_npro_i .GT. 1 ) THEN
  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox derivatives.')
  id = DNS_MPI_I_PARTIAL
  npage = kmax*jmax
  CALL DNS_MPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
       ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF

  IF ( ims_npro_k .GT. 1 ) THEN
  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Oz derivatives.')
  id = DNS_MPI_K_PARTIAL
  npage = imax*jmax
  CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
       ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

! -----------------------------------------------------------------------
  ! Immersed Boundary Method (IBM)
  IF (imode_ibm .EQ. 1) THEN
    id = DNS_MPI_J_PARTIAL
    npage = imax*kmax
    ims_size_j(id) = npage

    ! ------------------ !
    IF (ims_npro_i .GT. 1) THEN
    CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox IBM nobi.')
    id = DNS_MPI_I_IBM_NOB
    npage = g(2)%size * g(3)%size / ims_npro
    CALL DNS_MPI_TYPE_I(ims_npro_i, i1, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    ENDIF

    IF ( ims_npro_k .GT. 1 ) THEN
    CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Oz IBM nobk.')
    id = DNS_MPI_K_IBM_NOB
    npage = g(1)%size * g(2)%size / ims_npro
    CALL DNS_MPI_TYPE_K(ims_npro_k, i1, npage, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    ENDIF

    ! ------------------ !
    IF (ims_npro_i .GT. 1) THEN
    CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox IBM nobi_b and nobi_e.')
    id = DNS_MPI_I_IBM_NOB_BE
    npage = g(2)%size * g(3)%size / ims_npro
    CALL DNS_MPI_TYPE_I(ims_npro_i, xbars_geo(1), npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    ENDIF

    IF ( ims_npro_k .GT. 1 ) THEN
    CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Oz IBM nobk_b and nobk_e.')
    id = DNS_MPI_K_IBM_NOB_BE
    npage = g(1)%size * g(2)%size / ims_npro
    CALL DNS_MPI_TYPE_K(ims_npro_k, xbars_geo(1), npage, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    ENDIF

  ENDIF

! -----------------------------------------------------------------------
  IF ( ims_npro_i .GT. 1 .AND. ifourier .EQ. 1 ) THEN
  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox FFTW in Poisson solver.')
  id = DNS_MPI_I_POISSON1 
  npage = isize_txc_dimx ! isize_txc_field/imax
  CALL DNS_MPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
       ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Ox FFTW in Poisson solver.')
  id = DNS_MPI_I_POISSON2 ! isize_txc_field/(imax+2)
  npage = isize_txc_dimx
  CALL DNS_MPI_TYPE_I(ims_npro_i, imax+2, npage, i1, i1, i1, i1, &
       ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

  ENDIF

  IF ( ims_npro_k .GT. 1 .AND. ifourier .EQ. 1 ) THEN
  CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for Oz FFTW in Poisson solver.')
  id = DNS_MPI_K_POISSON
  npage = isize_txc_dimz ! isize_txc_field/kmax
  CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
       ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

! ######################################################################
! Work plans for circular transposes 
! ######################################################################
  DO ip=0,ims_npro_i-1
     ims_plan_trps_i(ip+1) = ip
     ims_plan_trpr_i(ip+1) = MOD(ims_npro_i-ip,ims_npro_i)
  ENDDO
  ims_plan_trps_i = CSHIFT(ims_plan_trps_i,  ims_pro_i  ) 
  ims_plan_trpr_i = CSHIFT(ims_plan_trpr_i,-(ims_pro_i) ) 


  DO ip=0,ims_npro_k-1 
     ims_plan_trps_k(ip+1) = ip 
     ims_plan_trpr_k(ip+1) = MOD(ims_npro_k-ip,ims_npro_k) 
  ENDDO
  ims_plan_trps_k = CSHIFT(ims_plan_trps_k,  ims_pro_k  ) 
  ims_plan_trpr_k = CSHIFT(ims_plan_trpr_k,-(ims_pro_k) ) 

  ! DO ip=0,ims_npro_i-1
  !    IF ( ims_pro .EQ. ip ) THEN 
  !       WRITE(*,*) ims_pro, ims_pro_i, 'SEND:', ims_plan_trps_i 
  !       WRITE(*,*) ims_pro, ims_pro_i, 'RECV:', ims_plan_trpr_i 
  !    ENDIF 
  !    CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  ! ENDDO

! #######################################################################
! Auxiliar depending on simmode
! #######################################################################
  ! IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
  ! CALL IO_WRITE_ASCII(lfile,'Initializing MPI types for spectra/correlations.')
  ! id = DNS_MPI_K_SHEAR
  ! CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, imax, i1, jmax, jmax, i1, &
  !      ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ! ENDIF

  CALL DNS_MPI_TAGRESET

  RETURN
END SUBROUTINE DNS_MPI_INITIALIZE

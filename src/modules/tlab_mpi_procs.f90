#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

MODULE TLAB_MPI_PROCS
  USE DNS_CONSTANTS, ONLY : lfile, efile
  USE TLAB_PROCS, ONLY : TLAB_WRITE_ASCII, TLAB_STOP
  USE TLAB_MPI_VARS
  IMPLICIT NONE
  SAVE
  PRIVATE

#include "mpif.h"

  INTEGER  :: ims_tag

  PUBLIC :: DNS_MPI_INITIALIZE
  PUBLIC :: DNS_MPI_TRPF_K
  PUBLIC :: DNS_MPI_TRPF_I
  PUBLIC :: DNS_MPI_TRPB_K
  PUBLIC :: DNS_MPI_TRPB_I
  PUBLIC :: DNS_MPI_TYPE_K
  PUBLIC :: DNS_MPI_TYPE_I
  PUBLIC :: DNS_MPI_WRITE_PE0_SINGLE

CONTAINS

  ! ######################################################################
  ! ######################################################################
  SUBROUTINE DNS_MPI_INITIALIZE
    USE DNS_GLOBAL, ONLY : imax,jmax,kmax
    USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_dimx
    USE DNS_GLOBAL, ONLY : imode_sim, ifourier
    IMPLICIT NONE

#include "integers.h"

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

    ALLOCATE(ims_size_p(ims_npro))            ! Particle information

#ifdef HLRS_HAWK
    ! On hawk, we tested that 192 yields optimum performace;
    ! Blocking will thus only take effect in very large cases
    ims_sizBlock_k=192
    ims_sizBlock_i=384
#else
    ! We assume that this will help to release some of the very heavy
    ! network load in transpositions on most systems
    ims_sizBlock_k=64
    ims_sizBlock_i=128
    ! ims_sizBlock_k=1e5   -- would essentially switch off the blocking
#endif

  ALLOCATE(ims_status (MPI_STATUS_SIZE,2*MAX(ims_sizBlock_i,ims_sizBlock_k)))
  ALLOCATE(ims_request(                2*MAX(ims_sizBlock_i,ims_sizBlock_k)))


    ! #######################################################################
    ims_pro_i = MOD(ims_pro,ims_npro_i) ! Starting at 0
    ims_pro_k =     ims_pro/ims_npro_i  ! Starting at 0

    ims_offset_i = ims_pro_i *imax
    ims_offset_j = 0
    ims_offset_k = ims_pro_k *kmax

    ims_map_i(1) = ims_pro_k*ims_npro_i
    DO ip = 2,ims_npro_i
      ims_map_i(ip) = ims_map_i(ip-1) + 1
    END DO

    ims_map_k(1) = ims_pro_i
    DO ip = 2,ims_npro_k
      ims_map_k(ip) = ims_map_k(ip-1) + ims_npro_i
    END DO

    ! #######################################################################
    CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI communicators.')

    ! the first index in the grid corresponds to k, the second to i
    dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .TRUE.; reorder = .FALSE.
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

    !  CALL MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
    !  coord(1) is ims_pro_k, and coord(2) is ims_pro_i

    remain_dims(1) = .FALSE.; remain_dims(2) = .TRUE.
    CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

    remain_dims(1) = .TRUE.;  remain_dims(2) = .FALSE.
    CALL MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

    ! ip = ims_pro
    ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_x, ims_err)
    ! print*, 'P:', ims_pro, 'Sum along X', id
    ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_z, ims_err)
    ! print*, 'P:', ims_pro, 'Sum along Z', id

    ! #######################################################################
    ! Main
    ! #######################################################################
    IF ( ims_npro_i > 1 ) THEN
      CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for Ox derivatives.')
      id = DNS_MPI_I_PARTIAL
      npage = kmax*jmax
      CALL DNS_MPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    END IF

    IF ( ims_npro_k > 1 ) THEN
      CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for Oz derivatives.')
      id = DNS_MPI_K_PARTIAL
      npage = imax*jmax
      CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    END IF

    ! -----------------------------------------------------------------------
    IF ( ims_npro_i > 1 .AND. ifourier == 1 ) THEN
      CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for Ox FFTW in Poisson solver.')
      id = DNS_MPI_I_POISSON1
      npage = isize_txc_dimx ! isize_txc_field/imax
      CALL DNS_MPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

      CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for Ox FFTW in Poisson solver.')
      id = DNS_MPI_I_POISSON2 ! isize_txc_field/(imax+2)
      npage = isize_txc_dimx
      CALL DNS_MPI_TYPE_I(ims_npro_i, imax+2, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

    END IF

    IF ( ims_npro_k > 1 .AND. ifourier == 1 ) THEN
      CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for Oz FFTW in Poisson solver.')
      id = DNS_MPI_K_POISSON
      npage = isize_txc_dimz ! isize_txc_field/kmax
      CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    END IF

    ! ######################################################################
    ! Work plans for circular transposes
    ! ######################################################################
    DO ip=0,ims_npro_i-1
      ims_plan_trps_i(ip+1) = ip
      ims_plan_trpr_i(ip+1) = MOD(ims_npro_i-ip,ims_npro_i)
    END DO
    ims_plan_trps_i = CSHIFT(ims_plan_trps_i,  ims_pro_i  )
    ims_plan_trpr_i = CSHIFT(ims_plan_trpr_i,-(ims_pro_i) )


    DO ip=0,ims_npro_k-1
      ims_plan_trps_k(ip+1) = ip
      ims_plan_trpr_k(ip+1) = MOD(ims_npro_k-ip,ims_npro_k)
    END DO
    ims_plan_trps_k = CSHIFT(ims_plan_trps_k,  ims_pro_k  )
    ims_plan_trpr_k = CSHIFT(ims_plan_trpr_k,-(ims_pro_k) )

    ! DO ip=0,ims_npro_i-1
    !    IF ( ims_pro == ip ) THEN
    !       WRITE(*,*) ims_pro, ims_pro_i, 'SEND:', ims_plan_trps_i
    !       WRITE(*,*) ims_pro, ims_pro_i, 'RECV:', ims_plan_trpr_i
    !    END IF
    !    CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
    ! END DO

    ! #######################################################################
    ! Auxiliar depending on simmode
    ! #######################################################################
    ! IF ( imode_sim == DNS_MODE_TEMPORAL ) THEN
    ! CALL TLAB_WRITE_ASCII(lfile,'Initializing MPI types for spectra/correlations.')
    ! id = DNS_MPI_K_SHEAR
    ! CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, imax, i1, jmax, jmax, i1, &
    !      ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    ! END IF

    CALL DNS_MPI_TAGRESET

    RETURN
  END SUBROUTINE DNS_MPI_INITIALIZE

  ! ######################################################################
  ! ######################################################################
  SUBROUTINE DNS_MPI_TAGUPDT
    IMPLICIT NONE

    ims_tag = ims_tag+1
    IF ( ims_tag > 32000 ) THEN
      CALL DNS_MPI_TAGRESET
    END IF

    RETURN
  END SUBROUTINE DNS_MPI_TAGUPDT

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TAGRESET
    IMPLICIT NONE

    ims_tag = 0

    RETURN
  END SUBROUTINE DNS_MPI_TAGRESET

  ! ######################################################################
  ! Pointers and types for transposition across ims_npro processors
  ! ######################################################################
  SUBROUTINE DNS_MPI_TYPE_I(npro_i, nmax, npage, nd, md, n1, n2, &
      nsize, sdisp, rdisp, stype, rtype)
    IMPLICIT NONE

    INTEGER npro_i
    TINTEGER npage, nmax, nsize
    TINTEGER nd, md, n1, n2
    TINTEGER sdisp(*), rdisp(*)
    INTEGER, INTENT(OUT) ::  stype, rtype

    ! -----------------------------------------------------------------------
    TINTEGER i
    INTEGER ims_ss, ims_rs
    INTEGER ims_tmp1, ims_tmp2, ims_tmp3
    CHARACTER*64 str, line

    ! #######################################################################
    IF ( MOD(npage,npro_i) == 0 ) THEN
      nsize = npage/npro_i
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_TYPE_I. Ratio npage/npro_i not an integer.')
      CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
    END IF

    ! Calculate Displacements in Forward Send/Receive
    sdisp(1) = 0
    rdisp(1) = 0
    DO i = 2,npro_i
      sdisp(i) = sdisp(i-1) + nmax *nd *nsize
      rdisp(i) = rdisp(i-1) + nmax *md
    END DO

    ! #######################################################################
    ims_tmp1 = nsize *n1 ! count
    ims_tmp2 = nmax  *n2 ! block
    ims_tmp3 = ims_tmp2  ! stride = block because things are together
    CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
    CALL MPI_TYPE_COMMIT(stype, ims_err)

    ims_tmp1 = nsize           *n1 ! count
    ims_tmp2 = nmax            *n2 ! block
    ims_tmp3 = nmax*npro_i *n2 ! stride is a multiple of nmax_total=nmax*npro_i
    CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
    CALL MPI_TYPE_COMMIT(rtype, ims_err)

    CALL MPI_TYPE_SIZE(stype, ims_ss, ims_err)
    CALL MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

    IF ( ims_ss /= ims_rs ) THEN
      WRITE(str, *) ims_ss; WRITE(line,*) ims_rs
      line='Send size '//TRIM(ADJUSTL(str))//'differs from recv size '//TRIM(ADJUSTL(line))
      WRITE(str, *) 1  ! i
      line=TRIM(ADJUSTL(line))//' in message '//TRIM(ADJUSTL(str))
      CALL TLAB_WRITE_ASCII(efile, line)
      CALL TLAB_STOP(DNS_ERROR_MPITYPECHECK)
    END IF

    RETURN
  END SUBROUTINE DNS_MPI_TYPE_I

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TYPE_K(npro_k, nmax, npage, nd, md, n1, n2, &
      nsize, sdisp, rdisp, stype, rtype)
    IMPLICIT NONE

    INTEGER npro_k
    TINTEGER npage, nmax, nsize
    TINTEGER nd, md, n1, n2
    TINTEGER sdisp(*), rdisp(*)
    INTEGER, INTENT(OUT) ::  stype, rtype

    ! -----------------------------------------------------------------------
    TINTEGER i
    INTEGER ims_ss, ims_rs
    INTEGER ims_tmp1, ims_tmp2, ims_tmp3

    ! #######################################################################
    IF ( MOD(npage,npro_k) == 0 ) THEN
      nsize = npage/npro_k
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_TYPE_K. Ratio npage/npro_k not an integer.')
      CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
    END IF

    ! Calculate Displacements in Forward Send/Receive
    sdisp(1) = 0
    rdisp(1) = 0
    DO i = 2,npro_k
      sdisp(i) = sdisp(i-1) + nsize *nd
      rdisp(i) = rdisp(i-1) + nsize *md *nmax
    END DO

    ! #######################################################################
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

    IF ( ims_ss /= ims_rs ) THEN
      PRINT *, 'Message   : ', 1, ' size is wrong' ! i
      PRINT *, 'Send size : ', ims_ss
      PRINT *, 'Recv size : ', ims_rs
      CALL TLAB_STOP(DNS_ERROR_MPITYPECHECK)
    END IF

    RETURN
  END SUBROUTINE DNS_MPI_TYPE_K

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TRPF_K(a, b, dsend, drecv, tsend, trecv)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
    USE DNS_GLOBAL, ONLY : itime
    IMPLICIT NONE

    INTERFACE
      FUNCTION DNS_USLEEP (useconds)  BIND ( C, name="usleep" )
        IMPORT
        INTEGER (c_int) :: nb3dfft_nbc_usleep
        INTEGER (c_int), INTENT (in), VALUE :: useconds
      END FUNCTION DNS_USLEEP
    END INTERFACE

    TREAL,    DIMENSION(*),          INTENT(IN)  :: a
    TREAL,    DIMENSION(*),          INTENT(OUT) :: b
    TINTEGER, DIMENSION(ims_npro_k), INTENT(IN)  :: dsend, drecv ! displacements
    ! INTEGER,  DIMENSION(ims_npro_k), INTENT(IN)  :: tsend, trecv ! types
    INTEGER, INTENT(IN) :: tsend, trecv

    ! -----------------------------------------------------------------------
    TINTEGER j,l,m,n,ns,nr,ips,ipr,idummy
    !INTEGER status(MPI_STATUS_SIZE,2*ims_npro_k)
    !INTEGER mpireq(                2*ims_npro_k)
    INTEGER ip
    TREAL rdum
    CHARACTER*256 line
#ifdef PROFILE_ON
    TREAL time_loc_1, time_loc_2
#endif

    ! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif
    DO j=1,ims_npro_k,ims_sizBlock_k
       l=0
       DO m=j,MIN(j+ims_sizBlock_k-1,ims_npro_k)
          ns=ims_plan_trps_k(m)+1; ips=ns-1
          nr=ims_plan_trpr_k(m)+1; ipr=nr-1
          IF ( ims_trp_mode_k == DNS_MPI_TRP_ASYNCHRONOUS) THEN
             l = l + 1
             CALL MPI_ISEND(a(dsend(ns)+1), 1, tsend, ips, ims_tag, ims_comm_z, ims_request(l), ims_err)
             l = l + 1
             CALL MPI_IRECV(b(drecv(nr)+1), 1, trecv, ipr, ims_tag, ims_comm_z, ims_request(l), ims_err)
          ELSEIF ( ims_trp_mode_k == DNS_MPI_TRP_SENDRECV) THEN
             CALL MPI_SENDRECV(&
                  a(dsend(ns)+1), 1, tsend, ips, ims_tag, &
                  b(drecv(nr)+1), 1, trecv, ipr, ims_tag, ims_comm_z, ims_status(1,1), ims_err)
          ELSE;  CONTINUE     ! No transpose
          END IF
       END DO

       IF ( ims_trp_mode_k == DNS_MPI_TRP_ASYNCHRONOUS ) &
            CALL MPI_WAITALL(l, ims_request(1), ims_status(1,1), ims_err)

       CALL DNS_MPI_TAGUPDT
    ENDDO
#ifdef PROFILE_ON
    time_loc_2 = MPI_WTIME()
    ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

    RETURN
  END SUBROUTINE DNS_MPI_TRPF_K

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TRPF_I(a, b, dsend, drecv, tsend, trecv)
    IMPLICIT NONE

    TREAL,    DIMENSION(*),          INTENT(IN)  :: a
    TREAL,    DIMENSION(*),          INTENT(OUT) :: b
    TINTEGER, DIMENSION(ims_npro_i), INTENT(IN)  :: dsend, drecv ! displacements
    !INTEGER,  DIMENSION(ims_npro_i), INTENT(IN)  :: tsend, trecv ! types
    INTEGER, INTENT(IN) :: tsend, trecv

    ! -----------------------------------------------------------------------
    TINTEGER j,l,m,ns,nr,ips,ipr
    !INTEGER status(MPI_STATUS_SIZE,2*ims_npro_i)
    !INTEGER mpireq(                2*ims_npro_i)
    INTEGER ip

    DO j=1,ims_npro_i,ims_sizBlock_i
       l = 0
       DO m=j,MIN(j+ims_sizBlock_i-1,ims_npro_i)
          ns=ims_plan_trps_i(m)+1;   ips=ns-1
          nr=ims_plan_trpr_i(m)+1;   ipr=nr-1
          IF ( ims_trp_mode_i == DNS_MPI_TRP_ASYNCHRONOUS ) THEN
             l = l + 1
             CALL MPI_ISEND(a(dsend(ns)+1), 1, tsend, ips, ims_tag, ims_comm_x, ims_request(l), ims_err)
             l = l + 1
             CALL MPI_IRECV(b(drecv(nr)+1), 1, trecv, ipr, ims_tag, ims_comm_x, ims_request(l), ims_err)
          ELSEIF ( ims_trp_mode_i == DNS_MPI_TRP_SENDRECV ) THEN
             CALL MPI_SENDRECV(&
                  a(dsend(ns)+1),1,tsend,ips, ims_tag, &
                  b(drecv(nr)+1),1,trecv,ipr, ims_tag,ims_comm_x,ims_status(1,1),ims_err)
          ELSE; CONTINUE ! No transpose
          END IF
       END DO

       IF ( ims_trp_mode_i == DNS_MPI_TRP_ASYNCHRONOUS ) &
            CALL MPI_WAITALL(l, ims_request(1), ims_status(1,1), ims_err)

       CALL DNS_MPI_TAGUPDT
    ENDDO

    RETURN
  END SUBROUTINE DNS_MPI_TRPF_I

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)
    IMPLICIT NONE

    TREAL,    DIMENSION(*),          INTENT(IN)  :: b
    TREAL,    DIMENSION(*),          INTENT(OUT) :: a
    TINTEGER, DIMENSION(ims_npro_k), INTENT(IN)  :: dsend, drecv
    INTEGER,  DIMENSION(ims_npro_k), INTENT(IN)  :: tsend, trecv

    ! -----------------------------------------------------------------------
    TINTEGER j,l,m,ns,nr,ips,ipr
    !INTEGER status(MPI_STATUS_SIZE,2*ims_npro_k)
    !INTEGER mpireq(                2*ims_npro_k)
    INTEGER ip

#ifdef PROFILE_ON
    TREAL time_loc_1, time_loc_2
#endif

    ! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif
    DO j=1,ims_npro_k,ims_sizBlock_k
       l = 0
       DO m=j,MIN(j+ims_sizBlock_k-1,ims_npro_k)
          ns=ims_plan_trps_k(m)+1; ips=ns-1
          nr=ims_plan_trpr_k(m)+1; ipr=nr-1
          IF ( ims_trp_mode_k == DNS_MPI_TRP_ASYNCHRONOUS ) THEN
             l = l + 1
             CALL MPI_ISEND(b(drecv(nr)+1), 1, trecv(1), ipr, ims_tag, ims_comm_z, ims_request(l), ims_err)
             l = l + 1
             CALL MPI_IRECV(a(dsend(ns)+1), 1, tsend(1), ips, ims_tag, ims_comm_z, ims_request(l), ims_err)
          ELSEIF ( ims_trp_mode_k == DNS_MPI_TRP_SENDRECV ) THEN
             CALL MPI_SENDRECV(&
                  b(drecv(nr)+1), 1, trecv(1), ipr, ims_tag,  &
                  a(dsend(ns)+1), 1, tsend(1), ips, ims_tag, ims_comm_z, ims_status(1,1), ims_err)
          ELSE; CONTINUE   ! No transpose
          END IF
       END DO

       IF ( ims_trp_mode_k == DNS_MPI_TRP_ASYNCHRONOUS ) &
            CALL MPI_WAITALL(l, ims_request(1), ims_status(1,1), ims_err)

       CALL DNS_MPI_TAGUPDT
    ENDDO

#ifdef PROFILE_ON
    time_loc_2 = MPI_WTIME()
    ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

    RETURN
  END SUBROUTINE DNS_MPI_TRPB_K

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_TRPB_I(b, a, dsend, drecv, tsend, trecv)
    IMPLICIT NONE

    TREAL,    DIMENSION(*),          INTENT(IN)  :: b
    TREAL,    DIMENSION(*),          INTENT(OUT) :: a
    TINTEGER, DIMENSION(ims_npro_i), INTENT(IN)  :: dsend, drecv ! displacements
    INTEGER,  DIMENSION(ims_npro_i), INTENT(IN)  :: tsend, trecv ! types

    ! -----------------------------------------------------------------------
    TINTEGER j,l,m,ns,nr,ips,ipr
    !INTEGER status(MPI_STATUS_SIZE,2*ims_npro_i)
    !INTEGER mpireq(                2*ims_npro_i)
    INTEGER ip

    DO j=1,ims_npro_i,ims_sizBlock_i
       l = 0
       DO m = j,MIN(j+ims_sizBlock_i-1,ims_npro_i)
          ns=ims_plan_trps_i(m)+1; ips=ns-1
          nr=ims_plan_trpr_i(m)+1; ipr=nr-1
          IF ( ims_trp_mode_i == DNS_MPI_TRP_ASYNCHRONOUS ) THEN
             l = l + 1
             CALL MPI_ISEND(b(drecv(nr)+1), 1, trecv(1), ipr, ims_tag, ims_comm_x, ims_request(l), ims_err)
             l = l + 1
             CALL MPI_IRECV(a(dsend(ns)+1), 1, tsend(1), ips, ims_tag, ims_comm_x, ims_request(l), ims_err)
          ELSEIF ( ims_trp_mode_i == DNS_MPI_TRP_SENDRECV ) THEN
             CALL MPI_SENDRECV(&
                  b(drecv(nr)+1), 1, trecv(1), ipr, ims_tag,&
                  a(dsend(ns)+1), 1, tsend(1), ips, ims_tag, ims_comm_x, ims_status(1,1), ims_err)
          ELSE; CONTINUE    ! No transpose
          END IF
       END DO

       IF ( ims_trp_mode_i == DNS_MPI_TRP_ASYNCHRONOUS ) &
            CALL MPI_WAITALL(l, ims_request(1), ims_status(1,1), ims_err)

       CALL DNS_MPI_TAGUPDT
    ENDDO

    RETURN
  END SUBROUTINE DNS_MPI_TRPB_I

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_WRITE_PE0_SINGLE(iunit, nx,ny,nz, subdomain, u, tmp1, tmp2)
    IMPLICIT NONE

    TINTEGER iunit, nx,ny,nz, subdomain(6)
    TREAL, DIMENSION(nx*ny*nz),        TARGET :: u, tmp1
    TREAL, DIMENSION(nx*ims_npro_i,*), TARGET :: tmp2

    ! -------------------------------------------------------------------
    TINTEGER nx_total, ny_total, nz_total
    TINTEGER nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
    TINTEGER nyz

    TINTEGER ip_i, ip_k, joffset_loc, koffset_loc, id
    TINTEGER i, jk, j_loc, k_loc
    INTEGER mpio_size, mpio_ip
    INTEGER status(MPI_STATUS_SIZE)

    TREAL, DIMENSION(:), POINTER :: p_org

    ! ###################################################################
    nx_total = nx*ims_npro_i
    ny_total = ny
    nz_total = nz*ims_npro_k

    nx_min = subdomain(1); nx_max = subdomain(2)
    ny_min = subdomain(3); ny_max = subdomain(4)
    nz_min = subdomain(5); nz_max = subdomain(6)

    koffset_loc = 0
    joffset_loc = 0

    id = DNS_MPI_I_PARTIAL

    ! -------------------------------------------------------------------
    ! Transposing along Ox
    ! -------------------------------------------------------------------
    IF ( ims_npro_i > 1 ) THEN
      CALL DNS_MPI_TRPF_I(u, tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
      p_org => tmp1
      nyz = ims_size_i(id)
    ELSE
      p_org => u
      nyz = ny*nz
    END IF
    mpio_size = nyz*nx_total

    CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    ! -------------------------------------------------------------------
    ! Passing all data through PE#0
    ! -------------------------------------------------------------------
    IF ( ims_pro == 0 ) THEN

      DO ip_k = 1,ims_npro_k
        koffset_loc = nz *(ip_k-1)

        DO ip_i = 1,ims_npro_i
          joffset_loc = nyz *(ip_i-1) ! Remember that data is Ox-transposed

          mpio_ip = ims_npro_i *(ip_k-1) + ip_i-1
          IF ( mpio_ip == 0 ) THEN
            tmp2(1:mpio_size,1) = p_org(1:mpio_size)
          ELSE
            CALL MPI_RECV(tmp2, mpio_size, MPI_REAL8, mpio_ip, ims_tag, MPI_COMM_WORLD, status, ims_err)
          END IF

          DO jk = 1,nyz
            j_loc = MOD( (jk-1+joffset_loc) , ny_total ) + 1
            k_loc =    ( (jk-1+joffset_loc) / ny_total ) + 1 + koffset_loc

            IF ( ( j_loc >= ny_min ) .AND. ( j_loc <= ny_max )  .AND. &
                ( k_loc >= nz_min ) .AND. ( k_loc <= nz_max ) ) THEN
              WRITE(iunit) (SNGL(tmp2(i,jk)),i=nx_min,nx_max)
            END IF

          END DO

        END DO
      END DO

    ELSE
      CALL MPI_SEND(p_org, mpio_size, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
    END IF

    RETURN
  END SUBROUTINE DNS_MPI_WRITE_PE0_SINGLE

  !########################################################################
  !# Moving plane information between adjacent PEs, circulant version.
  !# npl is smaller than 2*kmax
  !# The number of plane to move is given by npl
  !########################################################################
  SUBROUTINE DNS_MPI_COPYPLN_1(ijmax, kmax, npl, a, bl, br)
    IMPLICIT NONE

    TINTEGER ijmax, kmax, npl
    TREAL a(ijmax,*)
    TREAL bl(ijmax,*)
    TREAL br(ijmax,*)

    ! -----------------------------------------------------------------------
    INTEGER status(MPI_STATUS_SIZE,4)
    INTEGER mpireq(4)
    INTEGER ims_pro_l, ims_pro_r
    INTEGER icount

    ! #######################################################################
    IF ( ims_npro > 1 ) THEN

      ! Careful in case only 2 PEs
      IF ( ims_npro == 2 ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_COPYPLN_1. Undeveloped for 2 PEs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
      END IF

      ! left and right PEs
      ims_pro_l = MOD(ims_pro-1+ims_npro,ims_npro)
      ims_pro_r = MOD(ims_pro+1+ims_npro,ims_npro)

      icount   = ijmax*npl

      CALL MPI_IRECV(bl(1,1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
      CALL MPI_IRECV(br(1,1), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

      CALL MPI_ISEND(a(1,1),          icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
      CALL MPI_ISEND(a(1,kmax+1-npl), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

      CALL MPI_WAITALL(4, mpireq, status, ims_err)

      CALL DNS_MPI_TAGUPDT

    END IF

    RETURN
  END SUBROUTINE DNS_MPI_COPYPLN_1

  !########################################################################
  !########################################################################
  SUBROUTINE DNS_MPI_COPYPLN_2(ijmax, kmax, npl, a, bl, br)
    IMPLICIT NONE

    TINTEGER ijmax, kmax, npl
    TREAL a(ijmax,kmax)
    TREAL bl(ijmax,npl)
    TREAL br(ijmax,npl)

    ! -----------------------------------------------------------------------
    TINTEGER npl_loc
    INTEGER status(MPI_STATUS_SIZE,8)
    INTEGER mpireq(8)
    INTEGER ims_pro_l, ims_pro_r
    INTEGER icount

    ! #######################################################################
    IF ( ims_npro > 1 ) THEN

      ! Careful in case only 2 PEs
      IF ( ims_npro == 2 ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_COPYPLN_2. Undeveloped for 2 PEs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
      END IF

      ! -----------------------------------------------------------------------
      ! left and right PEs. Same as in routine DNS_MPI_COPYPLN_1
      ! -----------------------------------------------------------------------
      npl_loc = kmax

      ims_pro_l = MOD(ims_pro-1+ims_npro,ims_npro)
      ims_pro_r = MOD(ims_pro+1+ims_npro,ims_npro)

      icount   = ijmax*npl_loc

      CALL MPI_IRECV(bl(1,npl-kmax+1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
      CALL MPI_IRECV(br(1,1),          icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

      CALL MPI_ISEND(a(1,1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
      CALL MPI_ISEND(a(1,1), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

      CALL MPI_WAITALL(4, mpireq, status, ims_err)

      CALL DNS_MPI_TAGUPDT

      ! -----------------------------------------------------------------------
      ! second-left and second-right PEs.
      ! -----------------------------------------------------------------------
      npl_loc = npl-kmax

      ims_pro_l = MOD(ims_pro-2+ims_npro,ims_npro)
      ims_pro_r = MOD(ims_pro+2+ims_npro,ims_npro)

      icount   = ijmax*npl_loc

      CALL MPI_IRECV(bl(1,1),      icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(7), ims_err)
      CALL MPI_IRECV(br(1,1+kmax), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(8), ims_err)

      CALL MPI_ISEND(a(1,1),              icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(5), ims_err)
      CALL MPI_ISEND(a(1,kmax+1-npl_loc), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(6), ims_err)

      CALL MPI_WAITALL(4, mpireq(5), status, ims_err)

      CALL DNS_MPI_TAGUPDT

    END IF

    RETURN
  END SUBROUTINE DNS_MPI_COPYPLN_2

  ! ######################################################################
  ! Initialization of PSFFT Library for nonblocking communication
  ! ######################################################################
#ifdef USE_PSFFT
#include "nb3dfft_defines.inc"
#endif

  SUBROUTINE DNS_NB3DFFT_INITIALIZE
#ifdef USE_PSFFT
    USE DNS_GLOBAL, ONLY : imax,jmax,kmax
    USE DNS_GLOBAL, ONLY : g
    USE NB3DFFT, ONLY : nb3dfft_test_setup, nb3dfft_setup, get_dims
#endif

    IMPLICIT NONE

    ! #######################################################################
#ifdef USE_PSFFT
    CALL TLAB_WRITE_ASCII(lfile,'Initialize nonblocking communication.')

    ims_nb_proc_grid = (/ ims_npro_i, ims_npro_k /)
    CALL NB3DFFT_SETUP(ims_nb_proc_grid, g(1)%size,g(2)%size,g(3)%size, &
        ims_nb_msize)

    CALL GET_DIMS(ims_nb_xsrt,ims_nb_xend,ims_nb_xsiz,1,1)
    CALL GET_DIMS(ims_nb_ysrt,ims_nb_yend,ims_nb_ysiz,1,2)
    CALL GET_DIMS(ims_nb_zsrt,ims_nb_zend,ims_nb_zsiz,1,3)

    IF (       ims_nb_xsrt(1) == 1 .AND. ims_nb_xend(1) == g(1)%size &
        .AND. ims_nb_xsiz(2)*ims_nb_xsiz(3) == ims_size_i(DNS_MPI_I_PARTIAL) ) THEN
      ! Decomp standing in X okay
    ELSE
      CALL TLAB_WRITE_ASCII(efile,'Decomp standing in X-BAD')
      CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
    END IF

    IF (      ims_nb_ysrt(1) == ims_offset_i+1 &
        .AND.ims_nb_ysrt(2) == ims_offset_j+1 &
        .AND.ims_nb_ysrt(3) == ims_offset_k+1 &
        .AND.ims_nb_ysiz(1) == imax &
        .AND.ims_nb_ysiz(2) == jmax &
        .AND.ims_nb_ysiz(3) == kmax ) THEN
    ELSE
      CALL TLAB_WRITE_ASCII(efile,'Decomp standing in Y--BAD')
      CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
    END IF

    IF (      ims_nb_zsrt(3) == 1 .AND. ims_nb_zend(3) == g(3)%size &
        .AND.ims_nb_zsiz(1)*ims_nb_zsiz(2) == ims_size_k(DNS_MPI_K_PARTIAL) ) THEN
      ! Decomp standing in Z okay
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'Decomp standing in Z--BAD')
      CALL TLAB_STOP(DNS_ERROR_PARPARTITION)
    END IF

    CALL TLAB_WRITE_ASCII(lfile,'Checking that NB3DFFT and DNS domain decompositions agree.')

    CALL nb3dfft_test_setup()

#else
    CALL TLAB_WRITE_ASCII(efile,'Compiler flag USE_PSFFT needs to be used.')
    CALL TLAB_STOP(DNS_ERROR_PARPARTITION)

#endif

    RETURN
  END SUBROUTINE DNS_NB3DFFT_INITIALIZE

END MODULE TLAB_MPI_PROCS

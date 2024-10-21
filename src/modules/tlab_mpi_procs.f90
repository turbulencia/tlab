#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

module TLAB_MPI_PROCS
    use MPI
    use TLAB_CONSTANTS, only: lfile, efile, wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, isize_txc_dimx, isize_txc_dimz
    use TLAB_VARS, only: fourier_on
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use TLAB_MPI_VARS
    implicit none
    save
    private

    integer :: ims_tag

    public :: TLAB_MPI_INITIALIZE
    public :: TLAB_MPI_TRPF_K
    public :: TLAB_MPI_TRPF_I
    public :: TLAB_MPI_TRPB_K
    public :: TLAB_MPI_TRPB_I
    public :: TLAB_MPI_TYPE_K
    public :: TLAB_MPI_TYPE_I
    public :: TLAB_MPI_PANIC
    public :: TLAB_MPI_WRITE_PE0_SINGLE

contains

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_MPI_INITIALIZE

        ! -----------------------------------------------------------------------
        integer(wi) id, ip, npage
        integer(wi) dims(2)
        logical period(2), remain_dims(2), reorder

        ! #######################################################################
        allocate (ims_map_i(ims_npro_i))
        allocate (ims_size_i(TLAB_MPI_I_MAXTYPES))
        allocate (ims_ds_i(ims_npro_i, TLAB_MPI_I_MAXTYPES))
        allocate (ims_dr_i(ims_npro_i, TLAB_MPI_I_MAXTYPES))
        allocate (ims_ts_i(1, TLAB_MPI_I_MAXTYPES))     ! ims_npro_i, TLAB_MPI_I_MAXTYPES
        allocate (ims_tr_i(1, TLAB_MPI_I_MAXTYPES))     ! ims_npro_i, TLAB_MPI_I_MAXTYPES
        allocate (ims_plan_trps_i(ims_npro_i))
        allocate (ims_plan_trpr_i(ims_npro_i))

        allocate (ims_map_k(ims_npro_k))
        allocate (ims_size_k(TLAB_MPI_K_MAXTYPES))
        allocate (ims_ds_k(ims_npro_k, TLAB_MPI_K_MAXTYPES))
        allocate (ims_dr_k(ims_npro_k, TLAB_MPI_K_MAXTYPES))
        allocate (ims_ts_k(1, TLAB_MPI_K_MAXTYPES))     ! ims_npro_k,TLAB_MPI_K_MAXTYPES
        allocate (ims_tr_k(1, TLAB_MPI_K_MAXTYPES))     ! ims_npro_k,TLAB_MPI_K_MAXTYPES
        allocate (ims_plan_trps_k(ims_npro_k))
        allocate (ims_plan_trpr_k(ims_npro_k))

#ifdef HLRS_HAWK
        ! On hawk, we tested that 192 yields optimum performace;
        ! Blocking will thus only take effect in very large cases
        ims_sizBlock_k = 192
        ims_sizBlock_i = 384
#else
        ! We assume that this will help to release some of the very heavy
        ! network load in transpositions on most systems
        ims_sizBlock_k = 64
        ims_sizBlock_i = 128
        ! ims_sizBlock_k=1e5   -- would essentially switch off the blocking
#endif

        allocate (ims_status(MPI_STATUS_SIZE, 2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))
        allocate (ims_request(2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))

        ! #######################################################################
        ims_pro_i = mod(ims_pro, ims_npro_i) ! Starting at 0
        ! ims_pro_i = ims_pro / ims_npro_k
        ims_pro_k = ims_pro/ims_npro_i  ! Starting at 0
        ! ims_pro_k = mod(ims_pro,ims_npro_k)

        ims_offset_i = ims_pro_i*imax
        ims_offset_j = 0
        ims_offset_k = ims_pro_k*kmax

        ims_map_i(1) = ims_pro_k*ims_npro_i
        do ip = 2, ims_npro_i
            ims_map_i(ip) = ims_map_i(ip - 1) + 1
        end do

        ims_map_k(1) = ims_pro_i
        do ip = 2, ims_npro_k
            ims_map_k(ip) = ims_map_k(ip - 1) + ims_npro_i
        end do

        ! #######################################################################
        call TLAB_WRITE_ASCII(lfile, 'Initializing MPI communicators.')

        ! the first index in the grid corresponds to k, the second to i
        dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .true.; reorder = .false.
        ! dims(1) = ims_npro_i; dims(2) = ims_npro_k; period = .true.; reorder = .false.
        call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

        !  CALL MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
        !  coord(1) is ims_pro_k, and coord(2) is ims_pro_i

        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)
        !call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

        remain_dims(1) = .true.; remain_dims(2) = .false.
        call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)
        !call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

        ! ip = ims_pro
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_x, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along X', id
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_z, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along Z', id

        ! #######################################################################
        ! Main
        ! #######################################################################
        if (ims_npro_i > 1) then
            call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for Ox derivatives.')
            id = TLAB_MPI_I_PARTIAL
            npage = kmax*jmax
            call TLAB_MPI_TYPE_I(ims_npro_i, imax, npage, 1, 1, 1, 1, &
                                 ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        end if

        if (ims_npro_k > 1) then
            call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for Oz derivatives.')
            id = TLAB_MPI_K_PARTIAL
            npage = imax*jmax
            call TLAB_MPI_TYPE_K(ims_npro_k, kmax, npage, 1, 1, 1, 1, &
                                 ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        end if

        ! -----------------------------------------------------------------------
        if (ims_npro_i > 1 .and. fourier_on) then
            call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for Ox FFTW in Poisson solver.')
            id = TLAB_MPI_I_POISSON1
            npage = isize_txc_dimx ! isize_txc_field/imax
            call TLAB_MPI_TYPE_I(ims_npro_i, imax, npage, 1, 1, 1, 1, &
                                 ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for Ox FFTW in Poisson solver.')
            id = TLAB_MPI_I_POISSON2 ! isize_txc_field/(imax+2)
            npage = isize_txc_dimx
            call TLAB_MPI_TYPE_I(ims_npro_i, imax + 2, npage, 1, 1, 1, 1, &
                                 ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

        end if

        if (ims_npro_k > 1 .and. fourier_on) then
            call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for Oz FFTW in Poisson solver.')
            id = TLAB_MPI_K_POISSON
            npage = isize_txc_dimz ! isize_txc_field/kmax
            call TLAB_MPI_TYPE_K(ims_npro_k, kmax, npage, 1, 1, 1, 1, &
                                 ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        end if

        ! ######################################################################
        ! Work plans for circular transposes
        ! ######################################################################
        do ip = 0, ims_npro_i - 1
            ims_plan_trps_i(ip + 1) = ip
            ims_plan_trpr_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
        end do
        ims_plan_trps_i = cshift(ims_plan_trps_i, ims_pro_i)
        ims_plan_trpr_i = cshift(ims_plan_trpr_i, -(ims_pro_i))

        do ip = 0, ims_npro_k - 1
            ims_plan_trps_k(ip + 1) = ip
            ims_plan_trpr_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
        end do
        ims_plan_trps_k = cshift(ims_plan_trps_k, ims_pro_k)
        ims_plan_trpr_k = cshift(ims_plan_trpr_k, -(ims_pro_k))

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
        ! id = TLAB_MPI_K_SHEAR
        ! CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, imax, i1, jmax, jmax, i1, &
        !      ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
        ! END IF

        call TLAB_MPI_TAGRESET

        return
    end subroutine TLAB_MPI_INITIALIZE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_MPI_TAGUPDT

        ims_tag = ims_tag + 1
        if (ims_tag > 32000) then
            call TLAB_MPI_TAGRESET
        end if

        return
    end subroutine TLAB_MPI_TAGUPDT

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TAGRESET

        ims_tag = 0

        return
    end subroutine TLAB_MPI_TAGRESET

    ! ######################################################################
    ! Pointers and types for transposition across ims_npro processors
    ! ######################################################################
    subroutine TLAB_MPI_TYPE_I(npro_i, nmax, npage, nd, md, n1, n2, &
                               nsize, sdisp, rdisp, stype, rtype)

        integer npro_i
        integer(wi) npage, nmax, nsize
        integer(wi) nd, md, n1, n2
        integer(wi) sdisp(*), rdisp(*)
        integer, intent(out) :: stype, rtype

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_ss, ims_rs
        integer ims_tmp1, ims_tmp2, ims_tmp3
        character*64 str, line

        ! #######################################################################
        if (mod(npage, npro_i) == 0) then
            nsize = npage/npro_i
        else
            call TLAB_WRITE_ASCII(efile, 'TLAB_MPI_TYPE_I. Ratio npage/npro_i not an integer.')
            call TLAB_STOP(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate Displacements in Forward Send/Receive
        sdisp(1) = 0
        rdisp(1) = 0
        do i = 2, npro_i
            sdisp(i) = sdisp(i - 1) + nmax*nd*nsize
            rdisp(i) = rdisp(i - 1) + nmax*md
        end do

        ! #######################################################################
        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block because things are together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
        call MPI_TYPE_COMMIT(stype, ims_err)

        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = nmax*npro_i*n2 ! stride is a multiple of nmax_total=nmax*npro_i
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
        call MPI_TYPE_COMMIT(rtype, ims_err)

        call MPI_TYPE_SIZE(stype, ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            write (str, *) 1  ! i
            line = trim(adjustl(line))//' in message '//trim(adjustl(str))
            call TLAB_WRITE_ASCII(efile, line)
            call TLAB_STOP(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end subroutine TLAB_MPI_TYPE_I

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TYPE_K(npro_k, nmax, npage, nd, md, n1, n2, &
                               nsize, sdisp, rdisp, stype, rtype)

        integer npro_k
        integer(wi) npage, nmax, nsize
        integer(wi) nd, md, n1, n2
        integer(wi) sdisp(*), rdisp(*)
        integer, intent(out) :: stype, rtype

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_ss, ims_rs
        integer ims_tmp1, ims_tmp2, ims_tmp3

        ! #######################################################################
        if (mod(npage, npro_k) == 0) then
            nsize = npage/npro_k
        else
            call TLAB_WRITE_ASCII(efile, 'TLAB_MPI_TYPE_K. Ratio npage/npro_k not an integer.')
            call TLAB_STOP(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate Displacements in Forward Send/Receive
        sdisp(1) = 0
        rdisp(1) = 0
        do i = 2, npro_k
            sdisp(i) = sdisp(i - 1) + nsize*nd
            rdisp(i) = rdisp(i - 1) + nsize*md*nmax
        end do

        ! #######################################################################
        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = npage*n2 ! stride
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
        call MPI_TYPE_COMMIT(stype, ims_err)

        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block to put things together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
        call MPI_TYPE_COMMIT(rtype, ims_err)

        call MPI_TYPE_SIZE(stype, ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            print *, 'Message   : ', 1, ' size is wrong' ! i
            print *, 'Send size : ', ims_ss
            print *, 'Recv size : ', ims_rs
            call TLAB_STOP(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end subroutine TLAB_MPI_TYPE_K

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TRPF_K(a, b, dsend, drecv, tsend, trecv)
        use, intrinsic :: iso_c_binding, only: c_int

        interface
            function DNS_USLEEP(useconds) bind(C, name="usleep")
                import
                integer(c_int) :: nb3dfft_nbc_usleep
                integer(c_int), intent(in), value :: useconds
            end function DNS_USLEEP
        end interface

        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        integer(wi), dimension(ims_npro_k), intent(in) :: dsend, drecv ! displacements
        integer, intent(in) :: tsend, trecv

        ! -----------------------------------------------------------------------
        integer(wi) j, l, m, ns, nr, ips, ipr

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        do j = 1, ims_npro_k, ims_sizBlock_k
            l = 0
            do m = j, min(j + ims_sizBlock_k - 1, ims_npro_k)
                ns = ims_plan_trps_k(m) + 1; ips = ns - 1
                nr = ims_plan_trpr_k(m) + 1; ipr = nr - 1
                if (ims_trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    l = l + 1
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, ims_comm_z, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, ims_comm_z, ims_request(l), ims_err)
                elseif (ims_trp_mode_k == TLAB_MPI_TRP_SENDRECV) then
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, ims_comm_z, ims_status(:, 1), ims_err)
                else; continue     ! No transpose
                end if
            end do

            if (ims_trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) &
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)

            call TLAB_MPI_TAGUPDT
        end do
#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLAB_MPI_TRPF_K

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TRPF_I(a, b, dsend, drecv, tsend, trecv)

        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        integer(wi), dimension(ims_npro_i), intent(in) :: dsend, drecv ! displacements
        integer, intent(in) :: tsend, trecv

        ! -----------------------------------------------------------------------
        integer(wi) j, l, m, ns, nr, ips, ipr

        do j = 1, ims_npro_i, ims_sizBlock_i
            l = 0
            do m = j, min(j + ims_sizBlock_i - 1, ims_npro_i)
                ns = ims_plan_trps_i(m) + 1; ips = ns - 1
                nr = ims_plan_trpr_i(m) + 1; ipr = nr - 1
                if (ims_trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    l = l + 1
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, ims_comm_x, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, ims_comm_x, ims_request(l), ims_err)
                elseif (ims_trp_mode_i == TLAB_MPI_TRP_SENDRECV) then
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, ims_comm_x, ims_status(:, 1), ims_err)
                else; continue ! No transpose
                end if
            end do

            if (ims_trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) &
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)

            call TLAB_MPI_TAGUPDT
        end do

        return
    end subroutine TLAB_MPI_TRPF_I

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)

        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        integer(wi), dimension(ims_npro_k), intent(in) :: dsend, drecv
        integer, dimension(ims_npro_k), intent(in) :: tsend, trecv

        ! -----------------------------------------------------------------------
        integer(wi) j, l, m, ns, nr, ips, ipr

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        do j = 1, ims_npro_k, ims_sizBlock_k
            l = 0
            do m = j, min(j + ims_sizBlock_k - 1, ims_npro_k)
                ns = ims_plan_trps_k(m) + 1; ips = ns - 1
                nr = ims_plan_trpr_k(m) + 1; ipr = nr - 1
                if (ims_trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    l = l + 1
                    call MPI_ISEND(b(drecv(nr) + 1), 1, trecv(1), ipr, ims_tag, ims_comm_z, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(a(dsend(ns) + 1), 1, tsend(1), ips, ims_tag, ims_comm_z, ims_request(l), ims_err)
                elseif (ims_trp_mode_k == TLAB_MPI_TRP_SENDRECV) then
                    call MPI_SENDRECV(b(drecv(nr) + 1), 1, trecv(1), ipr, ims_tag, &
                                      a(dsend(ns) + 1), 1, tsend(1), ips, ims_tag, ims_comm_z, ims_status(:, 1), ims_err)
                else; continue   ! No transpose
                end if
            end do

            if (ims_trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) &
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)

            call TLAB_MPI_TAGUPDT
        end do

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLAB_MPI_TRPB_K

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_TRPB_I(b, a, dsend, drecv, tsend, trecv)

        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        integer(wi), dimension(ims_npro_i), intent(in) :: dsend, drecv ! displacements
        integer, dimension(ims_npro_i), intent(in) :: tsend, trecv ! types

        ! -----------------------------------------------------------------------
        integer(wi) j, l, m, ns, nr, ips, ipr

        do j = 1, ims_npro_i, ims_sizBlock_i
            l = 0
            do m = j, min(j + ims_sizBlock_i - 1, ims_npro_i)
                ns = ims_plan_trps_i(m) + 1; ips = ns - 1
                nr = ims_plan_trpr_i(m) + 1; ipr = nr - 1
                if (ims_trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    l = l + 1
                    call MPI_ISEND(b(drecv(nr) + 1), 1, trecv(1), ipr, ims_tag, ims_comm_x, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(a(dsend(ns) + 1), 1, tsend(1), ips, ims_tag, ims_comm_x, ims_request(l), ims_err)
                elseif (ims_trp_mode_i == TLAB_MPI_TRP_SENDRECV) then
                    call MPI_SENDRECV(b(drecv(nr) + 1), 1, trecv(1), ipr, ims_tag, &
                                      a(dsend(ns) + 1), 1, tsend(1), ips, ims_tag, ims_comm_x, ims_status(:, 1), ims_err)
                else; continue    ! No transpose
                end if
            end do

            if (ims_trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) &
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)

            call TLAB_MPI_TAGUPDT
        end do

        return
    end subroutine TLAB_MPI_TRPB_I

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_MPI_PANIC(location, mpi_error_code)
        character(len=*), intent(in) :: location
        integer, intent(in) :: mpi_error_code

        !##############################
        character error_string*1024
        integer error_local, error_len

        call MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
        call TLAB_WRITE_ASCII(efile, 'MPI-ERROR: Source file'//trim(adjustl(LOCATION)), .true.)
        call TLAB_WRITE_ASCII(efile, error_string, .true.)

        call TLAB_STOP(mpi_error_code)
        ! Not supposed to return from this subroutine

    end subroutine TLAB_MPI_PANIC

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_WRITE_PE0_SINGLE(iunit, nx, ny, nz, subdomain, u, tmp1, tmp2)

        integer(wi) iunit, nx, ny, nz, subdomain(6)
        real(wp), dimension(nx*ny*nz), target :: u, tmp1
        real(wp), dimension(nx*ims_npro_i, *), target :: tmp2

        ! -------------------------------------------------------------------
        integer(wi) nx_total, ny_total, nz_total
        integer(wi) nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
        integer(wi) nyz

        integer(wi) ip_i, ip_k, joffset_loc, koffset_loc, id
        integer(wi) i, jk, j_loc, k_loc
        integer mpio_size, mpio_ip
        integer status(MPI_STATUS_SIZE)

        real(wp), dimension(:), pointer :: p_org

        ! ###################################################################
        nx_total = nx*ims_npro_i
        ny_total = ny
        nz_total = nz*ims_npro_k

        nx_min = subdomain(1); nx_max = subdomain(2)
        ny_min = subdomain(3); ny_max = subdomain(4)
        nz_min = subdomain(5); nz_max = subdomain(6)

        koffset_loc = 0
        joffset_loc = 0

        id = TLAB_MPI_I_PARTIAL

        ! -------------------------------------------------------------------
        ! Transposing along Ox
        ! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            call TLAB_MPI_TRPF_I(u, tmp1, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
            p_org => tmp1
            nyz = ims_size_i(id)
        else
            p_org => u
            nyz = ny*nz
        end if
        mpio_size = nyz*nx_total

        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        ! -------------------------------------------------------------------
        ! Passing all data through PE#0
        ! -------------------------------------------------------------------
        if (ims_pro == 0) then

            do ip_k = 1, ims_npro_k
                koffset_loc = nz*(ip_k - 1)

                do ip_i = 1, ims_npro_i
                    joffset_loc = nyz*(ip_i - 1) ! Remember that data is Ox-transposed

                    mpio_ip = ims_npro_i*(ip_k - 1) + ip_i - 1
                    if (mpio_ip == 0) then
                        tmp2(1:mpio_size, 1) = p_org(1:mpio_size)
                    else
                        call MPI_RECV(tmp2, mpio_size, MPI_REAL8, mpio_ip, ims_tag, MPI_COMM_WORLD, status, ims_err)
                    end if

                    do jk = 1, nyz
                        j_loc = mod((jk - 1 + joffset_loc), ny_total) + 1
                        k_loc = ((jk - 1 + joffset_loc)/ny_total) + 1 + koffset_loc

                        if ((j_loc >= ny_min) .and. (j_loc <= ny_max) .and. &
                            (k_loc >= nz_min) .and. (k_loc <= nz_max)) then
                            write (iunit) (SNGL(tmp2(i, jk)), i=nx_min, nx_max)
                        end if

                    end do

                end do
            end do

        else
            call MPI_SEND(p_org, mpio_size, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
        end if

        return
    end subroutine TLAB_MPI_WRITE_PE0_SINGLE

    !########################################################################
    !# Moving plane information between adjacent PEs, circulant version.
    !# npl is smaller than 2*kmax
    !# The number of plane to move is given by npl
    !########################################################################
    subroutine TLAB_MPI_COPYPLN_1(ijmax, kmax, npl, a, bl, br)

        integer(wi) ijmax, kmax, npl
        real(wp) a(ijmax, *)
        real(wp) bl(ijmax, *)
        real(wp) br(ijmax, *)

        ! -----------------------------------------------------------------------
        integer status(MPI_STATUS_SIZE, 4)
        integer mpireq(4)
        integer ims_pro_l, ims_pro_r
        integer icount

        ! #######################################################################
        if (ims_npro > 1) then

            ! Careful in case only 2 PEs
            if (ims_npro == 2) then
                call TLAB_WRITE_ASCII(efile, 'TLAB_MPI_COPYPLN_1. Undeveloped for 2 PEs.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            ! left and right PEs
            ims_pro_l = mod(ims_pro - 1 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 1 + ims_npro, ims_npro)

            icount = ijmax*npl

            call MPI_IRECV(bl(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
            call MPI_IRECV(br(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
            call MPI_ISEND(a(1, kmax + 1 - npl), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

            call MPI_WAITALL(4, mpireq, status, ims_err)

            call TLAB_MPI_TAGUPDT

        end if

        return
    end subroutine TLAB_MPI_COPYPLN_1

    !########################################################################
    !########################################################################
    subroutine TLAB_MPI_COPYPLN_2(ijmax, kmax, npl, a, bl, br)

        integer(wi) ijmax, kmax, npl
        real(wp) a(ijmax, kmax)
        real(wp) bl(ijmax, npl)
        real(wp) br(ijmax, npl)

        ! -----------------------------------------------------------------------
        integer(wi) npl_loc
        integer status(MPI_STATUS_SIZE, 8)
        integer mpireq(8)
        integer ims_pro_l, ims_pro_r
        integer icount

        ! #######################################################################
        if (ims_npro > 1) then

            ! Careful in case only 2 PEs
            if (ims_npro == 2) then
                call TLAB_WRITE_ASCII(efile, 'TLAB_MPI_COPYPLN_2. Undeveloped for 2 PEs.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            ! -----------------------------------------------------------------------
            ! left and right PEs. Same as in routine TLAB_MPI_COPYPLN_1
            ! -----------------------------------------------------------------------
            npl_loc = kmax

            ims_pro_l = mod(ims_pro - 1 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 1 + ims_npro, ims_npro)

            icount = ijmax*npl_loc

            call MPI_IRECV(bl(1, npl - kmax + 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
            call MPI_IRECV(br(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

            call MPI_WAITALL(4, mpireq, status, ims_err)

            call TLAB_MPI_TAGUPDT

            ! -----------------------------------------------------------------------
            ! second-left and second-right PEs.
            ! -----------------------------------------------------------------------
            npl_loc = npl - kmax

            ims_pro_l = mod(ims_pro - 2 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 2 + ims_npro, ims_npro)

            icount = ijmax*npl_loc

            call MPI_IRECV(bl(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(7), ims_err)
            call MPI_IRECV(br(1, 1 + kmax), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(8), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(5), ims_err)
            call MPI_ISEND(a(1, kmax + 1 - npl_loc), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(6), ims_err)

            call MPI_WAITALL(4, mpireq(5:), status, ims_err)

            call TLAB_MPI_TAGUPDT

        end if

        return
    end subroutine TLAB_MPI_COPYPLN_2

    ! ######################################################################
    ! Initialization of PSFFT Library for nonblocking communication
    ! ######################################################################
#ifdef USE_PSFFT
#include "nb3dfft_defines.inc"
#endif

    subroutine DNS_NB3DFFT_INITIALIZE
#ifdef USE_PSFFT
        use NB3DFFT, only: nb3dfft_test_setup, nb3dfft_setup, get_dims
#endif

        implicit none

        ! #######################################################################
#ifdef USE_PSFFT
        call TLAB_WRITE_ASCII(lfile, 'Initialize nonblocking communication.')

        ims_nb_proc_grid = (/ims_npro_i, ims_npro_k/)
        call NB3DFFT_SETUP(ims_nb_proc_grid, g(1)%size, g(2)%size, g(3)%size, &
                           ims_nb_msize)

        call GET_DIMS(ims_nb_xsrt, ims_nb_xend, ims_nb_xsiz, 1, 1)
        call GET_DIMS(ims_nb_ysrt, ims_nb_yend, ims_nb_ysiz, 1, 2)
        call GET_DIMS(ims_nb_zsrt, ims_nb_zend, ims_nb_zsiz, 1, 3)

        if (ims_nb_xsrt(1) == 1 .and. ims_nb_xend(1) == g(1)%size &
            .and. ims_nb_xsiz(2)*ims_nb_xsiz(3) == ims_size_i(TLAB_MPI_I_PARTIAL)) then
            ! Decomp standing in X okay
        else
            call TLAB_WRITE_ASCII(efile, 'Decomp standing in X-BAD')
            call TLAB_STOP(DNS_ERROR_PARPARTITION)
        end if

        if (ims_nb_ysrt(1) == ims_offset_i + 1 &
            .and. ims_nb_ysrt(2) == ims_offset_j + 1 &
            .and. ims_nb_ysrt(3) == ims_offset_k + 1 &
            .and. ims_nb_ysiz(1) == imax &
            .and. ims_nb_ysiz(2) == jmax &
            .and. ims_nb_ysiz(3) == kmax) then
        else
            call TLAB_WRITE_ASCII(efile, 'Decomp standing in Y--BAD')
            call TLAB_STOP(DNS_ERROR_PARPARTITION)
        end if

        if (ims_nb_zsrt(3) == 1 .and. ims_nb_zend(3) == g(3)%size &
            .and. ims_nb_zsiz(1)*ims_nb_zsiz(2) == ims_size_k(TLAB_MPI_K_PARTIAL)) then
            ! Decomp standing in Z okay
        else
            call TLAB_WRITE_ASCII(efile, 'Decomp standing in Z--BAD')
            call TLAB_STOP(DNS_ERROR_PARPARTITION)
        end if

        call TLAB_WRITE_ASCII(lfile, 'Checking that NB3DFFT and DNS domain decompositions agree.')

        call nb3dfft_test_setup()

#else
        call TLAB_WRITE_ASCII(efile, 'Compiler flag USE_PSFFT needs to be used.')
        call TLAB_STOP(DNS_ERROR_PARPARTITION)

#endif

        return
    end subroutine DNS_NB3DFFT_INITIALIZE

end module TLAB_MPI_PROCS

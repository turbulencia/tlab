!mpif90 -fpp  -nbs -save-temps -xHost -simd -vec-threshold50 -unroll-aggressive    -axcommon-avx512,SSE4.2  -qopt-prefetch -O3 vmpi_transpose.f90
!mpif90 -cpp -ffree-form -ffree-line-length-2048 -fno-automatic -O3 -fconvert=little-endian -mtune=native -ffast-math -ffinite-math-only -funroll-loops
! from dns_const.h

! from dns_const_mpi.h
#define TLabMPI_K_PARTIAL   1 ! tags and sizes for MPI data
#define TLabMPI_I_PARTIAL   1

#define TLabMPI_K_MAXTYPES 10
#define TLabMPI_I_MAXTYPES  6

module DNS_MPI
    implicit none
    save

    logical :: ims_trp_blocking

    integer(KIND=4) :: imax, jmax, kmax
    integer(KIND=4) :: imax_total, jmax_total, kmax_total

    integer :: ims_npro
    integer :: ims_npro_i, ims_npro_j, ims_npro_k     ! number of tasks in Ox and Oz (no decomposition along Oy)

    integer(KIND=4) :: nmax   ! number of repetitions of operations

! Data below should not be changed
    integer :: ims_pro, ims_pro_i, ims_pro_j, ims_pro_k ! task positioning

    integer :: ims_comm_xz, ims_comm_x, ims_comm_z      ! communicators
    integer :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_err, ims_tag

    integer, dimension(:), allocatable :: ims_map_i
    integer(KIND=4), dimension(:), allocatable :: ims_size_i
    integer(KIND=4), dimension(:, :), allocatable :: ims_ds_i, ims_dr_i
    integer, dimension(:, :), allocatable :: ims_ts_i, ims_tr_i

    integer, dimension(:), allocatable :: ims_map_k
    integer(KIND=4), dimension(:), allocatable :: ims_size_k
    integer(KIND=4), dimension(:, :), allocatable :: ims_ds_k, ims_dr_k
    integer, dimension(:, :), allocatable :: ims_ts_k, ims_tr_k

    integer(KIND=4), dimension(:), allocatable :: ims_plan_trps_i, ims_plan_trpr_i
    integer(KIND=4), dimension(:), allocatable :: ims_plan_trps_k, ims_plan_trpr_k

end module DNS_MPI

!########################################################################
! Main program to test forwards and backwards transposition
!########################################################################
program VMPI_RTRANSPOSE

    use TLabMPI_VARS

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(:, :), allocatable :: a
    real(KIND=8), dimension(:), allocatable :: wrk3d

! -------------------------------------------------------------------
    real(KIND=8) residual                                      ! Control
    integer(KIND=4) t_srt, t_end, t_dif, PROC_CYCLES, MAX_CYCLES ! Time
    integer(KIND=4) n
    character*64 str
    character*256 line

    real(KIND=8) rdum
    integer(KIND=4) idum, id, narg
    character*32 cdum

    real(KIND=8) :: t_x, t_x2, t_z, t_z2

! ###################################################################
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

    narg = command_argument_count()
    if (narg < 4) then
        if (ims_pro == 0) write (*, *) 'Usage ./vmpi_rtranspose.x <nrun> <nx> <ny> <nz>'
        call MPI_FINALIZE(ims_err)
        stop
    end if

! Master rank processes input
    call GETARG(1, cdum); read (cdum, *) nmax
    call GETARG(2, cdum); read (cdum, *) imax_total
    call GETARG(3, cdum); read (cdum, *) jmax_total
    call GETARG(4, cdum); read (cdum, *) kmax_total

    ims_npro_i = 2**((exponent(real(ims_npro)))/2)
    ims_npro_j = 1
    ims_npro_k = ims_npro/ims_npro_i
    imax = imax_total/ims_npro_i
    jmax = jmax_total/ims_npro_j
    kmax = kmax_total/ims_npro_k

    if (ims_pro == 0) then
        write (*, *) '=== Initialization of Grid and Decomposition ==='
        write (*, *) 'GRID:        ', imax_total, ' x ', jmax_total, ' x ', kmax_total
        write (*, *) 'DECOMP: ranks', ims_npro_i, ' x ', ims_npro_j, ' x ', ims_npro_k
        write (*, *) '        grid ', imax, ' x ', jmax, ' x ', kmax
    end if

    if (ims_npro_i*ims_npro_k /= ims_npro .or. &
        ims_npro_i*imax /= imax_total .or. &
        ims_npro_k*kmax /= kmax_total .or. &
        ims_npro_j*jmax /= jmax_total) then ! check
        if (ims_pro == 0) write (*, '(a)') ims_pro, ': Inconsistency in Decomposition'
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
        call MPI_FINALIZE(ims_err)
        stop
    end if

    call TLabMPI_Initialize()

    allocate (a(imax*jmax*kmax, 18)) ! Number of 3d arrays commonly used in the code
    allocate (wrk3d(imax*jmax*kmax))

! ###################################################################
! ###################################################################
! Create random array
    call random_number(a(1:imax*jmax*kmax, 1))

    if (ims_pro == 0) &
        write (*, *) 'Executing everything once to get caches / stack / network in production state'

    if (ims_npro_k > 1) then
        id = TLabMPI_K_PARTIAL
        call TLabMPI_TRPF_K(a(1, 1), wrk3d, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        call TLabMPI_TRPB_K(wrk3d, a(1, 2), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
    end if
    if (ims_npro_i > 1) then
        id = TLabMPI_I_PARTIAL
        call TLabMPI_TRPF_I(a(1, 1), wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        call TLabMPI_TRPB_I(wrk3d, a(1, 2), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
    end if

    if (IMS_PRO == 0) then
        write (*, *) '======'
        write (*, *) '===== STARTING MEASUREMENT ====='
        write (*, *) '====='
    end if
    do n = 0, 2*nmax - 1

        if (n == 0) then
            ims_trp_blocking = .true.
            t_x = 0.; t_z = 0.; t_x2 = 0.; t_z2 = 0.; 
            if (ims_pro == 0) write (*, *) '======== BLOCKING TRANSPOSES ========'
        elseif (n == nmax) then
            ims_trp_blocking = .false.
            t_x = 0.; t_z = 0.; t_x2 = 0.; t_z2 = 0.; 
            if (ims_pro == 0) write (*, *) '====== NONBLOCKING TRANSPOSES========'
        end if
! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            id = TLabMPI_I_PARTIAL

            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)

            call TLabMPI_TRPF_I(a(1, 1), wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
            call TLabMPI_TRPB_I(wrk3d, a(1, 2), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

            idum = t_end - t_srt
            call MPI_REDUCE(idum, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            t_x = t_x + real(t_dif)/PROC_CYCLES

            rdum = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(rdum, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (line, '(E13.5E3)') residual

            line = ' transposition for Ox derivatives: Residual ' &
                   //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
            if (residual == 0) then
                line = 'PASSED'//' transposition for Ox derivatives '//trim(adjustl(str))//' sec.'
            else
                line = 'FAILED'//' transposition for Ox derivatives. Residual:'//trim(adjustl(line)) &
                       //' Duration '//trim(adjustl(str))//' sec.'
            end if

            if (ims_pro == 0) then
                write (*, '(a)') trim(adjustl(line))
            end if

        end if

! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
        if (ims_npro_k > 1) then
            id = TLabMPI_K_PARTIAL

            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)

            call TLabMPI_TRPF_K(a(1, 1), wrk3d, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            call TLabMPI_TRPB_K(wrk3d, a(1, 2), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))

            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

            idum = t_end - t_srt
            call MPI_REDUCE(idum, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            t_z = t_z + real(t_dif)/PROC_CYCLES

            rdum = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(rdum, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (line, '(E13.5E3)') residual

            if (residual == 0) then
                line = 'PASSED'//' transposition for Oz derivatives '//trim(adjustl(str))//' sec.'
            else
                line = 'FAILED'//' transposition for Oz derivatives. Residual:'//trim(adjustl(line)) &
                       //' Duration '//trim(adjustl(str))//' sec.'
            end if

            if (ims_pro == 0) then
                write (*, '(a)') trim(adjustl(line))
            end if

        end if

        if ((n == nmax - 1 .or. n == 2*nmax - 1) .and. ims_pro == 0) then
            write (*, *) 'X TRANSPOSES TIMING: ', t_x/nmax
            write (*, *) 'Z TRANSPOSES TIMING: ', t_z/nmax
        end if
    end do

    call MPI_FINALIZE(ims_err)

end program VMPI_RTRANSPOSE

! #######################################################################
! Rest of routines
! #######################################################################
subroutine TLabMPI_Initialize()

    use TLabMPI_VARS

    implicit none

#include "mpif.h"

! -----------------------------------------------------------------------
    integer(KIND=4) id, ip, npage
    integer(KIND=4) i1, dims(2)
    logical period(2), remain_dims(2), reorder

! #######################################################################
    allocate (ims_map_i(ims_npro_i))
    allocate (ims_size_i(TLabMPI_I_MAXTYPES))
    allocate (ims_ds_i(ims_npro_i, TLabMPI_I_MAXTYPES))
    allocate (ims_dr_i(ims_npro_i, TLabMPI_I_MAXTYPES))
    allocate (ims_ts_i(ims_npro_i, TLabMPI_I_MAXTYPES))
    allocate (ims_tr_i(ims_npro_i, TLabMPI_I_MAXTYPES))

    allocate (ims_map_k(ims_npro_k))
    allocate (ims_size_k(TLabMPI_K_MAXTYPES))
    allocate (ims_ds_k(ims_npro_k, TLabMPI_K_MAXTYPES))
    allocate (ims_dr_k(ims_npro_k, TLabMPI_K_MAXTYPES))
    allocate (ims_ts_k(ims_npro_k, TLabMPI_K_MAXTYPES))
    allocate (ims_tr_k(ims_npro_k, TLabMPI_K_MAXTYPES))

    allocate (ims_plan_trps_i(ims_npro_i))
    allocate (ims_plan_trpr_i(ims_npro_i))
    allocate (ims_plan_trps_k(ims_npro_k))
    allocate (ims_plan_trpr_k(ims_npro_k))

! #######################################################################
    ims_pro_i = mod(ims_pro, ims_npro_i) ! Starting at 0
    ims_pro_k = ims_pro/ims_npro_i  ! Starting at 0

    ims_map_i(1) = ims_pro_k*ims_npro_i
    do ip = 2, ims_npro_i
        ims_map_i(ip) = ims_map_i(ip - 1) + 1
    end do

    ims_map_k(1) = ims_pro_i
    do ip = 2, ims_npro_k
        ims_map_k(ip) = ims_map_k(ip - 1) + ims_npro_i
    end do

! #######################################################################
! Communicators
! #######################################################################
! the first index in the grid corresponds to k, the second to i
    dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .true.; reorder = .false.
    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

!  CALL MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
!  coord(1) is ims_pro_k, and coord(2) is ims_pro_i

    remain_dims(1) = .false.; remain_dims(2) = .true.
    call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

    remain_dims(1) = .true.; remain_dims(2) = .false.
    call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

! #######################################################################
! Derived MPI types to deal with the strides when tranposing data
! #######################################################################
    i1 = 1

    if (ims_npro_i > 1) then
!  CALL TLab_Write_ASCII(lfile,'Initializing MPI types for Ox derivatives.')
        id = TLabMPI_I_PARTIAL
        npage = kmax*jmax
        call TLabMPI_TYPE_I(ims_npro_i, imax, npage, i1, i1, i1, i1, &
                             ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
    end if

    if (ims_npro_k > 1) then
!  CALL TLab_Write_ASCII(lfile,'Initializing MPI types for Oz derivatives.')
        id = TLabMPI_K_PARTIAL
        npage = imax*jmax
        call TLabMPI_TYPE_K(ims_npro_k, kmax, npage, i1, i1, i1, i1, &
                             ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
    end if

! ######################################################################
! Work plans for circular transposes
! ######################################################################
    do ip = 0, ims_npro_i - 1
        ims_plan_trps_i(ip + 1) = ip
        ims_plan_trpr_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
    end do

    do ip = 0, ims_npro_k - 1
        ims_plan_trps_k(ip + 1) = ip
        ims_plan_trpr_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
    end do

    ims_plan_trps_i = cshift(ims_plan_trps_i, ims_pro_i)
    ims_plan_trpr_i = cshift(ims_plan_trpr_i, -(ims_pro_i))

    ims_plan_trps_k = cshift(ims_plan_trps_k, ims_pro_k)
    ims_plan_trpr_k = cshift(ims_plan_trpr_k, -(ims_pro_k))

    do ip = 0, ims_npro_i - 1
        if (ims_pro == ip) then
            write (*, *) ims_pro, ims_pro_i, 'SEND:', ims_plan_trps_i
            write (*, *) ims_pro, ims_pro_i, 'RECV:', ims_plan_trpr_i
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    end do

    call TLabMPI_TAGRESET

    return
end subroutine TLabMPI_Initialize

! ###################################################################
! ###################################################################
subroutine TLabMPI_TYPE_I(ims_npro, imax, npage, nd, md, n1, n2, &
                           nsize, sdisp, rdisp, stype, rtype)

    use TLabMPI_VARS, only: ims_pro

    implicit none

#include "mpif.h"

    integer ims_npro
    integer(KIND=4) npage, imax, nsize
    integer(KIND=4) nd, md, n1, n2
    integer(KIND=4) sdisp(*), rdisp(*)
    integer stype(*), rtype(*)

! -----------------------------------------------------------------------
    integer(KIND=4) i
    integer ims_ss, ims_rs, ims_err
    integer ims_tmp1, ims_tmp2, ims_tmp3
!  CHARACTER*64 str, line

! #######################################################################
    if (mod(npage, ims_npro) == 0) then
        nsize = npage/ims_npro
    else
        if (ims_pro == 0) then
            write (*, '(a)') 'Ratio npage/ims_npro_i not an integer'
        end if
        call MPI_FINALIZE(ims_err)
        stop
    end if

! Calculate Displacements in Forward Send/Receive
    sdisp(1) = 0
    rdisp(1) = 0
    do i = 2, ims_npro
        sdisp(i) = sdisp(i - 1) + imax*nd*nsize
        rdisp(i) = rdisp(i - 1) + imax*md
    end do

! #######################################################################
    do i = 1, ims_npro

        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = imax*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block because things are together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype(i), ims_err)
        call MPI_TYPE_COMMIT(stype(i), ims_err)

        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = imax*n2 ! block
        ims_tmp3 = imax*ims_npro*n2 ! stride is a multiple of imax_total=imax*ims_npro_i
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype(i), ims_err)
        call MPI_TYPE_COMMIT(rtype(i), ims_err)

        call MPI_TYPE_SIZE(stype(i), ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype(i), ims_rs, ims_err)

    end do

    return
end subroutine TLabMPI_TYPE_I

!########################################################################
!########################################################################
subroutine TLabMPI_TYPE_K(ims_npro, nmax, npage, nd, md, n1, n2, &
                           nsize, sdisp, rdisp, stype, rtype)

    use TLabMPI_VARS, only: ims_pro

    implicit none

#include "mpif.h"

    integer ims_npro
    integer(KIND=4) npage, nmax, nsize
    integer(KIND=4) nd, md, n1, n2
    integer(KIND=4) sdisp(*), rdisp(*)
    integer stype(*), rtype(*)

! -----------------------------------------------------------------------
    integer(KIND=4) i
    integer ims_ss, ims_rs, ims_err
    integer ims_tmp1, ims_tmp2, ims_tmp3

! #######################################################################
    if (mod(npage, ims_npro) == 0) then
        nsize = npage/ims_npro
    else
        if (ims_pro == 0) then
            write (*, '(a)') 'Ratio npage/ims_npro_k not an integer'
        end if
        call MPI_FINALIZE(ims_err)
        stop
    end if

! Calculate Displacements in Forward Send/Receive
    sdisp(1) = 0
    rdisp(1) = 0
    do i = 2, ims_npro
        sdisp(i) = sdisp(i - 1) + nsize*nd
        rdisp(i) = rdisp(i - 1) + nsize*md*nmax
    end do

! #######################################################################
    do i = 1, ims_npro

        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = npage*n2 ! stride
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype(i), ims_err)
        call MPI_TYPE_COMMIT(stype(i), ims_err)

        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block to put things together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype(i), ims_err)
        call MPI_TYPE_COMMIT(rtype(i), ims_err)

        call MPI_TYPE_SIZE(stype(i), ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype(i), ims_rs, ims_err)

    end do

    return
end subroutine TLabMPI_TYPE_K

! ###################################################################
! ###################################################################
subroutine TLabMPI_TRPF_K(a, b, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: ims_plan_trps_k, ims_plan_trpr_k, ims_trp_blocking

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: a
    real(KIND=8), dimension(*), intent(OUT) :: b
    integer(KIND=4), dimension(ims_npro_k), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_k), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) l, m, ns, nr
    integer status(MPI_STATUS_SIZE, 2*ims_npro_k)
    integer mpireq(2*ims_npro_k)
    integer ips, ipr

#ifdef PROFILE_ON
    real(KIND=8) time_loc_1, time_loc_2
#endif

! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif

    l = 0
    do m = 1, ims_npro_k
        ns = ims_plan_trps_k(m) + 1; ips = ns - 1
        nr = ims_plan_trpr_k(m) + 1; ipr = nr - 1
        if (.not. ims_trp_blocking) then
            l = l + 1
            call MPI_ISEND(a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_z, mpireq(l), ims_err)
        else
            call MPI_SENDRECV( &
                a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, &
                b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_z, status(1, 1), ims_err)
        end if

    end do

    if (.not. ims_trp_blocking) &
        call MPI_WAITALL(ims_npro_k*2, mpireq(1:), status(1, 1), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPF_K

!########################################################################
!########################################################################
subroutine TLabMPI_TRPF_I(a, b, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: ims_plan_trpr_i, ims_plan_trps_i
    use TLabMPI_VARS, only: ims_trp_blocking

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: a
    real(KIND=8), dimension(*), intent(OUT) :: b
    integer(KIND=4), dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_i), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) l, m
    integer status(MPI_STATUS_SIZE, 2*ims_npro_i)
    integer mpireq(2*ims_npro_i)
    integer ips, ipr, ns, nr

    l = 0

    do m = 1, ims_npro_i

        ns = ims_plan_trps_i(m) + 1; ips = ns - 1
        nr = ims_plan_trpr_i(m) + 1; ipr = nr - 1

        if (.not. ims_trp_blocking) then
            l = l + 1
            call MPI_ISEND(a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_x, mpireq(l), ims_err)
        else
            call MPI_SENDRECV( &
                a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, &
                b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_x, status(1, 1), ims_err)
        end if
    end do

    if (.not. ims_trp_blocking) &
        call MPI_WAITALL(ims_npro_i*2, mpireq(1:), status(1, 1), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPF_I

!########################################################################
!########################################################################
subroutine TLabMPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: ims_plan_trps_k, ims_plan_trpr_k, ims_trp_blocking

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: b
    real(KIND=8), dimension(*), intent(OUT) :: a
    integer(KIND=4), dimension(ims_npro_k), intent(IN) :: dsend, drecv
    integer, dimension(ims_npro_k), intent(IN) :: tsend, trecv

! -----------------------------------------------------------------------
    integer(KIND=4) l, m
    integer status(MPI_STATUS_SIZE, 2*ims_npro_k)
    integer mpireq(2*ims_npro_k)
    integer ips, ipr, ns, nr

#ifdef PROFILE_ON
    real(KIND=8) time_loc_1, time_loc_2
#endif

! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif

! #######################################################################
! Different processors
! #######################################################################
    l = 0
    !DO n = 1,ims_npro_k
    do m = 1, ims_npro_k
        ns = ims_plan_trps_k(m) + 1; ips = ns - 1
        nr = ims_plan_trpr_k(m) + 1; ipr = nr - 1
        if (.not. ims_trp_blocking) then
            l = l + 1
            call MPI_ISEND(b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_z, mpireq(l), ims_err)
        else
            call MPI_SENDRECV( &
                b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, &
                a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_z, status(1, m), ims_err)
        end if
    end do

    if (.not. ims_trp_blocking) &
        call MPI_WAITALL(ims_npro_k*2, mpireq(1:), status(1, 1), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPB_K

!########################################################################
!########################################################################
subroutine TLabMPI_TRPB_I(b, a, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: ims_plan_trpr_i, ims_plan_trps_i, ims_trp_blocking

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: b
    real(KIND=8), dimension(*), intent(OUT) :: a
    integer(KIND=4), dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_i), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) ns, nr, m, l
    integer status(MPI_STATUS_SIZE, 2*ims_npro_i)
    integer mpireq(2*ims_npro_i)
    integer ips, ipr

    l = 0
    do m = 1, ims_npro_i
        ns = ims_plan_trps_i(m) + 1; ips = ns - 1
        nr = ims_plan_trpr_i(m) + 1; ipr = nr - 1
        if (.not. ims_trp_blocking) then
            l = l + 1
            call MPI_ISEND(b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_x, mpireq(l), ims_err)
        else
            call MPI_SENDRECV( &
                b(drecv(nr) + 1), 1, trecv(nr), ipr, ims_tag, &
                a(dsend(ns) + 1), 1, tsend(ns), ips, ims_tag, ims_comm_x, status(1, m), ims_err)
        end if
    end do

    if (.not. ims_trp_blocking) &
        call MPI_WAITALL(ims_npro_i*2, mpireq(1:), status(1, 1), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPB_I

!########################################################################
!########################################################################
subroutine TLabMPI_TAGUPDT

    use TLabMPI_VARS, only: ims_tag

    implicit none

    ims_tag = ims_tag + 1

    if (ims_tag > 32000) then
        call TLabMPI_TAGRESET
    end if

    return
end subroutine TLabMPI_TAGUPDT

!########################################################################
!########################################################################
subroutine TLabMPI_TAGRESET

    use TLabMPI_VARS, only: ims_tag

    implicit none

    ims_tag = 0

    return
end subroutine TLabMPI_TAGRESET

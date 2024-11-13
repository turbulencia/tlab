!mpif90 -fpp  -nbs -save-temps -xHost -simd -vec-threshold50 -unroll-aggressive    -axcommon-avx512,SSE4.2  -qopt-prefetch -O3 vmpi_transpose.f90

! from dns_const.h
#define TREAL      REAL(8)    /* user-defined types */
#define TINTEGER   INTEGER(4)

! from dns_const_mpi.h
#define TLabMPI_K_PARTIAL   1 /* tags and sizes for MPI data*/
#define TLabMPI_I_PARTIAL   1

#define TLabMPI_K_MAXTYPES 10
#define TLabMPI_I_MAXTYPES  6

module DNS_MPI
    implicit none
    save

    integer(KIND=4) :: imax_total, jmax_total, kmax_total = 84   ! number of grid points per task
    integer(KIND=4) :: imax, jmax, kmax

    integer :: ims_npro_i, ims_npro_j, ims_npro_k ! number of tasks in Ox and Oz (no decomposition along Oy)

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

end module DNS_MPI

!########################################################################
! Main program to test forwards and backwards transposition
!########################################################################
program VMPI

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

    integer ims_npro
    real(KIND=8) dummy, t_x, t_z
    integer(KIND=4) idummy, id, narg, nmax
    character*32 cdum

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
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)  ! Make sure ims_pro gets to write the message above
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

    do n = 1, nmax

! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            id = TLabMPI_I_PARTIAL

            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)

            call TLabMPI_TRPF_I(a(1, 1), wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
            call TLabMPI_TRPB_I(wrk3d, a(1, 2), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

            idummy = t_end - t_srt
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            t_x = t_x + real(t_dif)/PROC_CYCLES

            dummy = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (line, '(E13.5E3)') residual

            line = 'Checking MPI transposition for Ox derivatives: Residual ' &
                   //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
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

            idummy = t_end - t_srt
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            t_z = t_z + real(t_dif)/PROC_CYCLES

            dummy = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (line, '(E13.5E3)') residual

            line = 'Checking MPI transposition for Oz derivatives: Residual ' &
                   //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
            if (ims_pro == 0) then
                write (*, '(a)') trim(adjustl(line))
            end if

        end if

    end do

    if (ims_pro == 0) then
        write (*, *) 'X TIMING', t_x/nmax, 'sec.'
        write (*, *) 'Z TIMING', t_z/nmax, 'sec.'
    end if

    call MPI_FINALIZE(ims_err)
    stop

end program VMPI

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

    use TLabMPI_VARS, only: ims_npro_k, ims_pro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: a
    real(KIND=8), dimension(*), intent(OUT) :: b
    integer(KIND=4), dimension(ims_npro_k), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_k), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) n, l
    integer status(MPI_STATUS_SIZE, 2*ims_npro_k)
    integer mpireq(2*ims_npro_k)
    integer ip

#ifdef PROFILE_ON
    real(KIND=8) time_loc_1, time_loc_2
#endif

! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_k; n = ip + 1
    call MPI_ISEND(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(1), ims_err)
    call MPI_IRECV(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_k
        ip = n - 1
        if (ip /= ims_pro_k) then
            l = l + 1
            call MPI_ISEND(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_k*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPF_K

!########################################################################
!########################################################################
subroutine TLabMPI_TRPF_I(a, b, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_i, ims_pro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: a
    real(KIND=8), dimension(*), intent(OUT) :: b
    integer(KIND=4), dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_i), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) n, l
    integer status(MPI_STATUS_SIZE, 2*ims_npro_i)
    integer mpireq(2*ims_npro_i)
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_i; n = ip + 1
    call MPI_ISEND(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(1), ims_err)
    call MPI_IRECV(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_i
        ip = n - 1
        if (ip /= ims_pro_i) then
            l = l + 1
            call MPI_ISEND(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_i*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPF_I

!########################################################################
!########################################################################
subroutine TLabMPI_TRPB_K(b, a, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_k, ims_pro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: b
    real(KIND=8), dimension(*), intent(OUT) :: a
    integer(KIND=4), dimension(ims_npro_k), intent(IN) :: dsend, drecv
    integer, dimension(ims_npro_k), intent(IN) :: tsend, trecv

! -----------------------------------------------------------------------
    integer(KIND=4) n, l
    integer status(MPI_STATUS_SIZE, 2*ims_npro_k)
    integer mpireq(2*ims_npro_k)
    integer ip

#ifdef PROFILE_ON
    real(KIND=8) time_loc_1, time_loc_2
#endif

! #######################################################################
#ifdef PROFILE_ON
    time_loc_1 = MPI_WTIME()
#endif

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_k; n = ip + 1
    call MPI_ISEND(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(1), ims_err)
    call MPI_IRECV(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_k
        ip = n - 1
        if (ip /= ims_pro_k) then
            l = l + 1
            call MPI_ISEND(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_k*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TAGUPDT

    return
end subroutine TLabMPI_TRPB_K

!########################################################################
!########################################################################
subroutine TLabMPI_TRPB_I(b, a, dsend, drecv, tsend, trecv)

    use TLabMPI_VARS, only: ims_npro_i, ims_pro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err

    implicit none

#include "mpif.h"

    real(KIND=8), dimension(*), intent(IN) :: b
    real(KIND=8), dimension(*), intent(OUT) :: a
    integer(KIND=4), dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    integer, dimension(ims_npro_i), intent(IN) :: tsend, trecv ! types

! -----------------------------------------------------------------------
    integer(KIND=4) n, l
    integer status(MPI_STATUS_SIZE, 2*ims_npro_i)
    integer mpireq(2*ims_npro_i)
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_i; n = ip + 1
    call MPI_ISEND(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(1), ims_err)
    call MPI_IRECV(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_i
        ip = n - 1
        if (ip /= ims_pro_i) then
            l = l + 1
            call MPI_ISEND(b(drecv(n) + 1), 1, trecv(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(n) + 1), 1, tsend(n), ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_i*2 - 2, mpireq(3:), status(1, 3), ims_err)

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

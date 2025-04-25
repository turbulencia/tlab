!mpif90 -fpp  -nbs -save-temps -xHost -simd -vec-threshold50 -unroll-aggressive    -axcommon-avx512,SSE4.2  -qopt-prefetch -O3 vmpi_transpose.f90

! from dns_const.h
#define TREAL      REAL(8)    ! user-defined types
#define TINTEGER   INTEGER(4)

! from dns_const_mpi.h
#define TLAB_MPI_TRP_K_PARTIAL   1 ! tags and sizes for MPI data
#define TLAB_MPI_TRP_I_PARTIAL   1

#define TLAB_MPI_TRP_K_MAXTYPES  1
#define TLAB_MPI_TRP_I_MAXTYPES  1

module DNS_MPI
    implicit none
    save

    TINTEGER, parameter :: imax = 84   ! number of grid points per task
    TINTEGER, parameter :: jmax = 480
    TINTEGER, parameter :: kmax = 56

    integer, parameter :: ims_npro_i = 64 ! number of tasks in Ox and Oz (no decomposition along Oy)
    integer, parameter :: ims_npro_k = 96

    TINTEGER, parameter :: nmax = 20000  ! number of repetitions of operations

! Data below should not be changed
    integer :: ims_pro, ims_pro_i, ims_pro_j, ims_pro_k ! task positioning

    integer :: ims_comm_xz, ims_comm_x, ims_comm_z      ! communicators
    integer :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_err, ims_tag

    integer, dimension(:), allocatable :: ims_map_i
    TINTEGER, dimension(:), allocatable :: ims_size_i
    TINTEGER, dimension(:, :), allocatable :: ims_ds_i, ims_dr_i

    integer, dimension(:), allocatable :: ims_map_k
    TINTEGER, dimension(:), allocatable :: ims_size_k
    TINTEGER, dimension(:, :), allocatable :: ims_ds_k, ims_dr_k

    integer, dimension(:, :), allocatable :: status
    integer, dimension(:), allocatable :: mpireq

end module DNS_MPI

!########################################################################
! Main program to test forwards and backwards transposition
!########################################################################
program VMPI

    use TLabMPI_VARS

    implicit none

#include "mpif.h"

    TREAL, dimension(:, :), allocatable :: a
    TREAL, dimension(:), allocatable :: wrk3d

! -------------------------------------------------------------------
    TREAL residual                                      ! Control
    TINTEGER t_srt, t_end, t_dif, PROC_CYCLES, MAX_CYCLES ! Time
    TINTEGER n
    character*64 str
    character*256 line

    integer ims_npro
    TREAL dummy
    TINTEGER idummy, id

! ###################################################################
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

    if (ims_npro_i*ims_npro_k /= ims_npro) then ! check
        if (ims_pro == 0) then
            write (*, '(a)') 'Inconsistency in total number of PEs'
        end if
        call MPI_FINALIZE(ims_err)
        stop
    end if

    call TLabMPI_Initialize(ifile)
call TLabMPI_Trp_Initialize(ifile)

    allocate (a(imax*jmax*kmax, 18)) ! Number of 3d arrays commonly used in the code
    allocate (wrk3d(imax*jmax*kmax))

! ###################################################################
! ###################################################################
! Create random array
    call random_number(a(1:imax*jmax*kmax, 1))

    do n = 1, nmax

! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            id = TLAB_MPI_TRP_I_PARTIAL

            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)

            call TLabMPI_Trp_ExecI_Forward(a(1, 1), wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_size_i(id))
            call TLabMPI_Trp_ExecI_Backward(wrk3d, a(1, 2), ims_ds_i(1, id), ims_dr_i(1, id), ims_size_i(id))

            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

            idummy = t_end - t_srt
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            line = '. Max. elapsed time '//trim(adjustl(str))//' sec.'

            dummy = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') residual
            line = 'Checking MPI transposition for Ox derivatives: Residual '//trim(adjustl(str))//trim(adjustl(line))

            write (str, '(I)') n
            line = 'It '//trim(adjustl(str))//'. '//trim(adjustl(line))

            if (ims_pro == 0) then
                write (*, '(a)') trim(adjustl(line))
            end if

        end if

! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
        if (ims_npro_k > 1) then
            id = TLAB_MPI_TRP_K_PARTIAL

            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)

            call TLabMPI_Trp_ExecK_Forward(a(1, 1), wrk3d, ims_ds_k(1, id), ims_dr_k(1, id), ims_size_k(id))
            call TLabMPI_Trp_ExecK_Backward(wrk3d, a(1, 2), ims_ds_k(1, id), ims_dr_k(1, id), ims_size_k(id))

            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

            idummy = t_end - t_srt
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') real(t_dif)/PROC_CYCLES
            line = '. Max. elapsed time '//trim(adjustl(str))//' sec.'

            dummy = maxval(abs(a(1:imax*jmax*kmax, 1) - a(1:imax*jmax*kmax, 2)))
            call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            write (str, '(E13.5E3)') residual
            line = 'Checking MPI transposition for Oz derivatives: Residual '//trim(adjustl(str))//trim(adjustl(line))

            write (str, '(I)') n
            line = 'It '//trim(adjustl(str))//'. '//trim(adjustl(line))

            if (ims_pro == 0) then
                write (*, '(a)') trim(adjustl(line))
            end if

        end if

    end do

    call MPI_FINALIZE(ims_err)

end program VMPI

! #######################################################################
! Rest of routines
! #######################################################################
subroutine TLabMPI_Initialize()

    use TLabMPI_VARS

    implicit none

#include "mpif.h"

! -----------------------------------------------------------------------
    TINTEGER id, ip, npage
    TINTEGER i1, dims(2)
    logical period(2), remain_dims(2), reorder

! #######################################################################
    allocate (ims_map_i(ims_npro_i))
    allocate (ims_size_i(TLAB_MPI_TRP_I_MAXTYPES))
    allocate (ims_ds_i(ims_npro_i, TLAB_MPI_TRP_I_MAXTYPES))
    allocate (ims_dr_i(ims_npro_i, TLAB_MPI_TRP_I_MAXTYPES))

    allocate (ims_map_k(ims_npro_k))
    allocate (ims_size_k(TLAB_MPI_TRP_K_MAXTYPES))
    allocate (ims_ds_k(ims_npro_k, TLAB_MPI_TRP_K_MAXTYPES))
    allocate (ims_dr_k(ims_npro_k, TLAB_MPI_TRP_K_MAXTYPES))

    allocate (status(MPI_STATUS_SIZE, 2*max(ims_npro_k, ims_npro_i)))
    allocate (mpireq(2*max(ims_npro_k, ims_npro_i)))

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
    if (ims_npro_i > 1) then
        id = TLAB_MPI_TRP_I_PARTIAL
        ! Calculate size
        ims_size_i(id) = imax*jmax*kmax/ims_npro_i
        ! Calculate Displacements in Forward Send/Receive
        ims_ds_i(1, id) = 0
        ims_dr_i(1, id) = 0
        do ip = 2, ims_npro_i
            ims_ds_i(ip, id) = ims_ds_i(ip - 1, id) + ims_size_i(id)
            ims_dr_i(ip, id) = ims_dr_i(ip - 1, id) + ims_size_i(id)
        end do
    end if

    if (ims_npro_k > 1) then
        id = TLAB_MPI_TRP_K_PARTIAL
        ! Calculate size
        ims_size_k(id) = imax*jmax*kmax/ims_npro_k
        ! Calculate Displacements in Forward Send/Receive
        ims_ds_k(1, id) = 0
        ims_dr_k(1, id) = 0
        do ip = 2, ims_npro_k
            ims_ds_k(ip, id) = ims_ds_k(ip - 1, id) + ims_size_k(id)
            ims_dr_k(ip, id) = ims_dr_k(ip - 1, id) + ims_size_k(id)
        end do
    end if

    call TLabMPI_TagReset

    return
end subroutine TLabMPI_Initialize

! ###################################################################
! ###################################################################
subroutine TLabMPI_Trp_ExecK_Forward(a, b, dsend, drecv, size)

    use TLabMPI_VARS, only: ims_npro_k, ims_pro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: status, mpireq

    implicit none

#include "mpif.h"

    TREAL, dimension(*), intent(IN) :: a
    TREAL, dimension(*), intent(OUT) :: b
    TINTEGER, dimension(ims_npro_k), intent(IN) :: dsend, drecv ! displacements
    TINTEGER, intent(IN) :: size

! -----------------------------------------------------------------------
    TINTEGER n, l
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_k; n = ip + 1
    call MPI_ISEND(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(1), ims_err)
    call MPI_IRECV(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_k
        ip = n - 1
        if (ip /= ims_pro_k) then
            l = l + 1
            call MPI_ISEND(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_k*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TagUpdate

    return
end subroutine TLabMPI_Trp_ExecK_Forward

!########################################################################
!########################################################################
subroutine TLabMPI_Trp_ExecI_Forward(a, b, dsend, drecv, size)

    use TLabMPI_VARS, only: ims_npro_i, ims_pro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: status, mpireq

    implicit none

#include "mpif.h"

    TREAL, dimension(*), intent(IN) :: a
    TREAL, dimension(*), intent(OUT) :: b
    TINTEGER, dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    TINTEGER, intent(IN) :: size

! -----------------------------------------------------------------------
    TINTEGER n, l
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_i; n = ip + 1
    call MPI_ISEND(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(1), ims_err)
    call MPI_IRECV(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_i
        ip = n - 1
        if (ip /= ims_pro_i) then
            l = l + 1
            call MPI_ISEND(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_i*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TagUpdate

    return
end subroutine TLabMPI_Trp_ExecI_Forward

!########################################################################
!########################################################################
subroutine TLabMPI_Trp_ExecK_Backward(b, a, dsend, drecv, size)

    use TLabMPI_VARS, only: ims_npro_k, ims_pro_k
    use TLabMPI_VARS, only: ims_comm_z
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: status, mpireq

    implicit none

#include "mpif.h"

    TREAL, dimension(*), intent(IN) :: b
    TREAL, dimension(*), intent(OUT) :: a
    TINTEGER, dimension(ims_npro_k), intent(IN) :: dsend, drecv
    TINTEGER, intent(IN) :: size

! -----------------------------------------------------------------------
    TINTEGER n, l
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_k; n = ip + 1
    call MPI_ISEND(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(1), ims_err)
    call MPI_IRECV(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_k
        ip = n - 1
        if (ip /= ims_pro_k) then
            l = l + 1
            call MPI_ISEND(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_z, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_k*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TagUpdate

    return
end subroutine TLabMPI_Trp_ExecK_Backward

!########################################################################
!########################################################################
subroutine TLabMPI_Trp_ExecI_Backward(b, a, dsend, drecv, size)

    use TLabMPI_VARS, only: ims_npro_i, ims_pro_i
    use TLabMPI_VARS, only: ims_comm_x
    use TLabMPI_VARS, only: ims_tag, ims_err
    use TLabMPI_VARS, only: status, mpireq

    implicit none

#include "mpif.h"

    TREAL, dimension(*), intent(IN) :: b
    TREAL, dimension(*), intent(OUT) :: a
    TINTEGER, dimension(ims_npro_i), intent(IN) :: dsend, drecv ! displacements
    TINTEGER, intent(IN) :: size

! -----------------------------------------------------------------------
    TINTEGER n, l
    integer ip

! #######################################################################
! Same processor
! #######################################################################
    ip = ims_pro_i; n = ip + 1
    call MPI_ISEND(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(1), ims_err)
    call MPI_IRECV(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(2), ims_err)

    call MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
    l = 2
    do n = 1, ims_npro_i
        ip = n - 1
        if (ip /= ims_pro_i) then
            l = l + 1
            call MPI_ISEND(b(drecv(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
            l = l + 1
            call MPI_IRECV(a(dsend(n) + 1), size, MPI_REAL8, ip, ims_tag, ims_comm_x, mpireq(l), ims_err)
        end if
    end do

    call MPI_WAITALL(ims_npro_i*2 - 2, mpireq(3:), status(1, 3), ims_err)

    call TLabMPI_TagUpdate

    return
end subroutine TLabMPI_Trp_ExecI_Backward

!########################################################################
!########################################################################
subroutine TLabMPI_TagUpdate

    use TLabMPI_VARS, only: ims_tag

    implicit none

    ims_tag = ims_tag + 1

    if (ims_tag > 32000) then
        call TLabMPI_TagReset
    end if

    return
end subroutine TLabMPI_TagUpdate

!########################################################################
!########################################################################
subroutine TLabMPI_TagReset

    use TLabMPI_VARS, only: ims_tag

    implicit none

    ims_tag = 0

    return
end subroutine TLabMPI_TagReset

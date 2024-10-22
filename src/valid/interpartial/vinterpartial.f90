#include "types.h"
#include "dns_const.h"
!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/01/21 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact interpolation schemes.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
program VINTERPARTIAL

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLab_WorkFlow
    use IO_FIELDS
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
    use TLabMPI_VARS
#endif
    use OPR_PARTIAL

    implicit none

#ifdef USE_MPI
    TREAL error2, dummy2
#else
    TINTEGER, parameter :: ims_pro = 0
#endif

    TREAL, dimension(:, :), allocatable, save, target :: x, y, z
    TREAL, dimension(:, :, :), allocatable :: a, a_int, a_dif
    TREAL, dimension(:, :, :), allocatable :: b, c
    TREAL, dimension(:, :), allocatable :: wrk1d, wrk2d
    TREAL, dimension(:), allocatable :: wrk3d, tmp1, d

    TINTEGER bcs(2, 2)
    TREAL dummy, error
! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL('tlab.ini')
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif

! -------------------------------------------------------------------
! Check input
! -------------------------------------------------------------------
    if (.not. stagger_on) then
        call TLAB_WRITE_ASCII(efile, 'VINTERPARTIAL. Set "StaggerGrid=yes" in tlab.ini!')
        call TLAB_STOP(0)
    end if

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))
    allocate (wrk1d(isize_wrk1d, inb_wrk1d + 1))
    allocate (wrk2d(isize_wrk2d, inb_wrk2d))
    allocate (a(imax, jmax, kmax), a_int(imax, jmax, kmax), a_dif(imax, jmax, kmax))
    allocate (b(imax, jmax, kmax), c(imax, jmax, kmax), d(imax*jmax*kmax))
    allocate (tmp1(isize_txc_field), wrk3d(isize_wrk3d))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    bcs = 0
! ###################################################################
! Define forcing term
! ###################################################################
    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)

! ###################################################################
! x-direction: Interpolation + interpolatory 1st derivative
! ###################################################################
! Interpolation: vel. <--> pre. grid
    call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), a, a_int)
    call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), a_int, b)
! Difference field and error
    a_dif = a - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) '=========================================================='
        write (*, *) '------ Check interpolatory operators in x-direction ------'
        write (*, *) '=========================================================='
        write (*, *) 'Interpolation, vel. <--> pre. grid '
        write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
    end if
    ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
! -------------------------------------------------------------------
! 1st interp. deriv + Interpolation: vel. <--> pre. grid
    call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), a, a_int)
    call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), a_int, b)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, c)
! Difference field and error
    a_dif = c - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Interpolation + interp. 1st derivative, vel. --> pre. grid'
        write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
    end if
    ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
! -------------------------------------------------------------------
! 1st interp. deriv + Interpolation: vel. <--> pre. grid
    call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), a, a_int)
    call OPR_PARTIAL_X(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(1), a_int, b)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, c)
! Difference field and error
    a_dif = c - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Interpolation + interp. 1st derivative, pre. --> vel. grid'
        write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
    end if
    ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
! ###################################################################
! z-direction: Interpolation + interpolatory 1st derivative
! ###################################################################
    if (g(3)%size > 1) then
        ! Interpolation: vel. <--> pre. grid
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), a, a_int)
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), a_int, b)
        ! Difference field and error
        a_dif = a - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) '  '
            write (*, *) '=========================================================='
            write (*, *) '------ Check interpolatory operators in z-direction ------'
            write (*, *) '=========================================================='
            write (*, *) 'Interpolation, vel. <--> pre. grid '
            write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
        end if
        ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
        ! -------------------------------------------------------------------
        ! 1st interp. deriv + Interpolation: vel. <--> pre. grid
        call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), a, a_int)
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), a_int, b)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, c)
        ! Difference field and error
        a_dif = c - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) 'Interpolation + interp. 1st derivative, vel. --> pre. grid'
            write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
        end if
        ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
        ! -------------------------------------------------------------------
        ! 1st interp. deriv + Interpolation: vel. <--> pre. grid
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), a, a_int)
        call OPR_PARTIAL_Z(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(3), a_int, b)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, c)
        ! Difference field and error
        a_dif = c - b; error = sum(a_dif**2); dummy = sum(a**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) 'Interpolation + interp. 1st derivative, pre. --> vel. grid'
            write (*, *) 'Relative error ...............: ', sqrt(error)/sqrt(dummy)
        end if
        ! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, 1, a_dif)
    end if

! ###################################################################

    call TLAB_STOP(0)
end program VINTERPARTIAL

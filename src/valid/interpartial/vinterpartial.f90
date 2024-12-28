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
    use TLab_Constants, only: wp, wi, ifile, gfile, efile
    use TLAB_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use IO_FIELDS
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: TLabMPI_Initialize
use TLabMPI_PROCS, only: TLabMPI_Transpose_Initialize
    use TLabMPI_VARS
#endif
    use FDM, only: g, FDM_Initialize
    use OPR_PARTIAL

    implicit none

#ifdef USE_MPI
    real(wp) error2, dummy2
#else
    integer(wi), parameter :: ims_pro = 0
#endif

    real(wp), dimension(:, :), allocatable, save, target :: x, y, z
    real(wp), dimension(:, :, :), allocatable :: a, a_int, a_dif
    real(wp), dimension(:, :, :), allocatable :: b, c
    real(wp), dimension(:, :), allocatable :: wrk1d, wrk2d
    real(wp), dimension(:), allocatable :: wrk3d, tmp1, d

    integer(wi) bcs(2, 2)
    real(wp) dummy, error
! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters('tlab.ini')
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
call TLabMPI_Transpose_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)

! -------------------------------------------------------------------
! Check input
! -------------------------------------------------------------------
    if (.not. stagger_on) then
        call TLab_Write_ASCII(efile, 'VINTERPARTIAL. Set "StaggerGrid=yes" in tlab.ini!')
        call TLab_Stop(0)
    end if

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d, inb_wrk1d + 1))
    allocate (wrk2d(isize_wrk2d, inb_wrk2d))
    allocate (a(imax, jmax, kmax), a_int(imax, jmax, kmax), a_dif(imax, jmax, kmax))
    allocate (b(imax, jmax, kmax), c(imax, jmax, kmax), d(imax*jmax*kmax))
    allocate (tmp1(isize_txc_field), wrk3d(isize_wrk3d))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_Initialize(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_Initialize(z, g(3), wrk1d(:,3), wrk1d(:,4))

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

    call TLab_Stop(0)
end program VINTERPARTIAL

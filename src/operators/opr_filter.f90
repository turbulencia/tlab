#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module OPR_FILTERS
    use TLAB_CONSTANTS, only: wp, wi, efile, MAX_PROF
    use TLAB_TYPES, only: grid_dt, filter_dt
    use TLAB_VARS, only: isize_txc_field, isize_txc_dimz, g
    use TLAB_PROCS
    use FLT_COMPACT
    use FLT_EXPLICIT
#ifdef USE_MPI
    use TLAB_MPI_VARS
    use TLAB_MPI_PROCS
#endif
    implicit none
    private

    public :: OPR_FILTER_INITIALIZE
    public :: OPR_FILTER
    public :: OPR_FILTER_1D

contains

    !###################################################################
    !###################################################################
    subroutine OPR_FILTER_INITIALIZE(g, f, wrk1d)
        type(grid_dt),   intent(in)    :: g
        type(filter_dt), intent(inout) :: f
        real(wp),        intent(inout) :: wrk1d(*)

        !###################################################################
        if (f%inb_filter > 0) allocate (f%coeffs(f%size, f%inb_filter))

        select case (f%type)

        case (DNS_FILTER_4E, DNS_FILTER_ADM)
            call FLT_E4_COEFFS(f%size, f%periodic, g%scale, g%nodes, f%coeffs)

        case (DNS_FILTER_TOPHAT)
            call FLT_T1_INI(g%scale, g%nodes, f, wrk1d)

        case (DNS_FILTER_COMPACT)
            call FLT_C4_LHS(f%size, f%bcsmin, f%bcsmax, f%parameters(1), f%coeffs(1, 6), f%coeffs(1, 7), f%coeffs(1, 8))
            if (f%periodic) then
                call TRIDPFS(f%size, f%coeffs(1, 6), f%coeffs(1, 7), f%coeffs(1, 8), f%coeffs(1, 9), f%coeffs(1, 10))
            else
                call TRIDFS(f%size, f%coeffs(1, 6), f%coeffs(1, 7), f%coeffs(1, 8))
            end if
            call FLT_C4_RHS_COEFFS(f%size, f%parameters(1), f%periodic, g%jac, f%coeffs(1, 1))

        case (DNS_FILTER_COMPACT_CUTOFF)
            if (f%periodic) then
                call FLT_C4P_CUTOFF_LHS(f%size, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), f%coeffs(1, 4), f%coeffs(1, 5))
                call PENTADPFS(f%size, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), &
                               f%coeffs(1, 4), f%coeffs(1, 5), f%coeffs(1, 6), f%coeffs(1, 7))
            else
                call FLT_C4_CUTOFF_LHS(f%size, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), f%coeffs(1, 4), f%coeffs(1, 5))
                call PENTADFS2(f%size, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), &
                               f%coeffs(1, 4), f%coeffs(1, 5))
            end if

        case (DNS_FILTER_HELMHOLTZ)
            f%parameters(2) = -1.0_wp/(f%parameters(1))**2

        end select

        return
    end subroutine OPR_FILTER_INITIALIZE

    !###################################################################
    !###################################################################
    subroutine OPR_FILTER(nx, ny, nz, f, u, wrk1d, wrk2d, txc)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

        integer(wi),     intent(in)    :: nx, ny, nz
        type(filter_dt), intent(in)    :: f(3)
        real(wp),        intent(inout) :: u(nx, ny, nz)    !Inplace operation
        real(wp),        intent(inout) :: wrk1d(ny, *)
        real(wp),        intent(inout) :: wrk2d(nx, nz, 2)
        real(wp),        intent(inout) :: txc(isize_txc_field, *)  !size 2 if ADM
        !size 4 if SPECTRAL, HELMHOLTZ
        !-------------------------------------------------------------------
        real(wp) dummy
        integer(wi) k, flag_bcs, n, bcs(2, 2), nxy, ip_b, ip_t, i2

        complex(wp), pointer :: s_wrk(:, :) => null()
        target txc

        !###################################################################
        nxy = nx*ny

        bcs = 0  !Boundary conditions for derivative operator set to biased, non-zero
        i2 = 2

        !Global filters
        select case (f(1)%type)

        case (DNS_FILTER_HELMHOLTZ)
            if (f(2)%BcsMin == DNS_FILTER_BCS_DIRICHLET) then
                do k = 1, nz
                    wrk2d(:, k, 1) = u(:, 1, k)
                    wrk2d(:, k, 2) = u(:, ny, k)
                end do
                flag_bcs = 0
            else if (f(2)%BcsMin == DNS_FILTER_BCS_NEUMANN) then
                call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, txc(1, 1), txc(1, 2), wrk2d, txc(1, 2))
                ip_b = 1
                ip_t = nx*(ny - 1) + 1
                do k = 1, nz
                    wrk2d(1:nx, k, 1) = txc(ip_b:ip_b + nx - 1, 1); ip_b = ip_b + nxy  !bottom
                    wrk2d(1:nx, k, 2) = txc(ip_t:ip_t + nx - 1, 1); ip_t = ip_t + nxy  !top
                end do
                flag_bcs = 3
            else if (f(2)%BcsMin == DNS_FILTER_BCS_SOLID) then
                wrk2d(:, :, 1) = 0.0_wp
                wrk2d(:, :, 2) = 0.0_wp
                flag_bcs = 3
            end if

            txc(1:nx*ny*nz, 1) = u(1:nx*ny*nz, 1, 1)*f(1)%parameters(2)  !I need extended arrays
            call OPR_HELMHOLTZ_FXZ(nx, ny, nz, g, flag_bcs, f(1)%parameters(2), &
                                   txc(1, 1), txc(1, 2), txc(1, 3), &
                                   wrk2d(1, 1, 1), wrk2d(1, 1, 2), wrk1d, wrk1d(1, 5), txc(1, 4))
            u(1:nx*ny*nz, 1, 1) = txc(1:nx*ny*nz, 1)

        case (DNS_FILTER_BAND)
            dummy = 1.0_wp/real(f(1)%size*f(3)%size, wp)
            txc(1:nx*ny*nz, 1) = u(1:nx*ny*nz, 1, 1)  !I need extended arrays
            call OPR_FOURIER_F(i2, nx, ny, nz, txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, txc(1, 4))
            call c_f_pointer(c_loc(txc(1, 2)), s_wrk, shape=[isize_txc_dimz/2, nz])
            call OPR_FILTER_BAND_2D(nx, ny, nz, f(1)%parameters, s_wrk)
            call OPR_FOURIER_B(i2, nx, ny, nz, txc(1, 2), txc(1, 3), txc(1, 4))
            u(1:nx*ny*nz, 1, 1) = txc(1:nx*ny*nz, 3)*dummy

        case (DNS_FILTER_ERF)
            dummy = 1.0_wp/real(f(1)%size*f(3)%size, wp)
            txc(1:nx*ny*nz, 1) = u(1:nx*ny*nz, 1, 1)  !I need extended arrays
            call OPR_FOURIER_F(i2, nx, ny, nz, txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, txc(1, 4))
            call c_f_pointer(c_loc(txc(1, 2)), s_wrk, shape=[isize_txc_dimz/2, nz])
            call OPR_FILTER_ERF_2D(nx, ny, nz, f(1)%parameters, s_wrk)
            call OPR_FOURIER_B(i2, nx, ny, nz, txc(1, 2), txc(1, 3), txc(1, 4))
            u(1:nx*ny*nz, 1, 1) = txc(1:nx*ny*nz, 3)*dummy

        case default

            !###################################################################
            !Directional filters
            if (f(1)%type /= DNS_FILTER_NONE) then
                do n = 1, f(1)%repeat
                    call OPR_FILTER_X(nx, ny, nz, f(1), u, txc(1, 2), wrk2d, txc(1, 1))
                end do
            end if

            if (f(2)%type /= DNS_FILTER_NONE) then
                do n = 1, f(2)%repeat
                    call OPR_FILTER_Y(nx, ny, nz, f(2), u, txc(1, 2), wrk2d, txc(1, 1))
                end do
            end if

            if (f(3)%type /= DNS_FILTER_NONE) then
                do n = 1, f(3)%repeat
                    call OPR_FILTER_Z(nx, ny, nz, f(3), u, txc(1, 2), wrk2d, txc(1, 1))
                end do
            end if

        end select

        return
    end subroutine OPR_FILTER

    !###################################################################
    !Filter kernel along one direction
    !###################################################################
    subroutine OPR_FILTER_1D(nlines, f, u, result, wrk2d, wrk3d)

        integer(wi),     intent(in)    :: nlines                  !# of lines to be solved
        type(filter_dt), intent(in)    :: f
        real(wp),        intent(in)    :: u(nlines, f%size)       !field to be filtered
        real(wp),        intent(out)   :: result(nlines, f%size)  !filtered filed
        real(wp),        intent(inout) :: wrk2d(*), wrk3d(*)

        !-------------------------------------------------------------------
        integer(wi) delta

        !###################################################################
        delta = int(f%parameters(1))

        select case (f%type)

        case (DNS_FILTER_COMPACT)
            call FLT_C4_RHS(f%size, nlines, f%periodic, f%bcsmin, f%bcsmax, f%coeffs, u, result)
            if (f%periodic) then
        call TRIDPSS(f%size, nlines, f%coeffs(1, 6), f%coeffs(1, 7), f%coeffs(1, 8), f%coeffs(1, 9), f%coeffs(1, 10), result, wrk2d)
            else
                call TRIDSS(f%size, nlines, f%coeffs(1, 6), f%coeffs(1, 7), f%coeffs(1, 8), result)
            end if

        case (DNS_FILTER_COMPACT_CUTOFF)
            if (f%periodic) then
                call FLT_C4P_CUTOFF_RHS(f%size, nlines, u, result)
                call PENTADPSS(f%size, nlines, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), &
                               f%coeffs(1, 4), f%coeffs(1, 5), f%coeffs(1, 6), f%coeffs(1, 7), result)
            else
                call FLT_C4_CUTOFF_RHS(f%size, nlines, u, result)
                call PENTADSS2(f%size, nlines, f%coeffs(1, 1), f%coeffs(1, 2), f%coeffs(1, 3), &
                               f%coeffs(1, 4), f%coeffs(1, 5), result)
            end if

        case (DNS_FILTER_6E)
            call FLT_E6(f%size, nlines, f%periodic, f%bcsmin, f%bcsmax, u, result)

        case (DNS_FILTER_4E)
            call FLT_E4(f%size, nlines, f%periodic, f%coeffs, u, result)

        case (DNS_FILTER_ADM)
            call FLT_ADM(f%size, nlines, f%periodic, f%coeffs, u, result, wrk3d)

        case (DNS_FILTER_TOPHAT)
            if (f%periodic) then
                if (f%uniform) then
                    if (delta == 2) then; call FLT_T1PD2(f%size, nlines, u, result)
                    else if (delta == 4) then; call FLT_T1PD4(f%size, nlines, u, result)
                    else; call FLT_T1P(f%size, nlines, delta, u, result)
                    end if
                else
                    call FLT_T1P_ND(f%size, nlines, delta, f%coeffs, u, result)
                end if
            else
                if (f%uniform) then
                    call FLT_T1(f%size, nlines, delta, f%coeffs, u, result)
                else
                    if (delta == 2) then; call FLT_T1NDD2(f%size, nlines, f%coeffs, u, result)
                    else if (delta == 4) then; call FLT_T1NDD4(f%size, nlines, f%coeffs, u, result)
                    else if (delta == 6) then; call FLT_T1NDD6(f%size, nlines, f%coeffs, u, result)
                    else; call FLT_T1ND(f%size, nlines, delta, f%coeffs, u, result)
                    end if
                end if
            end if

        end select

        return
    end subroutine OPR_FILTER_1D

    !###################################################################
    !Filter in Ox direction
    !###################################################################
    subroutine OPR_FILTER_X(nx, ny, nz, f, u, tmp, wrk2d, wrk3d)
        integer(wi),     intent(in)    :: nx, ny, nz
        type(filter_dt), intent(in)    :: f
        real(wp),        intent(inout) :: u(nx*ny*nz), wrk3d(nx*ny*nz)     !in-place operation
        real(wp),        intent(inout) :: wrk2d(ny*nz)        !Aux arrays
        real(wp),        intent(inout) :: tmp(nx*ny*nz)          !Aux array needed in ADM type

        !-------------------------------------------------------------------
        integer(wi) nyz
        real(wp), dimension(:), pointer :: p_a, p_b
        target u, wrk3d

#ifdef USE_MPI
        integer(wi) id
#endif

        !###################################################################
#ifdef USE_MPI
        id = f%mpitype  !TLAB_MPI_I_PARTIAL
#endif

        !-------------------------------------------------------------------
        !Transposition
        !-------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLAB_MPI_TRPF_I(u, wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
            p_a => wrk3d
            p_b => u
            nyz = ims_size_i(id)
        else
#endif
            p_a => u
            p_b => wrk3d
            nyz = ny*nz
#ifdef USE_MPI
        end if
#endif

        !-------------------------------------------------------------------
        !Make  x  direction the last one
        !-------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(p_a, f%size, f%size, nyz, p_b, nyz)
#else
        call DNS_TRANSPOSE(p_a, f%size, nyz, f%size, p_b, nyz)
#endif

        !###################################################################
        call OPR_FILTER_1D(nyz, f, p_b, p_a, wrk2d, tmp)

        !###################################################################
        !-------------------------------------------------------------------
        !Put arrays back in the order in which they came in
        !-------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(p_a, nyz, nyz, f%size, p_b, f%size)
#else
        call DNS_TRANSPOSE(p_a, nyz, f%size, nyz, p_b, f%size)
#endif

        !-------------------------------------------------------------------
        !Transposition
        !-------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLAB_MPI_TRPB_I(p_b, p_a, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        end if
#endif

        u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

        nullify (p_a, p_b)

        return
    end subroutine OPR_FILTER_X

    !###################################################################
    !Filter in Oy direction
    !###################################################################
    subroutine OPR_FILTER_Y(nx, ny, nz, f, u, tmp, wrk2d, wrk3d)

        integer(wi),     intent(in)    :: nx, ny, nz
        type(filter_dt), intent(in)    :: f
        real(wp),        intent(inout) :: u(nx*ny*nz), wrk3d(nx*ny*nz)  !in-place operation
        real(wp),        intent(inout) :: wrk2d(*)                      !Aux arrays
        real(wp),        intent(inout) :: tmp(nx*ny*nz)                 !Aux array needed in ADM type

        !-----------------------------------------------------------------------
        integer(wi) nxy, nxz
        real(wp), pointer :: p_org(:), p_dst(:)

        target u, wrk3d

        !#######################################################################
        nxy = nx*ny
        nxz = nx*nz

        !-------------------------------------------------------------------
        !Make y-direction the last one
        !-------------------------------------------------------------------
        if (nz == 1) then
            p_org => u
            p_dst => wrk3d
        else
#ifdef USE_ESSL
            call DGETMO(u, nxy, nxy, nz, wrk3d, nz)
#else
            call DNS_TRANSPOSE(u, nxy, nz, nxy, wrk3d, nz)
#endif
            p_org => wrk3d
            p_dst => u
        end if

        !###################################################################
        call OPR_FILTER_1D(nxz, f, p_org, p_dst, wrk2d, tmp)

        !-------------------------------------------------------------------
        !Put arrays back in the order in which they came in
        !-------------------------------------------------------------------
        if (nz > 1) then
#ifdef USE_ESSL
            call DGETMO(p_dst, nz, nz, nxy, p_org, nxy)
#else
            call DNS_TRANSPOSE(p_dst, nz, nxy, nz, p_org, nxy)
#endif
        end if

        u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

        nullify (p_org, p_dst)

        return
    end subroutine OPR_FILTER_Y

    !###################################################################
    !Filter in Oz direction
    !###################################################################
    subroutine OPR_FILTER_Z(nx, ny, nz, f, u, tmp, wrk2d, wrk3d)
        integer(wi),     intent(in)    :: nx, ny, nz
        type(filter_dt), intent(in)    :: f
        real(wp),        intent(inout) :: u(nx*ny*nz), wrk3d(nx*ny*nz)  !in-place operation
        real(wp),        intent(inout) :: wrk2d(nx*ny)                  !Aux arrays
        real(wp),        intent(inout) :: tmp(nx*ny*nz)                 !Aux array needed in ADM type

        !-------------------------------------------------------------------
        integer(wi) nxy
        real(wp), pointer :: p_a(:), p_b(:)
        target u, wrk3d

#ifdef USE_MPI
        integer(wi) id
#endif

        !###################################################################
#ifdef USE_MPI
        id = f%mpitype  !TLAB_MPI_K_PARTIAL
#endif

        !-------------------------------------------------------------------
        !Transposition
        !-------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPF_K(u, wrk3d, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            p_a => wrk3d
            p_b => u
            nxy = ims_size_k(id)
        else
#endif
            p_a => u
            p_b => wrk3d
            nxy = nx*ny
#ifdef USE_MPI
        end if
#endif

        !###################################################################
        call OPR_FILTER_1D(nxy, f, p_a, p_b, wrk2d, tmp)

        !###################################################################
        !-------------------------------------------------------------------
        !Transposition
        !-------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPB_K(p_b, p_a, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        end if
#endif

        u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

        nullify (p_a, p_b)

        return
    end subroutine OPR_FILTER_Z

    !########################################################################
    !Spectral band filter
    !########################################################################
    subroutine OPR_FILTER_BAND_2D(nx, ny, nz, spc_param, a)
        integer(wi), intent(in)    :: nx, ny, nz
        real(wp),    intent(in)    :: spc_param(*)
        complex(wp), intent(inout) :: a(isize_txc_dimz/2, nz)

        !-----------------------------------------------------------------------
        integer(wi) i, j, k, iglobal, kglobal, ip
        real(wp) fi, fk, f

        !#######################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= g(3)%size/2 + 1) then; fk = real(kglobal - 1, wp)/g(3)%scale
            else; fk = -real(g(3)%size + 1 - kglobal, wp)/g(3)%scale; end if

            do i = 1, nx/2 + 1
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                fi = real(iglobal - 1, wp)/g(1)%scale

                f = sqrt(fi**2 + fk**2)

                !apply spectral cutoff
                do j = 1, ny
                    ip = (j - 1)*(nx/2 + 1) + i
                    if ((f - spc_param(1))*(spc_param(2) - f) < 0.0_wp) a(ip, k) = 0.0_wp
                end do

            end do
        end do

        return
    end subroutine OPR_FILTER_BAND_2D

    !########################################################################
    !#
    !# Spectral filter with smooth (error-function) transition.
    !# The error function filter-response function is imposed
    !# in logarithmic wavenumber space.
    !#
    !########################################################################
    !# ARGUMENTS
    !#
    !#    spc_param(1) physical frequency for transition
    !#                 > 0: High-pass
    !#                 < 0: Low-pass
    !#    spc_param(2) width of transition in logarithmic wavenumber space
    !#    spc_param(3) normalise wavenumers in z-direction
    !#
    !########################################################################
    subroutine OPR_FILTER_ERF_2D(nx, ny, nz, spc_param, a)
        integer(wi), intent(in)    :: nx, ny, nz
        real(wp),    intent(in)    :: spc_param(*)
        complex(wp), intent(inout) :: a(isize_txc_dimz/2, nz)

        !-----------------------------------------------------------------------
        integer(wi) i, j, k, iglobal, kglobal, ip
        integer(wi) sign_pass, off_pass
        real(wp) fi, fk, f, fcut_log, damp

        !#######################################################################
        if (spc_param(1) > 0) then; 
            sign_pass = 1.    !HIGHPASS
            off_pass = 0.
        else                 !spc_param(1) <= 0
            sign_pass = -1.    !LOWPASS
            off_pass = 1.
        end if

        fcut_log = log(abs(spc_param(1)))
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= g(3)%size/2 + 1) then; fk = real(kglobal - 1, wp)/g(3)%scale/spc_param(3)
            else; fk = -real(g(3)%size + 1 - kglobal, wp)/g(3)%scale/spc_param(3); end if

            do i = 1, nx/2 + 1
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                fi = real(iglobal - 1, wp)/g(1)%scale

                f = sqrt(fi**2 + fk**2)
                if (f > 0) then; damp = (erf((log(f) - fcut_log)/spc_param(2)) + 1.)/2.
                else; damp = 0.0_wp; 
                end if

                !Set to high- or low-pass
                !high-pass: damp = 0.0 + damp
                !low-pass:  damp = 1.0 - damp
                damp = off_pass + sign_pass*damp

                !apply filter
                do j = 1, ny
                    ip = (j - 1)*(nx/2 + 1) + i
                    a(ip, k) = damp*a(ip, k)
                end do
            end do
        end do

        return
    end subroutine OPR_FILTER_ERF_2D

end module OPR_FILTERS

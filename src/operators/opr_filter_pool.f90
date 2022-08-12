#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

! ###################################################################
! ###################################################################
subroutine OPR_FILTER_INITIALIZE(g, f, wrk1d)
    use TLAB_TYPES, only: grid_dt, filter_dt
    use FLT_COMPACT
    use FLT_EXPLICIT

    implicit none

    type(grid_dt), intent(IN) :: g
    type(filter_dt), intent(INOUT) :: f
    TREAL, dimension(*), intent(INOUT) :: wrk1d

! -------------------------------------------------------------------

! ###################################################################
    if (f%inb_filter > 0) &
        allocate (f%coeffs(f%size, f%inb_filter))

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
        call FLT_C4_RHS_COEFFS(f%size, f%parameters(1), f%periodic, g%jac, f%coeffs(1,1))

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
        f%parameters(2) = -C_1_R/(f%parameters(1))**2

    end select

    return
end subroutine OPR_FILTER_INITIALIZE

! ###################################################################
! Filter kernel along one direction
! ###################################################################
subroutine OPR_FILTER_1D(nlines, f, u, result, wrk1d, wrk2d, wrk3d)
    use TLAB_TYPES, only: filter_dt
    use FLT_COMPACT
    use FLT_EXPLICIT

    implicit none

    TINTEGER, intent(IN) :: nlines     ! # of lines to be solved
    type(filter_dt), intent(IN) :: f
    TREAL, dimension(nlines, f%size), intent(IN) :: u          ! field to be filtered
    TREAL, dimension(nlines, f%size), intent(OUT) :: result     ! filtered filed
    TREAL, dimension(f%size, 5), intent(INOUT) :: wrk1d
    TREAL, dimension(*), intent(INOUT) :: wrk2d, wrk3d

! -------------------------------------------------------------------
    TINTEGER delta

! ###################################################################
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

! ###################################################################
! Filter in Ox direction
! ###################################################################
subroutine OPR_FILTER_X(nx, ny, nz, f, u, tmp, wrk1d, wrk2d, wrk3d)

    use TLAB_TYPES, only: filter_dt
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_i
    use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_PROCS
#endif

    implicit none

    TINTEGER, intent(IN) :: nx, ny, nz
    type(filter_dt), intent(IN) :: f
    TREAL, dimension(nx*ny*nz), intent(INOUT), target :: u, wrk3d    ! in-place operation
    TREAL, dimension(ny*nz), intent(INOUT) :: wrk2d       ! Aux arrays
    TREAL, dimension(f%size, *), intent(INOUT) :: wrk1d
    TREAL, dimension(nx*ny*nz), intent(INOUT) :: tmp         ! Aux array needed in ADM type

! -------------------------------------------------------------------
    TINTEGER nyz

    TREAL, dimension(:), pointer :: p_a, p_b

#ifdef USE_MPI
    TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
    id = f%mpitype ! TLAB_MPI_I_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
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

! -------------------------------------------------------------------
! Make  x  direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
    call DGETMO(p_a, f%size, f%size, nyz, p_b, nyz)
#else
    call DNS_TRANSPOSE(p_a, f%size, nyz, f%size, p_b, nyz)
#endif

! ###################################################################
    call OPR_FILTER_1D(nyz, f, p_b, p_a, wrk1d, wrk2d, tmp)

! ###################################################################
! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
#ifdef USE_ESSL
    call DGETMO(p_a, nyz, nyz, f%size, p_b, f%size)
#else
    call DNS_TRANSPOSE(p_a, nyz, f%size, nyz, p_b, f%size)
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLAB_MPI_TRPB_I(p_b, p_a, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
    end if
#endif

    u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

    nullify (p_a, p_b)

    return
end subroutine OPR_FILTER_X

! ###################################################################
! Filter in Oy direction
! ###################################################################
subroutine OPR_FILTER_Y(nx, ny, nz, f, u, tmp, wrk1d, wrk2d, wrk3d)

    use TLAB_TYPES, only: filter_dt

    implicit none

    TINTEGER, intent(IN) :: nx, ny, nz
    type(filter_dt), intent(IN) :: f
    TREAL, dimension(nx*ny*nz), intent(INOUT), target :: u, wrk3d    ! in-place operation
    TREAL, dimension(*), intent(INOUT) :: wrk2d       ! Aux arrays
    TREAL, dimension(f%size, *), intent(INOUT) :: wrk1d
    TREAL, dimension(nx*ny*nz), intent(INOUT) :: tmp         ! Aux array needed in ADM type

! -----------------------------------------------------------------------
    TINTEGER nxy, nxz

    TREAL, dimension(:), pointer :: p_org, p_dst

! #######################################################################
    nxy = nx*ny
    nxz = nx*nz

! -------------------------------------------------------------------
! Make y-direction the last one
! -------------------------------------------------------------------
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

! ###################################################################
    call OPR_FILTER_1D(nxz, f, p_org, p_dst, wrk1d, wrk2d, tmp)

! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
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

! ###################################################################
! Filter in Oz direction
! ###################################################################
subroutine OPR_FILTER_Z(nx, ny, nz, f, u, tmp, wrk1d, wrk2d, wrk3d)

    use TLAB_TYPES, only: filter_dt
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_k
    use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLAB_MPI_PROCS
#endif

    implicit none

    TINTEGER, intent(IN) :: nx, ny, nz
    type(filter_dt), intent(IN) :: f
    TREAL, dimension(nx*ny*nz), intent(INOUT), target :: u, wrk3d    ! in-place operation
    TREAL, dimension(nx*ny), intent(INOUT) :: wrk2d       ! Aux arrays
    TREAL, dimension(f%size, *), intent(INOUT) :: wrk1d
    TREAL, dimension(nx*ny*nz), intent(INOUT) :: tmp         ! Aux array needed in ADM type

! -------------------------------------------------------------------
    TINTEGER nxy

    TREAL, dimension(:), pointer :: p_a, p_b

#ifdef USE_MPI
    TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
    id = f%mpitype ! TLAB_MPI_K_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
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

! ###################################################################
    call OPR_FILTER_1D(nxy, f, p_a, p_b, wrk1d, wrk2d, tmp)

! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
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
! Spectral band filter
!########################################################################
subroutine OPR_FILTER_BAND_2D(nx, ny, nz, spc_param, a)

    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: g
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    TINTEGER nx, ny, nz
    TREAL, dimension(*) :: spc_param
    TCOMPLEX, dimension(isize_txc_dimz/2, nz), intent(INOUT) :: a

! -----------------------------------------------------------------------
    TINTEGER i, j, k, iglobal, kglobal, ip
    TREAL fi, fk, f

! #######################################################################
    do k = 1, nz
#ifdef USE_MPI
        kglobal = k + ims_offset_k
#else
        kglobal = k
#endif
        if (kglobal <= g(3)%size/2 + 1) then; fk = M_REAL(kglobal - 1)/g(3)%scale
        else; fk = -M_REAL(g(3)%size + 1 - kglobal)/g(3)%scale; end if

        do i = 1, nx/2 + 1
#ifdef USE_MPI
            iglobal = i + ims_offset_i/2
#else
            iglobal = i
#endif
            fi = M_REAL(iglobal - 1)/g(1)%scale

            f = sqrt(fi**2 + fk**2)

! apply spectral cutoff
            do j = 1, ny
                ip = (j - 1)*(nx/2 + 1) + i
                if ((f - spc_param(1))*(spc_param(2) - f) < C_0_R) a(ip, k) = C_0_R
            end do

        end do
    end do

    return
end subroutine OPR_FILTER_BAND_2D

!########################################################################
!# DESCRIPTION
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

    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: g
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    TINTEGER nx, ny, nz
    TREAL, dimension(*) :: spc_param
    TCOMPLEX, dimension(isize_txc_dimz/2, nz), intent(INOUT) :: a

! -----------------------------------------------------------------------
    TINTEGER i, j, k, iglobal, kglobal, ip
    TINTEGER sign_pass, off_pass
    TREAL fi, fk, f, fcut_log, damp

! #######################################################################
    if (spc_param(1) > 0) then; 
        sign_pass = 1.   ! HIGHPASS
        off_pass = 0.
    else                ! spc_param(1) <= 0
        sign_pass = -1.   ! LOWPASS
        off_pass = 1.
    end if

    fcut_log = log(abs(spc_param(1)))
    do k = 1, nz
#ifdef USE_MPI
        kglobal = k + ims_offset_k
#else
        kglobal = k
#endif
        if (kglobal <= g(3)%size/2 + 1) then; fk = M_REAL(kglobal - 1)/g(3)%scale/spc_param(3)
        else; fk = -M_REAL(g(3)%size + 1 - kglobal)/g(3)%scale/spc_param(3); end if

        do i = 1, nx/2 + 1
#ifdef USE_MPI
            iglobal = i + ims_offset_i/2
#else
            iglobal = i
#endif
            fi = M_REAL(iglobal - 1)/g(1)%scale

            f = sqrt(fi**2 + fk**2)
            if (f > 0) then; damp = (erf((log(f) - fcut_log)/spc_param(2)) + 1.)/2.
            else; damp = C_0_R; 
            end if

            ! Set to high- or low-pass
            ! high-pass: damp = 0.0 + damp
            ! low-pass:  damp = 1.0 - damp
            damp = off_pass + sign_pass*damp

            ! apply filter
            do j = 1, ny
                ip = (j - 1)*(nx/2 + 1) + i
                a(ip, k) = damp*a(ip, k)
            end do
        end do
    end do

    return
end subroutine OPR_FILTER_ERF_2D

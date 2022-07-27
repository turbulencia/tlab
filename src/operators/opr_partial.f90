#include "types.h"
#include "dns_const.h"

subroutine OPR_PARTIAL1(nlines, bcs, g, u, result, wrk2d)

    use TLAB_TYPES, only: grid_dt

    implicit none

    TINTEGER, intent(in) :: nlines              ! # of lines to be solved
    TINTEGER, dimension(2), intent(in) :: bcs   ! BCs at xmin (1) and xmax (2):
    !     0 biased, non-zero
    !     1 forced to zero
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines*g%size), intent(in) :: u
    TREAL, dimension(nlines*g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d

! -------------------------------------------------------------------
    TINTEGER ip

! ###################################################################
    if (g%periodic) then
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4P_RHS(g%size, nlines, u, result)

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)   ! Direct = Jacobian because uniform grid
            call FDM_C1N6P_RHS(g%size, nlines, u, result)

        case (FDM_COM6_JACPENTA)                    ! Direct = Jacobian because uniform grid
            call FDM_C1N6MP_RHS(g%size, nlines, u, result)

        case (FDM_COM8_JACOBIAN)
            call FDM_C1N8P_RHS(g%size, nlines, u, result)

        end select

        if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
            call TRIDPSS(g%size, nlines, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), result, wrk2d)
        else
            call PENTADPSS(g%size, nlines, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), &
                           g%lu1(1, 5), g%lu1(1, 6), g%lu1(1, 7), result)
        end if

! -------------------------------------------------------------------
    else
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_RHS(g%size, nlines, bcs(1), bcs(2), u, result)

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_RHS(g%size, nlines, bcs(1), bcs(2), u, result)

        case (FDM_COM6_JACPENTA)
            call FDM_C1N6M_RHS(g%size, nlines, bcs(1), bcs(2), u, result)

        case (FDM_COM8_JACOBIAN)
            call FDM_C1N8_RHS(g%size, nlines, bcs(1), bcs(2), u, result)

        case (FDM_COM6_DIRECT) ! Not yet implemented
            call FDM_C1N6_RHS(g%size, nlines, bcs(1), bcs(2), u, result)

        end select

        ip = (bcs(1) + bcs(2)*2)*5
        if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
            call TRIDSS(g%size, nlines, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3), result)
        else
            call PENTADSS2(g%size, nlines, g%lu1(1, ip + 1), g%lu1(1, ip + 2), g%lu1(1, ip + 3), g%lu1(1, ip + 4), &
                           g%lu1(1, ip + 5), result)
        end if

    end if

    return
end subroutine OPR_PARTIAL1

! ###################################################################
! ###################################################################
subroutine OPR_PARTIAL1_IBM(nlines, bcs, g, u, result, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: fld_ibm
    use IBM_VARS, only: nobi, nobj, nobk
    use IBM_VARS, only: nobi_b, nobj_b, nobk_b
    use IBM_VARS, only: nobi_e, nobj_e, nobk_e
    use IBM_VARS, only: isize_nobi, isize_nobj, isize_nobk
    use IBM_VARS, only: isize_nobi_be, isize_nobj_be, isize_nobk_be
    use IBM_VARS, only: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

    implicit none

#include "integers.h"

    TINTEGER, intent(in) :: nlines                  ! # of lines to be solved
    TINTEGER, dimension(2, *), intent(in) :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
    !     0 biased, non-zero
    !     1 forced to zero
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines*g%size), intent(in) :: u
    TREAL, dimension(nlines*g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d
    TREAL, dimension(nlines*g%size), intent(inout) :: wrk3d

    TINTEGER, parameter :: is = i0                  ! scalar index; if 0, then velocity

    ! -------------------------------------------------------------------
    ! modify incoming fields (fill solids with spline functions, depending on direction)

    select case (g%name)

    case ('x')
        if (ims_pro_ibm_x) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
            call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL1(nlines, bcs, g, u, result, wrk2d)  ! no splines needed
        end if

    case ('y')
        if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
            call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL1(nlines, bcs, g, u, result, wrk2d)  ! no splines needed
        end if

    case ('z')
        if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
            call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL1(nlines, bcs, g, u, result, wrk2d)  ! no splines needed
        end if

    end select

    return
end subroutine OPR_PARTIAL1_IBM

! ###################################################################
! ###################################################################
subroutine OPR_IBM(nlines, g, u, result, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: nobi, nobj, nobk
    use IBM_VARS, only: nobi_b, nobj_b, nobk_b
    use IBM_VARS, only: nobi_e, nobj_e, nobk_e
    use IBM_VARS, only: isize_nobi, isize_nobj, isize_nobk
    use IBM_VARS, only: isize_nobi_be, isize_nobj_be, isize_nobk_be

    implicit none

#include "integers.h"

    TINTEGER, intent(in) :: nlines
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines*g%size), intent(in) :: u
    TREAL, dimension(nlines*g%size), intent(out) :: result
    TREAL, dimension(nlines*g%size), intent(inout) :: wrk3d

    TINTEGER, parameter :: is = i0 ! scalar index; if 0, then velocity

    ! -------------------------------------------------------------------
    ! modify incoming fields (fill solids with spline functions, depending on direction)

    select case (g%name)
    case ('x')
        call IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
    case ('y')
        call IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
    case ('z')
        call IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
    end select

    return
end subroutine OPR_IBM

! ###################################################################################
! ###################################################################################
subroutine OPR_PARTIAL2(is, nlines, bcs, g, u, result, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt

    implicit none

    TINTEGER, intent(in) :: is                      ! premultiplying factor in second derivative
    ! -1            factor 1, pure derivative
    !  0            viscosity (velocity)
    !  1:inb_scal   diffusivity
    TINTEGER, intent(in) :: nlines                  ! # of lines to be solved
    TINTEGER, dimension(2, *), intent(in) :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
    !     0 biased, non-zero
    !     1 forced to zero
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines, g%size), intent(in) :: u
    TREAL, dimension(nlines, g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d
    TREAL, dimension(nlines, g%size), intent(inout) :: wrk3d  ! First derivative

    TREAL, dimension(:, :), pointer :: lu2_p

! -------------------------------------------------------------------
    TINTEGER ip

    ! ###################################################################
    ! Check whether to calculate 1. order derivative
    ! ###################################################################
    if (is >= 0) then   ! called from opr_burgers
        call OPR_PARTIAL1(nlines, bcs, g, u, wrk3d, wrk2d)
    else                ! needed anyhow to calculate 1. derivative
        if (.not. g%uniform) then
            if (g%mode_fdm == FDM_COM4_JACOBIAN .or. &
                g%mode_fdm == FDM_COM6_JACOBIAN .or. &
                g%mode_fdm == FDM_COM8_JACOBIAN) then
                call OPR_PARTIAL1(nlines, bcs, g, u, wrk3d, wrk2d)
            end if
        end if
    end if

    if (is >= 0) then
        if (g%periodic) then
            lu2_p => g%lu2d(:, is*5 + 1:)       ! periodic;     including diffusivity/viscosity
        else
            lu2_p => g%lu2d(:, is*3 + 1:)       ! non-periodic; including diffusivity/viscosity
        end if
    else
        if (g%periodic) then
            lu2_p => g%lu2(:, 1:)               ! periodic;     plain derivative
        else
            ip = (bcs(1, 2) + bcs(2, 2)*2)*3    ! non-periodic; plain derivative
            lu2_p => g%lu2(:, ip + 1:)
        end if
    end if

    ! ###################################################################
    if (g%periodic) then
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4P_RHS(g%size, nlines, u, result)

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA) ! Direct = Jacobian because uniform grid
            call FDM_C2N6HP_RHS(g%size, nlines, u, result)

        case (FDM_COM8_JACOBIAN)                  ! Not yet implemented
            call FDM_C2N6P_RHS(g%size, nlines, u, result)

        end select

        call TRIDPSS(g%size, nlines, lu2_p(1, 1), lu2_p(1, 2), lu2_p(1, 3), lu2_p(1, 4), lu2_p(1, 5), result, wrk2d)

        ! -------------------------------------------------------------------
    else
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN)
            if (g%uniform) then
                call FDM_C2N4_RHS(g%size, nlines, bcs(1, 2), bcs(2, 2), u, result)
            else ! Not yet implemented
            end if

        case (FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA)
            if (g%uniform) then
                call FDM_C2N6H_RHS(g%size, nlines, bcs(1, 2), bcs(2, 2), u, result)
            else        ! need first derivative from above
                call FDM_C2N6HNJ_RHS(g%size, nlines, bcs(1, 2), bcs(2, 2), g%jac, u, wrk3d, result)
            end if

        case (FDM_COM8_JACOBIAN) ! Not yet implemented; defaulting to 6. order
            if (g%uniform) then
                call FDM_C2N6_RHS(g%size, nlines, bcs(1, 2), bcs(2, 2), u, result)
            else        ! Need first derivative from above
                call FDM_C2N6NJ_RHS(g%size, nlines, bcs(1, 2), bcs(2, 2), g%jac, u, wrk3d, result)
            end if

        case (FDM_COM6_DIRECT)
            call FDM_C2N6ND_RHS(g%size, nlines, g%lu2(1, 4), u, result)

        end select

        call TRIDSS(g%size, nlines, lu2_p(1, 1), lu2_p(1, 2), lu2_p(1, 3), result)

    end if

    return
end subroutine OPR_PARTIAL2

! ###################################################################
! ###################################################################
subroutine OPR_PARTIAL2_IBM(is, nlines, bcs, g, u, result, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: fld_ibm
    use IBM_VARS, only: nobi, nobj, nobk
    use IBM_VARS, only: nobi_b, nobj_b, nobk_b
    use IBM_VARS, only: nobi_e, nobj_e, nobk_e
    use IBM_VARS, only: isize_nobi, isize_nobj, isize_nobk
    use IBM_VARS, only: isize_nobi_be, isize_nobj_be, isize_nobk_be
    use IBM_VARS, only: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

    implicit none

    TINTEGER, intent(in) :: is     ! scalar index; if 0, then velocity
    TINTEGER, intent(in) :: nlines ! # of lines to be solved
    TINTEGER, dimension(2, *), intent(in) :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
    !     0 biased, non-zero
    !     1 forced to zero
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines, g%size), intent(in), target :: u
    TREAL, dimension(nlines, g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d
    TREAL, dimension(nlines, g%size), intent(inout) :: wrk3d  ! First derivative

    TREAL, dimension(:, :), pointer :: p_fld
    TREAL, dimension(:), pointer :: p_fld_ibm

    ! -------------------------------------------------------------------

    ! pointer to field
    p_fld => u

    ! -------------------------------------------------------------------
    ! modify incoming fields (fill solids with spline functions, depending on direction)

    select case (g%name)

    case ('x')
        if (ims_pro_ibm_x) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
            p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, wrk2d, wrk3d)  ! no splines needed
        end if

    case ('y')
        if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
            p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, wrk2d, wrk3d)  ! no splines needed
        end if

    case ('z')
        if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
            call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
            p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
        else ! idle IBM-Tasks
            call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, wrk2d, wrk3d)  ! no splines needed
        end if

    end select

    ! -------------------------------------------------------------------

    nullify (p_fld, p_fld_ibm)

    return
end subroutine OPR_PARTIAL2_IBM
! ###################################################################
#include "dns_error.h"
! ###################################################################

subroutine OPR_PARTIAL0_INT(dir, nlines, g, u, result, wrk2d)

    use TLAB_TYPES, only: grid_dt
    use TLAB_PROCS, only: TLAB_STOP, TLAB_WRITE_ASCII
    use TLAB_CONSTANTS, only: efile

    implicit none

    TINTEGER, intent(in) :: dir     ! scalar direction flag
    !     0 'vp' --> vel. to pre. grid
    !     1 'pv' --> pre. to vel. grid
    TINTEGER, intent(in) :: nlines  ! number of lines to be solved
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines, g%size), intent(in) :: u
    TREAL, dimension(nlines, g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d

! ###################################################################
! Interpolation, direction 'vp': vel. --> pre. grid
    if (dir == 0) then
        if (g%periodic) then
            select case (g%mode_fdm)
            case DEFAULT
                call FDM_C0INTVP6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
        else
            call TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
            call TLAB_STOP(DNS_ERROR_NOTIMPL)
        end if
! Interpolation, direction 'pv': pre. --> vel. grid
    else if (dir == 1) then
        if (g%periodic) then
            select case (g%mode_fdm)
            case DEFAULT
                call FDM_C0INTPV6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
        else
            call TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
            call TLAB_STOP(DNS_ERROR_NOTIMPL)
        end if
    end if

    return
end subroutine OPR_PARTIAL0_INT

! ###################################################################
! ###################################################################

subroutine OPR_PARTIAL1_INT(dir, nlines, g, u, result, wrk2d)

    use TLAB_TYPES, only: grid_dt
    use TLAB_PROCS, only: TLAB_STOP, TLAB_WRITE_ASCII
    use TLAB_CONSTANTS, only: efile

    implicit none

    TINTEGER, intent(in) :: dir    ! scalar direction flag
    !     0 'vp' --> vel. to pre. grid
    !     1 'pv' --> pre. to vel. grid
    TINTEGER, intent(in) :: nlines ! number of lines to be solved
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nlines, g%size), intent(in) :: u
    TREAL, dimension(nlines, g%size), intent(out) :: result
    TREAL, dimension(nlines), intent(inout) :: wrk2d

! ###################################################################
! 1st interpolatory derivative, direction 'vp': vel. --> pre. grid
    if (dir == 0) then
        if (g%periodic) then
            select case (g%mode_fdm)
            case (FDM_COM6_JACOBIAN)
                call FDM_C1INTVP6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
        else
            call TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
            call TLAB_STOP(DNS_ERROR_NOTIMPL)
        end if
! 1st interpolatory derivative, direction 'pv': pre. --> vel. grid
    else if (dir == 1) then
        if (g%periodic) then
            select case (g%mode_fdm)
            case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM8_JACOBIAN)
                call FDM_C1INTPV6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
        else
            call TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
            call TLAB_STOP(DNS_ERROR_NOTIMPL)
        end if
    end if

    return
end subroutine OPR_PARTIAL1_INT

! ###################################################################
! ###################################################################

#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Routines for different specific directions
!########################################################################
subroutine OPR_PARTIAL_X(type, nx, ny, nz, bcs, g, u, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: ibm_partial
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_i
    use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_PROCS
#endif

    implicit none

#include "integers.h"

    TINTEGER, intent(in) :: type        ! OPR_P1           1.order derivative
    ! OPR_P2           2.order derivative
    ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
    ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
    ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
    TINTEGER, intent(in) :: nx, ny, nz  ! array sizes
    TINTEGER, dimension(2, *), intent(in) :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*ny*nz), intent(in) :: u
    TREAL, dimension(nx*ny*nz), intent(out) :: result
    TREAL, dimension(nx*ny*nz), intent(inout) :: tmp1, wrk3d
    TREAL, dimension(ny*nz), intent(inout) :: wrk2d

    target u, tmp1, result, wrk3d

! -------------------------------------------------------------------
    TINTEGER nyz
    TINTEGER, parameter :: is = -1 ! second derivative without viscosity/diffusivity

    TREAL, dimension(:), pointer :: p_a, p_b, p_c, p_d

#ifdef USE_MPI
    TINTEGER, parameter :: id = TLAB_MPI_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLAB_MPI_TRPF_I(u, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        p_a => result
        p_b => wrk3d
        p_c => result
        p_d => tmp1
        nyz = ims_size_i(id)
    else
#endif
        p_a => u
        p_b => result
        if (type == OPR_P2_P1) then
            p_c => tmp1
            p_d => wrk3d
        else
            p_c => wrk3d
            p_d => tmp1
        end if
        nyz = ny*nz
#ifdef USE_MPI
    end if
#endif

! -------------------------------------------------------------------
! Local transposition: make x-direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
    call DGETMO(p_a, g%size, g%size, nyz, p_b, nyz)
#else
    call DNS_TRANSPOSE(p_a, g%size, nyz, g%size, p_b, nyz)
#endif

! ###################################################################
    select case (type)

    case (OPR_P2)
        call OPR_PARTIAL2(is, nyz, bcs, g, p_b, p_c, wrk2d, p_d)

    case (OPR_P1)
        if (ibm_partial) then
            call OPR_PARTIAL1_IBM(nyz, bcs, g, p_b, p_c, wrk2d, p_d)
        else
            call OPR_PARTIAL1(nyz, bcs, g, p_b, p_c, wrk2d)
        end if

    case (OPR_P2_P1)
        call OPR_PARTIAL2(is, nyz, bcs, g, p_b, p_c, wrk2d, p_d)

! Check whether we need to calculate the 1. order derivative
        if (g%uniform .or. g%mode_fdm == FDM_COM6_DIRECT) then
            call OPR_PARTIAL1(nyz, bcs, g, p_b, p_d, wrk2d)
        end if

    case (OPR_P0_INT_VP)
        call OPR_PARTIAL0_INT(i0, nyz, g, p_b, p_c, wrk2d)

    case (OPR_P0_INT_PV)
        call OPR_PARTIAL0_INT(i1, nyz, g, p_b, p_c, wrk2d)

    case (OPR_P1_INT_VP)
        call OPR_PARTIAL1_INT(i0, nyz, g, p_b, p_c, wrk2d)

    case (OPR_P1_INT_PV)
        call OPR_PARTIAL1_INT(i1, nyz, g, p_b, p_c, wrk2d)

    case (OPR_P0_IBM)
        call OPR_IBM(nyz, g, p_b, p_c, p_d)

    end select

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
    call DGETMO(p_c, nyz, nyz, g%size, p_b, g%size)
#else
    call DNS_TRANSPOSE(p_c, nyz, g%size, nyz, p_b, g%size)
#endif

    if (type == OPR_P2_P1) then
#ifdef USE_ESSL
        call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
        call DNS_TRANSPOSE(p_d, nyz, g%size, nyz, p_c, g%size)
#endif
    end if

#ifdef USE_MPI
    if (ims_npro_i > 1) then
        if (type == OPR_P2_P1) then ! only if you really want first derivative back
            call TLAB_MPI_TRPB_I(p_c, tmp1, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        end if
        call TLAB_MPI_TRPB_I(p_b, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
    end if
#endif

    nullify (p_a, p_b, p_c, p_d)

    return
end subroutine OPR_PARTIAL_X

!########################################################################
!########################################################################
subroutine OPR_PARTIAL_Z(type, nx, ny, nz, bcs, g, u, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: ibm_partial
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_k
    use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLAB_MPI_PROCS
#endif

    implicit none

#include "integers.h"

    TINTEGER, intent(in) :: type        ! OPR_P1           1.order derivative
    ! OPR_P2           2.order derivative
    ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
    ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
    ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
    TINTEGER, intent(in) :: nx, ny, nz  ! array sizes
    TINTEGER, dimension(2, *), intent(in) :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*ny*nz), intent(in) :: u
    TREAL, dimension(nx*ny*nz), intent(out) :: result
    TREAL, dimension(nx*ny*nz), intent(inout) :: tmp1, wrk3d
    TREAL, dimension(nx*ny), intent(inout) :: wrk2d

    target u, tmp1, result, wrk3d

! -------------------------------------------------------------------
    TINTEGER nxy
    TINTEGER, parameter :: is = -1 ! second derivative without viscosity/diffusivity

    TREAL, dimension(:), pointer :: p_a, p_b, p_c

#ifdef USE_MPI
    TINTEGER, parameter :: id = TLAB_MPI_K_PARTIAL
#endif

! ###################################################################
    if (g%size == 1) then ! Set to zero in 2D case
        result = C_0_R
        if (type == OPR_P2_P1) tmp1 = C_0_R

    else
! ###################################################################
! -------------------------------------------------------------------
! MPI Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPF_K(u, result, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            p_a => result
            if (type == OPR_P2_P1) then
                p_b => tmp1
                p_c => wrk3d
            else
                p_b => wrk3d
                p_c => tmp1
            end if
            nxy = ims_size_k(id)
        else
#endif
            p_a => u
            p_b => result
            p_c => tmp1
            nxy = nx*ny
#ifdef USE_MPI
        end if
#endif

! ###################################################################
        select case (type)

        case (OPR_P2)
            call OPR_PARTIAL2(is, nxy, bcs, g, p_a, p_b, wrk2d, p_c)

        case (OPR_P1)
            if (ibm_partial) then
                call OPR_PARTIAL1_IBM(nxy, bcs, g, p_a, p_b, wrk2d, p_c)
            else
                call OPR_PARTIAL1(nxy, bcs, g, p_a, p_b, wrk2d)
            end if
        case (OPR_P2_P1)
            call OPR_PARTIAL2(is, nxy, bcs, g, p_a, p_b, wrk2d, p_c)

! Check whether we need to calculate the 1. order derivative
            if (g%uniform .or. g%mode_fdm == FDM_COM6_DIRECT) then
                call OPR_PARTIAL1(nxy, bcs, g, p_a, p_c, wrk2d)
            end if

        case (OPR_P0_INT_VP)
            call OPR_PARTIAL0_INT(i0, nxy, g, p_a, p_b, wrk2d)

        case (OPR_P0_INT_PV)
            call OPR_PARTIAL0_INT(i1, nxy, g, p_a, p_b, wrk2d)

        case (OPR_P1_INT_VP)
            call OPR_PARTIAL1_INT(i0, nxy, g, p_a, p_b, wrk2d)

        case (OPR_P1_INT_PV)
            call OPR_PARTIAL1_INT(i1, nxy, g, p_a, p_b, wrk2d)

        case (OPR_P0_IBM)
            call OPR_IBM(nxy, g, p_a, p_b, p_c)

        end select

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPB_K(p_b, result, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            if (type == OPR_P2_P1) then
                call TLAB_MPI_TRPB_K(p_c, tmp1, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            end if
        end if
#endif

        nullify (p_a, p_b, p_c)

    end if

    return
end subroutine OPR_PARTIAL_Z

!########################################################################
!########################################################################
subroutine OPR_PARTIAL_Y(type, nx, ny, nz, bcs, g, u, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use IBM_VARS, only: ibm_partial
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif

    implicit none

#include "integers.h"

    TINTEGER, intent(in) :: type        ! OPR_P1           1.order derivative
    ! OPR_P2           2.order derivative
    ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
    ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
    ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
    TINTEGER, intent(in) :: nx, ny, nz  ! array sizes
    TINTEGER, dimension(2, *), intent(in) :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*ny*nz), intent(in) :: u
    TREAL, dimension(nx*ny*nz), intent(out) :: result
    TREAL, dimension(nx*ny*nz), intent(inout) :: tmp1, wrk3d
    TREAL, dimension(nx*nz), intent(inout) :: wrk2d

    target u, tmp1, result, wrk3d

! -------------------------------------------------------------------
    TINTEGER nxy, nxz
    TINTEGER, parameter :: is = -1 ! second derivative without viscosity/diffusivity
    TREAL, dimension(:), pointer :: p_a, p_b, p_c

! ###################################################################
    if (g%size == 1) then ! Set to zero in 2D case
        result = C_0_R
        if (type == OPR_P2_P1) tmp1 = C_0_R

    else
! ###################################################################
        nxy = nx*ny
        nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y direction the last one
! -------------------------------------------------------------------
        if (nz == 1) then
            p_a => u
            p_b => result
            p_c => tmp1
        else
#ifdef USE_ESSL
            call DGETMO(u, nxy, nxy, nz, result, nz)
#else
            call DNS_TRANSPOSE(u, nxy, nz, nxy, result, nz)
#endif
            p_a => result
            if (type == OPR_P2_P1) then
                p_b => tmp1
                p_c => wrk3d
            else
                p_b => wrk3d
                p_c => tmp1
            end if
        end if

! ###################################################################
        select case (type)

        case (OPR_P2)
            call OPR_PARTIAL2(is, nxz, bcs, g, p_a, p_b, wrk2d, p_c)

        case (OPR_P1)
            if (ibm_partial) then
                call OPR_PARTIAL1_IBM(nxz, bcs, g, p_a, p_b, wrk2d, p_c)
            else
                call OPR_PARTIAL1(nxz, bcs, g, p_a, p_b, wrk2d)
            end if

        case (OPR_P2_P1)
            call OPR_PARTIAL2(is, nxz, bcs, g, p_a, p_b, wrk2d, p_c)

! Check whether we need to calculate the 1. order derivative
            if (g%uniform .or. g%mode_fdm == FDM_COM6_DIRECT) then
                call OPR_PARTIAL1(nxz, bcs, g, p_a, p_c, wrk2d)
            end if

        case (OPR_P0_INT_VP)
            call OPR_PARTIAL0_INT(i0, nxz, g, p_a, p_b, wrk2d)

        case (OPR_P0_INT_PV)
            call OPR_PARTIAL0_INT(i1, nxz, g, p_a, p_b, wrk2d)

        case (OPR_P1_INT_VP)
            call OPR_PARTIAL1_INT(i0, nxz, g, p_a, p_b, wrk2d)

        case (OPR_P1_INT_PV)
            call OPR_PARTIAL1_INT(i1, nxz, g, p_a, p_b, wrk2d)

        case (OPR_P0_IBM)
            call OPR_IBM(nxz, g, p_a, p_b, p_c)

        end select

! ###################################################################
! Put arrays back in the order in which they came in
        if (nz > 1) then
#ifdef USE_ESSL
            call DGETMO(p_b, nz, nz, nxy, result, nxy)
#else
            call DNS_TRANSPOSE(p_b, nz, nxy, nz, result, nxy)
#endif
            if (type == OPR_P2_P1) then
#ifdef USE_ESSL
                call DGETMO(p_c, nz, nz, nxy, tmp1, nxy)
#else
                call DNS_TRANSPOSE(p_c, nz, nxy, nz, tmp1, nxy)
#endif
            end if
        end if

        nullify (p_a, p_b, p_c)

    end if

    return
end subroutine OPR_PARTIAL_Y

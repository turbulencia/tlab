#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module OPR_PARTIAL
    use TLab_Constants, only: efile, wp, wi, BCS_DN, BCS_ND, BCS_NN
    use FDM, only: grid_dt
    use TLab_WorkFlow, only: TLab_Stop, TLab_Write_ASCII
    use IBM_VARS, only: ibm_partial
    use IBM_VARS, only: fld_ibm
    use IBM_VARS, only: nobi, nobj, nobk
    use IBM_VARS, only: nobi_b, nobj_b, nobk_b
    use IBM_VARS, only: nobi_e, nobj_e, nobk_e
    use IBM_VARS, only: isize_nobi, isize_nobj, isize_nobk
    use IBM_VARS, only: isize_nobi_be, isize_nobj_be, isize_nobk_be
    use IBM_VARS, only: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z
    use IBM_VARS, only: ibm_case_x, ibm_case_y, ibm_case_z
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_PROCS
#endif
    use FDM_MatMul
    implicit none
    private

    integer(wi) ip, ibc

    public :: OPR_PARTIAL_X
    public :: OPR_PARTIAL_Y
    public :: OPR_PARTIAL_Z
    public :: OPR_PARTIAL2, OPR_PARTIAL2_IBM

contains
! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL1(nlines, bcs, g, u, result)
        use TLab_Arrays, only: wrk2d
        integer(wi), intent(in) :: nlines   ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2)   ! BCs at xmin (1) and xmax (2):
        !                                   0 biased, non-zero
        !                                   1 forced to zero
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)

        integer(wi) nmin, nmax, nsize

! ###################################################################
        ibc = bcs(1) + bcs(2)*2
        ip = ibc*5

        nmin = 1; nmax = g%size
        if (any([BCS_ND, BCS_NN] == ibc)) then
            result(:, 1) = 0.0_wp      ! homogeneous bcs
            nmin = nmin + 1
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            result(:, g%size) = 0.0_wp
            nmax = nmax - 1
        end if
        nsize = nmax - nmin + 1

        select case (g%nb_diag_1(2))
        case (3)
            call MatMul_3d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, result, &
                                   g%periodic, ibc, g%rhs1_b, g%rhs1_t)
        case (5)
            call MatMul_5d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), u, result, &
                                   g%periodic, ibc, g%rhs1_b, g%rhs1_t)
        case (7)
            call MatMul_7d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), u, result, &
                                   g%periodic, ibc, g%rhs1_b, g%rhs1_t)
        end select

        if (g%periodic) then
            select case (g%nb_diag_1(1))
            case (3)
                call TRIDPSS(g%size, nlines, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), result, wrk2d)
            case (5)
                call PENTADPSS(g%size, nlines, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), g%lu1(1, 6), g%lu1(1, 7), result)
            end select

        else
            select case (g%nb_diag_1(1))
            case (3)
                call TRIDSS(nsize, nlines, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3), result(:, nmin:))
            case (5)
                call PENTADSS2(nsize, nlines, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3), g%lu1(nmin:, ip + 4), g%lu1(nmin:, ip + 5), result(:, nmin:))
            end select

        end if

        return
    end subroutine OPR_PARTIAL1

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL1_IBM(nlines, bcs, g, u, result)
        implicit none
        integer(wi), intent(in) :: nlines   ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2)   ! BCs at xmin (1,*) and xmax (2,*):
        !                                   0 biased, non-zero
        !                                   1 forced to zero
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines*g%size)
        real(wp), intent(out) :: result(nlines*g%size)

        integer(wi), parameter :: is = 0    ! scalar index; if 0, then velocity

        ! -------------------------------------------------------------------
        ! modify incoming fields (fill solids with spline functions, depending on direction)

        select case (g%name)

        case ('x')
            if (ims_pro_ibm_x) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
                call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL1(nlines, bcs, g, u, result)  ! no splines needed
            end if

        case ('y')
            if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
                call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL1(nlines, bcs, g, u, result)  ! no splines needed
            end if

        case ('z')
            if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
                call OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL1(nlines, bcs, g, u, result)  ! no splines needed
            end if

        end select

        return
    end subroutine OPR_PARTIAL1_IBM

! ###################################################################
! ###################################################################
    subroutine OPR_IBM(nlines, g, u, result)
        integer(wi), intent(in) :: nlines
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines*g%size)
        real(wp), intent(out) :: result(nlines*g%size)

        integer(wi), parameter :: is = 0 ! scalar index; if 0, then velocity

        ! -------------------------------------------------------------------
        ! modify incoming fields (fill solids with spline functions, depending on direction)

        select case (g%name)
        case ('x')
            call IBM_SPLINE_XYZ(is, u, result, g, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
        case ('y')
            call IBM_SPLINE_XYZ(is, u, result, g, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
        case ('z')
            call IBM_SPLINE_XYZ(is, u, result, g, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
        end select

        return
    end subroutine OPR_IBM

! ###################################################################################
! ###################################################################################
    subroutine OPR_PARTIAL2(is, nlines, bcs, g, u, result, du)
        ! bcs(*, 2) are not used, need to be updated
        use TLab_Arrays, only: wrk2d

        integer(wi), intent(in) :: is           ! premultiplying factor in second derivative
        !                                       -1            factor 1, pure derivative
        !                                       0            viscosity (velocity)
        !                                       1:inb_scal   diffusivity
        integer(wi), intent(in) :: nlines                  ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)    ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: du(nlines, g%size)  ! First derivative

        real(wp), pointer :: lu2_p(:, :)

        ! ###################################################################
        ! Check whether to calculate 1. order derivative
        ! ###################################################################
        if (is >= 0 .or. g%need_1der) then   ! called from opr_burgers or need for 1. order derivative
            call OPR_PARTIAL1(nlines, bcs(:, 1), g, u, du)
        end if

        if (is >= 0) then
            if (g%periodic) then
                lu2_p => g%lu2d(:, is*5 + 1:)       ! periodic;     including diffusivity/viscosity
            else
                lu2_p => g%lu2d(:, is*3 + 1:)       ! non-periodic; including diffusivity/viscosity
            end if
        else
            lu2_p => g%lu2(:, 1:)
        end if

        ! ###################################################################
        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm2)) then
            ! so far, only pentadiagonal cases
            call MatMul_5d(g%size, nlines, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), u, result)
        else
            select case (g%nb_diag_2(2))
            case (5)
                call MatMul_5d_sym(g%size, nlines, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), u, result, g%periodic)
            case (7)
     call MatMul_7d_sym(g%size, nlines, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%rhs2(:, 6), g%rhs2(:, 7), u, result, g%periodic)
            end select
            if (g%need_1der) then
                ip = g%nb_diag_2(2)      ! add Jacobian correction A_2 dx2 du
                ! so far, only tridiagonal systems
                call MatMul_3d_add(g%size, nlines, g%rhs2(:, ip + 1), g%rhs2(:, ip + 2), g%rhs2(:, ip + 3), du, result)
            end if
        end if

        if (g%periodic) then    ! so far, tridiagonal
            call TRIDPSS(g%size, nlines, lu2_p(1, 1), lu2_p(1, 2), lu2_p(1, 3), lu2_p(1, 4), lu2_p(1, 5), result, wrk2d)

        else
            call TRIDSS(g%size, nlines, lu2_p(:, 1), lu2_p(:, 2), lu2_p(:, 3), result)

        end if

        return
    end subroutine OPR_PARTIAL2

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL2_IBM(is, nlines, bcs, g, u, result, du)
        integer(wi), intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)     ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: du(nlines, g%size)  ! First derivative

        real(wp), pointer :: p_fld(:, :)
        real(wp), pointer :: p_fld_ibm(:)

        target u

        ! -------------------------------------------------------------------

        ! pointer to field
        p_fld => u

        ! -------------------------------------------------------------------
        ! modify incoming fields (fill solids with spline functions, depending on direction)

        select case (g%name)

        case ('x')
            if (ims_pro_ibm_x) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
                p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, du)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, du)  ! no splines needed
            end if

        case ('y')
            if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
                p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, du)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, du)  ! no splines needed
            end if

        case ('z')
            if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
                p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld_ibm, result, du)  ! now with modified u fields
            else ! idle IBM-Tasks
                call OPR_PARTIAL2(is, nlines, bcs, g, p_fld, result, du)  ! no splines needed
            end if

        end select

        ! -------------------------------------------------------------------

        nullify (p_fld, p_fld_ibm)

        return
    end subroutine OPR_PARTIAL2_IBM

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL0_INT(dir, nlines, g, u, result)
        use TLab_Arrays, only: wrk2d
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)

! ###################################################################
! Interpolation, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case DEFAULT
                    call FDM_C0INTVP6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
! Interpolation, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case DEFAULT
                    call FDM_C0INTPV6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
        end if

        return
    end subroutine OPR_PARTIAL0_INT

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL1_INT(dir, nlines, g, u, result)
        use TLab_Arrays, only: wrk2d
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)

! ###################################################################
! 1st interpolatory derivative, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case (FDM_COM6_JACOBIAN)
                    call FDM_C1INTVP6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
! 1st interpolatory derivative, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)
                    call FDM_C1INTPV6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
        end if

        return
    end subroutine OPR_PARTIAL1_INT

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL_X(type, nx, ny, nz, bcs, g, u, result, tmp1)
        use TLab_Arrays, only: wrk3d
        integer(wi), intent(in) :: type     ! OPR_P1         1.order derivative
        !                                   OPR_P2           2.order derivative
        !                                   OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
        !                                   OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
        !                                   OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(:, :)       ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nyz
        integer(wi), parameter :: is = -1 ! second derivative without viscosity/diffusivity

        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d

#ifdef USE_MPI
        integer(wi), parameter :: id = TLAB_MPI_TRP_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(u, result, id)
            p_a => result
            p_b => wrk3d
            p_c => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_d => tmp1
            end if
            nyz = ims_size_i(id)
        else
#endif
            p_a => u
            p_b => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_c => tmp1
                p_d => wrk3d
            else
                p_c => wrk3d
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
        call TLab_Transpose(p_a, g%size, nyz, g%size, p_b, nyz)
#endif

! ###################################################################
        select case (type)

        case (OPR_P2)
            call OPR_PARTIAL2(is, nyz, bcs, g, p_b, p_c, p_d)

        case (OPR_P2_P1)
            call OPR_PARTIAL2(is, nyz, bcs, g, p_b, p_c, p_d)

            if (.not. g%need_1der) then  ! Check whether we need to calculate the 1. order derivative
                call OPR_PARTIAL1(nyz, bcs(:, 1), g, p_b, p_d)
            end if

        case (OPR_P1)
            if (ibm_partial) then
                call OPR_PARTIAL1_IBM(nyz, bcs(:, 1), g, p_b, p_c)
            else
                call OPR_PARTIAL1(nyz, bcs(:, 1), g, p_b, p_c)
            end if

        case (OPR_P0_INT_VP)
            call OPR_PARTIAL0_INT(0, nyz, g, p_b, p_c)

        case (OPR_P0_INT_PV)
            call OPR_PARTIAL0_INT(1, nyz, g, p_b, p_c)

        case (OPR_P1_INT_VP)
            call OPR_PARTIAL1_INT(0, nyz, g, p_b, p_c)

        case (OPR_P1_INT_PV)
            call OPR_PARTIAL1_INT(1, nyz, g, p_b, p_c)

        case (OPR_P0_IBM)
            call OPR_IBM(nyz, g, p_b, p_c)

        end select

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_c, nyz, nyz, g%size, p_b, g%size)
#else
        call TLab_Transpose(p_c, nyz, g%size, nyz, p_b, g%size)
#endif

        if (type == OPR_P2_P1) then
#ifdef USE_ESSL
            call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
            call TLab_Transpose(p_d, nyz, g%size, nyz, p_c, g%size)
#endif
        end if

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            if (type == OPR_P2_P1) then ! only if you really want first derivative back
                call TLabMPI_TransposeI_Backward(p_c, tmp1, id)
            end if
            call TLabMPI_TransposeI_Backward(p_b, result, id)
        end if
#endif

        nullify (p_a, p_b, p_c)
        if (associated(p_c)) nullify (p_d)

        return
    end subroutine OPR_PARTIAL_X

!########################################################################
!########################################################################
    subroutine OPR_PARTIAL_Z(type, nx, ny, nz, bcs, g, u, result, tmp1)
#ifdef USE_MPI
        use TLab_Arrays, only: wrk3d
#endif
        integer(wi), intent(in) :: type     ! OPR_P1           1.order derivative
        !                                   OPR_P2           2.order derivative
        !                                   OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
        !                                   OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
        !                                   OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(:, :)       ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nxy
        integer(wi), parameter :: is = -1 ! second derivative without viscosity/diffusivity

        real(wp), dimension(:), pointer :: p_a, p_b, p_c

#ifdef USE_MPI
        integer(wi), parameter :: id = TLAB_MPI_TRP_K_PARTIAL
#endif

! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp

        else
! ###################################################################
! -------------------------------------------------------------------
! MPI Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TransposeK_Forward(u, result, id)
                p_a => result
                if (any([OPR_P2, OPR_P2_P1] == type)) then
                    p_b => tmp1
                    p_c => wrk3d
                else
                    p_b => wrk3d
                end if
                nxy = ims_size_k(id)
            else
#endif
                p_a => u
                p_b => result
                if (any([OPR_P2, OPR_P2_P1] == type)) then
                    p_c => tmp1
                end if
                nxy = nx*ny
#ifdef USE_MPI
            end if
#endif

! ###################################################################
            select case (type)

            case (OPR_P2)
                call OPR_PARTIAL2(is, nxy, bcs, g, p_a, p_b, p_c)

            case (OPR_P2_P1)
                call OPR_PARTIAL2(is, nxy, bcs, g, p_a, p_b, p_c)

                if (.not. g%need_1der) then  ! Check whether we need to calculate the 1. order derivative
                    call OPR_PARTIAL1(nxy, bcs(:, 1), g, p_a, p_c)
                end if

            case (OPR_P1)
                if (ibm_partial) then
                    call OPR_PARTIAL1_IBM(nxy, bcs(:, 1), g, p_a, p_b)
                else
                    call OPR_PARTIAL1(nxy, bcs(:, 1), g, p_a, p_b)
                end if

            case (OPR_P0_INT_VP)
                call OPR_PARTIAL0_INT(0, nxy, g, p_a, p_b)

            case (OPR_P0_INT_PV)
                call OPR_PARTIAL0_INT(1, nxy, g, p_a, p_b)

            case (OPR_P1_INT_VP)
                call OPR_PARTIAL1_INT(0, nxy, g, p_a, p_b)

            case (OPR_P1_INT_PV)
                call OPR_PARTIAL1_INT(1, nxy, g, p_a, p_b)

            case (OPR_P0_IBM)
                call OPR_IBM(nxy, g, p_a, p_b)

            end select

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TransposeK_Backward(p_b, result, id)
                if (type == OPR_P2_P1) then
                    call TLabMPI_TransposeK_Backward(p_c, tmp1, id)
                end if
            end if
#endif

            nullify (p_a, p_b)
            if (associated(p_c)) nullify (p_c)

        end if

        return
    end subroutine OPR_PARTIAL_Z

!########################################################################
!########################################################################
    subroutine OPR_PARTIAL_Y(type, nx, ny, nz, bcs, g, u, result, tmp1)
        use TLab_Arrays, only: wrk3d
        integer(wi), intent(in) :: type     ! OPR_P1           1.order derivative
        !                                   OPR_P2           2.order derivative
        !                                   OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
        !                                   OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
        !                                   OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(:, :)       ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nxy, nxz
        integer(wi), parameter :: is = -1 ! second derivative without viscosity/diffusivity
        real(wp), dimension(:), pointer :: p_a, p_b, p_c

! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp

        else
! ###################################################################
            nxy = nx*ny
            nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y direction the last one
! -------------------------------------------------------------------
            if (nz > 1) then
#ifdef USE_ESSL
                call DGETMO(u, nxy, nxy, nz, result, nz)
#else
                call TLab_Transpose(u, nxy, nz, nxy, result, nz)
#endif
                p_a => result
                if (any([OPR_P2, OPR_P2_P1] == type)) then
                    p_b => tmp1
                    p_c => wrk3d
                else
                    p_b => wrk3d
                end if

            else
                p_a => u
                p_b => result
                if (any([OPR_P2, OPR_P2_P1] == type)) then
                    p_c => tmp1
                end if

            end if

! ###################################################################
            select case (type)

            case (OPR_P2)
                call OPR_PARTIAL2(is, nxz, bcs, g, p_a, p_b, p_c)

            case (OPR_P2_P1)
                call OPR_PARTIAL2(is, nxz, bcs, g, p_a, p_b, p_c)

                if (.not. g%need_1der) then  ! Check whether we need to calculate the 1. order derivative
                    call OPR_PARTIAL1(nxz, bcs(:, 1), g, p_a, p_c)
                end if

            case (OPR_P1)
                if (ibm_partial) then
                    call OPR_PARTIAL1_IBM(nxz, bcs(:, 1), g, p_a, p_b)
                else
                    call OPR_PARTIAL1(nxz, bcs(:, 1), g, p_a, p_b)
                end if

            case (OPR_P0_INT_VP)
                call OPR_PARTIAL0_INT(0, nxz, g, p_a, p_b)

            case (OPR_P0_INT_PV)
                call OPR_PARTIAL0_INT(1, nxz, g, p_a, p_b)

            case (OPR_P1_INT_VP)
                call OPR_PARTIAL1_INT(0, nxz, g, p_a, p_b)

            case (OPR_P1_INT_PV)
                call OPR_PARTIAL1_INT(1, nxz, g, p_a, p_b)

            case (OPR_P0_IBM)
                call OPR_IBM(nxz, g, p_a, p_b)

            end select

! ###################################################################
! Put arrays back in the order in which they came in
            if (nz > 1) then
#ifdef USE_ESSL
                call DGETMO(p_b, nz, nz, nxy, result, nxy)
#else
                call TLab_Transpose(p_b, nz, nxy, nz, result, nxy)
#endif
                if (type == OPR_P2_P1) then
#ifdef USE_ESSL
                    call DGETMO(p_c, nz, nz, nxy, tmp1, nxy)
#else
                    call TLab_Transpose(p_c, nz, nxy, nz, tmp1, nxy)
#endif
                end if
            end if

            nullify (p_a, p_b)
            if (associated(p_c)) nullify (p_c)

        end if

        return
    end subroutine OPR_PARTIAL_Y

end module OPR_PARTIAL

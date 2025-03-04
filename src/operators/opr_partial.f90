#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI

#endif

module OPR_PARTIAL
    use TLab_Constants, only: efile, wp, wi, BCS_DN, BCS_ND, BCS_NN
    use FDM, only: fdm_dt, FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM4_DIRECT, FDM_COM6_DIRECT, FDM_COM6_JACOBIAN_PENTA
    use FDM, only: FDM_Der2_Solve
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
    use TLabMPI_Transpose
#endif
    use TLab_Arrays, only: wrk2d
    use FDM_MatMul
    implicit none
    private

    integer(wi) ip, ibc

    public :: OPR_PARTIAL_X
    public :: OPR_PARTIAL_Y
    public :: OPR_PARTIAL_Z
    public :: OPR_PARTIAL1 !, OPR_PARTIAL2

contains
! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL1(nlines, bcs, g, u, result)
        integer(wi), intent(in) :: nlines   ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2)   ! BCs at xmin (1) and xmax (2):
        !                                   0 biased, non-zero
        !                                   1 forced to zero
        type(fdm_dt), intent(in) :: g
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
        type(fdm_dt), intent(in) :: g
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
    subroutine OPR_IBM(is, nlines, g, u, result)
        integer(wi), intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines*g%size)
        real(wp), intent(out) :: result(nlines*g%size)

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

! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL0_INT(dir, nlines, g, u, result)
        use TLab_Arrays, only: wrk2d
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(fdm_dt), intent(in) :: g
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
        type(fdm_dt), intent(in) :: g
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
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nyz

        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d

! #ifdef USE_MPI
!         integer(wi), parameter :: id = TLAB_MPI_TRP_I_PARTIAL
! #endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(u, result, tmpi_plan_dx)
            p_a => result
            p_b => wrk3d
            p_c => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_d => tmp1
            end if
            ! nyz = ims_size_i(id)
            nyz = tmpi_plan_dx%nlines
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
            if (g%need_1der) call OPR_PARTIAL1(nyz, bcs(:, 1), g, p_b, p_d)
            call FDM_Der2_Solve(nyz, g, g%lu2, p_b, p_c, p_d, wrk2d)

        case (OPR_P2_P1)
            call OPR_PARTIAL1(nyz, bcs(:, 1), g, p_b, p_d)
            call FDM_Der2_Solve(nyz, g, g%lu2, p_b, p_c, p_d, wrk2d)

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
            call OPR_IBM(0, nyz, g, p_b, p_c)

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
                call TLabMPI_TransposeI_Backward(p_c, tmp1, tmpi_plan_dx)
            end if
            call TLabMPI_TransposeI_Backward(p_b, result, tmpi_plan_dx)
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
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nxy

        real(wp), dimension(:), pointer :: p_a, p_b, p_c

! #ifdef USE_MPI
!         integer(wi), parameter :: id = TLAB_MPI_TRP_K_PARTIAL
! #endif

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
                call TLabMPI_TransposeK_Forward(u, result, tmpi_plan_dz)
                p_a => result
                if (any([OPR_P2, OPR_P2_P1] == type)) then
                    p_b => tmp1
                    p_c => wrk3d
                else
                    p_b => wrk3d
                end if
                ! nxy = ims_size_k(id)
                nxy = tmpi_plan_dz%nlines
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
                if (g%need_1der) call OPR_PARTIAL1(nxy, bcs(:, 1), g, p_a, p_c)
                call FDM_Der2_Solve(nxy, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P2_P1)
                call OPR_PARTIAL1(nxy, bcs(:, 1), g, p_a, p_c)
                call FDM_Der2_Solve(nxy, g, g%lu2, p_a, p_b, p_c, wrk2d)

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
                call OPR_IBM(0, nxy, g, p_a, p_b)

            end select

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TransposeK_Backward(p_b, result, tmpi_plan_dz)
                if (type == OPR_P2_P1) then
                    call TLabMPI_TransposeK_Backward(p_c, tmp1, tmpi_plan_dz)
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
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)

        target u, tmp1, result

! -------------------------------------------------------------------
        integer(wi) nxy, nxz
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
                if (g%need_1der) call OPR_PARTIAL1(nxz, bcs(:, 1), g, p_a, p_c)
                call FDM_Der2_Solve(nxz, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P2_P1)
                call OPR_PARTIAL1(nxz, bcs(:, 1), g, p_a, p_c)
                call FDM_Der2_Solve(nxz, g, g%lu2, p_a, p_b, p_c, wrk2d)

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
                call OPR_IBM(0, nxz, g, p_a, p_b)

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

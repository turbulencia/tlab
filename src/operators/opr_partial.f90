#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI

#endif

module OPR_PARTIAL
    use TLab_Constants, only: efile, wp, wi, BCS_DN, BCS_ND, BCS_NN
    use FDM, only: fdm_dt
    use FDM, only: FDM_Der1_Solve, FDM_Der2_Solve
    use FDM, only: FDM_Interpol, FDM_Interpol_Der1
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
    implicit none
    private

    public :: OPR_PARTIAL_X
    public :: OPR_PARTIAL_Y
    public :: OPR_PARTIAL_Z

contains
! ###################################################################
! ###################################################################
    subroutine OPR_PARTIAL1_IBM(nlines, bcs, g, u, result)
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
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, fld_ibm, result, wrk2d)  ! now with modified u fields
            else ! idle IBM-Tasks
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, u, result, wrk2d)  ! no splines needed
            end if

        case ('y')
            if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, fld_ibm, result, wrk2d)  ! now with modified u fields
            else ! idle IBM-Tasks
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, u, result, wrk2d)  ! no splines needed
            end if

        case ('z')
            if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, fld_ibm, result, wrk2d)  ! now with modified u fields
            else ! idle IBM-Tasks
                call FDM_Der1_Solve(nlines, bcs, g, g%lu1, u, result, wrk2d)  ! no splines needed
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
            if (g%need_1der) call FDM_Der1_Solve(nyz, bcs(:, 1), g, g%lu1, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nyz, g, g%lu2, p_b, p_c, p_d, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nyz, bcs(:, 1), g, g%lu1, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nyz, g, g%lu2, p_b, p_c, p_d, wrk2d)

        case (OPR_P1)
            if (ibm_partial) then
                call OPR_PARTIAL1_IBM(nyz, bcs(:, 1), g, p_b, p_c)
            else
                call FDM_Der1_Solve(nyz, bcs(:, 1), g, g%lu1, p_b, p_c, wrk2d)
            end if

        case (OPR_P0_INT_VP)
            call FDM_Interpol(0, nyz, g, p_b, p_c, wrk2d)

        case (OPR_P0_INT_PV)
            call FDM_Interpol(1, nyz, g, p_b, p_c, wrk2d)

        case (OPR_P1_INT_VP)
            call FDM_Interpol_Der1(0, nyz, g, p_b, p_c, wrk2d)

        case (OPR_P1_INT_PV)
            call FDM_Interpol_Der1(1, nyz, g, p_b, p_c, wrk2d)

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
                if (g%need_1der) call FDM_Der1_Solve(nxy, bcs(:, 1), g, g%lu1, p_a, p_c, wrk2d)
                call FDM_Der2_Solve(nxy, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P2_P1)
                call FDM_Der1_Solve(nxy, bcs(:, 1), g, g%lu1, p_a, p_c, wrk2d)
                call FDM_Der2_Solve(nxy, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P1)
                if (ibm_partial) then
                    call OPR_PARTIAL1_IBM(nxy, bcs(:, 1), g, p_a, p_b)
                else
                    call FDM_Der1_Solve(nxy, bcs(:, 1), g, g%lu1, p_a, p_b, wrk2d)
                end if

            case (OPR_P0_INT_VP)
                call FDM_Interpol(0, nxy, g, p_a, p_b, wrk2d)

            case (OPR_P0_INT_PV)
                call FDM_Interpol(1, nxy, g, p_a, p_b, wrk2d)

            case (OPR_P1_INT_VP)
                call FDM_Interpol_Der1(0, nxy, g, p_a, p_b, wrk2d)

            case (OPR_P1_INT_PV)
                call FDM_Interpol_Der1(1, nxy, g, p_a, p_b, wrk2d)

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
                if (g%need_1der) call FDM_Der1_Solve(nxz, bcs(:, 1), g, g%lu1, p_a, p_c, wrk2d)
                call FDM_Der2_Solve(nxz, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P2_P1)
                call FDM_Der1_Solve(nxz, bcs(:, 1), g, g%lu1, p_a, p_c, wrk2d)
                call FDM_Der2_Solve(nxz, g, g%lu2, p_a, p_b, p_c, wrk2d)

            case (OPR_P1)
                if (ibm_partial) then
                    call OPR_PARTIAL1_IBM(nxz, bcs(:, 1), g, p_a, p_b)
                else
                    call FDM_Der1_Solve(nxz, bcs(:, 1), g, g%lu1, p_a, p_b, wrk2d)
                end if

            case (OPR_P0_INT_VP)
                call FDM_Interpol(0, nxz, g, p_a, p_b, wrk2d)

            case (OPR_P0_INT_PV)
                call FDM_Interpol(1, nxz, g, p_a, p_b, wrk2d)

            case (OPR_P1_INT_VP)
                call FDM_Interpol_Der1(0, nxz, g, p_a, p_b, wrk2d)

            case (OPR_P1_INT_PV)
                call FDM_Interpol_Der1(1, nxz, g, p_a, p_b, wrk2d)

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

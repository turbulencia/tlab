#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI

#endif

module OPR_Burgers
    use TLab_Constants, only: wp, wi, efile, lfile
    use FDM, only: fdm_dt, g
    use IBM_VARS, only: ibm_burgers
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_Transpose
#endif
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk2d, wrk3d
    use NavierStokes, only: visc, schmidt
    use OPR_FILTERS
    use OPR_PARTIAL
    use LargeScaleForcing, only: subsidenceProps, TYPE_SUB_CONSTANT_LOCAL
    implicit none
    private

    public :: OPR_Burgers_Initialize
    public :: OPR_Burgers_X
    public :: OPR_Burgers_Y
    public :: OPR_Burgers_Z

    ! Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
    ! the argument ivel indicates 2 options:
    integer, parameter, public :: OPR_B_SELF = 0        ! velocity component is the scalar itself, the transposed velocity is returned
    integer, parameter, public :: OPR_B_U_IN = 1        ! velocity component is passed through u, or u_t if transposed required

    type(filter_dt) :: Dealiasing(3)
    real(wp), allocatable, target :: wrkdea(:, :)       ! Work arrays for dealiasing (scratch space)

    type :: rho_anelastic                               ! 1/rho in diffusion term in anelastic formulation
        sequence
        logical :: active = .false.
        real(wp), allocatable :: values(:)
    end type rho_anelastic
    type(rho_anelastic) :: rhoinv(3)                    ! one for each direction

    type :: fdm_diffusion_dt
        sequence
        real(wp), allocatable :: lu(:, :, :)
    end type fdm_diffusion_dt
    type(fdm_diffusion_dt) :: fdmDiffusion(3)

contains
    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Initialize(inifile)
        use TLab_Memory, only: isize_field, imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        use NavierStokes, only: nse_eqns
        use TLab_Memory, only: inb_scal
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_k, ims_npro_k
#endif
        use THERMO_ANELASTIC, only: rbackground, ribackground

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) ig, is, j, ip, nlines, offset, idummy
        real(wp) dummy
        character*32 bakfile

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        call FILTER_READBLOCK(bakfile, inifile, 'Dealiasing', Dealiasing)

        ! ###################################################################
        ! Initialize LU factorization of the second-order derivative times the diffusivities
        do ig = 1, 3
            if (g(ig)%size == 1) cycle

            if (g(ig)%der2%nb_diag(1) /= 3) then
                call TLab_Write_ASCII(efile, __FILE__//'. Undeveloped for more than 3 LHS diagonals in 2. order derivatives.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            if (g(ig)%periodic) then
                idummy = g(ig)%der2%nb_diag(1) + 2
            else
                idummy = g(ig)%der2%nb_diag(1)
            end if
            allocate (fdmDiffusion(ig)%lu(g(ig)%size, idummy, 0:inb_scal))

            ip = 0
            do is = 0, inb_scal ! case 0 for the reynolds number
                if (is == 0) then
                    dummy = visc
                else
                    dummy = visc/schmidt(is)
                end if

                if (g(ig)%periodic) then                        ! Check routines TRIDPFS and TRIDPSS
                    fdmDiffusion(ig)%lu(:, 1, is) = g(ig)%der2%lu(:, 1)         ! matrix L; 1. subdiagonal
                    fdmDiffusion(ig)%lu(:, 2, is) = g(ig)%der2%lu(:, 2)*dummy   ! matrix L; 1/diagonal
                    fdmDiffusion(ig)%lu(:, 3, is) = g(ig)%der2%lu(:, 3)         ! matrix U is the same
                    fdmDiffusion(ig)%lu(:, 4, is) = g(ig)%der2%lu(:, 4)/dummy   ! matrix L; Additional row/column
                    fdmDiffusion(ig)%lu(:, 5, is) = g(ig)%der2%lu(:, 5)         ! matrix U is the same

                else                                            ! Check routines TRIDFS and TRIDSS
                    fdmDiffusion(ig)%lu(:, 1, is) = g(ig)%der2%lu(:, 1)         ! matrix L is the same
                    fdmDiffusion(ig)%lu(:, 2, is) = g(ig)%der2%lu(:, 2)*dummy   ! matrix U; 1/diagonal
                    fdmDiffusion(ig)%lu(:, 3, is) = g(ig)%der2%lu(:, 3)/dummy   ! matrix U; 1. superdiagonal

                end if

            end do
        end do

        ! ###################################################################
        ! Initialize dealiasing
        do ig = 1, 3
            if (Dealiasing(ig)%type /= DNS_FILTER_NONE) call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
        end do

        if (any(Dealiasing(:)%type /= DNS_FILTER_NONE)) then
            call TLab_Allocate_Real(__FILE__, wrkdea, [isize_field, 2], 'wrk-dealiasing')
        end if

        ! ###################################################################
        ! Initialize anelastic density correction
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
            rhoinv(1)%active = .true.
#ifdef USE_MPI
            if (ims_npro_i > 1) then
                ! nlines = ims_size_i(TLAB_MPI_TRP_I_PARTIAL)
                nlines = tmpi_plan_dx%nlines
                offset = nlines*ims_pro_i
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rhoinv(1)%values(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, g(2)%size) + 1
                rhoinv(1)%values(j) = ribackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y; see FDM_Initialize
            ! we implement it directly in the tridiagonal system
            do is = 0, inb_scal ! case 0 for the velocity
                fdmDiffusion(2)%lu(:, 2, is) = fdmDiffusion(2)%lu(:, 2, is)*ribackground(:)  ! matrix U; 1/diagonal
                fdmDiffusion(2)%lu(:g(2)%size - 1, 3, is) = fdmDiffusion(2)%lu(:g(2)%size - 1, 3, is)*rbackground(2:) ! matrix U; 1. superdiagonal
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Z
            ! g(3)%anelastic = .true.
            rhoinv(3)%active = .true.
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                ! nlines = ims_size_k(TLAB_MPI_TRP_K_PARTIAL)
                nlines = tmpi_plan_dz%nlines
                offset = nlines*ims_pro_k
            else
#endif
                nlines = imax*jmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rhoinv(3)%values(nlines))
            do j = 1, nlines
                ip = (offset + j - 1)/imax + 1
                rhoinv(3)%values(j) = ribackground(ip)
            end do

        end if

        return
    end subroutine OPR_Burgers_Initialize

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_X(ivel, is, nx, ny, nz, bcs, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2)                ! BCs at xmin (1,*) and xmax (2,*)
        real(wp), intent(in) :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)      ! transposed velocity
        real(wp), intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

        ! -------------------------------------------------------------------
        integer(wi) nyz
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d, p_vel

        ! ###################################################################
        if (g(1)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! ###################################################################
        ! -------------------------------------------------------------------
        ! MPI transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(s, result, tmpi_plan_dx)
            p_a => result
            p_b => tmp1
            p_c => wrk3d
            p_d => result
            nyz = tmpi_plan_dx%nlines
        else
#endif
            p_a => s
            p_b => tmp1
            p_c => result
            p_d => wrk3d
            nyz = ny*nz
#ifdef USE_MPI
        end if
#endif

        ! pointer to velocity
        if (ivel == OPR_B_SELF) then    ! velocity is the scalar itself
            p_vel => p_b
        else                            ! transposed velocity is passed through argument
            p_vel => u_t
        end if

        ! maybe check that nx is equal to g(1)%size

        ! -------------------------------------------------------------------
        ! Local transposition: make x-direction the last one
        ! -------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(p_a, g(1)%size, g(1)%size, nyz, p_b, nyz)
#else
        call TLab_Transpose(p_a, g(1)%size, nyz, g(1)%size, p_b, nyz)
#endif

        ! ###################################################################
        call OPR_Burgers_1D(is, nyz, bcs, g(1), fdmDiffusion(1)%lu(:, :, is), Dealiasing(1), rhoinv(1), p_b, p_vel, p_d, p_c)

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_d, nyz, nyz, g(1)%size, p_c, g(1)%size)
#else
        call TLab_Transpose(p_d, nyz, g(1)%size, nyz, p_c, g(1)%size)
#endif

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Backward(p_c, result, tmpi_plan_dx)
        end if
#endif

        nullify (p_a, p_b, p_c, p_d, p_vel)

        return
    end subroutine OPR_Burgers_X

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Y(ivel, is, nx, ny, nz, bcs, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        real(wp), intent(in) :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)      ! transposed velocity
        real(wp), intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

        ! -------------------------------------------------------------------
        integer(wi) nxy, nxz, j
        real(wp), dimension(:), pointer :: p_org, p_vel
        real(wp), dimension(:, :), pointer :: p_dst1, p_dst2    ! need (nx*nz,ny) shape

        ! ###################################################################
        if (g(2)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! ###################################################################
        nxy = nx*ny
        nxz = nx*nz

        ! -------------------------------------------------------------------
        ! Local transposition: Make y-direction the last one
        ! -------------------------------------------------------------------
        if (nz == 1) then
            p_org => s
            p_dst1(1:nx*nz, 1:ny) => wrk3d(1:nx*nz*ny)
            p_dst2(1:nx*nz, 1:ny) => result(1:nx*nz*ny)
        else
#ifdef USE_ESSL
            call DGETMO(s, nxy, nxy, nz, tmp1, nz)
#else
            call TLab_Transpose(s, nxy, nz, nxy, tmp1, nz)
#endif
            p_org => tmp1
            p_dst1(1:nx*nz, 1:ny) => result(1:nx*nz*ny)
            p_dst2(1:nx*nz, 1:ny) => wrk3d(1:nx*nz*ny)
        end if

        ! pointer to velocity
        if (ivel == OPR_B_SELF) then
            p_vel => p_org
        else
            if (nz == 1) then  ! I do not need the transposed
                p_vel => u
            else               ! I do need the transposed
                p_vel => u_t
            end if
        end if

        ! ###################################################################
        call OPR_Burgers_1D(is, nxz, bcs, g(2), fdmDiffusion(2)%lu(:, :, is), Dealiasing(2), rhoinv(2), p_org, p_vel, p_dst2, p_dst1)

        if (subsidenceProps%type == TYPE_SUB_CONSTANT_LOCAL) then
            do j = 1, ny
                p_dst2(:, j) = p_dst2(:, j) + g(2)%nodes(j)*subsidenceProps%parameters(1)*p_dst1(:, j)
            end do
        end if

        ! ###################################################################
        ! Put arrays back in the order in which they came in
        if (nz > 1) then
#ifdef USE_ESSL
            call DGETMO(p_dst2, nz, nz, nxy, result, nxy)
#else
            call TLab_Transpose(p_dst2, nz, nxy, nz, result, nxy)
#endif
        end if

        nullify (p_org, p_dst1, p_dst2, p_vel)

        return
    end subroutine OPR_Burgers_Y

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Z(ivel, is, nx, ny, nz, bcs, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        real(wp), intent(in) :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

        ! -------------------------------------------------------------------
        integer(wi) nxy
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_vel

        ! ###################################################################
        if (g(3)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! ###################################################################
        ! -------------------------------------------------------------------
        ! MPI transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_TransposeK_Forward(s, tmp1, tmpi_plan_dz)
            p_a => tmp1
            p_b => result
            p_c => wrk3d
            nxy = tmpi_plan_dz%nlines
        else
#endif
            p_a => s
            p_b => wrk3d
            p_c => result
            nxy = nx*ny
#ifdef USE_MPI
        end if
#endif

        ! pointer to velocity
        if (ivel == OPR_B_SELF) then
            p_vel => p_a
        else
#ifdef USE_MPI
            if (ims_npro_k > 1) then        ! I do need the transposed
                p_vel => u_t
            else
#endif
                p_vel => u                 ! I do not need the transposed
#ifdef USE_MPI
            end if
#endif
        end if

        ! ###################################################################
        call OPR_Burgers_1D(is, nxy, bcs, g(3), fdmDiffusion(3)%lu(:, :, is), Dealiasing(3), rhoinv(3), p_a, p_vel, p_c, p_b)

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_TransposeK_Backward(p_c, result, tmpi_plan_dz)
        end if
#endif

        nullify (p_a, p_b, p_c, p_vel)

        return
    end subroutine OPR_Burgers_Z

    !########################################################################
    !# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
    !# along generic direction x to nlines lines of data
    !#
    !# Second derivative uses LE decomposition including diffusivity coefficient
    !########################################################################
    subroutine OPR_Burgers_1D(is, nlines, bcs, g, lu2d, dealiasing, rhoinv, s, u, result, dsdx)
        use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)    ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: lu2d(:, :)      ! LU decomposition including the diffusion parameter for corresponding field is
        type(filter_dt), intent(in) :: dealiasing
        type(rho_anelastic), intent(in) :: rhoinv
        real(wp), intent(in) :: s(nlines, g%size), u(nlines, g%size)  ! argument field and velocity field
        real(wp), intent(out) :: result(nlines, g%size)                ! N(u) applied to s
        real(wp), intent(inout) :: dsdx(nlines, g%size)                  ! dsdx

        ! -------------------------------------------------------------------
        integer(wi) ij
        real(wp), pointer :: uf(:, :), dsf(:, :)

        ! ###################################################################
        if (bcs(1, 2) + bcs(2, 2) > 0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Only developed for biased BCs.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! dsdx: 1st derivative; result: 2nd derivative including diffusivity
        if (ibm_burgers) then
            call OPR_PARTIAL2_IBM(is, nlines, bcs, g, lu2d, s, result, dsdx)
        else
            call FDM_Der1_Solve(nlines, bcs(:, 1), g%der1, g%der1%lu, s, dsdx, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, lu2d, s, result, dsdx, wrk2d)
        end if

        ! ###################################################################
        ! Operation; diffusivity included in 2.-order derivativelu2_p
        ! ###################################################################
        if (dealiasing%type /= DNS_FILTER_NONE) then
            uf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 1)
            dsf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 2)
            call OPR_FILTER_1D(nlines, dealiasing, u, uf)
            call OPR_FILTER_1D(nlines, dealiasing, dsdx, dsf)

            ! We duplicate a few lines of code instead of using pointers becasue
            ! creating pointers take some running time
            if (rhoinv%active) then
                do ij = 1, g%size
                    result(:, ij) = result(:, ij)*rhoinv%values(:) - uf(:, ij)*dsf(:, ij)
                end do

            else
!$omp parallel default( shared ) private( ij )
!$omp do
                do ij = 1, nlines*g%size
                    result(ij, 1) = result(ij, 1) - uf(ij, 1)*dsf(ij, 1)
                end do
!$omp end do
!$omp end parallel
            end if

            nullify (uf, dsf)

        else
            if (rhoinv%active) then
                do ij = 1, g%size
                    result(:, ij) = result(:, ij)*rhoinv%values(:) - u(:, ij)*dsdx(:, ij)
                end do

            else
!$omp parallel default( shared ) private( ij )
!$omp do
                do ij = 1, nlines*g%size
                    result(ij, 1) = result(ij, 1) - u(ij, 1)*dsdx(ij, 1)
                end do
!$omp end do
!$omp end parallel
            end if
        end if

        return
    end subroutine OPR_Burgers_1D

! ###################################################################
! ###################################################################
    ! modify incoming fields (fill solids with spline functions, depending on direction)
    subroutine OPR_PARTIAL2_IBM(is, nlines, bcs, g, lu2, u, result, du)
        use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
        use IBM_VARS
        integer(wi), intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)     ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: lu2(:, :)
        real(wp), intent(in) :: u(nlines*g%size)
        real(wp), intent(out) :: result(nlines*g%size)
        real(wp), intent(inout) :: du(nlines*g%size)  ! First derivative

        target u

        ! -------------------------------------------------------------------
        real(wp), pointer :: p_fld(:)

        ! -------------------------------------------------------------------
        select case (g%name)
        case ('x')
            if (ims_pro_ibm_x) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
                p_fld => fld_ibm
            else ! idle IBM-Tasks
                p_fld => u
            end if

        case ('y')
            if (ims_pro_ibm_y) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
                p_fld => fld_ibm
            else ! idle IBM-Tasks
                p_fld => u
            end if

        case ('z')
            if (ims_pro_ibm_z) then ! only active IBM-Tasks (with objects in their subdomain) enter IBM-routines
                call IBM_SPLINE_XYZ(is, u, fld_ibm, g, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
                p_fld => fld_ibm
            else ! idle IBM-Tasks
                p_fld => u
            end if

        end select

        call FDM_Der1_Solve(nlines, bcs(:, 1), g%der1, g%der1%lu, p_fld, du, wrk2d)
        call FDM_Der2_Solve(nlines, g%der2, lu2, p_fld, result, du, wrk2d)  ! no splines needed

        nullify (p_fld)

        return
    end subroutine OPR_PARTIAL2_IBM

end module OPR_Burgers

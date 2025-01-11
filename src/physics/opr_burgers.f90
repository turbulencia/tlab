#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module OPR_Burgers
    use TLab_Constants, only: wp, wi, efile, lfile
    use FDM, only: grid_dt
    use IBM_VARS, only: ibm_burgers
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_Transpose
#endif
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk2d, wrk3d
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
    integer, parameter, public :: OPR_B_SELF = 0 ! velocity component is the scalar itself, the transposed velocity is returned
    integer, parameter, public :: OPR_B_U_IN = 1 ! velocity component is passed through u, or u_t if transposed required

    type(filter_dt) :: Dealiasing(3)

    real(wp), allocatable, target :: wrkdea(:, :)                   ! Work arrays for dealiasing (scratch space)

contains
    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Initialize(inifile)
        use TLAB_VARS, only: isize_field, imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        use TLAB_VARS, only: imode_eqns, inb_scal
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_k, ims_npro_k
        ! use TLabMPI_Transpose, only: ims_size_i, ims_size_k
        use TLabMPI_Transpose, only: ims_trp_plan_i, ims_trp_plan_k
#endif
        use FDM, only: g
        use THERMO_ANELASTIC, only: rbackground, ribackground

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) ig, is, j, ip, nlines, offset

        character*32 bakfile

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        call FILTER_READBLOCK(bakfile, inifile, 'Dealiasing', Dealiasing)

        ! ###################################################################
        ! Initialize dealiasing
        do ig = 1, 3
            if (Dealiasing(ig)%type /= DNS_FILTER_NONE) call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
        end do

        if (any(Dealiasing(:)%type /= DNS_FILTER_NONE)) then ! to be moved to OPR_Burgers_Initialize
            call TLab_Allocate_Real(__FILE__, wrkdea, [isize_field, 2], 'wrk-dealiasing')
        end if

        ! ###################################################################
        ! Initialize anelastic density correction
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
            g(1)%anelastic = .true.
#ifdef USE_MPI
            if (ims_npro_i > 1) then
                ! nlines = ims_size_i(TLAB_MPI_TRP_I_PARTIAL)
                nlines = ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL)%nlines
                offset = nlines*ims_pro_i
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (g(1)%rhoinv(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, g(2)%size) + 1
                g(1)%rhoinv(j) = ribackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y; see FDM_Initialize
            ! we implement it directly in the tridiagonal system
            ip = 0
            do is = 0, inb_scal ! case 0 for the velocity
                g(2)%lu2d(:, ip + 2) = g(2)%lu2d(:, ip + 2)*ribackground(:)  ! matrix U; 1/diagonal
                g(2)%lu2d(:g(2)%size - 1, ip + 3) = g(2)%lu2d(:, ip + 3)*rbackground(2:) ! matrix U; 1. superdiagonal
                ip = ip + 3
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Z
            g(3)%anelastic = .true.
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                ! nlines = ims_size_k(TLAB_MPI_TRP_K_PARTIAL)
                nlines = ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL)%nlines
                offset = nlines*ims_pro_k
            else
#endif
                nlines = imax*jmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (g(3)%rhoinv(nlines))
            do j = 1, nlines
                ip = (offset + j - 1)/imax + 1
                g(3)%rhoinv(j) = ribackground(ip)
            end do

        end if

        return
    end subroutine OPR_Burgers_Initialize

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_X(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2)                ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)      ! transposed velocity
        real(wp), intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

        ! -------------------------------------------------------------------
        integer(wi) nyz
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d, p_vel
! #ifdef USE_MPI
!         integer(wi), parameter :: id = TLAB_MPI_TRP_I_PARTIAL
! #endif

        ! ###################################################################
        ! -------------------------------------------------------------------
        ! MPI transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(s, result, ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL))
            p_a => result
            p_b => tmp1
            p_c => wrk3d
            p_d => result
            ! nyz = ims_size_i(id)
            nyz = ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL)%nlines
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

        ! -------------------------------------------------------------------
        ! Local transposition: make x-direction the last one
        ! -------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(p_a, g%size, g%size, nyz, p_b, nyz)
#else
        call TLab_Transpose(p_a, g%size, nyz, g%size, p_b, nyz)
#endif

        ! ###################################################################
        call OPR_Burgers_1D(is, nyz, bcs, g, Dealiasing(1), p_b, p_vel, p_d, p_c)

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
        call TLab_Transpose(p_d, nyz, g%size, nyz, p_c, g%size)
#endif

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Backward(p_c, result, ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL))
        end if
#endif

        nullify (p_a, p_b, p_c, p_d, p_vel)

        return
    end subroutine OPR_Burgers_X

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Y(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
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
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp

        else
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
            call OPR_Burgers_1D(is, nxz, bcs, g, Dealiasing(2), p_org, p_vel, p_dst2, p_dst1)

            if (subsidenceProps%type == TYPE_SUB_CONSTANT_LOCAL) then
                do j = 1, ny
                    p_dst2(:, j) = p_dst2(:, j) + g%nodes(j)*subsidenceProps%parameters(1)*p_dst1(:, j)
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

        end if

        return
    end subroutine OPR_Burgers_Y

    !########################################################################
    !########################################################################
    subroutine OPR_Burgers_Z(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer, intent(in) :: ivel
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

        ! -------------------------------------------------------------------
        integer(wi) nxy
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_vel
! #ifdef USE_MPI
!         integer(wi), parameter :: id = TLAB_MPI_TRP_K_PARTIAL
! #endif

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp

        else
            ! ###################################################################
            ! -------------------------------------------------------------------
            ! MPI transposition
            ! -------------------------------------------------------------------
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TransposeK_Forward(s, tmp1, ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL))
                p_a => tmp1
                p_b => result
                p_c => wrk3d
                ! nxy = ims_size_k(id)
                nxy = ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL)%nlines
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
            call OPR_Burgers_1D(is, nxy, bcs, g, Dealiasing(3), p_a, p_vel, p_c, p_b)

            ! ###################################################################
            ! Put arrays back in the order in which they came in
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TransposeK_Backward(p_c, result, ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL))
            end if
#endif

            nullify (p_a, p_b, p_c, p_vel)

        end if

        return
    end subroutine OPR_Burgers_Z

    !########################################################################
    !# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
    !# along generic direction x to nlines lines of data
    !#
    !# Second derivative uses LE decomposition including diffusivity coefficient
    !########################################################################
    subroutine OPR_Burgers_1D(is, nlines, bcs, g, dealiasing, s, u, result, dsdx)
        integer, intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)    ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(grid_dt), intent(in) :: g
        type(filter_dt), intent(in) :: dealiasing
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
            call OPR_PARTIAL2_IBM(is, nlines, bcs, g, s, result, dsdx)
        else
            call OPR_PARTIAL2(is, nlines, bcs, g, s, result, dsdx)
        end if

        ! ###################################################################
        ! Operation; diffusivity included in 2.-order derivative
        ! ###################################################################
        if (dealiasing%type /= DNS_FILTER_NONE) then
            uf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 1)
            dsf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 2)
            call OPR_FILTER_1D(nlines, dealiasing, u, uf)
            call OPR_FILTER_1D(nlines, dealiasing, dsdx, dsf)

            ! We duplicate a few lines of code instead of using pointers becasue
            ! creating pointers take some running time
            if (g%anelastic) then
                do ij = 1, g%size
                    result(:, ij) = result(:, ij)*g%rhoinv(:) - uf(:, ij)*dsf(:, ij)
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
            if (g%anelastic) then
                do ij = 1, g%size
                    result(:, ij) = result(:, ij)*g%rhoinv(:) - u(:, ij)*dsdx(:, ij)
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

end module OPR_Burgers

#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module OPR_BURGERS
    use TLab_Constants, only: efile, wp, wi
    use TLab_Types, only: grid_dt, filter_dt
    use IBM_VARS, only: ibm_burgers
    use TLAB_VARS, only: Dealiasing, subsidence
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i
    use TLabMPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLabMPI_VARS, only: ims_npro_k
    use TLabMPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLabMPI_PROCS
#endif
    use TLab_WorkFlow
    use TLab_Arrays, only: wrk2d, wrk3d
    use OPR_FILTERS
    use OPR_PARTIAL
    implicit none
    private

    ! Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
    ! the argument ivel indicates 2 options:
    integer, parameter :: OPR_B_SELF = 0 ! velocity component is the scalar itself, the transposed velocity is returned
    integer, parameter :: OPR_B_U_IN = 1 ! velocity component is passed through u, or u_t if transposed required

    public :: OPR_B_SELF, OPR_B_U_IN
    public :: OPR_BURGERS_X
    public :: OPR_BURGERS_Y
    public :: OPR_BURGERS_Z

contains
!########################################################################
!########################################################################
    subroutine OPR_BURGERS_X(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer,     intent(in) :: ivel
        integer,     intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2)                ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in)    :: g
        real(wp),      intent(in)    :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp),      intent(out)   :: result(nx*ny*nz)
        real(wp),      intent(inout) :: tmp1(nx*ny*nz)      ! transposed velocity
        real(wp),      intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

! -------------------------------------------------------------------
        integer(wi) nyz
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d, p_vel
#ifdef USE_MPI
        integer(wi), parameter :: id = TLabMPI_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TRPF_I(s, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
            p_a => result
            p_b => tmp1
            p_c => wrk3d
            p_d => result
            nyz = ims_size_i(id)
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
        call OPR_BURGERS_1D(is, nyz, bcs, g, Dealiasing(1), p_b, p_vel, p_d, p_c)

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
        call TLab_Transpose(p_d, nyz, g%size, nyz, p_c, g%size)
#endif

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_TRPB_I(p_c, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        end if
#endif

        nullify (p_a, p_b, p_c, p_d, p_vel)

        return
    end subroutine OPR_BURGERS_X

!########################################################################
!########################################################################
    subroutine OPR_BURGERS_Y(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer,     intent(in) :: ivel
        integer,     intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in)    :: g
        real(wp),      intent(in)    :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp),      intent(out)   :: result(nx*ny*nz)
        real(wp),      intent(inout) :: tmp1(nx*ny*nz)      ! transposed velocity
        real(wp),      intent(in), optional :: u_t(nx*ny*nz)

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
            call OPR_BURGERS_1D(is, nxz, bcs, g, Dealiasing(2), p_org, p_vel, p_dst2, p_dst1)

            if (subsidence%type == EQNS_SUB_CONSTANT_LOCAL) then
                do j = 1, ny
                    p_dst2(:, j) = p_dst2(:, j) + g%nodes(j)*subsidence%parameters(1)*p_dst1(:, j)
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
    end subroutine OPR_BURGERS_Y

!########################################################################
!########################################################################
    subroutine OPR_BURGERS_Z(ivel, is, nx, ny, nz, bcs, g, s, u, result, tmp1, u_t)
        integer,     intent(in) :: ivel
        integer,     intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: bcs(2, 2) ! BCs at xmin (1,*) and xmax (2,*)
        type(grid_dt), intent(in)    :: g
        real(wp),      intent(in)    :: s(nx*ny*nz), u(nx*ny*nz)
        real(wp),      intent(out)   :: result(nx*ny*nz)
        real(wp),      intent(inout) :: tmp1(nx*ny*nz)
        real(wp),      intent(in), optional :: u_t(nx*ny*nz)

        target s, u, result, tmp1, u_t

! -------------------------------------------------------------------
        integer(wi) nxy
        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_vel
#ifdef USE_MPI
        integer(wi), parameter :: id = TLabMPI_K_PARTIAL
#endif

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
                call TLabMPI_TRPF_K(s, tmp1, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
                p_a => tmp1
                p_b => result
                p_c => wrk3d
                nxy = ims_size_k(id)
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
            call OPR_BURGERS_1D(is, nxy, bcs, g, Dealiasing(3), p_a, p_vel, p_c, p_b)

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
            if (ims_npro_k > 1) then
                call TLabMPI_TRPB_K(p_c, result, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            end if
#endif

            nullify (p_a, p_b, p_c, p_vel)

        end if

        return
    end subroutine OPR_BURGERS_Z

!########################################################################
!# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
!# along generic direction x to nlines lines of data
!#
!# Second derivative uses LE decomposition including diffusivity coefficient
!########################################################################
    subroutine OPR_BURGERS_1D(is, nlines, bcs, g, dealiasing, s, u, result, dsdx)
        use TLab_Arrays, only : wrkdea
        integer,     intent(in) :: is           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2, 2)    ! BCs at xmin (1,*) and xmax (2,*):
        !                                       0 biased, non-zero
        !                                       1 forced to zero
        type(grid_dt),   intent(in)    :: g
        type(filter_dt), intent(in)    :: dealiasing
        real(wp),        intent(in)    :: s(nlines, g%size), u(nlines, g%size)  ! argument field and velocity field
        real(wp),        intent(out)   :: result(nlines, g%size)                ! N(u) applied to s
        real(wp),        intent(inout) :: dsdx(nlines, g%size)                  ! dsdx

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
            uf(1:nlines,1:g%size) => wrkdea(1:nlines*g%size,1)
            dsf(1:nlines,1:g%size) => wrkdea(1:nlines*g%size,2)
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

            nullify(uf, dsf)

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
    end subroutine OPR_BURGERS_1D

end module OPR_BURGERS

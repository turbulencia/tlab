#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# DESCRIPTION
!#
!# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
!# along generic direction x to nlines lines of data
!#
!# Second derivative uses LE decomposition including diffusivity coefficient
!#
!########################################################################
subroutine OPR_BURGERS(is, nlines, bcs, g, dealiasing, s, u, result, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt, filter_dt
    use TLAB_CONSTANTS, only: efile
    use IBM_VARS, only: ibm_burgers
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_PROCS
    implicit none

    TINTEGER, intent(in) :: is     ! scalar index; if 0, then velocity
    TINTEGER, intent(in) :: nlines ! # of lines to be solved
    TINTEGER, dimension(2, *), intent(in) :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
    !     0 biased, non-zero
    !     1 forced to zero
    type(grid_dt), intent(in) :: g
    type(filter_dt), intent(in) :: dealiasing
    TREAL, dimension(nlines, g%size), intent(in) :: s, u    ! argument field and velocity field
    TREAL, dimension(nlines, g%size), intent(out) :: result ! N(u) applied to s
    TREAL, dimension(nlines, g%size), intent(inout) :: wrk3d  ! dsdx
    TREAL, dimension(*), intent(inout) :: wrk2d

! -------------------------------------------------------------------
    TINTEGER ij
    TREAL, dimension(:, :), allocatable :: uf, dsf

! ###################################################################
    if (bcs(1, 2) + bcs(2, 2) > 0) then
        call TLAB_WRITE_ASCII(efile, 'OPR_BURGERS. Only developed for biased BCs.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if

    ! wrk3d: 1st derivative; result: 2nd derivative including diffusivity
    if (ibm_burgers) then
        call OPR_PARTIAL2_IBM(is, nlines, bcs, g, s, result, wrk2d, wrk3d)
    else
        call OPR_PARTIAL2(is, nlines, bcs, g, s, result, wrk2d, wrk3d)
    end if

! ###################################################################
! Operation; diffusivity included in 2.-order derivative
! ###################################################################
    if (dealiasing%type /= DNS_FILTER_NONE) then
        allocate (uf(nlines, g%size), dsf(nlines, g%size))
        call OPR_FILTER_1D(nlines, dealiasing, u, uf, wrk1d, wrk2d, wrk3d) ! wrk3d is not used in compact filter
        call OPR_FILTER_1D(nlines, dealiasing, wrk3d, dsf, wrk1d, wrk2d, wrk3d)

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

        deallocate(uf,dsf)

    else
        if (g%anelastic) then
            do ij = 1, g%size
                result(:, ij) = result(:, ij)*g%rhoinv(:) - u(:, ij)*wrk3d(:, ij)
            end do

        else
!$omp parallel default( shared ) private( ij )
!$omp do
            do ij = 1, nlines*g%size
                result(ij, 1) = result(ij, 1) - u(ij, 1)*wrk3d(ij, 1)
            end do
!$omp end do
!$omp end parallel
        end if
    end if

    return
end subroutine OPR_BURGERS

!########################################################################
!# Routines for different specific directions:
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# tmp1       Out  Transpose velocity
!########################################################################
subroutine OPR_BURGERS_X(ivel, is, nx, ny, nz, bcs, g, s, u1, u2, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: Dealiasing
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_i
    use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_PROCS
#endif

    implicit none

    TINTEGER ivel, is, nx, ny, nz
    TINTEGER, dimension(2, *), intent(in) :: bcs ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*ny*nz), intent(in), target :: s, u1, u2
    TREAL, dimension(nx*ny*nz), intent(out), target :: result
    TREAL, dimension(nx*ny*nz), intent(inout), target :: tmp1, wrk3d
    TREAL, dimension(ny*nz), intent(inout) :: wrk2d

! -------------------------------------------------------------------
    TINTEGER nyz
    TREAL, dimension(:), pointer :: p_a, p_b, p_c, p_d, p_vel
#ifdef USE_MPI
    TINTEGER, parameter :: id = TLAB_MPI_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLAB_MPI_TRPF_I(s, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
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
    if (ivel == 0) then; p_vel => p_b
    else; p_vel => u2; end if ! always transposed needed

! -------------------------------------------------------------------
! Local transposition: make x-direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
    call DGETMO(p_a, g%size, g%size, nyz, p_b, nyz)
#else
    call DNS_TRANSPOSE(p_a, g%size, nyz, g%size, p_b, nyz)
#endif

! ###################################################################
    call OPR_BURGERS(is, nyz, bcs, g, Dealiasing(1), p_b, p_vel, p_d, wrk2d, p_c)

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
    call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
    call DNS_TRANSPOSE(p_d, nyz, g%size, nyz, p_c, g%size)
#endif

#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLAB_MPI_TRPB_I(p_c, result, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
    end if
#endif

    nullify (p_a, p_b, p_c, p_d, p_vel)

    return
end subroutine OPR_BURGERS_X

!########################################################################
!########################################################################
subroutine OPR_BURGERS_Y(ivel, is, nx, ny, nz, bcs, g, s, u1, u2, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: subsidence
    use TLAB_VARS, only: Dealiasing
    implicit none

    TINTEGER ivel, is, nx, ny, nz
    TINTEGER, dimension(2, *), intent(in) :: bcs ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*nz, ny), intent(in), target :: s, u1, u2
    TREAL, dimension(nx*nz, ny), intent(out), target :: result
    TREAL, dimension(nx*nz, ny), intent(inout), target :: tmp1, wrk3d
    TREAL, dimension(nx*nz), intent(inout) :: wrk2d

! -------------------------------------------------------------------
    TINTEGER nxy, nxz, j
    TREAL, dimension(:, :), pointer :: p_org, p_dst1, p_dst2, p_vel

! ###################################################################
    if (g%size == 1) then ! Set to zero in 2D case
        result = C_0_R

    else
! ###################################################################
        nxy = nx*ny
        nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y-direction the last one
! -------------------------------------------------------------------
        if (nz == 1) then
            p_org => s
            p_dst1 => tmp1
            p_dst2 => result
        else
#ifdef USE_ESSL
            call DGETMO(s, nxy, nxy, nz, tmp1, nz)
#else
            call DNS_TRANSPOSE(s, nxy, nz, nxy, tmp1, nz)
#endif
            p_org => tmp1
            p_dst1 => result
            p_dst2 => wrk3d
        end if

! pointer to velocity
        if (ivel == 0) then; p_vel => p_org
        else; 
            if (nz == 1) then; p_vel => u1         ! I do not need the transposed
            else; p_vel => u2; end if  ! I do     need the transposed
        end if

! ###################################################################
        call OPR_BURGERS(is, nxz, bcs, g, Dealiasing(2), p_org, p_vel, p_dst2, wrk2d, p_dst1)

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
            call DNS_TRANSPOSE(p_dst2, nz, nxy, nz, result, nxy)
#endif
        end if

        nullify (p_org, p_dst1, p_dst2, p_vel)

    end if

    return
end subroutine OPR_BURGERS_Y

!########################################################################
!########################################################################
subroutine OPR_BURGERS_Z(ivel, is, nx, ny, nz, bcs, g, s, u1, u2, result, tmp1, wrk2d, wrk3d)

    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: Dealiasing
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_k
    use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLAB_MPI_PROCS
#endif
    implicit none

    TINTEGER ivel, is, nx, ny, nz
    TINTEGER, dimension(2, *), intent(in) :: bcs ! BCs at xmin (1,*) and xmax (2,*)
    type(grid_dt), intent(in) :: g
    TREAL, dimension(nx*ny*nz), intent(in), target :: s, u1, u2
    TREAL, dimension(nx*ny*nz), intent(out), target :: result
    TREAL, dimension(nx*ny*nz), intent(inout), target :: tmp1, wrk3d
    TREAL, dimension(nx*ny), intent(inout) :: wrk2d

! -------------------------------------------------------------------
    TINTEGER nxy
    TREAL, dimension(:), pointer :: p_a, p_b, p_c, p_vel
#ifdef USE_MPI
    TINTEGER, parameter :: id = TLAB_MPI_K_PARTIAL
#endif

! ###################################################################
    if (g%size == 1) then ! Set to zero in 2D case
        result = C_0_R

    else
! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPF_K(s, tmp1, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            p_a => tmp1
            p_b => result
            p_c => wrk3d
            nxy = ims_size_k(id)
        else
#endif
            p_a => s
            p_b => tmp1
            p_c => result
            nxy = nx*ny
#ifdef USE_MPI
        end if
#endif

! pointer to velocity
        if (ivel == 0) then; p_vel => p_a
        else
#ifdef USE_MPI
            if (ims_npro_k > 1) then ! I do     need the transposed
                p_vel => u2
            else
#endif
                p_vel => u1                ! I do not need the transposed
#ifdef USE_MPI
            end if
#endif
        end if

! ###################################################################
        call OPR_BURGERS(is, nxy, bcs, g, Dealiasing(3), p_a, p_vel, p_c, wrk2d, p_b)

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPB_K(p_c, result, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        end if
#endif

        nullify (p_a, p_b, p_c, p_vel)

    end if

    return
end subroutine OPR_BURGERS_Z

module Averages
    use TLab_Constants, only: wp, wi
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_comm_z, ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_err
#endif
    implicit none
    private

#ifdef USE_MPI
    real(wp) sum_mpi
#endif

    public :: AVG1V1D, SUM1V1D_V
    public :: AVG1V2D, AVG1V2D_V, AVG1V2D1G
    public :: INTER1V2D
    public :: COV2V1D
    public :: COV2V2D
    public :: AVG_IK, AVG_IK_V

contains
    ! this routine should initialize area and the arrays for the weights to calculate the
    ! integrals. We assume for the time being that it is uniform grid in the directions along which
    ! the integral is calculated
    ! this routine should probably allocate as well the wrk space that is needed for the mpi case...
    ! ! subroutine Averages_Initialize()
    !     if (present(area)) then
    !         area = scalex
    !         if (kmax > 1) area = area*scalez ! 3D case
    !     end if
    ! ! end subroutine Averages_Initialize
!########################################################################
! Average along line ij
!########################################################################
    function AVG1V1D(nx, ny, nz, i, j, imom, a) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, i, j, imom   ! Moment order
        real(wp),    intent(in) :: a(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            avg = avg + a(i, j, k)**imom
        end do

        avg = avg/real(nz, wp)
#ifdef USE_MPI
        sum_mpi = avg/real(ims_npro_k, wp)
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
#endif

        return
    end function AVG1V1D

!########################################################################
! Adding in k of the matrix a for all the elements in j
!########################################################################
    subroutine SUM1V1D_V(ny, nz, a, avg, wrk)
        integer(wi), intent(in)    :: ny, nz
        real(wp),    intent(in)    :: a(ny, nz)
        real(wp),    intent(out)   :: avg(ny)
        real(wp),    intent(inout) :: wrk(ny)

        ! -------------------------------------------------------------------
        integer(wi) k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            avg(:) = avg(:) + a(:, k)
        end do

#ifdef USE_MPI
        wrk = avg
        call MPI_ALLREDUCE(wrk, avg, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

        return
    end subroutine SUM1V1D_V

!########################################################################
! Covariance along line ij
!########################################################################
    function COV2V1D(nx, ny, nz, i, j, a, b) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, i, j
        real(wp),    intent(in) :: a(nx, ny, nz), b(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            avg = avg + a(i, j, k)*b(i, j, k)
        end do

        avg = avg/real(nz, wp)
#ifdef USE_MPI
        sum_mpi = avg/real(ims_npro_k, wp)
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
#endif

        return
    end function COV2V1D

!########################################################################
! Average within the plane j
!########################################################################
    function AVG1V2D(nx, ny, nz, j, imom, a) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, j, imom      ! Moment order
        real(wp),    intent(in) :: a(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) i, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do i = 1, nx
                avg = avg + a(i, j, k)**imom
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        sum_mpi = avg/real(ims_npro_i*ims_npro_k, wp)
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

        return
    end function AVG1V2D

! #######################################################################
! #######################################################################
! Vector form
    subroutine AVG1V2D_V(nx, ny, nz, imom, a, avg, wrk)
        integer(wi), intent(in)    :: nx, ny, nz, imom      ! Moment order
        real(wp),    intent(in)    :: a(nx, ny, nz)
        real(wp),    intent(out)   :: avg(ny)
        real(wp),    intent(inout) :: wrk(ny)

        ! -------------------------------------------------------------------
        integer(wi) i, j, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    avg(j) = avg(j) + a(i, j, k)**imom
                end do
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        wrk = avg/real(ims_npro_i*ims_npro_k, wp)
        call MPI_ALLREDUCE(wrk, avg, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

    end subroutine AVG1V2D_V

!########################################################################
! Average within the plane j conditioned on the intermittency signal given by array gate
!########################################################################
    function AVG1V2D1G(nx, ny, nz, j, igate, imom, a, gate) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, j, imom      ! Moment order
        integer(1),  intent(in) :: igate                    ! Gate level to use
        real(wp),    intent(in) :: a(nx, ny, nz)
        integer(1),  intent(in) :: gate(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) i, k, nsample

#ifdef USE_MPI
        integer(wi) nsample_mpi
#endif

        ! ###################################################################
        avg = 0.0_wp
        nsample = 0
        do k = 1, nz
            do i = 1, nx
                if (gate(i, j, k) == igate) then
                    avg = avg + a(i, j, k)**imom
                    nsample = nsample + 1
                end if
            end do
        end do

#ifdef USE_MPI
        nsample_mpi = nsample
        call MPI_ALLREDUCE(nsample_mpi, nsample, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)

        sum_mpi = avg
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

        if (nsample > 0) avg = avg/real(nsample, wp)

        return
    end function AVG1V2D1G

!########################################################################
! Intermittency factor within the plane j
!########################################################################
    function INTER1V2D(nx, ny, nz, j, igate, gate) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, j
        integer(1),  intent(in) :: gate(nx, ny, nz), igate
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) i, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do i = 1, nx
                if (gate(i, j, k) == igate) then
                    avg = avg + 1.0_wp
                end if
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        sum_mpi = avg/real(ims_npro_i*ims_npro_k, wp)
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

        return
    end function INTER1V2D

! ###################################################################
! Covariance within the plane j
! ###################################################################
    function COV2V2D(nx, ny, nz, j, a, b) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, j
        real(wp),    intent(in) :: a(nx, ny, nz), b(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) i, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do i = 1, nx
                avg = avg + a(i, j, k)*b(i, j, k)
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        sum_mpi = avg/real(ims_npro_i*ims_npro_k, wp)
        call MPI_ALLREDUCE(sum_mpi, avg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

        return
    end function COV2V2D

!########################################################################
!# Calculate the average of the plane j in array a over nonuniform grids.
!# Same as avg1v2d but specific for the common case imom = 1
!########################################################################
    function AVG_IK(nx, ny, nz, j, a) result(avg)
        integer(wi), intent(in) :: nx, ny, nz, j
        real(wp),    intent(in) :: a(nx, ny, nz)
        real(wp) avg

        ! -------------------------------------------------------------------
        integer(wi) i, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do i = 1, nx
                avg = avg + a(i, j, k)
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        call MPI_ALLREDUCE(avg, sum_mpi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        avg = sum_mpi/real(ims_npro_i*ims_npro_k, wp)
#endif

        return
    end function AVG_IK

!########################################################################
!########################################################################
! Vector form
    subroutine AVG_IK_V(nx, ny, nz, jm, a, avg, wrk)
        integer(wi), intent(in)    :: nx, ny, nz, jm
        real(wp),    intent(in)    :: a(nx, ny, nz)
        real(wp),    intent(out)   :: avg(jm)
        real(wp),    intent(inout) :: wrk(jm)

        ! -------------------------------------------------------------------
        integer(wi) i, j, k

        ! ###################################################################
        avg = 0.0_wp
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    avg(j) = avg(j) + a(i, j, k)
                end do
            end do
        end do

        avg = avg/real(nx*nz, wp)
#ifdef USE_MPI
        call MPI_ALLREDUCE(avg, wrk, ny, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        avg = wrk/real(ims_npro_i*ims_npro_k, wp)
#endif

        return
    end subroutine AVG_IK_V

! To be rewritten in terms of weights precalculated in Averages_Initialize
!########################################################################
! !# Calculate the average of the plane j in array a over nonuniform grids.
! !########################################################################
!     function AVG_XZ(nx, ny, nz, j, a, dx, dz) result(avg)
!         integer(wi), intent(in) :: nx, ny, nz, j
!         real(wp),    intent(in) :: dx(*), dz(*), area
!         real(wp),    intent(in) :: a(nx, ny, nz)
!         real(wp) avg

!         ! -------------------------------------------------------------------
!         integer(wi) i, k, idsp, kdsp
!         real(wp) sum

!         ! ###################################################################
! #ifdef USE_MPI
!         idsp = ims_offset_i
!         kdsp = ims_offset_k
! #else
!         idsp = 0
!         kdsp = 0
! #endif

!         ! number of + ops: nx*nz*1 + nz*1
!         ! number of * ops: nx*nz*1 + nz*1
!         avg = 0.0_wp
!         do k = 1, nz
!             sum = 0.0_wp
!             do i = 1, nx
!                 sum = sum + a(i, j, k)*dx(idsp + i)
!             end do
!             avg = avg + sum*dz(kdsp + k)
!         end do

!         avg = avg/area
! #ifdef USE_MPI
!         sum = avg
!         call MPI_ALLREDUCE(sum, avg, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
! #else
! #endif

!         return
!     end function AVG_XZ

! !########################################################################
! !########################################################################
! ! Vector form
!     subroutine AVG_XZ_V(nx, ny, nz, jm, a, dx, dz, avg, wrk)
!         integer(wi), intent(in)    :: nx, ny, nz, jm
!         real(wp),    intent(in)    :: dx(*), dz(*), area
!         real(wp),    intent(in)    :: a(nx, ny, nz)
!         real(wp),    intent(out)   :: avg(jm)
!         real(wp),    intent(inout) :: wrk(jm)

!         ! -------------------------------------------------------------------
!         integer(wi) i, j, idsp, kdsp
!         real(wp) sum

!         ! ###################################################################
! #ifdef USE_MPI
!         idsp = ims_offset_i
!         kdsp = ims_offset_k
! #else
!         idsp = 0
!         kdsp = 0
! #endif

!         avg = 0.0_wp
!         do k = 1, nz
!             do j = 1, jm
!                 sum = 0.0_wp
!                 do i = 1, nx
!                     sum = sum + a(i, j, k)*dx(idsp + i)
!                 end do
!                 avg(j) = avg(j) + sum*dz(k + kdsp)
!             end do
!         end do

!         avg = avg/area
! #ifdef USE_MPI
!         wrk = avg
!         call MPI_ALLREDUCE(wrk, avg, jm, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
! #endif

!         return
!     end subroutine AVG_XZ_V

end module Averages

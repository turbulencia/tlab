#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!# 2023/12/07 - S. Deshpande
!#              Modified
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   generates relevant geometry fields for IBM routines
!#     nobi,    nobj,   nobk   : number of objects in i/j/k
!#     nobi_b,  nobj_b, nobk_b : beginn of objects in i/j/k
!#     nobi_e,  nobj_e, nobk_e : end    of objects in i/j/k
!#
!#
!########################################################################
!# ARGUMENTS
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!#
!########################################################################

subroutine IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)

    use IBM_VARS
    use FDM, only: g
    use TLAB_VARS, only: isize_field, imax, jmax, kmax
    use TLab_Constants, only: wi, wp
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k, ims_err
    ! use TLabMPI_Transpose, only: ims_size_i, ims_size_k
    use TLabMPI_Transpose, only: ims_trp_plan_i, ims_trp_plan_k
#ifdef IBM_DEBUG
    use TLabMPI_VARS, only: ims_pro
#endif
#endif

    implicit none

    real(wp), dimension(isize_field), intent(in) :: epsi, epsj, epsk

! #ifdef USE_MPI
!     integer(wi), parameter :: idi = TLAB_MPI_TRP_I_PARTIAL
!     integer(wi), parameter :: idk = TLAB_MPI_TRP_K_PARTIAL
! #endif
    integer(wi) :: i, j, k, ij, ik, jk, ip, inum, rse
    integer(wi) :: nyz, nxz, nxy

#ifdef USE_MPI
    integer(wi) :: dummy
#else
#ifdef IBM_DEBUG
    integer(wi), parameter :: ims_pro = 0
#endif
#endif

    ! ================================================================== !
    ! npages
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        ! nyz = ims_size_i(idi)
        nyz = ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL)%nlines
    else
#endif
        nyz = jmax*kmax
#ifdef USE_MPI
    end if
#endif

    nxz = imax*kmax

#ifdef USE_MPI
    if (ims_npro_k > 1) then
        ! nxy = ims_size_k(idk)
        nxy = ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL)%nlines
    else
#endif
        nxy = imax*jmax
#ifdef USE_MPI
    end if
#endif

    ! initialize
    nobi(:) = 0; nobj(:) = 0; nobk(:) = 0
    nobi_b(:) = 0; nobj_b(:) = 0; nobk_b(:) = 0
    nobi_e(:) = 0; nobj_e(:) = 0; nobk_e(:) = 0
    nobi_max = 0; nobj_max = 0; nobk_max = 0

    ! ================================================================== !
    ! number, beginning and end of objects in x-direction
    do jk = 1, nyz            ! pages of   i-lines
        ip = 1
        do i = 1, g(1)%size - 1     ! contiguous i-lines
            if ((ip == 1) .and. (epsi(jk) == 1.0_wp)) then ! exception: check first plane for objects
                nobi(jk) = 1
            end if
            if ((epsi(ip + jk - 1) == 0.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then ! check for interface
                if (nobi_b(jk) /= i + 1) then ! check if the object is new
                    nobi(jk) = nobi(jk) + 1
                end if
            end if
            if ((i == 1) .and. (epsi(jk) == 1.0_wp)) then ! exception: check first plane for interface
                if (g(1)%periodic .eqv. .false.) then
                    nobi_b(jk) = 1 ! nobi_b
                else if ((g(1)%periodic .eqv. .true.) .and. ((epsi(nyz*(g(1)%size - 1) + jk) == 0.0_wp))) then ! check if object is not at the end of boundary
                    nobi_b(jk) = 1 ! nobi_b
                else if ((g(1)%periodic .eqv. .true.) .and. ((epsi(nyz*(g(1)%size - 1) + jk) == 1.0_wp))) then ! check if object spans across boundary
                    inum = 0
                    rse = 0
                    ! reverse loop to check the start of object on the boundary
                    do while ((nyz*(g(1)%size - 1) + jk - inum) >= jk)
                        if (epsi(nyz*(g(1)%size - 1) + jk - inum) /= 0) then
                            inum = inum + nyz
                            rse = rse + 1
                        else
                            exit
                        end if
                    end do
                    nobi_b(jk) = g(1)%size - rse + 1
                end if
            end if
            if ((epsi(ip + jk - 1) == 0.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then     ! nobi_b check for interface
                inum = 0
                do while (nobi_b(inum + jk) /= 0)
                    inum = inum + nyz
                end do
                if (nobi_b(jk) /= (i + 1)) then
                    nobi_b(inum + jk) = i + 1
                end if
            elseif ((epsi(ip + jk - 1) == 1.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 0.0_wp)) then ! nobi_e check for interface
                inum = 0
                do while (nobi_e(inum + jk) /= 0)
                    inum = inum + nyz
                end do
                nobi_e(inum + jk) = i
            end if
            if ((i == (g(1)%size - 1)) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (nobi_e(inum + jk) /= 0)
                    inum = inum + nyz
                end do
                if (g(1)%periodic .eqv. .false.) then
                    nobi_e(inum + jk) = g(1)%size
                else if ((epsi(jk) /= 1.0_wp) .and. (g(1)%periodic .eqv. .true.)) then ! check for objects in 1st plane
                    nobi_e(inum + jk) = g(1)%size
                else if ((nobi(jk) == 1) .and. (nobi_e(jk) == 0) .and. &
                         (g(1)%periodic .eqv. .true.)) then ! Number of object found == 1 and its end is not detected yet
                    nobi_e(inum + jk) = g(1)%size
                else if ((nobi(jk) == 1) .and. (nobi_e(jk) /= 0) .and. &
                         (g(1)%periodic .eqv. .true.)) then ! Object extends across the borders of domain
                    ! Do nothing. This is not the end of the object
                end if
            end if
            ip = ip + nyz
        end do
    end do

    ! ================================================================== !
    ! number, begin and end of objects in y-direction
    do ik = 1, nxz            ! pages of   j-lines
        ip = 1
        do j = 1, g(2)%size - 1     ! contiguous j-lines
            if ((ip == 1) .and. (epsj(ik) == 1.0_wp)) then ! exception: check first plane for objects
                nobj(ik) = 1
            end if
            if ((epsj(ip + ik - 1) == 0.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then ! check for interface
                if (nobj_b(ik) /= j + 1) then  ! check if the object is new
                    nobj(ik) = nobj(ik) + 1
                end if
            end if
            if ((j == 1) .and. (epsj(ik) == 1.0_wp)) then ! exception: check first plane for interface
                if (g(2)%periodic .eqv. .false.) then
                    nobj_b(ik) = j ! nobj_b
                else if ((g(2)%periodic .eqv. .true.) .and. ((epsj(nxz*(g(2)%size - 1) + ik) == 0.0_wp))) then
                    nobj_b(ik) = j ! nobj_b
                else if ((g(2)%periodic .eqv. .true.) .and. ((epsj(nxz*(g(2)%size - 1) + ik) == 1.0_wp))) then
                    inum = 0
                    rse = 0
                    ! reverse loop to check the start of object on the boundary
                    do while ((nxz*(g(2)%size - 1) + ik - inum) >= ik)
                        if (epsj(nxz*(g(2)%size - 1) + ik - inum) /= 0.0_wp) then
                            inum = inum + nxz
                            rse = rse + 1
                        else
                            exit
                        end if
                    end do
                    nobj_b(ik) = g(1)%size - rse + 1
                end if
            end if
            if ((epsj(ip + ik - 1) == 0.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then     ! nobj_b check for interface
                inum = 0
                do while (nobj_b(inum + ik) /= 0)
                    inum = inum + nxz
                end do
                if (nobj_b(ik) /= (j + 1)) then
                    nobj_b(inum + ik) = j + 1
                end if
            elseif ((epsj(ip + ik - 1) == 1.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 0.0_wp)) then ! nobj_e check for interface
                inum = 0
                do while (nobj_e(inum + ik) /= 0)
                    inum = inum + nxz
                end do
                nobj_e(inum + ik) = j
            end if
            if ((j == (g(2)%size - 1)) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (nobj_e(inum + ik) /= 0)
                    inum = inum + nxz
                end do
                if (g(2)%periodic .eqv. .false.) then
                    nobj_e(inum + ik) = g(2)%size
                else if ((epsj(ik) /= 1.0_wp) .and. (g(2)%periodic .eqv. .true.)) then
                    nobj_e(inum + ik) = g(2)%size
                else if ((nobj(ik) == 1) .and. (nobj_e(ik) == 0) .and. &
                         (g(2)%periodic .eqv. .true.)) then
                    nobj_e(inum + ik) = g(2)%size
                else if ((nobj(ik) == 1) .and. (nobj_e(ik) /= 0) .and. &
                         (g(2)%periodic .eqv. .true.)) then ! Object extends across the borders of domain
                    ! Do nothing. This is not the end of the object
                end if
            end if
            ip = ip + nxz
        end do
    end do

    ! ================================================================== !
    ! number, begin and end of objects in z-direction
    do ij = 1, nxy            ! pages of   k-lines
        ip = 1
        do k = 1, g(3)%size - 1     ! contiguous k-lines
            if ((ip == 1) .and. (epsk(ij) == 1.0_wp)) then ! exception: check first plane for objects
                nobk(ij) = 1
            end if
            if ((epsk(ip + ij - 1) == 0.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then ! check for interface
                if (nobk_b(ij) /= k) then  ! check if the object is new
                    nobk(ij) = nobk(ij) + 1
                end if
            end if
            if ((k == 1) .and. (epsk(ij) == 1.0_wp)) then ! exception: check first plane for interface
                if (g(2)%periodic .eqv. .false.) then
                    nobk_b(ij) = k ! nobj_b
                else if ((g(2)%periodic .eqv. .true.) .and. ((epsk(nxy*(g(3)%size - 1) + ij) == 0.0_wp))) then
                    nobk_b(ij) = k
                else if ((g(3)%periodic .eqv. .true.) .and. ((epsk(nxy*(g(3)%size - 1) + ij) == 1.0_wp))) then ! check if object spans across boundary
                    inum = 0
                    rse = 0
                    ! reverse loop to check the start of object on the boundary
                    do while ((nxy*(g(3)%size - 1) + ij - inum) >= ij)
                        if (epsk(nxy*(g(3)%size - 1) + ij - inum) /= 0.0_wp) then
                            inum = inum + nxy
                            rse = rse + 1
                        else
                            exit
                        end if
                    end do
                    nobk_b(ij) = g(3)%size - rse + 1
                end if
            end if
            if ((epsk(ip + ij - 1) == 0.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then     ! nobk_b check for interface
                inum = 0
                do while (nobk_b(inum + ij) /= 0)
                    inum = inum + nxy
                end do
                if (nobk_b(ij) /= (k + 1)) then
                    nobk_b(inum + ij) = k + 1
                end if
            elseif ((epsk(ip + ij - 1) == 1.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 0.0_wp)) then ! nobk_e check for interface
                inum = 0
                do while (nobk_e(inum + ij) /= 0)
                    inum = inum + nxy
                end do
                nobk_e(inum + ij) = k
            end if
            if ((k == (g(3)%size - 1)) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (nobk_e(inum + ij) /= 0)
                    inum = inum + nxy
                end do
                if (g(1)%periodic .eqv. .false.) then
                    nobk_e(inum + ij) = g(3)%size
                else if ((epsk(ij) /= 1.0_wp) .and. (g(3)%periodic .eqv. .true.)) then ! check for objects in 1st plane
                    nobk_e(inum + ij) = g(3)%size
                else if ((nobk(ij) == 1) .and. (nobk_e(ij) == 0) .and. &
                         (g(3)%periodic .eqv. .true.)) then ! Number of object found == 1 and its end is not detected yet
                    nobk_e(inum + ij) = g(3)%size
                else if ((nobi(jk) == 1) .and. (nobi_e(jk) /= 0) .and. &
                         (g(1)%periodic .eqv. .true.)) then ! Object extends across the borders of domain
                    ! Do nothing. This is not the end of the object
                end if
            end if
            ip = ip + nxy
        end do
    end do

    ! ================================================================== !
    nobi_max = maxval(nobi)
    nobj_max = maxval(nobj)
    nobk_max = maxval(nobk)
#ifdef USE_MPI
    dummy = nobi_max
    call MPI_ALLREDUCE(dummy, nobi_max, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
    dummy = nobj_max
    call MPI_ALLREDUCE(dummy, nobj_max, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
    dummy = nobk_max
    call MPI_ALLREDUCE(dummy, nobk_max, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

    ! ================================================================== !
#ifdef IBM_DEBUG
    if (ims_pro == 0) then
        write (*, *) '======== Max number of objects in each direction ========'
        write (*, *) 'max number of objects in x = ', nobi_max
        write (*, *) 'max number of objects in y = ', nobj_max
        write (*, *) 'max number of objects in z = ', nobk_max
    end if
#endif

    return
end subroutine IBM_GENERATE_GEOMETRY

!########################################################################

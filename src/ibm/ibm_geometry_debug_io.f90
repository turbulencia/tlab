#ifdef IBM_DEBUG
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   debug IO for all relevant geometry fields for IBM routines
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

subroutine IBM_GEOMETRY_DEBUG_IO(epsi, epsj, epsk, tmp1, tmp2, tmp3)

    use IBM_VARS
    use IO_FIELDS
    use FDM, only: g
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLab_Constants, only: wi, wp
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
    use TLabMPI_VARS, only: ims_pro
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k, ims_err
#endif

    implicit none

    real(wp), dimension(isize_field), intent(in) :: epsi, epsj, epsk
    real(wp), dimension(isize_field), intent(inout) :: tmp1, tmp2, tmp3

#ifdef USE_MPI
    integer(wi), parameter :: idi = TLAB_MPI_TRP_I_PARTIAL
    integer(wi), parameter :: idk = TLAB_MPI_TRP_K_PARTIAL
#endif
    integer(wi) :: i, j, k, ij, ik, jk, ip, inum
    integer(wi) :: nyz, nxz, nxy
#ifdef USE_MPI
#else
    integer(wi), parameter :: ims_pro = 0
#endif

    ! ================================================================== !
    ! npages
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        nyz = ims_size_i(idi)
    else
#endif
        nyz = jmax*kmax
#ifdef USE_MPI
    end if
#endif

    nxz = imax*kmax

#ifdef USE_MPI
    if (ims_npro_k > 1) then
        nxy = ims_size_k(idk)
    else
#endif
        nxy = imax*jmax
#ifdef USE_MPI
    end if
#endif

    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp; tmp3(:) = 0.0_wp

    ! ================================================================== !
    ! number of objects in x-direction
    ip = 1
    do i = 1, g(1)%size - 1     ! contiguous i-lines
        do jk = 1, nyz            ! pages of   i-lines
            if ((ip == 1) .and. (epsi(jk) == 1.0_wp)) then ! exception: check first plane for objects
                tmp1(jk) = 1.0_wp
            end if
            if ((epsi(ip + jk - 1) == 0.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then ! check for interface
                tmp1(jk) = tmp1(jk) + 1.0_wp
            end if
        end do
        ip = ip + nyz
    end do

    call TLab_Transpose(tmp1, nyz, g(1)%size, nyz, tmp2, g(1)%size)
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLabMPI_TransposeI_Backward(tmp2, tmp1, ims_ds_i(1, idi), ims_dr_i(1, idi), ims_ts_i(1, idi), ims_tr_i(1, idi))
    end if
    call IO_WRITE_FIELDS('nobi3d', IO_FLOW, imax, jmax, kmax, 1, tmp1)
#else
    call IO_WRITE_FIELDS('nobi3d', IO_FLOW, imax, jmax, kmax, 1, tmp2)
#endif
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp

    ! ================================================================== !
    ! number of objects in y-direction
    ip = 1
    do j = 1, g(2)%size - 1     ! contiguous j-lines
        do ik = 1, nxz            ! pages of   j-lines
            if ((ip == 1) .and. (epsj(ik) == 1.0_wp)) then ! exception: check first plane for objects
                tmp1(ik) = 1.0_wp
            end if
            if ((epsj(ip + ik - 1) == 0.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then ! check for interface
                tmp1(ik) = tmp1(ik) + 1.0_wp
            end if
        end do
        ip = ip + nxz
    end do

    call TLab_Transpose(tmp1, kmax, imax*jmax, kmax, tmp2, imax*jmax)
    call IO_WRITE_FIELDS('nobj3d', IO_FLOW, imax, jmax, kmax, 1, tmp2)
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp

    ! ================================================================== !
    ! number of objects in z-direction
    ip = 1
    do k = 1, g(3)%size - 1     ! contiguous k-lines
        do ij = 1, nxy            ! pages of   k-lines
            if ((ip == 1) .and. (epsk(ij) == 1.0_wp)) then ! exception: check first plane for objects
                tmp1(ij) = 1.0_wp
            end if
            if ((epsk(ip + ij - 1) == 0.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then ! check for interface
                tmp1(ij) = tmp1(ij) + 1.0_wp
            end if
        end do
        ip = ip + nxy
    end do

#ifdef USE_MPI
    if (ims_npro_k > 1) then
        call TLabMPI_TransposeK_Backward(tmp1, tmp2, ims_ds_k(1, idk), ims_dr_k(1, idk), ims_ts_k(1, idk), ims_tr_k(1, idk))
    end if
    call IO_WRITE_FIELDS('nobk3d', IO_FLOW, imax, jmax, kmax, 1, tmp2)
#else
    call IO_WRITE_FIELDS('nobk3d', IO_FLOW, imax, jmax, kmax, 1, tmp1)
#endif
    if (ims_pro == 0) write (*, *) '============= Writing all geometry fields ==============='
    if (ims_pro == 0) write (*, *) 'done writing files: nobi3d,   nobj3d,   nobk3d'
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp

    ! ================================================================== !
    ! begin and end of objects in x-direction
    ip = 1
    do i = 1, g(1)%size - 1     ! contiguous i-lines
        do jk = 1, nyz            ! pages of   i-lines
            if ((i == 1) .and. (epsi(jk) == 1.0_wp)) then ! exception: check first plane for interface
                tmp1(jk) = real(i, wp) ! nobi_b
            end if
            if ((epsi(ip + jk - 1) == 0.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then     ! nobi_b check for interface
                inum = 0
                do while (tmp1(inum + jk) /= 0.0_wp)
                    inum = inum + nyz
                end do
                tmp1(inum + jk) = real(i + 1, wp)
            elseif ((epsi(ip + jk - 1) == 1.0_wp) .and. (epsi(ip + jk - 1 + nyz) == 0.0_wp)) then ! nobi_e check for interface
                inum = 0
                do while (tmp2(inum + jk) /= 0.0_wp)
                    inum = inum + nyz
                end do
                tmp2(inum + jk) = real(i, wp)
            end if
            if ((i == (g(1)%size - 1)) .and. (epsi(ip + jk - 1 + nyz) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (tmp2(inum + jk) /= 0.0_wp)
                    inum = inum + nyz
                end do
                tmp2(inum + jk) = real(g(1)%size, wp)
            end if
        end do
        ip = ip + nyz
    end do

    call TLab_Transpose(tmp1, nyz, g(1)%size, nyz, tmp3, g(1)%size)
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLabMPI_TransposeI_Backward(tmp3, tmp1, ims_ds_i(1, idi), ims_dr_i(1, idi), ims_ts_i(1, idi), ims_tr_i(1, idi))
    end if
    call IO_WRITE_FIELDS('nobi3d_b', IO_FLOW, imax, jmax, kmax, 1, tmp1)
#else
    call IO_WRITE_FIELDS('nobi3d_b', IO_FLOW, imax, jmax, kmax, 1, tmp3)
#endif

    call TLab_Transpose(tmp2, nyz, g(1)%size, nyz, tmp3, g(1)%size)
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLabMPI_TransposeI_Backward(tmp3, tmp2, ims_ds_i(1, idi), ims_dr_i(1, idi), ims_ts_i(1, idi), ims_tr_i(1, idi))
    end if
    call IO_WRITE_FIELDS('nobi3d_e', IO_FLOW, imax, jmax, kmax, 1, tmp2)
#else
    call IO_WRITE_FIELDS('nobi3d_e', IO_FLOW, imax, jmax, kmax, 1, tmp3)
#endif
    if (ims_pro == 0) write (*, *) 'done writing files: nobi3d_b, nobi3d_e'
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp; tmp3(:) = 0.0_wp

    ! ================================================================== !
    ! begin and end of objects in y-direction
    ip = 1
    do j = 1, g(2)%size - 1     ! contiguous j-lines
        do ik = 1, nxz            ! pages of   j-lines
            if ((j == 1) .and. (epsj(ik) == 1.0_wp)) then ! exception: check first plane for interface
                tmp1(ik) = real(j, wp) ! nobj_b
            end if
            if ((epsj(ip + ik - 1) == 0.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then     ! nobj_b check for interface
                inum = 0
                do while (tmp1(inum + ik) /= 0.0_wp)
                    inum = inum + nxz
                end do
                tmp1(inum + ik) = real(j + 1, wp)
            elseif ((epsj(ip + ik - 1) == 1.0_wp) .and. (epsj(ip + ik - 1 + nxz) == 0.0_wp)) then ! nobj_e check for interface
                inum = 0
                do while (tmp2(inum + ik) /= 0.0_wp)
                    inum = inum + nxz
                end do
                tmp2(inum + ik) = real(j, wp)
            end if
            if ((j == (g(2)%size - 1)) .and. (epsj(ip + ik - 1 + nxz) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (tmp2(inum + ik) /= 0.0_wp)
                    inum = inum + nxz
                end do
                tmp2(inum + ik) = real(g(2)%size, wp)
            end if
        end do
        ip = ip + nxz
    end do

    call TLab_Transpose(tmp1, kmax, imax*jmax, kmax, tmp3, imax*jmax)
    call IO_WRITE_FIELDS('nobj3d_b', IO_FLOW, imax, jmax, kmax, 1, tmp3)
    call TLab_Transpose(tmp2, kmax, imax*jmax, kmax, tmp3, imax*jmax)
    call IO_WRITE_FIELDS('nobj3d_e', IO_FLOW, imax, jmax, kmax, 1, tmp3)
    if (ims_pro == 0) write (*, *) 'done writing files: nobj3d_b, nobj3d_e'
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp; tmp3(:) = 0.0_wp

    ! ================================================================== !
    ! begin and end of objects in z-direction
    ip = 1
    do k = 1, g(3)%size - 1     ! contiguous k-lines
        do ij = 1, nxy            ! pages of   k-lines
            if ((k == 1) .and. (epsk(ij) == 1.0_wp)) then ! exception: check first plane for interface
                tmp1(ij) = real(k, wp)    ! nobk_b
            end if
            if ((epsk(ip + ij - 1) == 0.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then     ! nobk_b check for interface
                inum = 0
                do while (tmp1(ij + inum) /= 0.0_wp)
                    inum = inum + nxy
                end do
                tmp1(ij + inum) = real(k + 1, wp)
            elseif ((epsk(ip + ij - 1) == 1.0_wp) .and. (epsk(ip + ij - 1 + nxy) == 0.0_wp)) then ! nobk_e check for interface
                inum = 0
                do while (tmp2(ij + inum) /= 0.0_wp)
                    inum = inum + nxy
                end do
                tmp2(ij + inum) = real(k, wp)
            end if
            if ((k == (g(3)%size - 1)) .and. (epsk(ip + ij - 1 + nxy) == 1.0_wp)) then ! exception: check last plane for interface
                inum = 0
                do while (tmp2(inum + ij) /= 0.0_wp)
                    inum = inum + nxy
                end do
                tmp2(inum + ij) = real(g(3)%size, wp)
            end if
        end do
        ip = ip + nxy
    end do

#ifdef USE_MPI
    if (ims_npro_k > 1) then
        call TLabMPI_TransposeK_Backward(tmp1, tmp3, ims_ds_k(1, idk), ims_dr_k(1, idk), ims_ts_k(1, idk), ims_tr_k(1, idk))
    end if
    call IO_WRITE_FIELDS('nobk3d_b', IO_FLOW, imax, jmax, kmax, 1, tmp3)
    if (ims_npro_k > 1) then
        call TLabMPI_TransposeK_Backward(tmp2, tmp3, ims_ds_k(1, idk), ims_dr_k(1, idk), ims_ts_k(1, idk), ims_tr_k(1, idk))
    end if
    call IO_WRITE_FIELDS('nobk3d_e', IO_FLOW, imax, jmax, kmax, 1, tmp3)
#else
    call IO_WRITE_FIELDS('nobk3d_b', IO_FLOW, imax, jmax, kmax, 1, tmp1)
    call IO_WRITE_FIELDS('nobk3d_e', IO_FLOW, imax, jmax, kmax, 1, tmp2)
#endif
    if (ims_pro == 0) write (*, *) 'done writing files: nobk3d_b, nobk3d_e'
    if (ims_pro == 0) write (*, *) '========================================================='
    tmp1(:) = 0.0_wp; tmp2(:) = 0.0_wp; tmp3(:) = 0.0_wp

    return
end subroutine IBM_GEOMETRY_DEBUG_IO

!########################################################################

#endif

#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module PARTICLE_INTERPOLATE
    use TLAB_CONSTANTS, only: wp, wi, efile, lfile
    use TLAB_TYPES, only: pointers_dt, pointers3d_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk3d
    use TLAB_VARS, only: g
    use TLAB_PROCS
    use PARTICLE_VARS
    use PARTICLE_ARRAYS, only: l_comm
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
#endif
    implicit none
    private

    public :: FIELD_TO_PARTICLE

contains
!#######################################################################
!#######################################################################
    subroutine FIELD_TO_PARTICLE(nvar, data_in, data_out, l_g, l_q, wrk3d)
        integer(wi), intent(in) :: nvar
        type(pointers3d_dt), intent(in) :: data_in(nvar)
        type(pointers_dt), intent(out) :: data_out(nvar)
        type(particle_dt), intent(inout) :: l_g
        real(wp), intent(inout) :: l_q(isize_part, inb_part_array)
        real(wp), intent(inout) :: wrk3d(isize_wrk3d)

! -------------------------------------------------------------------
        integer(wi) grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal
        integer(wi) npar_start, npar_end
        integer(wi) ip_i, ip_k, ip_ik, np_i, np_k, np_ik, iv

        type(pointers3d_dt), dimension(nvar) :: data_halo_i, data_halo_k, data_halo_ik

!#######################################################################
        if (nvar > inb_part_interp) then
            call TLAB_WRITE_ASCII(efile, 'FIELD_TO_PARTICLE. Not enough memory.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        np_i = 2*jmax*kmax; ip_i = 1
        np_k = imax*jmax*2; ip_k = ip_i + np_i*nvar
        np_ik = 2*jmax*2; ip_ik = ip_k + np_k*nvar
        do iv = 1, nvar
            data_halo_i(iv)%field(1:2, 1:jmax, 1:kmax) => l_comm(ip_i:ip_i + np_i - 1); ip_i = ip_i + np_i
            data_halo_k(iv)%field(1:imax, 1:jmax, 1:2) => l_comm(ip_k:ip_k + np_k - 1); ip_k = ip_k + np_k
            data_halo_ik(iv)%field(1:2, 1:jmax, 1:2) => l_comm(ip_ik:ip_ik + np_ik - 1); ip_ik = ip_ik + np_ik
        end do

! -------------------------------------------------------------------
! Setting fields in halo regions
        call Create_Halo_K(nvar, data_in, data_halo_k(1)%field, wrk3d(1), wrk3d(imax*jmax*nvar + 1))
        call Create_Halo_I_IK(nvar, data_in, data_halo_i(1)%field, data_halo_k(1)%field, data_halo_ik(1)%field, wrk3d(1), &
                              wrk3d(jmax*(kmax + 1)*nvar + 1))

! -------------------------------------------------------------------
! Sorting and counting particles for each zone
        call Sort_Into_Zones(l_g, l_q, data_out, grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal)

#ifdef USE_MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Interpolating
        npar_start = 1
        npar_end = grid_zone
        call Interpolate_Inside_Zones('grid_zone', nvar, data_in, data_out, l_g, l_q, npar_start, npar_end)

        if (halo_zone_x /= 0) then
            npar_start = npar_end + 1
            npar_end = npar_end + halo_zone_x
            call Interpolate_Inside_Zones('halo_zone_x', nvar, data_halo_i, data_out, l_g, l_q, npar_start, npar_end)
        end if

        if (halo_zone_z /= 0) then
            npar_start = npar_end + 1
            npar_end = npar_end + halo_zone_z
            call Interpolate_Inside_Zones('halo_zone_z', nvar, data_halo_k, data_out, l_g, l_q, npar_start, npar_end)
        end if

        if (halo_zone_diagonal /= 0) then
            npar_start = npar_end + 1
            npar_end = npar_end + halo_zone_diagonal
            call Interpolate_Inside_Zones('halo_zone_diagonal', nvar, data_halo_ik, data_out, l_g, l_q, npar_start, npar_end)
        end if

! -------------------------------------------------------------------
        do iv = 1, nvar
            nullify (data_halo_i(iv)%field, data_halo_k(iv)%field, data_halo_ik(iv)%field)
        end do

        return
    end subroutine FIELD_TO_PARTICLE

!#######################################################################
!#######################################################################
    subroutine Create_Halo_K(nvar, data, halo_field_k, buffer_send, buffer_recv)
        integer(wi), intent(in) :: nvar
        type(pointers3d_dt), intent(in) :: data(nvar)
        real(wp), intent(out) :: halo_field_k(imax, jmax, 2, nvar)
        real(wp), intent(inout) :: buffer_send(imax, jmax, 1, nvar), buffer_recv(imax, jmax, 1, nvar)

! -------------------------------------------------------------------
        integer(wi) i

#ifdef USE_MPI
        integer source, dest, l, size
        integer mpireq(ims_npro*2 + 2)
        integer status(MPI_STATUS_SIZE, ims_npro*2)
#endif

! ######################################################################
#ifdef USE_MPI
        if (ims_npro_k == 1) then
#endif
            do i = 1, nvar
                halo_field_k(1:imax, 1:jmax, 1, i) = data(i)%field(1:imax, 1:jmax, kmax)
                halo_field_k(1:imax, 1:jmax, 2, i) = data(i)%field(1:imax, 1:jmax, 1)
            end do

#ifdef USE_MPI
        else
            do i = 1, nvar
                halo_field_k(1:imax, 1:jmax, 1, i) = data(i)%field(1:imax, 1:jmax, kmax)
                buffer_send(1:imax, 1:jmax, 1, i) = data(i)%field(1:imax, 1:jmax, 1) ! data to be transfered
            end do

            mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
            l = 2*ims_pro + 1
            dest = ims_map_k(mod(ims_pro_k - 1 + ims_npro_k, ims_npro_k) + 1)
            source = ims_map_k(mod(ims_pro_k + 1 + ims_npro_k, ims_npro_k) + 1)
            size = imax*jmax*nvar
            call MPI_ISEND(buffer_send, size, MPI_REAL8, dest, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
            call MPI_IRECV(buffer_recv, size, MPI_REAL8, source, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)
            call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

            halo_field_k(1:imax, 1:jmax, 2, 1:nvar) = buffer_recv(1:imax, 1:jmax, 1, 1:nvar)

        end if
#endif

        return
    end subroutine Create_Halo_K

!#######################################################################
!#######################################################################
    subroutine Create_Halo_I_IK(nvar, data, halo_field_i, halo_field_k, halo_field_ik, buffer_send, buffer_recv)
        integer(wi), intent(in) :: nvar
        type(pointers3d_dt), intent(in) :: data(nvar)
        real(wp), intent(out) :: halo_field_i(2, jmax, kmax, nvar)
        real(wp), intent(in) :: halo_field_k(imax, jmax, 2, nvar)
        real(wp), intent(out) :: halo_field_ik(2, jmax, 2, nvar)
        real(wp), intent(inout) :: buffer_send(1, jmax, kmax + 1, nvar), buffer_recv(1, jmax, kmax + 1, nvar)

! -------------------------------------------------------------------
        integer(wi) i

#ifdef USE_MPI
        integer source, dest, l, size
        integer mpireq(ims_npro*2 + 2)
        integer status(MPI_STATUS_SIZE, ims_npro*2)
#endif

! ######################################################################
#ifdef USE_MPI
        if (ims_npro_i == 1) then
#endif
            do i = 1, nvar
                halo_field_i(1, 1:jmax, 1:kmax, i) = data(i)%field(imax, 1:jmax, 1:kmax)
                halo_field_i(2, 1:jmax, 1:kmax, i) = data(i)%field(1, 1:jmax, 1:kmax)
                halo_field_ik(2, 1:jmax, 2, i) = halo_field_k(1, 1:jmax, 2, i) ! top-right corner
            end do

#ifdef USE_MPI
        else
            do i = 1, nvar
                halo_field_i(1, 1:jmax, 1:kmax, i) = data(i)%field(imax, 1:jmax, 1:kmax)
                buffer_send(1, 1:jmax, 1:kmax, i) = data(i)%field(1, 1:jmax, 1:kmax)   ! data to be transfered
                buffer_send(1, 1:jmax, kmax + 1, i) = halo_field_k(1, 1:jmax, 2, i)
            end do

            mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
            l = 2*ims_pro + 1
            dest = ims_map_i(mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i) + 1)
            source = ims_map_i(mod(ims_pro_i + 1 + ims_npro_i, ims_npro_i) + 1)
            size = jmax*(kmax + 1)*nvar
            call MPI_ISEND(buffer_send, size, MPI_REAL8, dest, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
            call MPI_IRECV(buffer_recv, size, MPI_REAL8, source, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)
            call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

            halo_field_i(2, 1:jmax, 1:kmax, 1:nvar) = buffer_recv(1, 1:jmax, 1:kmax, 1:nvar)
            halo_field_ik(2, 1:jmax, 2, 1:nvar) = buffer_recv(1, 1:jmax, kmax + 1, 1:nvar) ! top-right corner

        end if
#endif

        halo_field_ik(1, 1:jmax, 1, 1:nvar) = halo_field_i(1, 1:jmax, kmax, 1:nvar)
        halo_field_ik(2, 1:jmax, 1, 1:nvar) = halo_field_i(2, 1:jmax, kmax, 1:nvar)
        halo_field_ik(1, 1:jmax, 2, 1:nvar) = halo_field_k(imax, 1:jmax, 2, 1:nvar)

        return
    end subroutine Create_Halo_I_IK

!########################################################################
!########################################################################
    subroutine Interpolate_Inside_Zones(zone, nvar, data_in, data_out, l_g, l_q, grid_start, grid_end)
        character(len=*), intent(in) :: zone
        integer(wi), intent(in) :: nvar, grid_start, grid_end
        type(pointers3d_dt), intent(in) :: data_in(nvar)
        type(pointers_dt), intent(out) :: data_out(nvar)
        type(particle_dt), intent(in) :: l_g
        real(wp), intent(in) :: l_q(isize_part, 3)

! -------------------------------------------------------------------
        real(wp) length_g_p(6), cube_g_p(4)
        integer(wi) g_p(10), g1loc, g2loc, g5loc, g6loc
        integer(wi) i, j
        real(wp) dx_loc_inv, dz_loc_inv

! ######################################################################
        dx_loc_inv = real(g(1)%size, wp)/g(1)%scale
        dz_loc_inv = real(g(3)%size, wp)/g(3)%scale

! Managing the zone option outside the loop
        g_p(7) = 1
        g_p(8) = 2
        g_p(9) = 1
        g_p(10) = 2

        g1loc = 1
        g2loc = 2
        g5loc = 5
        g6loc = 6

        select case (trim(adjustl(zone)))
        case ('halo_zone_x')
            g1loc = 7
            g2loc = 8

        case ('halo_zone_z')
            g5loc = 9
            g6loc = 10

        case ('halo_zone_diagonal')
            g1loc = 7
            g2loc = 8
            g5loc = 9
            g6loc = 10

        end select

! ######################################################################
        if (g(3)%size /= 1) then

            do i = grid_start, grid_end ! loop over all particles

                length_g_p(1) = l_q(i, 1)*dx_loc_inv            ! Local X position
                g_p(1) = floor(length_g_p(1))
                length_g_p(1) = length_g_p(1) - real(g_p(1), wp)
#ifdef USE_MPI
                g_p(1) = g_p(1) + 1 - ims_offset_i
#else
                g_p(1) = g_p(1) + 1
#endif
                g_p(2) = g_p(1) + 1
                length_g_p(2) = 1.0_wp - length_g_p(1)

                length_g_p(5) = l_q(i, 3)*dz_loc_inv            ! Local Z position
                g_p(5) = floor(length_g_p(5))
                length_g_p(5) = length_g_p(5) - real(g_p(5), wp)
#ifdef USE_MPI
                g_p(5) = g_p(5) + 1 - ims_offset_k
#else
                g_p(5) = g_p(5) + 1
#endif
                g_p(6) = g_p(5) + 1
                length_g_p(6) = 1.0_wp - length_g_p(5)

                g_p(3) = l_g%nodes(i)                    ! Local Y position
                g_p(4) = g_p(3) + 1
                length_g_p(3) = (l_q(i, 2) - g(2)%nodes(g_p(3)))/(g(2)%nodes(g_p(4)) - g(2)%nodes(g_p(3)))
                length_g_p(4) = 1.0_wp - length_g_p(3)

                cube_g_p(1) = length_g_p(1)*length_g_p(3) ! bilear cubes for X and Y
                cube_g_p(2) = length_g_p(1)*length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
                cube_g_p(3) = length_g_p(4)*length_g_p(2)
                cube_g_p(4) = length_g_p(2)*length_g_p(3)

! -------------------------------------------------------------------
! Trilinear interpolation
! Two bilinear interpolations for each k plane (g_p(5) and g_p(6)
! Then multipled by (1-length) for trilinear aspect
! -------------------------------------------------------------------
                do j = 1, nvar
                    data_out(j)%field(i) = data_out(j)%field(i) + &
                                           ((cube_g_p(3)*data_in(j)%field(g_p(g1loc), g_p(3), g_p(g5loc)) &
                                             + cube_g_p(4)*data_in(j)%field(g_p(g1loc), g_p(4), g_p(g5loc)) &
                                             + cube_g_p(1)*data_in(j)%field(g_p(g2loc), g_p(4), g_p(g5loc)) &
                                             + cube_g_p(2)*data_in(j)%field(g_p(g2loc), g_p(3), g_p(g5loc)))*length_g_p(6)) &
                                           + ((cube_g_p(3)*data_in(j)%field(g_p(g1loc), g_p(3), g_p(g6loc)) &
                                               + cube_g_p(4)*data_in(j)%field(g_p(g1loc), g_p(4), g_p(g6loc)) &
                                               + cube_g_p(1)*data_in(j)%field(g_p(g2loc), g_p(4), g_p(g6loc)) &
                                               + cube_g_p(2)*data_in(j)%field(g_p(g2loc), g_p(3), g_p(g6loc)))*length_g_p(5))
                end do

            end do

! ######################################################################
        else !2D case

            do i = grid_start, grid_end

                length_g_p(1) = l_q(i, 1)*dx_loc_inv
                g_p(1) = floor(length_g_p(1))
                length_g_p(1) = length_g_p(1) - real(g_p(1), wp)
#ifdef USE_MPI
                g_p(1) = g_p(1) + 1 - ims_offset_i
#else
                g_p(1) = g_p(1) + 1
#endif
                g_p(2) = g_p(1) + 1
                length_g_p(2) = 1.0_wp - length_g_p(1)

                g_p(3) = l_g%nodes(i)
                g_p(4) = g_p(3) + 1
                length_g_p(3) = (l_q(i, 2) - g(2)%nodes(g_p(3)))/(g(2)%nodes(g_p(4)) - g(2)%nodes(g_p(3)))
                length_g_p(4) = 1.0_wp - length_g_p(3)

                cube_g_p(1) = length_g_p(1)*length_g_p(3)
                cube_g_p(2) = length_g_p(1)*length_g_p(4)
                cube_g_p(3) = length_g_p(4)*length_g_p(2)
                cube_g_p(4) = length_g_p(2)*length_g_p(3)

! -------------------------------------------------------------------
! Bilinear interpolation
! -------------------------------------------------------------------
                do j = 1, nvar
                    data_out(j)%field(i) = data_out(j)%field(i) + &
                                           (cube_g_p(3)*data_in(j)%field(g_p(g1loc), g_p(3), 1) &
                                            + cube_g_p(4)*data_in(j)%field(g_p(g1loc), g_p(4), 1) &
                                            + cube_g_p(1)*data_in(j)%field(g_p(g2loc), g_p(4), 1) &
                                            + cube_g_p(2)*data_in(j)%field(g_p(g2loc), g_p(3), 1))
                end do

            end do

        end if

        return
    end subroutine Interpolate_Inside_Zones

!########################################################################
! Sorting structure grid-halo_x-halo_z-halo_diagonal
!########################################################################
    subroutine Sort_Into_Zones(l_g, l_q, data, grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal)
        type(pointers_dt), intent(inout) :: data(:)
        type(particle_dt), intent(inout) :: l_g
        real(wp), intent(inout) :: l_q(:,:)
        integer(wi), intent(out) :: grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal

        ! -------------------------------------------------------------------
        real(wp) dummy, right_limit, upper_limit
        integer(wi) idummy, nvar
        integer(longi) idummy8
        integer(wi) i, j, k

        !#######################################################################
#ifdef USE_MPI
        right_limit = g(1)%nodes(ims_offset_i + imax)  ! right_limit is east
        upper_limit = g(3)%nodes(ims_offset_k + kmax)  ! upper_limit is north
#else
        right_limit = g(1)%nodes(g(1)%size)
        upper_limit = g(3)%nodes(g(3)%size)
#endif

        nvar = size(data)
        
        !#######################################################################
        ! Group together particles outside of the grid zone at the end of the arrays
        grid_zone = 0
        i = 1       ! starting point of sorting algorithm
        j = l_g%np  ! end point of sorting algorithm

        do while (i <= j)

            if (l_q(i, 1) > right_limit .or. l_q(i, 3) > upper_limit) then          ! particle i is in halo
                do
                    if (l_q(j, 1) > right_limit .or. l_q(j, 3) > upper_limit) then  ! partcile j is in halo, leave it there
                        j = j - 1                                                   ! go to next particle
                        if (i >= j) exit                                            ! finished your upwards loop

                    else                                                            ! found a particle in grid, so swap
                        idummy = l_g%nodes(i)
                        l_g%nodes(i) = l_g%nodes(j)
                        l_g%nodes(j) = idummy

                        idummy8 = l_g%tags(i)
                        l_g%tags(i) = l_g%tags(j)
                        l_g%tags(j) = idummy8

                        do k = 1, inb_part_array
                            dummy = l_q(i, k)
                            l_q(i, k) = l_q(j, k)
                            l_q(j, k) = dummy
                        end do

                        do k = 1, nvar      ! swap also this data because of the substages in the Runge-Kutta
                            dummy = data(k)%field(i)
                            data(k)%field(i) = data(k)%field(j)
                            data(k)%field(j) = dummy
                        end do

                        j = j - 1
                        grid_zone = grid_zone + 1
                        exit

                    end if
                end do

            else            ! particle i is in the grid, the right place
                grid_zone = grid_zone + 1

            end if
            i = i + 1       ! go to next particle

        end do

        ! !Last particle might not be properly checked. We check it here.
        ! if (i == j) then !Probably not needed
        !     if ((l_q(i, 1) > right_limit) .or. (l_q(i, 3) > upper_limit)) then !Particle out the grid
        !         !Do nothing
        !     else
        !         grid_zone = grid_zone + 1 !The last particle was not checked but it is in the grid
        !     end if
        ! end if
        !Possible optimization to avoid if i. EQ. j. Not completed.
        ! ELSE IF (i-1 .GT. j+1) THEN
        !   j=j+1
        !   i=i-1
        !   grid_zone=grid_zone-1
        !   IF ( (l_q(j,1) .GT. right_limit) .OR. l_q(j,3) .GT. upper_limit) THEN
        !          DO k=1,3
        !           dummy=l_q(i,k)
        !           l_q(i,k)=l_q(j,k)
        !           l_q(j,k)=dummy

        !           idummy8=l_g%tags(i)
        !           l_g%tags(i)=l_g%tags(j)
        !           l_g%tags(j)=idummy8
        !         END DO
        !   END IF

        ! -------------------------------------------------------------------
        ! From the remaning particles, group together particles outside the East zone at the end of the array
        halo_zone_x = 0
        i = grid_zone + 1
        j = l_g%np

        do while (i <= j)
            if (l_q(i, 3) > upper_limit) then           ! particle i is in North
                do
                    if (l_q(j, 3) > upper_limit) then   ! particle j is in North, leave it here
                        j = j - 1
                        if (i >= j) exit

                    else                                ! found a particle in East, so swap
                        idummy = l_g%nodes(i)
                        l_g%nodes(i) = l_g%nodes(j)
                        l_g%nodes(j) = idummy

                        idummy8 = l_g%tags(i)
                        l_g%tags(i) = l_g%tags(j)
                        l_g%tags(j) = idummy8

                        do k = 1, inb_part_array
                            dummy = l_q(i, k)
                            l_q(i, k) = l_q(j, k)
                            l_q(j, k) = dummy
                        end do

                        do k = 1, nvar
                            dummy = data(k)%field(i)
                            data(k)%field(i) = data(k)%field(j)
                            data(k)%field(j) = dummy
                        end do

                        j = j - 1
                        halo_zone_x = halo_zone_x + 1
                        exit

                    end if
                end do

            else
                halo_zone_x = halo_zone_x + 1

            end if
            i = i + 1

        end do

        ! !Last particle might not be properly checked. We check it here.
        ! if (i == j) then !Probably not needed
        !     if (l_q(i, 3) > upper_limit) then !Particle is out North
        !         !Do nothing
        !     else
        !         halo_zone_x = halo_zone_x + 1
        !     end if
        ! end if

        ! -------------------------------------------------------------------
        ! From the remaning particles, group together particles outside the North zone at the end of the array
        halo_zone_z = 0
        i = grid_zone + halo_zone_x + 1
        j = l_g%np

        do while (i <= j)
            if (l_q(i, 1) > right_limit) then           ! particle i is in North-east
                do 
                    if (l_q(j, 1) > right_limit) then   ! particle j is in North-east, leave it here
                        j = j - 1
                        if (i >= j) exit

                    else                                ! found a particle in North, so swap
                        idummy = l_g%nodes(i)
                        l_g%nodes(i) = l_g%nodes(j)
                        l_g%nodes(j) = idummy

                        idummy8 = l_g%tags(i)
                        l_g%tags(i) = l_g%tags(j)
                        l_g%tags(j) = idummy8

                        do k = 1, inb_part_array
                            dummy = l_q(i, k)
                            l_q(i, k) = l_q(j, k)
                            l_q(j, k) = dummy
                        end do

                        do k = 1, nvar
                            dummy = data(k)%field(i)
                            data(k)%field(i) = data(k)%field(j)
                            data(k)%field(j) = dummy
                        end do

                        j = j - 1
                        halo_zone_z = halo_zone_z + 1
                        exit

                    end if
                end do

            else
                halo_zone_z = halo_zone_z + 1

            end if
            i = i + 1

        end do

        ! !Last particle might not be properly checked. We check it here.
        ! if (i == j) then !Probably not needed
        !     if (l_q(i, 1) > right_limit) then !Particle is out East
        !         !Do nothing
        !     else
        !         halo_zone_z = halo_zone_z + 1
        !     end if
        ! end if

        ! -------------------------------------------------------------------
        ! What remains at the end of the array is in the North-East zone
        halo_zone_diagonal = l_g%np - grid_zone - halo_zone_x - halo_zone_z

        return
    end subroutine Sort_Into_Zones

end module PARTICLE_INTERPOLATE

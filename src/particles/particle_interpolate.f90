#include "dns_error.h"

module PARTICLE_INTERPOLATE
    use TLAB_CONSTANTS, only: wp, wi, efile
    use TLAB_VARS, only: lfile
    use TLAB_TYPES, only: pointers_dt, pointers3d_dt
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use TLAB_PROCS
    use PARTICLE_VARS
    use PARTICLE_ARRAYS, only: halo_i, halo_k, halo_ik
    use PARTICLE_ARRAYS, only: p_halo_i, p_halo_k, p_halo_ik
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
    use PARTICLE_ARRAYS, only: halo_mpi_recv_i, halo_mpi_recv_k, halo_mpi_send_i, halo_mpi_send_k
#endif
    implicit none
    private

    integer iv, nvar
#ifdef USE_MPI
    integer source, dest, l, buff_size
#endif

    public :: FIELD_TO_PARTICLE

contains
!#######################################################################
!#######################################################################
    subroutine FIELD_TO_PARTICLE(data_in, data_out, l_g, l_q)
        type(pointers3d_dt), intent(in)    :: data_in(:)
        type(pointers_dt),   intent(out)   :: data_out(:)
        type(particle_dt),   intent(inout) :: l_g
        real(wp),            intent(inout) :: l_q(isize_part, inb_part_array)

! -------------------------------------------------------------------
        integer(wi) np_in_grid, np_in_halo_x, np_in_halo_z, np_in_halo_xz
        integer(wi) ip_start, ip_end

!#######################################################################
        if (size(data_in) > inb_part_interp) then
            call TLAB_WRITE_ASCII(efile, 'FIELD_TO_PARTICLE. Not enough memory.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

! -------------------------------------------------------------------
! Setting fields in halo regions
        call Create_Halo_K(data_in)
        call Create_Halo_I_IK(data_in)

! -------------------------------------------------------------------
! Sorting and counting particles for each zone
        ! call Sort_Into_Zones(l_g, l_q, data_out, np_in_grid, np_in_halo_x, np_in_halo_z, np_in_halo_xz)

        ! Group particles inside the grid at the beginning of the array, outside of the grid zone at the end of the arrays
        ip_start = 1
        call Sort_Into_Grid(l_g, l_q, data_out, ip_start, l_g%np, np_in_grid)
        ! From the remaning particles, group particles inside the East halo at the beginning of the array, outside the East halo at the end of the array
        ip_start = ip_start + np_in_grid
        call Sort_Into_Halos(3, l_g, l_q, data_out, ip_start, l_g%np, np_in_halo_x)
        ! From the remaning particles, group particles inside the North halo at the beginning of the array, outside the North zone at the end of the array
        ip_start = ip_start + np_in_halo_x
        call Sort_Into_Halos(1, l_g, l_q, data_out, ip_start, l_g%np, np_in_halo_z)
        ! What remains at the end of the array is in the North-East zone
        np_in_halo_xz = l_g%np - np_in_grid - np_in_halo_x - np_in_halo_z

#ifdef USE_MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Interpolating
        ip_start = 1
        ip_end = np_in_grid
        call Interpolate_Inside_Zones('grid', data_in, data_out, l_g, l_q, ip_start, ip_end)

        if (np_in_halo_x /= 0) then
            ip_start = ip_end + 1
            ip_end = ip_end + np_in_halo_x
            call Interpolate_Inside_Zones('halo_x', p_halo_i, data_out, l_g, l_q, ip_start, ip_end)
        end if

        if (np_in_halo_z /= 0) then
            ip_start = ip_end + 1
            ip_end = ip_end + np_in_halo_z
            call Interpolate_Inside_Zones('halo_z', p_halo_k, data_out, l_g, l_q, ip_start, ip_end)
        end if

        if (np_in_halo_xz /= 0) then
            ip_start = ip_end + 1
            ip_end = ip_end + np_in_halo_xz
            call Interpolate_Inside_Zones('halo_xz', p_halo_ik, data_out, l_g, l_q, ip_start, ip_end)
        end if

        return
    end subroutine FIELD_TO_PARTICLE

!#######################################################################
!#######################################################################
    subroutine Create_Halo_K(data)
        type(pointers3d_dt), intent(in) :: data(:)

        nvar = size(data)

#ifdef USE_MPI
        if (ims_npro_k == 1) then
#endif
            do iv = 1, nvar
                halo_k(1:imax, 1:jmax, 1, iv) = data(iv)%field(1:imax, 1:jmax, kmax)
                halo_k(1:imax, 1:jmax, 2, iv) = data(iv)%field(1:imax, 1:jmax, 1)
            end do

#ifdef USE_MPI
        else
            do iv = 1, nvar
                halo_k(1:imax, 1:jmax, 1, iv) = data(iv)%field(1:imax, 1:jmax, kmax)
                halo_mpi_send_k(1:imax, 1:jmax, iv) = data(iv)%field(1:imax, 1:jmax, 1) ! data to be transfered
            end do

            ims_request(1:ims_npro_k*2) = MPI_REQUEST_NULL
            l = 2*ims_pro_k + 1
            source = mod(ims_pro_k + 1, ims_npro_k)
            dest = mod(ims_pro_k - 1 + ims_npro_k, ims_npro_k)
            buff_size = imax*jmax*nvar
            call MPI_ISEND(halo_mpi_send_k, buff_size, MPI_REAL8, dest, 0, ims_comm_z, ims_request(l), ims_err)
            call MPI_IRECV(halo_mpi_recv_k, buff_size, MPI_REAL8, source, MPI_ANY_TAG, ims_comm_z, ims_request(l + 1), ims_err)
            call MPI_Waitall(ims_npro_k*2, ims_request, ims_status, ims_err)

            halo_k(1:imax, 1:jmax, 2, 1:nvar) = halo_mpi_recv_k(1:imax, 1:jmax, 1:nvar)

        end if
#endif

        return
    end subroutine Create_Halo_K

!#######################################################################
!#######################################################################
    subroutine Create_Halo_I_IK(data)
        type(pointers3d_dt), intent(in) :: data(:)

        nvar = size(data)

#ifdef USE_MPI
        if (ims_npro_i == 1) then
#endif
            do iv = 1, nvar
                halo_i(1, 1:jmax, 1:kmax, iv) = data(iv)%field(imax, 1:jmax, 1:kmax)
                halo_i(2, 1:jmax, 1:kmax, iv) = data(iv)%field(1, 1:jmax, 1:kmax)
                halo_ik(2, 1:jmax, 2, iv) = halo_k(1, 1:jmax, 2, iv) ! top-right corner
            end do

#ifdef USE_MPI
        else
            do iv = 1, nvar
                halo_i(1, 1:jmax, 1:kmax, iv)         = data(iv)%field(imax, 1:jmax, 1:kmax)
                halo_mpi_send_i(1:jmax, 1:kmax, iv)   = data(iv)%field(1, 1:jmax, 1:kmax)   ! data to be transfered
                halo_mpi_send_i(1:jmax, kmax + 1, iv) = halo_k(1, 1:jmax, 2, iv)
            end do

            ims_request(1:ims_npro_i*2) = MPI_REQUEST_NULL
            l = 2*ims_pro_i + 1
            source = mod(ims_pro_i + 1, ims_npro_i)
            dest = mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i)
            buff_size = jmax*(kmax + 1)*nvar
            call MPI_ISEND(halo_mpi_send_i, buff_size, MPI_REAL8, dest, 0, ims_comm_x, ims_request(l), ims_err)
            call MPI_IRECV(halo_mpi_recv_i, buff_size, MPI_REAL8, source, MPI_ANY_TAG, ims_comm_x, ims_request(l + 1), ims_err)
            call MPI_Waitall(ims_npro_i*2, ims_request, ims_status, ims_err)

            halo_i(2, 1:jmax, 1:kmax, 1:nvar) = halo_mpi_recv_i(1:jmax, 1:kmax, 1:nvar)
            halo_ik(2, 1:jmax, 2, 1:nvar) = halo_mpi_recv_i(1:jmax, kmax + 1, 1:nvar) ! top-right corner

        end if
#endif

        halo_ik(1, 1:jmax, 1, 1:nvar) = halo_i(1, 1:jmax, kmax, 1:nvar)
        halo_ik(2, 1:jmax, 1, 1:nvar) = halo_i(2, 1:jmax, kmax, 1:nvar)
        halo_ik(1, 1:jmax, 2, 1:nvar) = halo_k(imax, 1:jmax, 2, 1:nvar)

        return
    end subroutine Create_Halo_I_IK

!########################################################################
!########################################################################
    subroutine Interpolate_Inside_Zones(zone, data_in, data_out, l_g, l_q, ip_start, ip_end)
        character(len=*),    intent(in)  :: zone
        type(pointers3d_dt), intent(in)  :: data_in(:)
        type(pointers_dt),   intent(out) :: data_out(:)
        type(particle_dt),   intent(in)  :: l_g
        real(wp),            intent(in)  :: l_q(isize_part, 3)
        integer(wi),         intent(in)  :: ip_start, ip_end

! -------------------------------------------------------------------
        real(wp) length_g_p(6), cube_g_p(4)
        integer(wi) g_p(10), g1loc, g2loc, g5loc, g6loc
        integer(wi) i
        real(wp) dx_loc_inv, dz_loc_inv

! ######################################################################
        nvar = size(data_out)

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
        case ('halo_x')
            g1loc = 7
            g2loc = 8

        case ('halo_z')
            g5loc = 9
            g6loc = 10

        case ('halo_xz')
            g1loc = 7
            g2loc = 8
            g5loc = 9
            g6loc = 10

        end select

! ######################################################################
        if (g(3)%size /= 1) then

            do i = ip_start, ip_end ! loop over all particles

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

                g_p(3) = l_g%nodes(i)                           ! Local Y position
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
                do iv = 1, nvar
                    data_out(iv)%field(i) = data_out(iv)%field(i) + &
                                           ((cube_g_p(3)*data_in(iv)%field(g_p(g1loc), g_p(3), g_p(g5loc)) &
                                             + cube_g_p(4)*data_in(iv)%field(g_p(g1loc), g_p(4), g_p(g5loc)) &
                                             + cube_g_p(1)*data_in(iv)%field(g_p(g2loc), g_p(4), g_p(g5loc)) &
                                             + cube_g_p(2)*data_in(iv)%field(g_p(g2loc), g_p(3), g_p(g5loc)))*length_g_p(6)) &
                                           + ((cube_g_p(3)*data_in(iv)%field(g_p(g1loc), g_p(3), g_p(g6loc)) &
                                               + cube_g_p(4)*data_in(iv)%field(g_p(g1loc), g_p(4), g_p(g6loc)) &
                                               + cube_g_p(1)*data_in(iv)%field(g_p(g2loc), g_p(4), g_p(g6loc)) &
                                               + cube_g_p(2)*data_in(iv)%field(g_p(g2loc), g_p(3), g_p(g6loc)))*length_g_p(5))
                end do

            end do

! ######################################################################
        else !2D case

            do i = ip_start, ip_end

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
                do iv = 1, nvar
                    data_out(iv)%field(i) = data_out(iv)%field(i) + &
                                           (cube_g_p(3)*data_in(iv)%field(g_p(g1loc), g_p(3), 1) &
                                            + cube_g_p(4)*data_in(iv)%field(g_p(g1loc), g_p(4), 1) &
                                            + cube_g_p(1)*data_in(iv)%field(g_p(g2loc), g_p(4), 1) &
                                            + cube_g_p(2)*data_in(iv)%field(g_p(g2loc), g_p(3), 1))
                end do

            end do

        end if

        return
    end subroutine Interpolate_Inside_Zones

!########################################################################
! Group particles outside the grid zone at the end of the arrays
!########################################################################
    subroutine Sort_Into_Grid(l_g, l_q, data, ip_start, ip_end, counter)
        type(particle_dt), intent(inout) :: l_g
        real(wp),          intent(inout) :: l_q(:, :)
        type(pointers_dt), intent(inout) :: data(:)
        integer(wi),       intent(in)    :: ip_start, ip_end
        integer(wi),       intent(out)   :: counter

        ! -------------------------------------------------------------------
        real(wp) dummy, right_limit, upper_limit
        integer(longi) idummy8
        integer(wi) i, j, idummy

        !#######################################################################
#ifdef USE_MPI
        right_limit = g(1)%nodes(ims_offset_i + imax)  ! right_limit is east
        upper_limit = g(3)%nodes(ims_offset_k + kmax)  ! upper_limit is north
#else
        right_limit = g(1)%nodes(g(1)%size)
        upper_limit = g(3)%nodes(g(3)%size)
#endif

        nvar = size(data)

        ! -------------------------------------------------------------------
        i = ip_start        ! starting particle of sorting algorithm, upward loop
        j = ip_end          ! end particle of sorting algorithm, downward loop

        counter = 0

        do while (i <= j)

            if (l_q(i, 1) > right_limit .or. l_q(i, 3) > upper_limit) then          ! particle i is outside, look for particles inside to swap with
                do
                    if (l_q(j, 1) > right_limit .or. l_q(j, 3) > upper_limit) then  ! particle j is outside, leave it here
                        j = j - 1 
                        if (i >= j) exit                                            ! finished your upwards loop

                    else                                                            ! particle j is inside, so swap
                        idummy = l_g%nodes(i)
                        l_g%nodes(i) = l_g%nodes(j)
                        l_g%nodes(j) = idummy

                        idummy8 = l_g%tags(i)
                        l_g%tags(i) = l_g%tags(j)
                        l_g%tags(j) = idummy8

                        do iv = 1, inb_part_array
                            dummy = l_q(i, iv)
                            l_q(i, iv) = l_q(j, iv)
                            l_q(j, iv) = dummy
                        end do

                        do iv = 1, nvar      ! swap also this data because of the substages in the Runge-Kutta
                            dummy = data(iv)%field(i)
                            data(iv)%field(i) = data(iv)%field(j)
                            data(iv)%field(j) = dummy
                        end do

                        j = j - 1
                        counter = counter + 1
                        exit

                    end if
                end do

            else            ! particle i is inside, the right place
                counter = counter + 1

            end if
            i = i + 1       ! go up to the next particle

        end do

        return
    end subroutine Sort_Into_Grid

!########################################################################
!########################################################################
    subroutine Sort_Into_Halos(x_or_z, l_g, l_q, data, ip_start, ip_end, counter)
        integer,           intent(in)    :: x_or_z
        type(particle_dt), intent(inout) :: l_g
        real(wp),          intent(inout) :: l_q(:, :)
        type(pointers_dt), intent(inout) :: data(:)
        integer(wi),       intent(in)    :: ip_start, ip_end
        integer(wi),       intent(out)   :: counter

        ! -------------------------------------------------------------------
        real(wp) dummy, limit
        integer(longi) idummy8
        integer(wi) i, j, idummy

        !#######################################################################
        select case (x_or_z)
        case (1)
#ifdef USE_MPI
            limit = g(1)%nodes(ims_offset_i + imax)  ! right_limit is east
#else
            limit = g(1)%nodes(g(1)%size)
#endif

        case (3)
#ifdef USE_MPI
            limit = g(3)%nodes(ims_offset_k + kmax)  ! upper_limit is north
#else
            limit = g(3)%nodes(g(3)%size)
#endif

        end select

        nvar = size(data)

        ! -------------------------------------------------------------------
        i = ip_start
        j = ip_end

        counter = 0

        do while (i <= j)
            if (l_q(i, x_or_z) > limit) then            ! particle i is outside, look for particles inside to swap with
                do
                    if (l_q(j, x_or_z) > limit) then    ! particle j is outside, leave it here
                        j = j - 1
                        if (i >= j) exit

                    else                                ! particle j is inside, so swap
                        idummy = l_g%nodes(i)
                        l_g%nodes(i) = l_g%nodes(j)
                        l_g%nodes(j) = idummy

                        idummy8 = l_g%tags(i)
                        l_g%tags(i) = l_g%tags(j)
                        l_g%tags(j) = idummy8

                        do iv = 1, inb_part_array
                            dummy = l_q(i, iv)
                            l_q(i, iv) = l_q(j, iv)
                            l_q(j, iv) = dummy
                        end do

                        do iv = 1, nvar
                            dummy = data(iv)%field(i)
                            data(iv)%field(i) = data(iv)%field(j)
                            data(iv)%field(j) = dummy
                        end do

                        j = j - 1
                        counter = counter + 1
                        exit

                    end if
                end do

            else                                        ! particle i is inside
                counter = counter + 1

            end if
            i = i + 1

        end do

        return
    end subroutine Sort_Into_Halos

! !########################################################################
! ! Sorting structure grid-halo_x-halo_z-halo_diagonal
! !########################################################################
!     subroutine Sort_Into_Zones(l_g, l_q, data, np_in_grid, np_in_halo_x, np_in_halo_z, np_in_halo_xz)
!         type(pointers_dt), intent(inout) :: data(:)
!         type(particle_dt), intent(inout) :: l_g
!         real(wp), intent(inout) :: l_q(:, :)
!         integer(wi), intent(out) :: np_in_grid, np_in_halo_x, np_in_halo_z, np_in_halo_xz

!         ! -------------------------------------------------------------------
!         real(wp) dummy, right_limit, upper_limit
!         integer(wi) idummy, nvar
!         integer(longi) idummy8
!         integer(wi) i, j, k

!         !#######################################################################
! #ifdef USE_MPI
!         right_limit = g(1)%nodes(ims_offset_i + imax)  ! right_limit is east
!         upper_limit = g(3)%nodes(ims_offset_k + kmax)  ! upper_limit is north
! #else
!         right_limit = g(1)%nodes(g(1)%size)
!         upper_limit = g(3)%nodes(g(3)%size)
! #endif

!         nvar = size(data)

!         !#######################################################################
!         ! Group together particles outside of the grid zone at the end of the arrays
!         i = 1       ! starting point of sorting algorithm
!         j = l_g%np  ! end point of sorting algorithm

!         np_in_grid = 0

!         do while (i <= j)

!             if (l_q(i, 1) > right_limit .or. l_q(i, 3) > upper_limit) then          ! particle i is in halo
!                 do
!                     if (l_q(j, 1) > right_limit .or. l_q(j, 3) > upper_limit) then  ! partcile j is in halo, leave it there
!                         j = j - 1                                                   ! go to next particle
!                         if (i >= j) exit                                            ! finished your upwards loop

!                     else                                                            ! found a particle in grid, so swap
!                         idummy = l_g%nodes(i)
!                         l_g%nodes(i) = l_g%nodes(j)
!                         l_g%nodes(j) = idummy

!                         idummy8 = l_g%tags(i)
!                         l_g%tags(i) = l_g%tags(j)
!                         l_g%tags(j) = idummy8

!                         do k = 1, inb_part_array
!                             dummy = l_q(i, k)
!                             l_q(i, k) = l_q(j, k)
!                             l_q(j, k) = dummy
!                         end do

!                         do k = 1, nvar      ! swap also this data because of the substages in the Runge-Kutta
!                             dummy = data(k)%field(i)
!                             data(k)%field(i) = data(k)%field(j)
!                             data(k)%field(j) = dummy
!                         end do

!                         j = j - 1
!                         np_in_grid = np_in_grid + 1
!                         exit

!                     end if
!                 end do

!             else            ! particle i is in the grid, the right place
!                 np_in_grid = np_in_grid + 1

!             end if
!             i = i + 1       ! go to next particle

!         end do

!         ! !Last particle might not be properly checked. We check it here.
!         ! if (i == j) then !Probably not needed
!         !     if ((l_q(i, 1) > right_limit) .or. (l_q(i, 3) > upper_limit)) then !Particle out the grid
!         !         !Do nothing
!         !     else
!         !         np_in_grid = np_in_grid + 1 !The last particle was not checked but it is in the grid
!         !     end if
!         ! end if
!         !Possible optimization to avoid if i. EQ. j. Not completed.
!         ! ELSE IF (i-1 .GT. j+1) THEN
!         !   j=j+1
!         !   i=i-1
!         !   np_in_grid=np_in_grid-1
!         !   IF ( (l_q(j,1) .GT. right_limit) .OR. l_q(j,3) .GT. upper_limit) THEN
!         !          DO k=1,3
!         !           dummy=l_q(i,k)
!         !           l_q(i,k)=l_q(j,k)
!         !           l_q(j,k)=dummy

!         !           idummy8=l_g%tags(i)
!         !           l_g%tags(i)=l_g%tags(j)
!         !           l_g%tags(j)=idummy8
!         !         END DO
!         !   END IF

!         ! -------------------------------------------------------------------
!         ! From the remaning particles, group together particles outside the East zone at the end of the array
!         i = np_in_grid + 1
!         j = l_g%np

!         np_in_halo_x = 0

!         do while (i <= j)
!             if (l_q(i, 3) > upper_limit) then           ! particle i is in North
!                 do
!                     if (l_q(j, 3) > upper_limit) then   ! particle j is in North, leave it here
!                         j = j - 1
!                         if (i >= j) exit

!                     else                                ! found a particle in East, so swap
!                         idummy = l_g%nodes(i)
!                         l_g%nodes(i) = l_g%nodes(j)
!                         l_g%nodes(j) = idummy

!                         idummy8 = l_g%tags(i)
!                         l_g%tags(i) = l_g%tags(j)
!                         l_g%tags(j) = idummy8

!                         do k = 1, inb_part_array
!                             dummy = l_q(i, k)
!                             l_q(i, k) = l_q(j, k)
!                             l_q(j, k) = dummy
!                         end do

!                         do k = 1, nvar
!                             dummy = data(k)%field(i)
!                             data(k)%field(i) = data(k)%field(j)
!                             data(k)%field(j) = dummy
!                         end do

!                         j = j - 1
!                         np_in_halo_x = np_in_halo_x + 1
!                         exit

!                     end if
!                 end do

!             else
!                 np_in_halo_x = np_in_halo_x + 1

!             end if
!             i = i + 1

!         end do

!         ! !Last particle might not be properly checked. We check it here.
!         ! if (i == j) then !Probably not needed
!         !     if (l_q(i, 3) > upper_limit) then !Particle is out North
!         !         !Do nothing
!         !     else
!         !         np_in_halo_x = np_in_halo_x + 1
!         !     end if
!         ! end if

!         ! -------------------------------------------------------------------
!         ! From the remaning particles, group together particles outside the North zone at the end of the array
!         i = np_in_grid + np_in_halo_x + 1
!         j = l_g%np

!         np_in_halo_z = 0

!         do while (i <= j)
!             if (l_q(i, 1) > right_limit) then           ! particle i is in North-east
!                 do
!                     if (l_q(j, 1) > right_limit) then   ! particle j is in North-east, leave it here
!                         j = j - 1
!                         if (i >= j) exit

!                     else                                ! found a particle in North, so swap
!                         idummy = l_g%nodes(i)
!                         l_g%nodes(i) = l_g%nodes(j)
!                         l_g%nodes(j) = idummy

!                         idummy8 = l_g%tags(i)
!                         l_g%tags(i) = l_g%tags(j)
!                         l_g%tags(j) = idummy8

!                         do k = 1, inb_part_array
!                             dummy = l_q(i, k)
!                             l_q(i, k) = l_q(j, k)
!                             l_q(j, k) = dummy
!                         end do

!                         do k = 1, nvar
!                             dummy = data(k)%field(i)
!                             data(k)%field(i) = data(k)%field(j)
!                             data(k)%field(j) = dummy
!                         end do

!                         j = j - 1
!                         np_in_halo_z = np_in_halo_z + 1
!                         exit

!                     end if
!                 end do

!             else
!                 np_in_halo_z = np_in_halo_z + 1

!             end if
!             i = i + 1

!         end do

!         ! !Last particle might not be properly checked. We check it here.
!         ! if (i == j) then !Probably not needed
!         !     if (l_q(i, 1) > right_limit) then !Particle is out East
!         !         !Do nothing
!         !     else
!         !         np_in_halo_z = np_in_halo_z + 1
!         !     end if
!         ! end if

!         ! -------------------------------------------------------------------
!         ! What remains at the end of the array is in the North-East zone
!         np_in_halo_xz = l_g%np - np_in_grid - np_in_halo_x - np_in_halo_z

!         return
!     end subroutine Sort_Into_Zones

end module PARTICLE_INTERPOLATE

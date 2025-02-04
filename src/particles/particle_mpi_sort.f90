#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!#  Sorting of particles into grid-west-east
!#  Sending to processors in west and east
!#  Sorting of particles into grid-south-north
!#  Sending to processors in south and north
!#
!########################################################################
subroutine PARTICLE_MPI_SORT(x_or_z, l_g, l_q, l_hq, &
                             nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)
    use TLab_Constants, only: wp, wi, longi
    use TLab_Memory, only: imax, kmax
    use PARTICLE_VARS, only: isize_part, inb_part_array, inb_part
    use FDM, only: g
    use PARTICLE_TYPES, only: particle_dt
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    integer(wi) nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north, x_or_z
    type(particle_dt) :: l_g
    real(wp), dimension(isize_part, *) :: l_q
    real(wp), dimension(isize_part, *) :: l_hq

! -------------------------------------------------------------------
    real(wp) dx_grid, dz_grid
    real(wp) dummy, lower_limit, upper_limit
    integer(wi) nzone_west_south, nzone_east_north
    integer(wi) counter_swap, idummy
    integer(longi) idummy8
    integer(wi) i, j, k

!#######################################################################
    nzone_west_south = 0
    nzone_east_north = 0

    nzone_west = 0
    nzone_east = 0
    nzone_grid = 0

    dx_grid = g(1)%nodes(2) - g(1)%nodes(1)  ! Distance between gridpoints, to deal with periodicity
    dz_grid = g(3)%nodes(2) - g(3)%nodes(1)

    select case (x_or_z)

    case (1) !Sort in West-East direction
        lower_limit = g(1)%nodes(ims_offset_i + 1)              !lower_limit is West
        upper_limit = g(1)%nodes(ims_offset_i + imax) + dx_grid !upper_limit is East
    case (3) !Sort in South-North direction
        lower_limit = g(3)%nodes(ims_offset_k + 1)              !lower_limit is south
        upper_limit = g(3)%nodes(ims_offset_k + kmax) + dz_grid !upper_limit is north
    end select

!#######################################################################
!Sorting structure grid-west-east or grid-south-north
!First Algorythm sorts all grid-particle into first part of particle
!#######################################################################
    i = 1       !Starting point of sorting algorithm
    j = l_g%np  !End point of sorting algorithm

    do while (i < j)
        if (l_q(i, x_or_z) < lower_limit) then      !If particle is out to West

            counter_swap = 0                        !Particle must be swapped
            do while (counter_swap == 0)            !Enter in the swaping region
                if ((l_q(j, x_or_z) < lower_limit) .or. (l_q(j, x_or_z) > upper_limit)) then !Partcile here in right area
                    j = j - 1                       !Go to next particle in the swaping region
                    if (i == j) then                !You finished your upwards loop
                        counter_swap = 1
                    end if
                else                                ! Found a particle which belongs to grid so SWAP
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

                    do k = 1, inb_part
                        dummy = l_hq(i, k)
                        l_hq(i, k) = l_hq(j, k)
                        l_hq(j, k) = dummy
                    end do

                    j = j - 1
                    nzone_grid = nzone_grid + 1
                    counter_swap = 1
                end if
            end do
            i = i + 1                               !Go to next particle

        elseif (l_q(i, x_or_z) > upper_limit) then  ! If particle is out to the east
            counter_swap = 0
            do while (counter_swap == 0)            !Same procedure as above for east
!  IF (i .EQ. j ) THEN
!  counter_swap=1
                if ((l_q(j, x_or_z) < lower_limit) .or. (l_q(j, x_or_z) > upper_limit)) then
                    j = j - 1
                    if (i == j) then                !You finished your upwards loop
                        counter_swap = 1
                    end if
                else
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

                    do k = 1, inb_part
                        dummy = l_hq(i, k)
                        l_hq(i, k) = l_hq(j, k)
                        l_hq(j, k) = dummy
                    end do

                    j = j - 1
                    nzone_grid = nzone_grid + 1
                    counter_swap = 1
                end if
            end do
            i = i + 1

        else  ! Particle is in the grid
            i = i + 1
            nzone_grid = nzone_grid + 1

        end if

    end do

!Last particle might not be properly checked. We check it here.
    if (i == j) then !Probably not needed
        if ((l_q(i, x_or_z) < lower_limit) .or. (l_q(i, x_or_z) > upper_limit)) then !Particle out the grid
!Do nothing
        else
            nzone_grid = nzone_grid + 1 !The last particle was not checked but it is in the grid
        end if

    end if

! -------------------------------------------------------------------
    j = l_g%np
    i = nzone_grid + 1

    do while (i < j)
        if (l_q(i, x_or_z) < lower_limit) then
            i = i + 1
            nzone_west_south = nzone_west_south + 1
        else  !particle is out to the east
            counter_swap = 0
            do while (counter_swap == 0)
                if (l_q(j, x_or_z) < lower_limit) then !if particle is out to west
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

                    do k = 1, inb_part
                        dummy = l_hq(i, k)
                        l_hq(i, k) = l_hq(j, k)
                        l_hq(j, k) = dummy
                    end do

                    counter_swap = 1
                    j = j - 1
                    nzone_west_south = nzone_west_south + 1

                else
                    j = j - 1
                    if (i == j) then !You finished your upwards loop
                        counter_swap = 1
                    end if
                end if
            end do
            i = i + 1
        end if
    end do

!Last particle might not be properly checked. We check it here.
    if (i == j) then !Probably not needed
        if (l_q(i, x_or_z) > upper_limit) then !Particle is out east
!Do nothing
        else
            nzone_west_south = nzone_west_south + 1
        end if
    end if

!Calculating the number of particles send to east or north
    nzone_east_north = l_g%np - nzone_grid - nzone_west_south

    if (x_or_z == 1) then
        nzone_west = nzone_west_south
        nzone_east = nzone_east_north

    elseif (x_or_z == 3) then
        nzone_south = nzone_west_south
        nzone_north = nzone_east_north
    end if

    return
end subroutine PARTICLE_MPI_SORT

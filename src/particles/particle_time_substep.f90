#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!########################################################################
subroutine PARTICLE_TIME_SUBSTEP(dte, l_hq, l_comm)

    use TLAB_VARS,     only: g
    use TLAB_ARRAYS
    use PARTICLE_VARS, only: isize_part, inb_part
    use PARTICLE_VARS, only: isize_l_comm
    use PARTICLE_ARRAYS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
    use PARTICLE_VARS, only: isize_pbuffer
#endif

    implicit none

    real(wp) dte
    real(wp), dimension(isize_part, *) :: l_hq
    real(wp), dimension(isize_l_comm), target :: l_comm

! -------------------------------------------------------------------
    integer(wi) is, i

#ifdef USE_MPI
    real(wp), dimension(:), pointer :: p_buffer_1, p_buffer_2
    integer(wi) nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
#else
    real(wp) x_right, z_right
#endif

!#####################################################################
    call RHS_PARTICLE_GLOBAL(q, s, txc, l_q, l_hq, l_txc, l_comm, wrk1d, wrk2d, wrk3d)

!#######################################################################
! Update particle properties
!#######################################################################
    do is = 1, inb_part
        l_q(1:l_g%np, is) = l_q(1:l_g%np, is) + dte*l_hq(1:l_g%np, is)
    end do

!#####################################################################
! Boundary control to see if particles leave processor
!#####################################################################
#ifdef USE_MPI
    p_buffer_1(1:isize_pbuffer) => l_comm(1:isize_pbuffer)
    p_buffer_2(1:isize_pbuffer) => l_comm(isize_pbuffer + 1:isize_pbuffer*2)

! -------------------------------------------------------------------
!Particle sorting for Send/Recv X-Direction
! -------------------------------------------------------------------
    call PARTICLE_SORT(1, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    if (ims_pro_i == 0) then !Take care of periodic boundary conditions west
        if (nzone_west /= 0) then
            l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) = &
                l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) + g(1)%scale
        end if
    end if

    if (ims_pro_i == (ims_npro_i - 1)) then !Take care of periodic boundary conditions east
        if (nzone_east /= 0) then
            l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) = &
                l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) - g(1)%scale
        end if
    end if

    call PARTICLE_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, &
                              p_buffer_1, p_buffer_2, l_q, l_hq, l_g%tags, l_g%np)

! -------------------------------------------------------------------
!Particle sorting for Send/Recv Z-Direction
! -------------------------------------------------------------------
    call PARTICLE_SORT(3, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    if (ims_pro_k == 0) then !Take care of periodic boundary conditions south
        if (nzone_south /= 0) then
            l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) = &
                l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) + g(3)%scale
        end if
    end if

    if (ims_pro_k == (ims_npro_k - 1)) then !Take care of periodic boundary conditions north
        if (nzone_north /= 0) then
            l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) = &
                l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) - g(3)%scale
        end if
    end if

    call PARTICLE_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, &
                              p_buffer_1, p_buffer_2, l_q, l_hq, l_g%tags, l_g%np)

    nullify (p_buffer_1, p_buffer_2)

#else
!#######################################################################
! Serial
    x_right = g(1)%nodes(1) + g(1)%scale
    z_right = g(3)%nodes(1) + g(3)%scale

    do i = 1, l_g%np
        if (l_q(i, 1) > x_right) then
            l_q(i, 1) = l_q(i, 1) - g(1)%scale

        elseif (l_q(i, 1) < g(1)%nodes(1)) then
            l_q(i, 1) = l_q(i, 1) + g(1)%scale

        end if

        if (l_q(i, 3) > z_right) then
            l_q(i, 3) = l_q(i, 3) - g(3)%scale

        elseif (l_q(i, 3) < g(3)%nodes(1)) then
            l_q(i, 3) = l_q(i, 3) + g(3)%scale

        end if
    end do

#endif

!#######################################################################
! Recalculating closest node below in Y direction
!#######################################################################
    call PARTICLE_LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, g(2)%size, g(2)%nodes)

    return
end subroutine PARTICLE_TIME_SUBSTEP

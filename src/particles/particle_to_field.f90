#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!#
!# Interpolate particle information to a field
!#
!########################################################################
subroutine PARTICLE_TO_FIELD(l_q, particle_property, field_out, wrk3d)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use PARTICLE_VARS, only: isize_part
#ifdef USE_MPI
    use PARTICLE_ARRAYS, only: l_work
    use MPI
    use TLabMPI_VARS, only: ims_err
#endif

    implicit none

    real(wp), dimension(isize_part, 3) :: l_q
    real(wp), dimension(isize_part) :: particle_property
    real(wp), dimension(imax, jmax, kmax) :: field_out
    real(wp), dimension(imax + 1, jmax, kmax + 1) :: wrk3d

!#######################################################################
!Write data to fielf
!Field is an extended field with grid, halo1, halo2 and halo3
!#######################################################################
    wrk3d = 0.0_wp
    call PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, wrk3d)

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!#######################################################################
!SEND to the East and RECV from the West
!Sum the most left colum of field
!SEND to the North and RECV from the South
!Sum the lowest row of field
!#######################################################################
    call PARTICLE_TO_FIELD_SEND_RECV_EAST(l_work(1), l_work((jmax*(kmax + 1)) + 1), wrk3d)
    call PARTICLE_TO_FIELD_SEND_RECV_NORTH(l_work(1), l_work(((imax + 1)*jmax) + 1), wrk3d)

#endif

!#######################################################################
!Put wrk3d into continuous memory field_out
!#######################################################################
    field_out(1:imax, 1:jmax, 1:kmax) = wrk3d(1:imax, 1:jmax, 1:kmax)

    return
end subroutine PARTICLE_TO_FIELD

!########################################################################
!########################################################################
subroutine PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, field)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use PARTICLE_VARS, only: isize_part
    use PARTICLE_ARRAYS, only: l_g
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    real(wp), dimension(imax + 1, jmax, kmax + 1) :: field
    real(wp), dimension(isize_part, 3) :: l_q
    real(wp), dimension(isize_part) :: particle_property

    real(wp) length_g_p(6), cube_g_p(4)
    integer(wi) g_p(6)
    integer(wi) i
    real(wp) dx_loc_inv, dz_loc_inv

! ######################################################################
    dx_loc_inv = real(g(1)%size,wp)/g(1)%scale
    dz_loc_inv = real(g(3)%size,wp)/g(3)%scale

    g_p(5) = 1 ! Default is 2D
    g_p(6) = 1

! ######################################################################
    do i = 1, l_g%np

        length_g_p(1) = l_q(i, 1)*dx_loc_inv            ! Local X position
        g_p(1) = floor(length_g_p(1))
        length_g_p(1) = length_g_p(1) - real(g_p(1),wp)
#ifdef USE_MPI
        g_p(1) = g_p(1) + 1 - ims_offset_i
#else
        g_p(1) = g_p(1) + 1
#endif
        g_p(2) = g_p(1) + 1
        length_g_p(2) = 1.0_wp - length_g_p(1)

        if (g(3)%size /= 1) then
            length_g_p(5) = l_q(i, 3)*dz_loc_inv            ! Local Z position
            g_p(5) = floor(length_g_p(5))
            length_g_p(5) = length_g_p(5) - real(g_p(5),wp)
#ifdef USE_MPI
            g_p(5) = g_p(5) + 1 - ims_offset_k
#else
            g_p(5) = g_p(5) + 1
#endif
            g_p(6) = g_p(5) + 1
            length_g_p(6) = 1.0_wp - length_g_p(5)
        end if

        g_p(3) = l_g%nodes(i)                    ! Local Y position
        g_p(4) = g_p(3) + 1
        length_g_p(3) = (l_q(i, 2) - g(2)%nodes(g_p(3)))/(g(2)%nodes(g_p(4)) - g(2)%nodes(g_p(3)))
        length_g_p(4) = 1.0_wp - length_g_p(3)

        cube_g_p(1) = length_g_p(1)*length_g_p(3) ! bilear cubes for X and Y
        cube_g_p(2) = length_g_p(1)*length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4)*length_g_p(2)
        cube_g_p(4) = length_g_p(2)*length_g_p(3)

!###################################################################
!Two bilinear calculation for each k direction (g_p(5) and g_p(6))
!Then multipled by (1-length) for Trilinear aspect
!###################################################################
        if (g(3)%size == 1) then
!U
            field(g_p(1), g_p(3), g_p(5)) = field(g_p(1), g_p(3), g_p(5)) + particle_property(i)*cube_g_p(3)
            field(g_p(1), g_p(4), g_p(5)) = field(g_p(1), g_p(4), g_p(5)) + particle_property(i)*cube_g_p(4)
            field(g_p(2), g_p(4), g_p(5)) = field(g_p(2), g_p(4), g_p(5)) + particle_property(i)*cube_g_p(1)
            field(g_p(2), g_p(3), g_p(5)) = field(g_p(2), g_p(3), g_p(5)) + particle_property(i)*cube_g_p(2)

        else
!Z 1
            field(g_p(1), g_p(3), g_p(5)) = field(g_p(1), g_p(3), g_p(5)) + particle_property(i)*cube_g_p(3)*length_g_p(6)
            field(g_p(1), g_p(4), g_p(5)) = field(g_p(1), g_p(4), g_p(5)) + particle_property(i)*cube_g_p(4)*length_g_p(6)
            field(g_p(2), g_p(4), g_p(5)) = field(g_p(2), g_p(4), g_p(5)) + particle_property(i)*cube_g_p(1)*length_g_p(6)
            field(g_p(2), g_p(3), g_p(5)) = field(g_p(2), g_p(3), g_p(5)) + particle_property(i)*cube_g_p(2)*length_g_p(6)

!Z 2
            field(g_p(1), g_p(3), g_p(6)) = field(g_p(1), g_p(3), g_p(6)) + particle_property(i)*cube_g_p(3)*length_g_p(5)
            field(g_p(1), g_p(4), g_p(6)) = field(g_p(1), g_p(4), g_p(6)) + particle_property(i)*cube_g_p(4)*length_g_p(5)
            field(g_p(2), g_p(4), g_p(6)) = field(g_p(2), g_p(4), g_p(6)) + particle_property(i)*cube_g_p(1)*length_g_p(5)
            field(g_p(2), g_p(3), g_p(6)) = field(g_p(2), g_p(3), g_p(6)) + particle_property(i)*cube_g_p(2)*length_g_p(5)

        end if

    end do

    return
end subroutine PARTICLE_TO_FIELD_INTERPOLATE

#ifdef USE_MPI

!########################################################################
!#
!# Sends field_halo to neighbouring processors
!# Field has size imax+1 and kmax+1
!# Only halo columns are needed for interpolation to field
!#
!########################################################################
subroutine PARTICLE_TO_FIELD_SEND_RECV_EAST(f_buffer_1, f_buffer_2, field)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use MPI
    use TLabMPI_VARS

    implicit none

    integer(wi) source_east
    integer(wi) dest_east
    integer(wi) l
    integer(wi) mpireq(ims_npro*2)
    integer(wi) status(MPI_STATUS_SIZE, ims_npro*2)

    real(wp), dimension(imax + 1, jmax, kmax + 1) :: field
    real(wp), dimension(jmax, kmax + 1) :: f_buffer_1, f_buffer_2

    f_buffer_1 = 0.0_wp
    f_buffer_2 = 0.0_wp

!###################################################################
!Calculating the source/dest for east and north
!###################################################################
    if (ims_pro_i == 0) then ! Fisrt row
        dest_east = ims_pro_i + 1 + ims_npro_i*ims_pro_k !Destination of the message - to East
        source_east = ims_npro_i - 1 + ims_npro_i*ims_pro_k !Source of the message - from West

    else if (ims_pro_i == (ims_npro_i - 1)) then !Last row
        dest_east = ims_npro_i*ims_pro_k
        source_east = ims_pro_i - 1 + ims_npro_i*ims_pro_k

    else !Any case
        dest_east = ims_pro_i + 1 + ims_npro_i*ims_pro_k !Dest of the message
        source_east = ims_pro_i - 1 + ims_npro_i*ims_pro_k !source of the message
    end if

!#################################################################
!Setting up the plane which needs to be send EAST
!Send this plane EAST and receive the plane from WEST
!#################################################################
!The buffer is one point larger because the diagonal point is first send EAST
!and then NORTH
    f_buffer_1(1:jmax, 1:(kmax + 1)) = field(imax + 1, 1:jmax, 1:(kmax + 1))

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(f_buffer_1, jmax*(kmax + 1), MPI_REAL8, dest_east, 0, MPI_COMM_WORLD, mpireq(l), ims_err) ! Send particles to the west in buffer2
    call MPI_IRECV(f_buffer_2, jmax*(kmax + 1), MPI_REAL8, source_east, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err) !Get particles from the east in buffer1

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

    field(1, 1:jmax, 1:(kmax + 1)) = field(1, 1:jmax, 1:(kmax + 1)) + f_buffer_2(1:jmax, 1:(kmax + 1))

    return
end subroutine PARTICLE_TO_FIELD_SEND_RECV_EAST

!########################################################################
!########################################################################
subroutine PARTICLE_TO_FIELD_SEND_RECV_NORTH(f_buffer_1, f_buffer_2, field)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use MPI
    use TLabMPI_VARS

    implicit none

    integer(wi) source_north
    integer(wi) dest_north
    integer(wi) l
    integer(wi) mpireq(ims_npro*2)
    integer(wi) status(MPI_STATUS_SIZE, ims_npro*2)

    real(wp), dimension(imax + 1, jmax, kmax + 1) :: field
    real(wp), dimension(imax + 1, jmax) :: f_buffer_1, f_buffer_2

    f_buffer_1 = 0.0_wp
    f_buffer_2 = 0.0_wp

!###################################################################
!Calculating the source/dest for east and north
!###################################################################
    if (ims_pro_k == 0) then ! Fisrt row
        dest_north = ims_pro_i + ims_npro_i*(ims_pro_k + 1) !dest of the message - to east
        source_north = ims_pro_i + ims_npro_i*(ims_npro_k - 1) !source of the message - from wesst

    else if (ims_pro_k == (ims_npro_k - 1)) then !Last row
        dest_north = ims_pro_i
        source_north = ims_pro_i + ims_npro_i*(ims_pro_k - 1)

    else !Any case
        dest_north = ims_pro_i + ims_npro_i*(ims_pro_k + 1) !Dest of the message
        source_north = ims_pro_i + ims_npro_i*(ims_pro_k - 1)  !source of the message

    end if

!#################################################################
!Setting up the plane which needs to be send NORTH
!Send this plane NORTH and receive the plane from SOUTH
!#################################################################
    f_buffer_1(1:(imax + 1), 1:jmax) = field(1:(imax + 1), 1:jmax, kmax + 1)

!CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(f_buffer_1, (imax + 1)*jmax, MPI_REAL8, dest_north, 0, MPI_COMM_WORLD, mpireq(l), ims_err) !Send p_buffer_2 to the east
    call MPI_IRECV(f_buffer_2, (imax + 1)*jmax, MPI_REAL8, source_north, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err) !Get particles from the west into p_buffer_1

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

    field(1:(imax + 1), 1:jmax, 1) = field(1:(imax + 1), 1:jmax, 1) + f_buffer_2(1:(imax + 1), 1:jmax)

    return
end subroutine PARTICLE_TO_FIELD_SEND_RECV_NORTH

#endif

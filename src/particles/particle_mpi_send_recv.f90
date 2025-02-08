!########################################################################
!#
!# Send particles to neighbouring processors
!# Subroutine is divided into 2 optional parts
!# One for East - West direction or South - North
!# First amount of to be sended particles is communicated
!# Them the actual send is realized
!#
!########################################################################

!#######################################################################
!#######################################################################
subroutine PARTICLE_MPI_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, l_q, l_hq, l_tags, particle_number)
    use TLab_Constants, only: wp, wi
    use PARTICLE_VARS, only: isize_part, inb_part_array, inb_part
    use PARTICLE_ARRAYS, only: p_buffer_1, p_buffer_2
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_k, ims_npro_k, ims_pro, ims_npro, ims_err

    implicit none

    integer(wi) nzone_grid, nzone_west, nzone_east

    real(wp), dimension(isize_part, *) :: l_q
    real(wp), dimension(isize_part, *) :: l_hq
    real(wp), dimension(isize_part) :: l_tags !Attention. Chosen real(wp) on purpose.
    integer(wi) particle_number

! -------------------------------------------------------------------
    integer(wi) nzone_send_west, nzone_send_east
    integer(wi) source_west, source_east
    integer(wi) dest_west, dest_east
    integer(wi) l
    type(MPI_Request) mpireq(ims_npro*4)
    type(MPI_Status) status(ims_npro*4)
    integer(wi) i, k, k_loc, m

!#######################################################################
    m = (inb_part_array*2) + 1 !Sending size of the buffer_parts

    particle_number = nzone_grid

!###################################################################
!Setting up destination and source + send/recv number of particles which will be send
!###################################################################
    if (ims_pro_i == 0) then ! Fisrt row
        dest_west = ims_npro_i - 1 + ims_npro_i*ims_pro_k !Destination of the message - to West
        source_west = ims_pro_i + 1 + ims_npro_i*ims_pro_k !Source of the message - from East

    else if (ims_pro_i == (ims_npro_i - 1)) then !Last row
        dest_west = ims_pro_i - 1 + ims_npro_i*ims_pro_k
        source_west = ims_npro_i*ims_pro_k

    else !Any case
        dest_west = ims_pro_i - 1 + ims_npro_i*ims_pro_k
        source_west = ims_pro_i + 1 + ims_npro_i*ims_pro_k
    end if

    mpireq(1:ims_npro*4) = MPI_REQUEST_NULL
    l = 4*ims_pro + 1
!Amount of particles which will be sent/received
    call MPI_ISEND(nzone_west, 1, MPI_INTEGER4, dest_west, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
    call MPI_IRECV(nzone_send_west, 1, MPI_INTEGER4, source_west, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)

!    CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err) !ONLY ONE MPI_WAITALL FOR THE 4 COMMUNICATIONS

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

!    mpireq(1:ims_npro*2)=MPI_REQUEST_NULL  !ONLY ONE MPI_WAITALL NOW
!    l = 2*ims_pro +1
!Amount of particles which will be sent/received
    call MPI_ISEND(nzone_east, 1, MPI_INTEGER4, dest_east, 0, MPI_COMM_WORLD, mpireq(l + 2), ims_err)
    call MPI_IRECV(nzone_send_east, 1, MPI_INTEGER4, source_east, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 3), ims_err)

    call MPI_Waitall(ims_npro*4, mpireq, status, ims_err)

!#################################################################
! Send nzone_west particles to west and get nzone_send_west particles from east
!#################################################################
!I get nzone_send_west particles form the east
    particle_number = particle_number + nzone_send_west

!Construct sending buffer to west in p_buffer_2
    if (nzone_west /= 0) then !I have something to send

!Send particles and co. west
        k = 0  ! Array index 1,3 for positions, 4,6 for rhs velocities and 7 for id
        do k_loc = 1, inb_part_array !Position
            k = k + 1
            do i = 1, nzone_west
                p_buffer_2(i + ((k - 1)*nzone_west)) = l_q(nzone_grid + i, k_loc)
            end do
        end do

        do k_loc = 1, inb_part !Right hand side
            k = k + 1
            do i = 1, nzone_west
                p_buffer_2(i + ((k - 1)*nzone_west)) = l_hq(nzone_grid + i, k_loc)
            end do

        end do

        k = k + 1
        do i = 1, nzone_west
            p_buffer_2(i + ((k - 1)*nzone_west)) = l_tags(nzone_grid + i)
        end do

    end if

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(p_buffer_2, nzone_west*m, MPI_REAL8, dest_west, 0, MPI_COMM_WORLD, mpireq(l), ims_err) ! Send particles to the west in buffer2
    call MPI_IRECV(p_buffer_1, nzone_send_west*m, MPI_REAL8, source_west, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err) !Get particles from the east in buffer1

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

!Construct sending buffer to east in p_buffer_2
    if (nzone_east /= 0) then  ! I have something to send to the east

!Send particles and co. east
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_east
                p_buffer_2(i + ((k - 1)*nzone_east)) = l_q(nzone_grid + nzone_west + i, k_loc)
            end do

        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_east
                p_buffer_2(i + ((k - 1)*nzone_east)) = l_hq(nzone_grid + nzone_west + i, k_loc)
            end do

        end do

        k = k + 1
        do i = 1, nzone_east
            p_buffer_2(i + ((k - 1)*nzone_east)) = l_tags(nzone_grid + nzone_west + i)
        end do

    end if

!Write nzone_send_west recieved particles from east (step one) that are in p_buffer_1 to the particle array
!Afterwards, p_buffer_1 is free
    if (nzone_send_west /= 0) then
!Received particles and co. from east
        k = 0
        do k_loc = 1, inb_part_array
            k = k + 1
            do i = 1, nzone_send_west
                l_q(nzone_grid + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_west))
            end do
        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_send_west
                l_hq(nzone_grid + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_west))
            end do

        end do

        k = k + 1
        do i = 1, nzone_send_west
            l_tags(nzone_grid + i) = p_buffer_1(i + ((k - 1)*nzone_send_west))
        end do

    end if

!#################################################################
! Send to east and get from west
!#################################################################
! Increase the particle number by the particles recieved from west
    particle_number = particle_number + nzone_send_east

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(p_buffer_2, nzone_east*m, MPI_REAL8, dest_east, 0, MPI_COMM_WORLD, mpireq(l), ims_err) !Send p_buffer_2 to the east
    call MPI_IRECV(p_buffer_1, nzone_send_east*m, MPI_REAL8, source_east, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err) !Get particles from the west into p_buffer_1

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

!Write nzone_send_east recieved particles from west (step one) that are in p_buffer_1 to the particle array
    if (nzone_send_east /= 0) then
!Received particles and co. from west
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_send_east
                l_q(nzone_grid + nzone_send_west + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_east))
            end do

        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_send_east
                l_hq(nzone_grid + nzone_send_west + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_east))
            end do

        end do

        k = k + 1
        do i = 1, nzone_send_east
            l_tags(nzone_grid + nzone_send_west + i) = p_buffer_1(i + ((k - 1)*nzone_send_east))

        end do

    end if

    return
end subroutine PARTICLE_MPI_SEND_RECV_I

!#######################################################################
!#######################################################################
subroutine PARTICLE_MPI_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, l_q, l_hq, l_tags, particle_number)
    use TLab_Constants, only: wp, wi
    use PARTICLE_VARS, only: isize_part, inb_part_array, inb_part
    use PARTICLE_ARRAYS, only: p_buffer_1, p_buffer_2
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_k, ims_npro_k, ims_pro, ims_npro, ims_err

    implicit none

    integer(wi) nzone_grid, nzone_south, nzone_north

    real(wp), dimension(isize_part, *) :: l_q
    real(wp), dimension(isize_part, *) :: l_hq
    real(wp), dimension(isize_part) :: l_tags !Attention. Chosen real(wp) on purpose.
    integer(wi) particle_number

! -------------------------------------------------------------------
    integer(wi) nzone_send_south, nzone_send_north
    integer(wi) source_south, source_north
    integer(wi) dest_south, dest_north
    integer(wi) l
    type(MPI_Request) mpireq(ims_npro*4)
    type(MPI_Status) status(ims_npro*4)
    integer(wi) i, k, k_loc, m

!#######################################################################
    m = (inb_part_array*2) + 1 !Sending size of the buffer_parts

    particle_number = nzone_grid

!###################################################################
!Setting up destination and source + send/recv number of particles which will be send
!###################################################################
    if (ims_pro_k == 0) then ! Fisrt row
        dest_south = ims_pro_i + ims_npro_i*(ims_npro_k - 1)!Destination of the message
        source_south = ims_pro_i + ims_npro_i !Source of the message

    else if (ims_pro_k == (ims_npro_k - 1)) then !Last row
        dest_south = ims_pro_i + ims_npro_i*(ims_npro_k - 2)
        source_south = ims_pro_i

    else !Any case
        dest_south = ims_pro_i + ims_npro_i*(ims_pro_k - 1)  !Destination of the message
        source_south = ims_pro_i + ims_npro_i*(ims_pro_k + 1) !Source of the message
    end if

    mpireq(1:ims_npro*4) = MPI_REQUEST_NULL
    l = 4*ims_pro + 1
!Amount of particles which will be sent/received
    call MPI_ISEND(nzone_south, 1, MPI_INTEGER4, dest_south, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
    call MPI_IRECV(nzone_send_south, 1, MPI_INTEGER4, source_south, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)

!    CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  ! ALBERTO CONSIDER DELETING THIS ONE

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

!    mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
!    l = 2*ims_pro +1
!Amount of particles which will be sent/received
    call MPI_ISEND(nzone_north, 1, MPI_INTEGER4, dest_north, 0, MPI_COMM_WORLD, mpireq(l + 2), ims_err)
    call MPI_IRECV(nzone_send_north, 1, MPI_INTEGER4, source_north, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 3), ims_err)

    call MPI_Waitall(ims_npro*4, mpireq, status, ims_err)

!#################################################################
! Send to west and get from east
!#################################################################
    particle_number = particle_number + nzone_send_south

    if (nzone_south /= 0) then

!Send particles and co. to south
!Fill p_buffer_2
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_south
                p_buffer_2(i + ((k - 1)*nzone_south)) = l_q(nzone_grid + i, k_loc)
            end do

        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_south
                p_buffer_2(i + ((k - 1)*nzone_south)) = l_hq(nzone_grid + i, k_loc)
            end do

        end do

        k = k + 1
        do i = 1, nzone_south
            p_buffer_2(i + ((k - 1)*nzone_south)) = l_tags(nzone_grid + i)
        end do

    end if

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(p_buffer_2, nzone_south*m, MPI_REAL8, dest_south, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
    call MPI_IRECV(p_buffer_1, nzone_send_south*m, MPI_REAL8, source_south, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

    if (nzone_north /= 0) then

!Send particles and co. to north
!Fill p_buffer_2 again
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_north
                p_buffer_2(i + ((k - 1)*nzone_north)) = l_q(nzone_grid + nzone_south + i, k_loc)
            end do

        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_north
                p_buffer_2(i + ((k - 1)*nzone_north)) = l_hq(nzone_grid + nzone_south + i, k_loc)
            end do

        end do

        k = k + 1
        do i = 1, nzone_north
            p_buffer_2(i + ((k - 1)*nzone_north)) = l_tags(nzone_grid + nzone_south + i)
        end do

    end if

    if (nzone_send_south /= 0) then

!Received particles and co. from north
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_send_south
                l_q(nzone_grid + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_south))
            end do
        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_send_south
                l_hq(nzone_grid + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_south))
            end do

        end do

        k = k + 1
        do i = 1, nzone_send_south
            l_tags(nzone_grid + i) = p_buffer_1(i + ((k - 1)*nzone_send_south))
        end do

    end if

!#################################################################
! Send to east and get from west
!#################################################################
    particle_number = particle_number + nzone_send_north

    mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
    l = 2*ims_pro + 1
!Actual sending/receiving of particles
    call MPI_ISEND(p_buffer_2, nzone_north*m, MPI_REAL8, dest_north, 0, MPI_COMM_WORLD, mpireq(l), ims_err)
    call MPI_IRECV(p_buffer_1, nzone_send_north*m, MPI_REAL8, source_north, MPI_ANY_TAG, MPI_COMM_WORLD, mpireq(l + 1), ims_err)

    call MPI_Waitall(ims_npro*2, mpireq, status, ims_err)

    if (nzone_send_north /= 0) then
!Received particles and co. from south
        k = 0
        do k_loc = 1, inb_part_array

            k = k + 1
            do i = 1, nzone_send_north
                l_q(nzone_grid + nzone_send_south + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_north))
            end do

        end do

        do k_loc = 1, inb_part

            k = k + 1
            do i = 1, nzone_send_north
                l_hq(nzone_grid + nzone_send_south + i, k_loc) = p_buffer_1(i + ((k - 1)*nzone_send_north))
            end do

        end do

        k = k + 1
        do i = 1, nzone_send_north
            l_tags(nzone_grid + nzone_send_south + i) = p_buffer_1(i + ((k - 1)*nzone_send_north))

        end do

    end if

    return
end subroutine PARTICLE_MPI_SEND_RECV_K

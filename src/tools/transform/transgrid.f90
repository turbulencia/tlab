program TRANSGRID
    use TLab_Constants, only: wp, wi
    use TLab_Types, only: grid_dt
    implicit none

    type(grid_dt), dimension(3) :: g, g_ref

    integer(wi) option, direction, n, isize_wrk1d
    character*32 ifile, ffile, sfile, file_ref
    logical flag_exit

    real(wp), dimension(:, :), allocatable :: wrk1d
    real(wp) offset, factor1, factor2, dummy

    ! ###################################################################
    ! Initialize and read reference data
    ! ###################################################################
    g(1)%name = 'x'
    g(2)%name = 'y'
    g(3)%name = 'z'

    write (*, '("- Reference grid file ? ", $)')
    read (*, *) ifile

    ffile = trim(adjustl(ifile))//'.trn'
    sfile = trim(adjustl(ffile))//'.sts'

    call RD_GRIDSIZE(ifile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale)

    isize_wrk1d = max(g(1)%size, max(g(2)%size, g(3)%size))

    allocate (g(1)%nodes(2*g(1)%size)) ! Allocation of memory is doubled to allow introduction of planes
    allocate (g(2)%nodes(2*g(2)%size))
    allocate (g(3)%nodes(2*g(3)%size))
    allocate (wrk1d(isize_wrk1d, 3))

  CALL IO_READ_GRID(ifile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes, dummy)

    ! #######################################################################
    ! Main loop
    ! #######################################################################
    flag_exit = .false.

    do while (.not. flag_exit)

        write (*, '(20(/),"TransGrid Main Menu",/)')
        write (*, '(a)') '0. Dump ASCII values to file'
        write (*, '(a)') '1. Offsets'
        write (*, '(a)') '2. Scaling'
        write (*, '(a)') '3. Drop planes'
        write (*, '(a)') '4. Introduce planes'
        write (*, '(a)') '5. Transfer grid between files'
        write (*, '(a)') '6. Stretching'
        write (*, '(a,/)') '9. Exit'
        read (*, *) option

        if (option < 9) then
            write (*, '(/,"Direction (1 for x, 2 for y, 3 for z) ? ", $)')
            read (*, *) direction
        end if

        select case (option)

        case (0)
            open (23, file=trim(adjustl(g(direction)%name))//'.dat')
            do n = 1, g(direction)%size
                write (23, *) g(direction)%nodes(n)
            end do
            close (23)

        case (1)
            write (*, '(/,"Offset ? ", $)')
            read (*, *) offset
            g(direction)%nodes(:) = g(direction)%nodes(:) + offset

        case (2)
            write (*, '(/,"Scaling factor ? ", $)')
            read (*, *) factor1
            g(direction)%nodes(:) = g(direction)%nodes(1) + (g(direction)%nodes(:) - g(direction)%nodes(1))*factor1
            g(direction)%scale = g(direction)%scale*factor1

        case (3) ! Dropping planes
            call TRANS_DROP_PLANES(g(direction)%size, g(direction)%nodes, g(direction)%scale)

        case (4) ! Introducing planes
            call TRANS_ADD_PLANES(g(direction)%size, g(direction)%nodes, g(direction)%scale)
            flag_exit = .true. ! Exit directly to avoid memory allocation problems

        case (5) ! Transfering
            write (*, '(/,"Reference file to copy the data from ? ", $)')
            read (*, *) file_ref

            call RD_GRIDSIZE(file_ref, g_ref(1)%size, g_ref(2)%size, g_ref(3)%size, g_ref(1)%scale, g_ref(2)%scale, g_ref(3)%scale)
            do n = 1, 3
                if (g_ref(n)%size > 2*g(n)%size) then
                    write (*, *) 'Error. Reference grid too big.'
                    stop
                else
                    allocate (g_ref(n)%nodes(g_ref(n)%size))
                end if
            end do
        CALL IO_READ_GRID(ifile, g_ref(1)%size,g_ref(2)%size,g_ref(3)%size, g_ref(1)%scale,g_ref(2)%scale,g_ref(3)%scale, g_ref(1)%nodes,g_ref(2)%nodes,g_ref(3)%nodes, dummy)

            g(direction)%size = g_ref(direction)%size
            g(direction)%scale = g_ref(direction)%scale
            g(direction)%nodes(1:g(direction)%size) = g_ref(direction)%nodes(1:g(direction)%size)

            flag_exit = .true. ! Exit directly to avoid memory allocation problems

        case (6) ! Stretching
            write (*, '(/,"Stretching parameters? ", $)')
            read (*, *) factor1, factor2

            g(direction)%nodes = g(direction)%nodes*(1.0_wp + factor1*exp(-factor2*g(direction)%nodes))

        case (9)
            flag_exit = .true.

        end select

    end do

    ! #######################################################################
    ! Statistics and writing the data
    ! #######################################################################
    if (option > 0) then
        do direction = 1, 3
            call TRANS_DATA(sfile, g(direction), wrk1d(:, 1), wrk1d(:, 2))
        end do
  call IO_WRITE_GRID(ffile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, g(1)%nodes, g(2)%nodes, g(3)%nodes)
    end if

    stop

contains

    ! ###################################################################
    ! ###################################################################
    subroutine RD_GRIDSIZE(name, imax, jmax, kmax, scalex, scaley, scalez)
        implicit none

        character*(*) name
        integer(wi) imax, jmax, kmax
        real(wp) scalex, scaley, scalez

        real(wp) scale(3)

        open (50, file=name, status='old', form='unformatted')
        rewind (50)
        read (50) imax, jmax, kmax
        read (50) scale
        scalex = scale(1)
        scaley = scale(2)
        scalez = scale(3)
        close (50)

        return
    end subroutine RD_GRIDSIZE

    ! ###################################################################
    ! ###################################################################
    subroutine TRANS_DROP_PLANES(nmax, a, scale)
        implicit none

        integer(wi) nmax
        real(wp) a(nmax), scale

        integer(wi) nplanes, option, n
        real(wp) correction, tolerance

        tolerance = 1.0e-10_wp

        write (*, '(20(/),"Drop planes",/)')
        write (*, '(a)') '1. Drop planes symmetrically'
        write (*, '(a)') '2. Drop planes at the beginning'
        write (*, '(a)') '3. Drop planes at the end'
        write (*, '(a)') '4. Drop intermediate planes'
        write (*, '(a,/)') '9. Exit'
        read (*, *) option

        if (option < 4) then
            write (*, '(/,"Total number of planes to drop ? ", $)')
            read (*, *) nplanes
            if (nplanes >= nmax) then
                write (*, *) 'Error. Trying to drop equal/more planes than exist.'
                stop
            end if
        end if

        select case (option)

        case (1) ! Drop planes symmetrically
            nplanes = nplanes/2
            correction = scale - (a(nmax) - a(1)) ! The variable correction takes care of the periodic case
            scale = a(nmax - nplanes) - a(nplanes + 1) + correction
            nmax = nmax - 2*nplanes
            do n = 1, nmax
                a(n) = a(n + nplanes)
            end do

        case (2) ! Drop planes at the beginning
            correction = scale - (a(nmax) - a(1))
            scale = a(nmax) - a(nplanes + 1) + correction
            do n = nplanes + 1, nmax
                a(n - nplanes) = a(n)
            end do
            nmax = nmax - nplanes
            if (nmax == 1) scale = 1.0_wp

        case (3) ! Drop planes at the end
            correction = scale - (a(nmax) - a(1))
            scale = a(nmax - nplanes) - a(1) + correction
            nmax = nmax - nplanes
            if (nmax == 1) scale = 1.0_wp

        case (4) ! Drop intermediate planes
            correction = scale - (a(nmax) - a(1))
            do n = 1, nmax/2
                a(n) = a(2*n - 1)
            end do
            nmax = nmax/2
            if (correction < tolerance) then
                do n = 1, nmax
                    a(n) = a(1) + (a(n) - a(1))/(a(nmax) - a(1))*scale
                end do
            end if

        end select

        return
    end subroutine TRANS_DROP_PLANES

    ! ###################################################################
    ! ###################################################################
    subroutine TRANS_ADD_PLANES(nmax, a, scale)
        implicit none

        integer(wi) nmax
        real(wp) a(2*nmax), scale

        integer(wi) nplanes, option, n, inipos
        real(wp) correction, tolerance, deltaup, deltadown, posnew

        tolerance = 1.0e-10_wp

        write (*, '(20(/),"Introduce planes and exit",/)')
        write (*, '(a)') '1. Introduce planes symmetrically'
        write (*, '(a)') '2. Introduce planes at the beginning'
        write (*, '(a)') '3. Introduce planes at the end'
        write (*, '(a)') '4. Introduce planes at midpositions'
        write (*, '(a)') '5. Introduce zone of planes'
        write (*, '(a)') '6. Introduce particular planes by stdin'
        write (*, '(a,/)') '9. Exit'
        read (*, *) option

        if (option < 4 .or. option == 5) then
            write (*, '(/,"Total number of planes to introduce? ", $)')
            read (*, *) nplanes
            if (nplanes > nmax) then
                write (*, *) 'Error: Not enough space in grid arrays.'
                stop
            end if
        end if

        select case (option)

        case (1) ! Introduce planes symmetrically
            nplanes = nplanes/2
            write (*, '(/,"Initial position ? ", $)')
            read (*, *) deltadown
            write (*, '(/,"Final position ? ", $)')
            read (*, *) deltaup

            correction = scale - a(nmax) ! The variable correction takes care of the periodic case
            if (deltaup > scale) then; deltaup = (deltaup - scale)/real(nplanes, wp)
            else; deltaup = a(nmax) - a(nmax - 1)
            end if
            if (deltadown > a(1)) then; deltadown = a(2) - a(1)
            else; deltadown = (a(1) - deltadown)/real(nplanes, wp)
            end if
            do n = nmax, 1, -1
                a(n + nplanes) = a(n)
            end do
            do n = 1, nplanes
                a(nmax + nplanes + n) = a(nmax + nplanes + n - 1) + deltaup
                a(nplanes + 1 - n) = a(nplanes + 1 - n + 1) - deltadown
            end do
            nmax = nmax + 2*nplanes
            scale = a(nmax) - a(1) + correction

        case (2) ! Introduce planes at the beginning
            write (*, '(/,"Initial position ? ", $)')
            read (*, *) deltadown

            correction = scale - a(nmax)
            if (deltadown > a(1)) then; deltadown = a(2) - a(1)
            else; deltadown = (a(1) - deltadown)/real(nplanes, wp)
            end if
            do n = nmax, 1, -1
                a(n + nplanes) = a(n)
            end do
            do n = 1, nplanes
                a(nplanes + 1 - n) = a(nplanes + 1 - n + 1) - deltadown
            end do
            nmax = nmax + nplanes
            scale = a(nmax) - a(1) + correction

        case (3) ! Introduce planes at the end
            write (*, '(/,"Final position ? ", $)')
            read (*, *) deltaup

            correction = scale - a(nmax)
            if (deltaup > scale) then; deltaup = (deltaup - scale)/real(nplanes, wp)
            else; deltaup = a(nmax) - a(nmax - 1)
            end if
            do n = 1, nplanes
                a(nmax + n) = a(nmax + n - 1) + deltaup
            end do
            nmax = nmax + nplanes
            scale = a(nmax) - a(1) + correction

        case (4) ! Introduce planes at midpositions
            correction = scale - a(nmax)
            do n = nmax, 1, -1
                a(2*n - 1) = a(n)
            end do
            do n = 2, 2*nmax - 2, 2
                a(n) = (a(n + 1) + a(n - 1))*0.5_wp
            end do
            a(2*nmax) = a(2*nmax - 1) + (a(2*nmax - 1) - a(2*nmax - 2))
            nmax = 2*nmax
            if (correction < tolerance) then
                do n = 1, nmax
                    a(n) = a(1) + (a(n) - a(1))/(a(nmax) - a(1))*scale
                end do
            end if

        case (5) ! Introduce zone of planes
            write (*, '(/,"Initial position ? ", $)')
            read (*, *) inipos

            correction = scale - a(nmax)
            deltadown = a(inipos) - a(inipos - 1)
            do n = nmax, inipos, -1
                a(n + nplanes) = a(n) + deltadown*real(nplanes, wp)
            end do
            do n = 1, nplanes
                a(inipos + n - 1) = a(inipos + n - 2) + deltadown
            end do
            nmax = nmax + nplanes
            scale = a(nmax) - a(1) + correction

        case (6) ! Introduce particular planes by stdin
            write (*, '(/,"Plane position ? ", $)') ! to be written in terms of a list_real
            read (*, *) posnew

            correction = scale - a(nmax)
            if (correction > tolerance) then
                write (*, *) "ERROR 2: Periodic direction. Use other option."
                stop
            end if
            do n = nmax, 1, -1
                if (a(n) > posnew) then
                    a(n + 1) = a(n)
                else
                    a(n + 1) = posnew
                    exit
                end if
            end do
            if (a(1) > posnew) then
                a(1) = posnew
            end if
            nmax = nmax + 1
            scale = a(nmax) - a(1)

        end select

        return

    end subroutine TRANS_ADD_PLANES

    ! ###################################################################
    ! ###################################################################
    subroutine TRANS_DATA(name, grid, work1, work2)
        implicit none

        type(grid_dt) grid

        character*(*) name
        real(wp) work1(*), work2(*)

        real(wp) dxmx, dxmn, axmx, axmn
        integer(wi) n

        open (20, file=name, status='unknown', position='append')

        write (20, '(a13)') '['//trim(adjustl(grid%name))//'-direction]'

        if (grid%size > 1) then

            if (mod(grid%size, 2) /= 0) then
                write (20, *) 'Warning. Not an even number of points.'
            end if

            work1(2) = grid%nodes(2) - grid%nodes(1)
            do n = 3, grid%size
                work1(n) = grid%nodes(n) - grid%nodes(n - 1)
                work2(n) = work1(n)/work1(n - 1)
            end do
            dxmx = maxval(work1(2:grid%size)); dxmn = minval(work1(2:grid%size))
            axmx = maxval(work2(3:grid%size)); axmn = minval(work2(3:grid%size))

            write (20, '(a25,i5)') 'number of points .......: ', grid%size
            write (20, '(a25,e12.5)') 'origin .................: ', grid%nodes(1)
            write (20, '(a25,e12.5)') 'end point ..............: ', grid%nodes(grid%size)
            write (20, '(a25,e12.5)') 'scale ..................: ', grid%scale
            write (20, '(a25,e12.5)') 'minimum step ...........: ', dxmn
            write (20, '(a25,e12.5)') 'maximum step ...........: ', dxmx
            write (20, '(a25,e12.5)') 'minimum stretching .....: ', axmn
            write (20, '(a25,e12.5)') 'maximum stretching .....: ', axmx

        else
            write (20, '(a7)') '2D case'

        end if

        close (20)

        return
    end subroutine TRANS_DATA

end program TRANSGRID

!# Grid generation tool. Origin is set always to (0,0,0)
#include "dns_error.h"

#define C_FILE_LOC "INIGRID"

program INIGRID
    use FDM, only: fdm_dt
    use TLab_Constants, only: wp, gfile, ifile, lfile, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Arrays, only: wrk1d, wrk2d
    use TLab_Grid
    use GRID_LOCAL
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif
    implicit none

    character*32 sfile, bakfile
    type(fdm_dt) :: g(3)
    ! real(wp), allocatable :: wrk1d(:, :)
    integer(wi) idir, iseg, isize_wrk1d, n, nmax, iloc
    real(wp) scale_old, scale_new, ds
    character(len=16), parameter :: block(3) = ['IniGridOx', 'IniGridOy', 'IniGridOz']

    ! #######################################################################
    ! Initialize
    ! #######################################################################
    sfile = trim(adjustl(gfile))//'.sts'
    bakfile = trim(adjustl(ifile))//'.bak'

    g(1)%name = 'x'
    g(2)%name = 'y'
    g(3)%name = 'z'

    call TLab_Start()

    do idir = 1, 3
        call GRID_READBLOCK(bakfile, ifile, block(idir), g_build(idir), g(idir)%periodic)

        g(idir)%size = g_build(idir)%size(1)                    ! Calculate total number of points
        do iseg = 2, g_build(idir)%nseg
            g(idir)%size = g(idir)%size + g_build(idir)%size(iseg) - 1
        end do
        if (g_build(idir)%mirrored) g(idir)%size = 2*g(idir)%size - 2

        allocate (g(idir)%nodes(g(idir)%size))                   ! memory space for the grid nodes

    end do

    isize_wrk1d = max(g(1)%size, max(g(2)%size, g(3)%size))
    allocate (wrk1d(isize_wrk1d, 8))
    allocate (wrk2d(isize_wrk1d, 2))

    ! #######################################################################
    ! Construct grid
    ! #######################################################################
    do idir = 1, 3

        iloc = 1
        if (g_build(idir)%mirrored) iloc = g(idir)%size/2   ! mirrored case; first point in array is imax/2

        g(idir)%nodes(iloc) = 0.0_wp                        ! set origin at zero

        do iseg = 1, g_build(idir)%nseg                     ! Loop over the segments that build the grid
            nmax = g_build(idir)%size(iseg)                 ! for readability below
            ! create uniform reference grid s starting at zero
            if (nmax > 1) then
                ds = (g_build(idir)%end(iseg) - g(idir)%nodes(iloc))/real(nmax - 1, wp)
                g(idir)%nodes(iloc:) = [(real(n - 1, wp), n=1, nmax)]*ds + g(idir)%nodes(iloc)
                ! do n = 1, nmax - 1
                !     g(idir)%nodes(iloc + n) = g(idir)%nodes(iloc + n - 1) + ds
                ! end do

                select case (g_build(idir)%opts(1, iseg))

                case (GTYPE_UNIFORM)
                    ! already done

                case (GTYPE_TANH)
                    call BLD_TANH(idir, iseg, g(idir)%nodes(iloc:), nmax, wrk1d)

                case (GTYPE_EXP)
                    call BLD_EXP(idir, iseg, g(idir)%nodes(iloc:), nmax, wrk1d)

                case DEFAULT
                    call BLD_THEREST(idir, iseg, g(idir)%nodes(iloc:), nmax)

                end select

                iloc = iloc + nmax - 1

            end if
        end do

        if (g_build(idir)%mirrored) call GRID_MIRROR(g(idir)%size, g(idir)%nodes)

        if (g(idir)%size > 1) then
            g(idir)%scale = g(idir)%nodes(g(idir)%size) - g(idir)%nodes(1)
        else
            g(idir)%scale = 1.0_wp
        end if

        if (g_build(idir)%fixed_scale > 0.0_wp) then                ! rescale on exact fixed value
            scale_new = g_build(idir)%fixed_scale
            scale_old = g(idir)%scale
            g(idir)%nodes = (g(idir)%nodes/scale_old)*scale_new     ! rescale nodes
            g(idir)%nodes(g(idir)%size) = scale_new                 ! avoid rounding error
            g(idir)%scale = g(idir)%nodes(g(idir)%size)             ! update scale
        end if

        if (g(idir)%periodic) g(idir)%size = g(idir)%size - 1

    end do

    ! #######################################################################
    ! Statistics
    ! #######################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        open (20, file=sfile)

        do idir = 1, 3
            write (20, 3000) '['//trim(adjustl(g(idir)%name))//'-direction]'

            if (g(idir)%size > 1) then
                wrk1d(2, 1) = g(idir)%nodes(2) - g(idir)%nodes(1)
                do n = 3, g(idir)%size
                    wrk1d(n, 1) = g(idir)%nodes(n) - g(idir)%nodes(n - 1)
                    wrk1d(n, 2) = wrk1d(n, 1)/wrk1d(n - 1, 1)
                end do

                write (20, 2000) 'number of points .......: ', g(idir)%size
                write (20, 1000) 'origin .................: ', g(idir)%nodes(1)
                write (20, 1000) 'end point ..............: ', g(idir)%nodes(g(idir)%size)
                write (20, 1000) 'scale ..................: ', g(idir)%scale
                write (20, 1000) 'minimum step ...........: ', minval(wrk1d(2:g(idir)%size, 1))
                write (20, 1000) 'maximum step ...........: ', maxval(wrk1d(2:g(idir)%size, 1))
                write (20, 1000) 'minimum stretching .....: ', minval(wrk1d(3:g(idir)%size, 2))
                write (20, 1000) 'maximum stretching .....: ', maxval(wrk1d(3:g(idir)%size, 2))

            else
                write (20, '(a7)') '2D case'

            end if

        end do

        close (20)

        ! #######################################################################
        ! Writing data
        ! #######################################################################
        call TLab_Write_ASCII(lfile, 'Writing grid.')
        call TLab_Grid_Write(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, g(1)%nodes, g(2)%nodes, g(3)%nodes)

#ifdef USE_MPI
    end if
#endif

    call TLab_Stop(0)

1000 format(a25, e12.5)
2000 format(a25, i5)
3000 format(a13)

contains
    ! #######################################################################
    ! #######################################################################
    subroutine GRID_READBLOCK(bakfile, inifile, block, var, periodic)

        character(LEN=*), intent(in) :: bakfile, inifile, block
        type(grid_build_dt), intent(out) :: var
        logical, intent(OUT) :: periodic

! -------------------------------------------------------------------
        integer(wi) idummy
        character(LEN=512) sRes, str

! #######################################################################
        call TLab_Write_ASCII(bakfile, '['//block//']')
        call TLab_Write_ASCII(bakfile, 'segments=<number of segments>')
        call TLab_Write_ASCII(bakfile, 'periodic=<yes/no>')
        call TLab_Write_ASCII(bakfile, 'mirrored=<yes/no>')
        call TLab_Write_ASCII(bakfile, 'fixed_scale=<value>')

        call ScanFile_Int(bakfile, inifile, block, 'segments', '1', var%nseg)

        periodic = .false.
        call ScanFile_Char(bakfile, inifile, block, 'periodic', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') periodic = .true.

        var%mirrored = .false.
        call ScanFile_Char(bakfile, inifile, block, 'mirrored', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') var%mirrored = .true.

        call ScanFile_Real(bakfile, inifile, block, 'fixed_scale', '-1.0', var%fixed_scale)

        if (periodic .and. var%mirrored) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Periodicity with mirroring is not supported.')
            call TLab_Stop(DNS_ERROR_GRID_SCALE)
        end if

! -------------------------------------------------------------------
! Loop over the segments
! -------------------------------------------------------------------
        do iseg = 1, var%nseg
            write (str, *) iseg

            call TLab_Write_ASCII(bakfile, 'Segment number '//trim(adjustl(str)))
            call TLab_Write_ASCII(bakfile, 'scales_'//trim(adjustl(str))//'=<physical end of the segment>')
            call TLab_Write_ASCII(bakfile, 'points_'//trim(adjustl(str))//'=<points in the segment>')
            call TLab_Write_ASCII(bakfile, 'opts_'//trim(adjustl(str))//'=<option>')
            call TLab_Write_ASCII(bakfile, 'vals_'//trim(adjustl(str))//'=<values>')

            call ScanFile_Int(bakfile, inifile, block, 'points_'//trim(adjustl(str)), '1', var%size(iseg))
            call ScanFile_Real(bakfile, inifile, block, 'scales_'//trim(adjustl(str)), '-1.0', var%end(iseg))

            var%opts(:, iseg) = 0
            call ScanFile_Char(bakfile, inifile, block, 'opts_'//trim(adjustl(str)), '1', sRes)
            if (trim(adjustl(sRes)) == 'uniform') then; var%opts(1, iseg) = GTYPE_UNIFORM
            else if (trim(adjustl(sRes)) == 'tanh') then; var%opts(1, iseg) = GTYPE_TANH
            else if (trim(adjustl(sRes)) == 'exp') then; var%opts(1, iseg) = GTYPE_EXP
            else
                idummy = MAX_PARAMES
                call LIST_INTEGER(sRes, idummy, var%opts(1, iseg))
            end if

            var%vals(:, iseg) = 0
            call ScanFile_Char(bakfile, inifile, block, 'vals_'//trim(adjustl(str)), '1.0', sRes)
            idummy = MAX_PARAMES
            call LIST_REAL(sRes, idummy, var%vals(1, iseg))

        end do

        return
    end subroutine GRID_READBLOCK

    ! #######################################################################
    ! #######################################################################
    subroutine GRID_MIRROR(imax, x)
        implicit none

        integer(wi), intent(IN) :: imax
        real(wp), intent(INOUT) :: x(imax)

        ! -----------------------------------------------------------------------
        integer(wi) i
        real(wp) offset

        ! #######################################################################
        ! Offset for even number of points
        offset = (x(imax/2 + 1) - x(imax/2))/2.0_wp
        do i = imax/2, imax
            x(i) = x(i) - offset
        end do

        ! Mirroring
        do i = 1, imax/2 - 1
            x(i) = -x(imax + 1 - i)
        end do

        ! Global translation to set origin at zero
        offset = x(1)
        x = x - offset

        return
    end subroutine GRID_MIRROR

end program INIGRID

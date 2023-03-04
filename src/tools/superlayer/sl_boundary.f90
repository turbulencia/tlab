#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
program SL_BOUNDARY

    use TLAB_VARS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN

    implicit none

! -------------------------------------------------------------------
! Grid and associated arrays
    TREAL, dimension(:, :), allocatable, save, target :: x, y, z

! Flow variables
    TREAL, dimension(:, :), allocatable, target :: q
    TREAL, dimension(:), allocatable :: s, field

! Work arrays
    TREAL, dimension(:), allocatable :: wrk1d, wrk2d, wrk3d

! Surface arrays
    TREAL, dimension(:, :), allocatable :: sl, samples

    TREAL txc(:)
    allocatable txc

    TREAL pdf(:)
    allocatable pdf

! Pointers to existing allocated space
    TREAL, dimension(:), pointer :: u, v, w

    TINTEGER iopt, iint, isl, ith, itxc_size, nfield, np
    TINTEGER iread_flow, iread_scal, jmin_loc, jmax_loc, idummy
    TREAL threshold, vmin, vmax
    TINTEGER buff_nps_u_jmin, buff_nps_u_jmax
    character*64 str
    character*32 fname, bakfile

    TINTEGER itime_size_max, itime_size, i
    parameter(itime_size_max=128)
    TINTEGER itime_vec(itime_size_max)
    TINTEGER iopt_size_max, iopt_size
    parameter(iopt_size_max=10)
    TREAL opt_vec(iopt_size_max)
    character*512 sRes
#ifdef USE_MPI
    integer icount
#endif

    TREAL, dimension(:, :), pointer :: dx, dy, dz

! ###################################################################
    bakfile = trim(adjustl(ifile))//'.bak'

    call DNS_START

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    call SCANINIINT(bakfile, ifile, 'BufferZone', 'PointsUJmin', '0', buff_nps_u_jmin)
    call SCANINIINT(bakfile, ifile, 'BufferZone', 'PointsUJmax', '0', buff_nps_u_jmax)

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

    allocate (sl(imax*kmax, 6))
    allocate (wrk1d(isize_wrk1d*5))
    allocate (wrk2d(isize_wrk2d*10))
    allocate (wrk3d(isize_field))

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        call SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'Files', '-1', sRes)
        if (sRes == '-1') then
            write (*, *) 'Integral Iterations ?'
            read (*, '(A512)') sRes
        end if
        itime_size = itime_size_max
        call LIST_INTEGER(sRes, itime_size, itime_vec)
#ifdef USE_MPI
    end if
    call MPI_BCAST(itime_size, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    icount = itime_size
    call MPI_BCAST(itime_vec, icount, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        call SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'ParamSuperlayer', '-1', sRes)
        iopt_size = iopt_size_max
        call LIST_REAL(sRes, iopt_size, opt_vec)

        if (sRes == '-1') then
            write (*, *) 'Option ?'
            write (*, *) '1. Extract scalar gradient envelope surfaces to file'
            write (*, *) '2. PDFs conditioned on the envelope surface'
            write (*, *) '3. RQ JPDF conditioned on the envelope surface'
            write (*, *) '4. WS JPDF conditioned on the envelope surface'
            read (*, *) iopt

            write (*, *) 'Intermittency conditioning ?'
            write (*, *) ' 1. Based on scalar'
            write (*, *) ' 2. Based on vorticity'
            write (*, *) ' 3. Based on scalar gradient'
            read (*, *) iint

            write (*, *) 'Threshold based on relative (1) or absolute (2) values?'
            read (*, *) ith
            write (*, *) 'Threshold value ?'
            read (*, *) threshold
            if (iopt > 1) then
                write (*, *) 'Upper (1), lower (2) or both (3) envelope surfaces ?'
                read (*, *) isl
                write (*, *) 'Number of PDF bins ?'
                read (*, *) np
            end if
        else
            iopt = int(opt_vec(1))
            iint = int(opt_vec(2))
            ith = int(opt_vec(3))
            threshold = opt_vec(4)
            isl = int(opt_vec(5))
            np = int(opt_vec(6))
        end if

#ifdef USE_MPI
    end if

    call MPI_BCAST(iopt, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(iint, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(isl, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(np, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(ith, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(threshold, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
    if (iopt == 1) then; itxc_size = isize_field*2; nfield = 1; 
    else if (iopt == 2) then; itxc_size = isize_field*6; nfield = 5; iread_flow = 1; iread_scal = 1
    else if (iopt >= 3) then; itxc_size = isize_field*6; nfield = 4; iread_flow = 1; iread_scal = 0
    end if

    if (iint == 1) then; iread_scal = 1
    else if (iint == 2) then; iread_flow = 1
    else if (iint == 3) then; iread_scal = 1
    end if

    if (iread_flow == 1) allocate (q(isize_field, 3))
    if (iread_scal == 1) allocate (s(isize_field))

    allocate (pdf(nfield*np))
    allocate (samples(imax*kmax, nfield*2))
    allocate (txc(itxc_size))

    if (iopt <= 2) allocate (field(isize_field))

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

! ###################################################################
! Define pointers
! ###################################################################
    dx => x(:, 2:) ! to be removed
    dy => y(:, 2:)
    dz => z(:, 2:)

    u => q(:, 1)
    v => q(:, 2)
    w => q(:, 3)

! ###################################################################
! Postprocess given list of files
! ###################################################################
    do i = 1, itime_size

        itime = itime_vec(i)

        if (iread_flow == 1) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 3, 0, q)
        end if

        if (iread_scal == 1) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, inb_scal, s)
        end if

! -------------------------------------------------------------------
! Extract scalar gradient surfaces to file
! -------------------------------------------------------------------
        if (iopt == 1) then
            jmin_loc = max(1, buff_nps_u_jmin)                 ! remove buffers
            jmax_loc = min(jmax, jmax - buff_nps_u_jmax + 1)

! Based on scalar
            if (iint == 1) then
                if (ith == 1) then ! relative to max
                    call MINMAX(imax, jmax, kmax, s, vmin, vmax)
                    vmin = threshold*threshold*vmax
                else if (ith == 2) then ! absolute
                    vmin = threshold
                end if
                call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, y, s, txc, sl(1, 1), wrk2d)
                call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, y, s, txc, sl(1, 2), wrk2d)

! Based on vorticity
            else if (iint == 2) then
                call TLAB_WRITE_ASCII(lfile, 'Calculating vorticity...')
                call FI_VORTICITY(imax, jmax, kmax, u, v, w, field, txc(1), txc(1 + isize_field))
                if (ith == 1) then ! relative to max
                    call MINMAX(imax, jmax, kmax, field, vmin, vmax)
                    vmin = threshold*threshold*vmax
                else if (ith == 2) then ! absolute
                    vmin = threshold
                end if

                call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, y, field, txc, sl(1, 1), wrk2d)
                call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, y, field, txc, sl(1, 2), wrk2d)

! Based on scalar gradient
            else if (iint == 3) then
                call TLAB_WRITE_ASCII(lfile, 'Calculating scalar gradient...')
                call FI_GRADIENT(imax, jmax, kmax, s, field, txc)
                call MINMAX(imax, jmax, kmax, field, vmin, vmax)
                write (str, '(E22.15E3,E22.15E3)') vmin, vmax; str = 'Bounds '//trim(adjustl(str))
                call TLAB_WRITE_ASCII(lfile, str)
                if (ith == 1) then ! relative to max
                    call MINMAX(imax, jmax, kmax, field, vmin, vmax)
                    vmin = threshold*threshold*vmax
                else if (ith == 2) then ! absolute
                    vmin = threshold
                end if

                call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, y, field, txc, sl(1, 1), wrk2d)
                call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, y, field, txc, sl(1, 2), wrk2d)

            end if

! write threshold
            write (str, '(E22.15E3)') vmin; str = 'Threshold '//trim(adjustl(str))
            call TLAB_WRITE_ASCII(lfile, str)

! save surfaces w/o header
            write (fname, *) itime; fname = 'sl'//trim(adjustl(fname))
            ! idummy = g(2)%size; g(2)%size = 1
            ! CALL DNS_WRITE_FIELDS(fname, i0, imax,i1,kmax, i2, isize_field, sl, wrk3d)
            ! g(2)%size = idummy
            call TLAB_WRITE_ASCII(efile, 'SL_BOUNDARY. To be written in terms of IO_SUBARRAY as in averages.x')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)

! save scalar dissipation as scalar field
            ! WRITE(fname,*) itime; fname = 'chi'//TRIM(ADJUSTL(fname))
            ! CALL IO_WRITE_FIELDS(fname, IO_SCAL, imax,jmax,kmax, i1, field)

! -------------------------------------------------------------------
! Surface PDFs
! -------------------------------------------------------------------
        else if (iopt == 2) then
            call SL_BOUNDARY_VORTICITY_PDF(isl, ith, np, nfield, itxc_size, threshold, buff_nps_u_jmax, &
                                           u, v, w, s, field, sl, samples, pdf, txc, wrk1d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! Surface JPDFs
! -------------------------------------------------------------------
        else if (iopt >= 3) then
            call SL_BOUNDARY_VORTICITY_JPDF(iopt, isl, ith, np, nfield, itxc_size, threshold, buff_nps_u_jmax, &
                                            u, v, w, sl, samples, txc, wrk1d, wrk2d, wrk3d)

        end if

    end do

    call TLAB_STOP(0)
end program SL_BOUNDARY

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
    use TLab_Constants, only: wp, wi

    use TLAB_VARS
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    use IO_Grid

    implicit none

! -------------------------------------------------------------------
! Grid and associated arrays
    real(wp), dimension(:, :), allocatable, save, target :: x, y, z

! Flow variables
    real(wp), dimension(:, :), allocatable, target :: q
    real(wp), dimension(:), allocatable :: s, field

! Work arrays
    real(wp), dimension(:), allocatable :: wrk1d, wrk2d, wrk3d

! Surface arrays
    real(wp), dimension(:, :), allocatable :: sl, samples

    real(wp) txc(:)
    allocatable txc

    real(wp) pdf(:)
    allocatable pdf

! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w

    integer(wi) iopt, iint, isl, ith, itxc_size, nfield, np
    integer(wi) jmin_loc, jmax_loc, idummy
    logical iread_flow, iread_scal
    real(wp) threshold, vmin, vmax
    integer(wi) buff_nps_u_jmin, buff_nps_u_jmax
    character*64 str
    character*32 fname, bakfile

    integer(wi) itime_size_max, itime_size, i
    parameter(itime_size_max=128)
    integer(wi) itime_vec(itime_size_max)
    integer(wi) iopt_size_max, iopt_size
    parameter(iopt_size_max=10)
    real(wp) opt_vec(iopt_size_max)
    character*512 sRes
#ifdef USE_MPI
    integer icount
#endif

    real(wp), dimension(:, :), pointer :: dx, dy, dz
    real(wp) params(1)
    
! ###################################################################
    bakfile = trim(adjustl(ifile))//'.bak'

    call DNS_START

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
call Rotation_Initialize(ifile)

    call ScanFile_Int(bakfile, ifile, 'BufferZone', 'PointsUJmin', '0', buff_nps_u_jmin)
    call ScanFile_Int(bakfile, ifile, 'BufferZone', 'PointsUJmax', '0', buff_nps_u_jmax)

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
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
        call ScanFile_Char(lfile, 'tlab.ini', 'PostProcessing', 'Files', '-1', sRes)
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
        call ScanFile_Char(lfile, 'tlab.ini', 'PostProcessing', 'ParamSuperlayer', '-1', sRes)
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
    else if (iopt == 2) then; itxc_size = isize_field*6; nfield = 5; iread_flow = .true.; iread_scal = .true.
    else if (iopt >= 3) then; itxc_size = isize_field*6; nfield = 4; iread_flow = .true.; iread_scal = .false.
    end if

    if (iint == 1) then; iread_scal = .true.
    else if (iint == 2) then; iread_flow = .true.
    else if (iint == 3) then; iread_scal = .true.
    end if

    if (iread_flow) allocate (q(isize_field, 3))
    if (iread_scal) allocate (s(isize_field))

    allocate (pdf(nfield*np))
    allocate (samples(imax*kmax, nfield*2))
    allocate (txc(itxc_size))

    if (iopt <= 2) allocate (field(isize_field))

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, x, y, z, [g(1)%size, g(2)%size, g(3)%size])
    call FDM_Initialize(g(1), x)
    call FDM_Initialize(g(2), y)
    call FDM_Initialize(g(3), z)

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

        if (iread_flow) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, 3, 0, q, params)
            rtime = params(1)
        end if

        if (iread_scal) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, inb_scal, s, params)
            rtime = params(1)
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
                call TLab_Write_ASCII(lfile, 'Calculating vorticity...')
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
                call TLab_Write_ASCII(lfile, 'Calculating scalar gradient...')
                call FI_GRADIENT(imax, jmax, kmax, s, field, txc)
                call MINMAX(imax, jmax, kmax, field, vmin, vmax)
                write (str, '(E22.15E3,E22.15E3)') vmin, vmax; str = 'Bounds '//trim(adjustl(str))
                call TLab_Write_ASCII(lfile, str)
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
            call TLab_Write_ASCII(lfile, str)

! save surfaces w/o header
            write (fname, *) itime; fname = 'sl'//trim(adjustl(fname))
            ! idummy = g(2)%size; g(2)%size = 1
            ! CALL DNS_WRITE_FIELDS(fname, i0, imax,i1,kmax, i2, isize_field, sl, wrk3d)
            ! g(2)%size = idummy
            call TLab_Write_ASCII(efile, 'SL_BOUNDARY. To be written in terms of IO_SUBARRAY as in averages.x')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)

! save scalar dissipation as scalar field
            ! WRITE(fname,*) itime; fname = 'chi'//TRIM(ADJUSTL(fname))
            ! CALL IO_Write_Fields(fname, IO_SCAL, imax,jmax,kmax, i1, field)

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

    call TLab_Stop(0)
end program SL_BOUNDARY

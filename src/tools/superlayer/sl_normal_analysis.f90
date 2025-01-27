#include "dns_const.h"
#include "dns_error.h"

program SL_NORMAL_ANALYSIS
    use TLab_Constants, only: wp, wi

    use TLAB_VARS
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use IO_FIELDS

    implicit none

! -------------------------------------------------------------------
! Grid and associated arrays
    real(wp), dimension(:, :), allocatable, save, target :: x, y, z

! Flow variables
    real(wp), dimension(:, :), pointer :: q

! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w, p

    real(wp) z1(:)
    allocatable z1
    real(wp) field(:)
    allocatable field
    real(wp) sl(:)
    allocatable sl
    real(wp) txc(:)
    allocatable txc
    real(wp) profiles(:)
    allocatable profiles
    real(wp) mean(:)
    allocatable mean
    real(wp) wrk1d(:)
    allocatable wrk1d
    real(wp) wrk2d(:)
    allocatable wrk2d
    real(wp) wrk3d(:)
    allocatable wrk3d

    integer(wi) iopt, isl, ith, itxc_size, iavg
    real(wp) threshold
    integer(wi) ibuffer_npy
    integer(wi) nmax, istep, kstep, nprof_size, nfield
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

    call ScanFile_Int(bakfile, ifile, 'BufferZone', 'NumPointsY', '0', ibuffer_npy)

    itxc_size = imax*jmax*kmax*7

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (u(imax*jmax*kmax))
    allocate (v(imax*jmax*kmax))
    allocate (w(imax*jmax*kmax))
    allocate (p(imax*jmax*kmax))
    allocate (z1(imax*jmax*kmax))
    allocate (field(imax*jmax*kmax))
    allocate (sl(imax*kmax))
    allocate (txc(itxc_size))
    allocate (wrk1d(isize_wrk1d*5))
    allocate (wrk2d(isize_wrk2d))
    allocate (wrk3d(isize_wrk3d))

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
        call ScanFile_Char(lfile, 'tlab.ini', 'PostProcessing', 'Superlayer', '-1', sRes)
        iopt_size = iopt_size_max
        call LIST_REAL(sRes, iopt_size, opt_vec)

        if (sRes == '-1') then
            write (*, *) 'Option ?'
            write (*, *) '1. Based on vorticity magnitude w_i w_i'
            write (*, *) '2. Based on scalar gradient G_i G_i'
            write (*, *) '3. Based on strain rate s_ij s_ij'
            read (*, *) iopt
            write (*, *) 'Upper (1) or lower (2) envelope surface ?'
            read (*, *) isl
            write (*, *) 'Number of normal points ?'
            read (*, *) nmax
            write (*, *) 'Averages (1) or instantaneous profiles (2) ?'
            read (*, *) iavg
            write (*, *) 'Step along OX ?'
            read (*, *) istep
            write (*, *) 'Step along OZ ?'
            read (*, *) kstep
            write (*, *) 'Threshold based on maximum (1) or mean (2) ?'
            read (*, *) ith
            write (*, *) 'Threshold value ?'
            read (*, *) threshold
        else
            iopt = int(opt_vec(1))
            isl = int(opt_vec(2))
            nmax = int(opt_vec(3))
            iavg = int(opt_vec(4))
            istep = int(opt_vec(5))
            kstep = int(opt_vec(6))
            ith = int(opt_vec(7))
            threshold = opt_vec(8)
        end if

#ifdef USE_MPI
    end if

    call MPI_BCAST(iopt, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(isl, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(nmax, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(iavg, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(istep, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(kstep, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(ith, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(threshold, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
    nfield = 13
    nprof_size = (imax/istep)*(kmax/kstep)*nmax*nfield

    allocate (profiles(nprof_size))
    allocate (mean(nmax*nfield*2))

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1), wrk1d(:, 4))
    call FDM_Initialize(y, g(2), wrk1d(:, 2), wrk1d(:, 4))
    call FDM_Initialize(z, g(3), wrk1d(:, 3), wrk1d(:, 4))

! ###################################################################
! Define pointers
! ###################################################################
    dx => x(:, 2:) ! to be removed
    dy => y(:, 2:)
    dz => z(:, 2:)

    u => q(:, 1)
    v => q(:, 2)
    w => q(:, 3)
    p => q(:, 4)

! ###################################################################
! Postprocess given list of files
! ###################################################################
    do i = 1, itime_size

        itime = itime_vec(i)

! -------------------------------------------------------------------
! Binary data
! -------------------------------------------------------------------
        write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 4, 0, q)

        write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, inb_scal, z1)

        call THERMO_CALORIC_TEMPERATURE(imax*jmax*kmax, z1, p, field, txc, wrk3d)
        call THERMO_THERMAL_PRESSURE(imax*jmax*kmax, z1, field, txc, p)

! -------------------------------------------------------------------
! Vorticity analysis
! -------------------------------------------------------------------
        if (iopt == 1) then
            call SL_NORMAL_VORTICITY(isl, ith, iavg, nmax, istep, kstep, nfield, itxc_size, &
                                     threshold, ibuffer_npy, u, v, w, p, z1, field, sl, profiles, txc, mean, wrk1d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! Scalar gradient analysis
! -------------------------------------------------------------------
        else if (iopt == 2) then
            call SL_NORMAL_GRADIENT(isl, nmax, istep, kstep, ibuffer_npy, &
                                    u, v, w, z1, field, sl, profiles, txc, wrk1d, wrk2d, wrk3d)
        end if

    end do

    call TLab_Stop(0)
end program SL_NORMAL_ANALYSIS

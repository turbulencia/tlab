#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI

#endif

#define SPEC_SINGLE  0
#define SPEC_AVERAGE 1

#define C_FILE_LOC "SPECTRA"

!########################################################################
!# Tool/Library SPECTRA
!#
!########################################################################
!# HISTORY
!#
!# 2012/04/26 - C. Ansorge
!#              Created
!# 2014/01/01 - J.P. Mellado
!#              Adding correlation, and cross terms
!# 2015/01/01 - J.P. Mellado
!#              Parallelizing the 1D spectra; radial spectre not yet
!# 2022/02/23 - J. Kostelecky
!#              adding IBM
!#
!########################################################################
!# DESCRIPTION
!#
!# Postprocessing tool to compute spectra in temporal mode
!#
!########################################################################
program SPECTRA

    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal, tag_part
    use TLab_Pointers, only: pointers_dt
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, isize_wrk3d, inb_flow, inb_scal, isize_txc_field, isize_field, isize_wrk2d, inb_wrk2d
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, fourier_on, flow_on, scal_on
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_pro, ims_npro_k
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose
#endif
    use FDM, only: g, FDM_Initialize
    use Thermodynamics, only: imixture, Thermodynamics_Initialize_Parameters
    use NavierStokes, only: nse_eqns, froude
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use Gravity, only: Gravity_Initialize, buoyancy, Gravity_Buoyancy
    use Rotation, only: Rotation_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use Thermo_Anelastic
    use Radiation
    use Microphysics
    use Chemistry
    use IBM_VARS
    use IO_FIELDS
    use OPR_FOURIER
    use OPR_FILTERS
    use OPR_Burgers, only: OPR_Burgers_Initialize
    use OPR_ELLIPTIC
    use Averages, only: AVG1V2D, COV2V2D
#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use SpectraMod

    implicit none

! Parameter definitions
    integer(wi), parameter :: itime_size_max = 512
    integer(wi), parameter :: iopt_size_max = 20

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), dimension(:), allocatable :: p_aux, y_aux, samplesize
    real(wp), dimension(:, :), allocatable :: out2d, outx, outz, outr

    type(pointers_dt), dimension(16) :: vars

    target p_aux

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
    character*32 fname, bakfile
    character*32 varname(16)
    character*64 str, line
    character*8 tag_file, tag_name, tag_var(16)
    integer(wi) p_pairs(16, 2)

    integer(wi) opt_main, opt_ffmt, opt_time, opt_block, flag_buoyancy
    integer(wi) flag_mode, ierr
    logical iread_flow, iread_scal
    integer(wi) isize_out2d, isize_aux, sizes(5)
    integer(wi) nfield, nfield_ref
    integer(wi) is, iv, iv_offset, iv1, iv2, ip, j, ig
    integer(wi) jmax_aux, kxmax, kymax, kzmax
    integer(wi) icalc_radial
    real(wp) norm, dummy

    integer(wi) kx_total, ky_total, kz_total, kr_total, isize_spec2dr

    integer(wi) inb_scal_min, inb_scal_max ! Iterval of scalars to calculate, to be able reduce memory constraints (hard coded)

! Reading variables
    character*512 sRes

    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi) iopt_size
    real(wp) opt_vec(iopt_size_max)

    integer, parameter :: i0 = 0 !, i1 = 1, i2 = 2, i3 = 3

! #ifdef USE_MPI
!     integer(wi) id
! #endif

    real(wp) params(1)

!########################################################################
!########################################################################
    bakfile = trim(adjustl(ifile))//'.bak'

    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)
    call Chemistry_Initialize(ifile)

    call TLab_Consistency_Check()

    ! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, ifile, 'IBMParameter', 'Status', 'off', sRes)
    if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
    else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
    else
        call TLab_Write_ASCII(efile, 'SPECTRA. Wrong IBM Status option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (y_aux(g(2)%size)) ! Reduced vertical grid

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Additional options
! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_block = 1
    opt_ffmt = 0
    opt_time = 0

    call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamSpectra', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_REAL(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Option ?'
        write (*, *) '1. Main variables 2D spectra'
        write (*, *) '2. Main variables 2D cross-spectra'
        write (*, *) '3. Main variables 2D correlation'
        write (*, *) '4. Main variables 2D cross-correlation'
        write (*, *) '5. Main variables 3D spectra'
        read (*, *) opt_main

        if (opt_main < 5) then
            write (*, *) 'Planes block size ?'
            read (*, *) opt_block
        end if
        write (*, *) 'Save full spectra fields to disk (1-yes/0-no) ?'
        read (*, *) opt_ffmt
        write (*, *) 'Average over time (1-yes/0-no) ?'
        read (*, *) opt_time
#endif
    else
        opt_main = int(opt_vec(1))
        if (iopt_size > 1) opt_block = int(opt_vec(2))
        if (iopt_size > 2) opt_ffmt = int(opt_vec(3))
        if (iopt_size > 3) opt_time = int(opt_vec(4))
    end if

    if (opt_main < 0) then ! Check
        call TLab_Write_ASCII(efile, 'SPECTRA. Missing input [ParamSpectra] in tlab.ini.')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    if (opt_block < 1) then
        call TLab_Write_ASCII(efile, 'SPECTRA. Invalid value of opt_block.')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    if (opt_time /= SPEC_SINGLE .and. opt_time /= SPEC_AVERAGE) then
        call TLab_Write_ASCII(efile, 'SPECTRA. Invalid value of opt_time.')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    ! -------------------------------------------------------------------
    ! Read local options - IBM parameters and geometry
    ! -------------------------------------------------------------------
    if (imode_ibm == 1) then
        call IBM_READ_INI(ifile)
    end if
! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
! in case g(2)%size is not divisible by opt_block, drop the upper most planes
    jmax_aux = g(2)%size/opt_block

    flag_buoyancy = 0 ! default

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
! in case we need the buoyancy statistics
        if (buoyancy%type == EQNS_BOD_QUADRATIC .or. &
            buoyancy%type == EQNS_BOD_BILINEAR .or. &
            imixture == MIXT_TYPE_AIRWATER .or. &
            imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            flag_buoyancy = 1
            inb_scal_array = inb_scal_array + 1             ! space for the buoyancy field
        end if
    end if

    if (opt_main == 1) then; flag_mode = 1 ! spectra
    else if (opt_main == 2) then; flag_mode = 1
    else if (opt_main == 3) then; flag_mode = 2 ! correlations
    else if (opt_main == 4) then; flag_mode = 2
    else if (opt_main == 5) then; flag_mode = 1 ! spectra
    end if

    if (flag_mode == 1) then                 ! spectra
        kxmax = imax/2; kymax = jmax/2; kzmax = kmax/2
    else                                              ! correlation
        kxmax = imax; kymax = jmax; kzmax = kmax
    end if
    isize_out2d = imax*jmax_aux*kmax                  ! accumulation of 2D data

! -------------------------------------------------------------------
!  maximum wavenumber & length lag; radial data is not really parallelized yet
    kx_total = max(g(1)%size/2, 1); ky_total = max(g(2)%size/2, 1); kz_total = max(g(3)%size/2, 1)

    if (opt_main == 4) then    ! Cross-correlations need the full length
        kx_total = max(g(1)%size, 1); ky_total = max(g(2)%size, 1); kz_total = max(g(3)%size, 1)
    end if

    if (opt_main >= 5) then ! 3D spectrum
!     kr_total =  INT(SQRT(real( (kx_total-1)**2 + (kz_total-1)**2 + (ky_total-1)**2))) + 1 ! Use if need to check Parseval's in output data
        kr_total = min(kx_total, min(ky_total, kz_total))
    else
!     kr_total =  INT(SQRT(real( (kx_total-1)**2 + (kz_total-1)**2))) + 1 ! Use if need to check Parseval's in output data
        kr_total = min(kx_total, kz_total)
    end if

    if (opt_main >= 5) then; isize_spec2dr = kr_total            ! 3D spectrum
    else; isize_spec2dr = kr_total*jmax_aux; end if

! -------------------------------------------------------------------
! Define MPI-type for writing spectra
! -------------------------------------------------------------------
#ifdef USE_MPI
    call SPECTRA_MPIO_AUX(opt_main, opt_block)
#else
    io_subarrays(:)%offset = 0
    io_subarrays(:)%precision = IO_TYPE_SINGLE
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
    iread_flow = flow_on
    iread_scal = scal_on

    inb_scal_min = 1              ! Change this values if you want to reduce the number of scalars to process
    inb_scal_max = inb_scal_array ! and thereby reduced memory requirements
    ! inb_scal_min = 4
    ! inb_scal_max = 4

    nfield_ref = 0     ! defining the number of accesible fields
    if (flow_on) nfield_ref = nfield_ref + inb_flow + 1 ! pressure
    if (scal_on) nfield_ref = nfield_ref + inb_scal_array

    nfield = 0          ! defining the number of accessed fields
    if (opt_main == 1 .or. opt_main == 3) then ! Auto-spectra & correlations
        if (flow_on) nfield = nfield + inb_flow + 1 ! pressure
        if (scal_on) nfield = nfield + (inb_scal_max - inb_scal_min + 1)

    else if (opt_main == 2 .or. opt_main == 4) then ! cross-spectra and cross-correlations
        if (flow_on) nfield = nfield + 3
        if (scal_on) nfield = nfield + 3*(inb_scal_max - inb_scal_min + 1)

    else
        nfield = nfield_ref

    end if

#ifdef IBM_DEBUG
    inb_txc = 6
#else
    inb_txc = 5 ! default
#endif

    isize_aux = jmax_aux
#ifdef USE_MPI
    if (ims_npro_k > 1) then
        if (mod(jmax_aux, ims_npro_k) /= 0) then
            isize_aux = ims_npro_k*(jmax_aux/ims_npro_k + 1)
        end if

        tmpi_plan_z = TLabMPI_Trp_TypeK_Create(kmax, isize_aux, 1, 1, 1, 1, 'type-2 Oz spectra integration.')

    end if
#endif

    isize_wrk2d = max(isize_wrk2d, isize_aux*kmax); inb_wrk2d = max(inb_wrk2d, 6)

    allocate (outx(kxmax*jmax_aux, nfield))
    allocate (outz(kzmax*jmax_aux, nfield))

    write (str, *) nfield; line = 'Allocating array outr  of size '//trim(adjustl(str))//'x'
    write (str, *) isize_spec2dr; line = trim(adjustl(line))//trim(adjustl(str))
    call TLab_Write_ASCII(lfile, line)
    allocate (outr(isize_spec2dr, nfield), stat=ierr)
    if (ierr /= 0) then
        call TLab_Write_ASCII(efile, 'SPECTRA. Not enough memory for spectral data.')
        call TLab_Stop(DNS_ERROR_ALLOC)
    end if

    if (flag_mode == 2) then
        allocate (samplesize(kr_total))
    end if

    if (opt_ffmt == 1) then ! need additional space for 2d spectra
        write (str, *) nfield; line = 'Allocating array out2d of size '//trim(adjustl(str))//'x'
        write (str, *) isize_out2d; line = trim(adjustl(line))//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, line)
        allocate (out2d(isize_out2d, nfield), stat=ierr)
        if (ierr /= 0) then
            call TLab_Write_ASCII(efile, 'SPECTRA. Not enough memory for spectral data.')
            call TLab_Stop(DNS_ERROR_ALLOC)
        end if
    end if

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        write (str, *) isize_txc_field; line = 'Allocating array p_aux of size '//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, line)
        allocate (p_aux(isize_txc_field), stat=ierr)
        if (ierr /= 0) then
            call TLab_Write_ASCII(efile, 'SPECTRA. Not enough memory for p_aux.')
            call TLab_Stop(DNS_ERROR_ALLOC)
        end if
    end if

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(C_FILE_LOC)
    end if

! extend array by complex nyquist frequency in x (+1 complex(wp) = +2 real(wp))
!              by boundary conditions in y       (+1 complex(wp) = +2 real(wp))

    isize_wrk3d = max(isize_wrk3d, isize_spec2dr) ! space needed in INTEGRATE_SPECTRUM

    call TLab_Initialize_Memory(C_FILE_LOC)

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1))
    call FDM_Initialize(y, g(2), wrk1d(:, 2))
    call FDM_Initialize(z, g(3), wrk1d(:, 3))

    call TLab_Initialize_Background(ifile)

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    end do

    call OPR_Burgers_Initialize(ifile)

    call OPR_Elliptic_Initialize(ifile)

    icalc_radial = 0
    if (flag_mode == 1 .and. g(1)%size == g(3)%size) icalc_radial = 1 ! Calculate radial spectra
    if (flag_mode == 2 .and. g(1)%jac(1, 1) == g(3)%jac(1, 1)) icalc_radial = 1 ! Calculate radial correlations

! ------------------------------------------------------------------------
! Define size of blocks
! ------------------------------------------------------------------------
    y_aux(:) = 0
    do j = 1, jmax
        is = (j - 1)/opt_block + 1
        y_aux(is) = y_aux(is) + y(j, 1)/real(opt_block, wp)
    end do

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
    if (fourier_on) then
        call OPR_FOURIER_INITIALIZE()
    end if

    call OPR_CHECK()

! -------------------------------------------------------------------
! Initialize IBM geometry
! -------------------------------------------------------------------
    if (imode_ibm == 1) then
        call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
    end if

! -------------------------------------------------------------------
! Initialize
! -------------------------------------------------------------------
    outx = 0.0_wp; outz = 0.0_wp; outr = 0.0_wp
    if (opt_ffmt == 1) out2d = 0.0_wp

! Normalization
    if (opt_main >= 5) then ! 3D spectra
        norm = 1.0_wp/real(g(1)%size*g(3)%size*g(2)%size, wp)
    else
        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)
    end if

! Define tags
    if (flag_mode == 1) then; tag_file = 'sp'; tag_name = 'E' ! spectra
    else if (flag_mode == 2) then; tag_file = 'cr'; tag_name = 'C' ! correlations
    end if

! Define reference pointers and tags
    iv = 0
    if (flow_on) then
        iv = iv + 1; vars(iv)%field => q(:, 1); tag_var(iv) = 'u'
        iv = iv + 1; vars(iv)%field => q(:, 2); tag_var(iv) = 'v'
        iv = iv + 1; vars(iv)%field => q(:, 3); tag_var(iv) = 'w'
        if (nse_eqns == DNS_EQNS_INTERNAL .or. nse_eqns == DNS_EQNS_TOTAL) then
            iv = iv + 1; vars(iv)%field => q(:, 6); tag_var(iv) = 'p'
            iv = iv + 1; vars(iv)%field => q(:, 5); tag_var(iv) = 'r'
            iv = iv + 1; vars(iv)%field => q(:, 7); tag_var(iv) = 't'
        else
            iv = iv + 1; vars(iv)%field => p_aux; tag_var(iv) = 'p'
        end if
    end if
    iv_offset = iv

    if (scal_on) then
        do is = 1, inb_scal_array
            write (sRes, *) is
            iv = iv + 1; vars(iv)%field => s(:, is); tag_var(iv) = trim(adjustl(sRes))
        end do
    end if

    if (nfield_ref /= iv) then ! Check
        call TLab_Write_ASCII(efile, 'SPECTRA. Array space nfield_ref incorrect.')
        call TLab_Stop(DNS_ERROR_WRKSIZE)
    end if

! Define pairs
    iv = 0
    if (opt_main == 1 .or. opt_main == 3) then ! Auto-spectra & correlations
        do ip = 1, iv_offset
            iv = iv + 1; p_pairs(iv, 1) = iv; p_pairs(iv, 2) = iv
        end do
        if (scal_on) then
            do is = inb_scal_min, inb_scal_max
                ip = is + iv_offset
                iv = iv + 1; p_pairs(iv, 1) = ip; p_pairs(iv, 2) = ip
            end do
        end if

    else if (opt_main == 2 .or. opt_main == 4) then ! Cross-spectra & correlations
        if (flow_on) then
            iv = iv + 1; p_pairs(iv, 1) = 1; p_pairs(iv, 2) = 2
            iv = iv + 1; p_pairs(iv, 1) = 1; p_pairs(iv, 2) = 3
            iv = iv + 1; p_pairs(iv, 1) = 2; p_pairs(iv, 2) = 3
            if (scal_on) then
                do is = inb_scal_min, inb_scal_max
                    ip = is + iv_offset
                    iv = iv + 1; p_pairs(iv, 1) = 1; p_pairs(iv, 2) = ip
                    iv = iv + 1; p_pairs(iv, 1) = 2; p_pairs(iv, 2) = ip
                    iv = iv + 1; p_pairs(iv, 1) = 3; p_pairs(iv, 2) = ip
                end do
            end if
! block to calculate the pressure-velocity and triple-velocity correlation terms in turbulent transport
            ! iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 4
            ! iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = 4
            ! iv = iv+1; p_pairs(iv,1) = 3; p_pairs(iv,2) = 4
            ! IF ( scal_on ) THEN ! aux array for u_iu_i/2
            !    s(:,1) = 0.5_wp*( q(:,1)*q(:,1) + q(:,2)*q(:,2) + q(:,3)*q(:,3) ); tag_var(5) = 'q'
            !    iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 5
            !    iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = 5
            !    iv = iv+1; p_pairs(iv,1) = 3; p_pairs(iv,2) = 5
            ! ENDIF
        else
            call TLab_Write_ASCII(efile, 'SPECTRA. Cross-spectra needs flow fields.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

    end if

    if (nfield /= iv) then ! Check
        call TLab_Write_ASCII(efile, 'SPECTRA. Array space nfield incorrect.')
        call TLab_Stop(DNS_ERROR_WRKSIZE)
    end if

    do iv = 1, nfield ! define variable names
        varname(iv) = tag_name(1:1)//trim(adjustl(tag_var(p_pairs(iv, 1))))//trim(adjustl(tag_var(p_pairs(iv, 2))))
    end do

! ###################################################################
! Calculating statistics
! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLab_Write_ASCII(lfile, sRes)

        if (iread_flow) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        if (iread_scal) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        end if

        if (imode_ibm == 1) then
            call IBM_BCS_FIELD_COMBINED(i0, q)
            if (scal_on) call IBM_INITIALIZE_SCAL(i0, s)
        end if

        call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

! Calculate additional diagnostic quantities to be processed
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            call FI_PRESSURE_BOUSSINESQ(q, s, p_aux, txc(1, 1), txc(1, 2), txc(1, 3), DCMP_TOTAL)
            if (flag_buoyancy == 1) then
                if (buoyancy%type == EQNS_EXPLICIT) then
                    call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, s(1, inb_scal_array))
                else
                    wrk1d(1:jmax, 1) = 0.0_wp
                    call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, s(1, inb_scal_array), wrk1d)
                end if
                dummy = 1.0_wp/froude
                s(:, inb_scal_array) = s(:, inb_scal_array)*dummy
            end if
        end if

! remove mean -- fluctuation only for spectrum
        if (opt_main >= 5) then ! 3D spectra
            do iv = 1, nfield_ref
                dummy = AVG1V2D(imax*jmax, 1, kmax, 1, 1, vars(iv)%field)  ! 3D average
                vars(iv)%field = vars(iv)%field - dummy
            end do
        else
            do iv = 1, nfield_ref
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, vars(iv)%field)
            end do
        end if

! If IBM is active: remove mean values in solid regions from fluctuations (except pressure)
        if (imode_ibm == 1) then
            if (flow_on) then
                do iv = 1, 3 ! u,v,w fields - skip pressure
                    call IBM_BCS_FIELD(vars(iv)%field)
                end do
                if (nfield_ref > 4) then
                    do iv = 5, nfield_ref ! r,t,(s) fields
                        call IBM_BCS_FIELD(vars(iv)%field)
                    end do
                end if
            else if ((.not. flow_on) .and. (scal_on)) then
                do iv = 1, nfield_ref    ! s fields (no pressure)
                    call IBM_BCS_FIELD(vars(iv)%field)
                end do
            end if
        end if

! reset if needed
        if (opt_time == SPEC_SINGLE) then
            outx = 0.0_wp; outz = 0.0_wp; outr = 0.0_wp
            if (opt_ffmt == 1) out2d = 0.0_wp
        end if

! ###################################################################
! 2D Spectra & Correlations
! ###################################################################
        if (opt_main <= 4) then

! -------------------------------------------------------------------
! Calculate 2d spectra into array out2d and 1d spectra into arrays outX
! -------------------------------------------------------------------
            do iv = 1, nfield
                iv1 = p_pairs(iv, 1); iv2 = p_pairs(iv, 2)

                wrk1d(:, 1:3) = 0.0_wp ! variance to normalize and check Parseval's relation
                do j = 1, jmax
                    wrk1d(j, 1) = COV2V2D(imax, jmax, kmax, j, vars(iv1)%field, vars(iv2)%field)
                    wrk1d(j, 2) = COV2V2D(imax, jmax, kmax, j, vars(iv1)%field, vars(iv1)%field)
                    wrk1d(j, 3) = COV2V2D(imax, jmax, kmax, j, vars(iv2)%field, vars(iv2)%field)
                end do

                txc(1:isize_field, 1) = vars(iv1)%field(1:isize_field)
                if (iv2 == iv1) then
                    call OPR_FOURIER_CONVOLUTION_FXZ('auto', flag_mode, imax, jmax, kmax, &
                                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                else
                    txc(1:isize_field, 2) = vars(iv2)%field(1:isize_field)
                    call OPR_FOURIER_CONVOLUTION_FXZ('cross', flag_mode, imax, jmax, kmax, &
                                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                end if

                if (flag_mode == 1) then ! Spectra
                    txc(:, 1) = txc(:, 1)*norm*norm

! Reduce 2D spectra into array wrk3d
                    wrk3d = 0.0_wp
                    call REDUCE_SPECTRUM(imax, jmax, kmax, opt_block, &
                                         (txc(:, 1)), wrk3d, txc(:, 3), wrk1d(1, 4))

! Calculate and accumulate 1D spectra; only the half of wrk3d with the power data is necessary
                    call INTEGRATE_SPECTRUM(imax/2, jmax_aux, kmax, kr_total, isize_aux, &
                                            wrk3d, outx(1, iv), outz(1, iv), outr(1, iv), wrk2d(1, 1), wrk2d(1, 3), wrk2d(1, 5))

                else if (flag_mode == 2) then  ! Correlations
                    txc(:, 2) = txc(:, 2)*norm*norm

! Reduce 2D correlation into array wrk3d and accumulate 1D correlation
                    wrk3d = 0.0_wp
                    call REDUCE_CORRELATION(imax, jmax, kmax, opt_block, kr_total, &
                                            txc(1, 2), wrk3d, outx(1, iv), outz(1, iv), outr(1, iv), wrk1d(1, 2), wrk1d(1, 4), icalc_radial)
                end if

! Check Parseval's relation
                ip = g(2)%size - mod(g(2)%size, opt_block)  ! Drop the uppermost ny%nblock
                write (line, 100) maxval(abs(wrk1d(1:ip, 4) - wrk1d(1:ip, 1)))
                write (str, *) maxloc(abs(wrk1d(1:ip, 4) - wrk1d(1:ip, 1)))
                line = 'Checking Parseval: Maximum residual '//trim(adjustl(line))//' at level '//trim(adjustl(str))//'.'
                call TLab_Write_ASCII(lfile, line)

! Accumulate 2D information, if needed
                if (opt_ffmt == 1) out2d(1:isize_out2d, iv) = out2d(1:isize_out2d, iv) + wrk3d(1:isize_out2d)

            end do

            if (flag_mode == 2 .and. icalc_radial == 1) then  ! Calculate sampling size for radial correlation
                samplesize = 0.0_wp
                call RADIAL_SAMPLESIZE(imax, kmax, kr_total, samplesize)
            end if

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
            if (opt_time == SPEC_SINGLE .or. it == itime_size) then

! Normalizing accumulated spectra
                ip = opt_block
                if (opt_time == SPEC_AVERAGE) ip = ip*itime_size
                dummy = 1.0_wp/real(ip, wp)
                if (ip > 1) then
                    outx = outx*dummy; outz = outz*dummy; outr = outr*dummy
                    if (opt_ffmt == 1) out2d = out2d*dummy
                end if

! Reducing radial data
#ifdef USE_MPI
                do iv = 1, nfield
                    call MPI_Reduce(outr(1, iv), wrk3d, isize_spec2dr, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
                    if (ims_pro == 0) outr(1:isize_spec2dr, iv) = wrk3d(1:isize_spec2dr)
                end do

                if (flag_mode == 2 .and. icalc_radial == 1) then ! Calculate sampling size for radial correlation
                    call MPI_Reduce(samplesize, wrk3d, kr_total, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
                    if (ims_pro == 0) samplesize(1:kr_total) = wrk3d(1:kr_total)
                end if
#endif

! Normalize radial correlation
#ifdef USE_MPI
                if (ims_pro == 0) then
#endif
                    if (flag_mode == 2 .and. icalc_radial == 1) then
                        do iv1 = 1, kr_total
                            if (samplesize(iv1) > 0.0_wp) samplesize(iv1) = 1.0_wp/samplesize(iv1)
                        end do

                        do iv = 1, nfield
                            do j = 1, jmax_aux
                                do iv1 = 1, kr_total
                                    ip = iv1 + (j - 1)*kr_total
                                    outr(ip, iv) = outr(ip, iv)*samplesize(iv1)
                                end do
                            end do
                        end do

                    end if
#ifdef USE_MPI
                end if
#endif

! Saving 1D fields
                if (opt_time == SPEC_AVERAGE) then
                    write (str, *) itime; write (fname, *) itime_vec(1); str = trim(adjustl(fname))//'-'//trim(adjustl(str))
                else
                    write (str, *) itime; 
                end if
                fname = 'x'//trim(adjustl(tag_file))//trim(adjustl(str))
                sizes(1) = kxmax*jmax_aux; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                call IO_Write_Subarray(io_subarrays(1), fname, varname, outx, sizes)

                if (g(3)%size > 1) then
                    fname = 'z'//trim(adjustl(tag_file))//trim(adjustl(str))
                    sizes(1) = kzmax*jmax_aux; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                    call IO_Write_Subarray(io_subarrays(2), fname, varname, outz, sizes)
                end if

                if (icalc_radial == 1) then
                    fname = 'r'//trim(adjustl(tag_file))//trim(adjustl(str))
                    call WRITE_SPECTRUM1D(fname, varname, kr_total*jmax_aux, nfield, outr)
                end if

! Saving 2D fields
                if (opt_ffmt == 1) then
                    if (flag_mode == 2) then ! correlations
                        fname = 'cor'//trim(adjustl(str))
                        sizes(1) = isize_out2d; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                        call IO_Write_Subarray(io_subarrays(3), fname, varname, out2d, sizes)

                    else                         ! spectra
                        fname = 'pow'//trim(adjustl(str))
                        sizes(1) = isize_out2d; sizes(2) = 1; sizes(3) = sizes(1)/2; sizes(4) = 1; sizes(5) = nfield
                        call IO_Write_Subarray(io_subarrays(3), fname, varname, out2d, sizes)

                        fname = 'pha'//trim(adjustl(str))
                        sizes(1) = isize_out2d; sizes(2) = 1 + sizes(1)/2; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                        call IO_Write_Subarray(io_subarrays(3), fname, varname, out2d, sizes)

                    end if

                end if

            end if

! ###################################################################
! 3D Spectra
! ###################################################################
        else if (opt_main == 5) then

            do iv = 1, nfield
                txc(1:isize_field, 1) = vars(iv)%field(1:isize_field)
                call OPR_FOURIER_F(3, imax, jmax, kmax, txc(1, 1), txc(1, 2), txc(1, 3))

                call OPR_FOURIER_SPECTRA_3D(imax, jmax, kmax, isize_spec2dr, txc(1, 2), outr(1, iv))
            end do

            outr = outr*norm*norm

            write (fname, *) itime; fname = 'rsp'//trim(adjustl(fname))
            call WRITE_SPECTRUM1D(fname, varname, kr_total, nfield, outr)

        end if

    end do ! Loop in itime

100 format(G_FORMAT_R)
    call TLab_Stop(0)
end program SPECTRA

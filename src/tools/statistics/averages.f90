#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

#define C_FILE_LOC "AVERAGES"

program AVERAGES

    use TLab_Pointers, only: pointers_dt
    use TLab_Constants, only: wp, wi, small_wp, MAX_AVG_TEMPORAL
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal, tag_part
    use TLAB_VARS
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use FDM, only: g, FDM_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize, buoyancy, Gravity_Buoyancy, Gravity_Buoyancy_Source
    use Thermodynamics, only: imixture, Thermodynamics_Initialize_Parameters
    use THERMO_ANELASTIC
    use Radiation
    use Microphysics
    use Chemistry
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS
    use IBM_VARS
    use IO_FIELDS
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    use OPR_PARTIAL
    use OPR_FOURIER
    use OPR_FILTERS
    use OPR_Burgers, only: OPR_Burgers_Initialize
    use OPR_ELLIPTIC
    use AVG_PHASE

    implicit none

    ! Parameter definitions
    integer, parameter :: itime_size_max = 512
    integer, parameter :: iopt_size_max = 20
    integer, parameter :: igate_size_max = 8
    integer, parameter :: params_size_max = 2

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), allocatable, save :: mean(:), y_aux(:)
    integer(1), allocatable, save :: gate(:)
    type(pointers_dt) :: vars(16)
    real(wp), allocatable, save :: surface(:, :, :)       ! Gate envelopes

    integer, parameter :: i1 = 1  ! to be removed, still needed for reduce procedures

    ! -------------------------------------------------------------------
    ! Local variables
    ! -------------------------------------------------------------------
    character*512 sRes
    character*32 fname, bakfile
    character*32 varname(16)

    integer opt_main, opt_block, opt_order
    integer opt_cond, opt_cond_scal, opt_cond_relative
    integer nfield, ifield, is, k, bcs(2, 2), ig
    real(wp) eloc1, eloc2, eloc3, cos1, cos2, cos3, dummy
    integer jmax_aux
    logical iread_flow, iread_scal
    integer(wi) ij, idummy

    ! Gates for the definition of the intermittency function (partition of the fields)
    integer(wi) igate_size
    real(wp) gate_threshold(igate_size_max)
    integer(1) gate_level

    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi) iopt_size
    real(wp) opt_vec(iopt_size_max)
    real(wp) opt_vec2(iopt_size_max)

    integer(wi) params_size
    real(wp) params(params_size_max)

    integer(wi) io_sizes(5), id

    ! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w, p

    !########################################################################
    !########################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif
    call Particle_Initialize_Parameters(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)
    call Chemistry_Initialize(ifile)

    ! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, ifile, 'IBMParameter', 'Status', 'off', sRes)
    if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
    else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
    else
        call TLab_Write_ASCII(efile, 'AVERAGES. Wrong IBM Status option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    ! -------------------------------------------------------------------
    ! File names
    ! -------------------------------------------------------------------
#include "dns_read_times.h"

    ! -------------------------------------------------------------------
    ! Additional options
    ! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_block = 1
    gate_level = 0
    opt_order = 1

    call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamAverages', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_REAL(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Option ?'
        write (*, *) ' 1. Conventional averages'
        write (*, *) ' 2. Intermittency or gate function'
        write (*, *) ' 3. Momentum equation'
        write (*, *) ' 4. Main variables'
        write (*, *) ' 5. Enstrophy W_iW_i/2 equation'
        write (*, *) ' 6. Strain 2S_ijS_ij/2 equation'
        write (*, *) ' 7. Scalar gradient G_iG_i/2 equation'
        write (*, *) ' 8. Velocity gradient invariants'
        write (*, *) ' 9. Scalar gradient components'
        write (*, *) '10. Eigenvalues of rate-of-strain tensor'
        write (*, *) '11. Eigenframe of rate-of-strain tensor'
        write (*, *) '12. Longitudinal velocity derivatives'
        write (*, *) '13. Vertical fluxes'
        write (*, *) '14. Pressure partition'
        write (*, *) '15. Dissipation'
        write (*, *) '16. Third-order scalar covariances'
        write (*, *) '17. Potential vorticity'
        write (*, *) '18. Phase Average'
        read (*, *) opt_main

        write (*, *) 'Planes block size ?'
        read (*, *) opt_block

        if (opt_main > 2) then
            write (*, *) 'Gate level to be used ?'
            read (*, *) gate_level
            write (*, *) 'Number of moments ?'
            read (*, *) opt_order
        end if

#endif
    else
        opt_main = int(opt_vec(1))
        if (iopt_size >= 2) opt_block = int(opt_vec(2))
        if (iopt_size >= 3) gate_level = int(opt_vec(3), KIND=1)
        if (iopt_size >= 3) opt_order = int(opt_vec(4))

    end if

    if (opt_main < 0) then ! Check
        call TLab_Write_ASCII(efile, 'AVERAGES. Missing input [ParamAverages] in tlab.ini.')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    if (opt_block < 1) then
        call TLab_Write_ASCII(efile, 'AVERAGES. Invalid value of opt_block.')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    ! -------------------------------------------------------------------
    iread_flow = .false.
    iread_scal = .false.
    inb_txc = 0
    nfield = 2

    select case (opt_main)
    case (1)
        inb_txc = max(inb_txc, 9)
        iread_flow = flow_on; iread_scal = scal_on
    case (2)
        fourier_on = .false.
    case (3)
        nfield = 14
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 12)
    case (4)
        nfield = 6 + inb_scal
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 3)
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) inb_txc = max(inb_txc, 6)
    case (5) ! enstrophy
        nfield = 7
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 8)
    case (6)
        nfield = 5
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 8)
    case (7) ! scalar gradient
        nfield = 5
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6)
    case (8)
        nfield = 3
        iread_flow = .true.; inb_txc = max(inb_txc, 6)
    case (9)
        nfield = 5
        iread_scal = .true.; inb_txc = max(inb_txc, 4)
    case (10) ! eigenvalues
        nfield = 3
        iread_flow = .true.; inb_txc = max(inb_txc, 9)
    case (11) ! eigenframe
        nfield = 6
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 10)
    case (12) ! longitudinal velocity derivatives
        nfield = 3
        iread_flow = .true.; iread_scal = .false.; inb_txc = max(inb_txc, 3)
    case (13) ! Vertical flux
        nfield = 2*(3 + inb_scal_array)
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(max(inb_txc, 3 + inb_scal_array), 4)
    case (14) ! pressure partition
        nfield = 3
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
    case (15) ! dissipation partition
        nfield = 1
        iread_flow = .true.; iread_scal = .false.; inb_txc = max(inb_txc, 6)
    case (16) ! third-order scalar covariances
        nfield = 3
        iread_flow = .false.; iread_scal = .true.; inb_txc = max(inb_txc, 3)
    case (17) ! potential vorticity
        nfield = 2
        iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6)
    case (18) ! Phase average
        PhAvg%active = .true.
        PhAvg%stride = 1
        nfield = 3
        inb_txc = max(inb_txc, 9)
        iread_flow = flow_on; iread_scal = scal_on
    end select

    if (imode_ibm == 1) then ! check if enough memory is provided for the IBM
#ifdef IBM_DEBUG
        inb_txc = max(inb_txc, 6)
#else
        inb_txc = max(inb_txc, 3)
#endif
    end if

    ! -------------------------------------------------------------------
    ! Read local options - IBM parameters and geometry
    ! -------------------------------------------------------------------
    if (imode_ibm == 1) then
        call IBM_READ_INI(ifile)
    end if

    ! -------------------------------------------------------------------
    ! Defining gate levels for conditioning
    ! -------------------------------------------------------------------
    opt_cond = 0 ! default values
    opt_cond_relative = 0
    igate_size = 0

    if (opt_main == 2 .or. gate_level /= 0) then
#include "dns_read_partition.h"
        if (opt_cond > 1) inb_txc = max(inb_txc, 5)
    end if

    ! -------------------------------------------------------------------
    ! Allocating memory space
    ! -------------------------------------------------------------------
    allocate (gate(isize_field))

    ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
    jmax_aux = g(2)%size/opt_block
    allocate (y_aux(jmax_aux)) ! Reduced vertical grid

    ! Subarray information to read and write envelopes
    if (opt_main == 2 .and. opt_cond > 1 .and. igate_size /= 0) then
        ! Info for IO routines: total size, lower bound, upper bound, stride, # variables
        idummy = imax*2*igate_size*kmax; io_sizes = [idummy, 1, idummy, 1, 1]
        id = IO_SUBARRAY_ENVELOPES
        io_aux(id)%offset = 0
        io_aux(id)%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
        io_aux(id)%active = .true.
        io_aux(id)%communicator = MPI_COMM_WORLD
        io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ(imax, igate_size*2, kmax, MPI_REAL4)
#endif

        allocate (surface(imax, 2*igate_size, kmax))

    end if

    if (opt_main == 1) then
        allocate (mean(jmax*MAX_AVG_TEMPORAL))
    else if (opt_main == 2) then
        allocate (mean(igate_size*(jmax_aux + 1)))
    else
        allocate (mean(opt_order*nfield*(jmax_aux + 1)))
    end if

    isize_wrk3d = max(isize_wrk3d, opt_order*nfield*jmax)

    call TLab_Initialize_Memory(C_FILE_LOC)

    call Particle_Initialize_Memory(C_FILE_LOC)

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(C_FILE_LOC)
    end if

    if (opt_main == 18) then
        call AvgPhaseInitializeMemory(__FILE__, -1)
    end if
    ! -------------------------------------------------------------------
    ! Initialize
    ! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1), wrk1d(:, 4))
    call FDM_Initialize(y, g(2), wrk1d(:, 2), wrk1d(:, 4))
    call FDM_Initialize(z, g(3), wrk1d(:, 3), wrk1d(:, 4))

    call TLab_Initialize_Background(ifile)  ! Initialize thermodynamic quantities

    call OPR_Burgers_Initialize(ifile)

    call OPR_Elliptic_Initialize(ifile)

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    end do

    if (fourier_on) then         ! For Poisson solver
        call OPR_FOURIER_INITIALIZE()
    end if

    if (iread_flow) then       ! We need array space
        call OPR_CHECK()
    end if

    y_aux(:) = 0                        ! Reduced vertical grid
    do ij = 1, jmax_aux*opt_block
        is = (ij - 1)/opt_block + 1
        y_aux(is) = y_aux(is) + y(ij, 1)/real(opt_block, wp)
    end do

    if (iread_flow) then
        u => q(:, 1)
        v => q(:, 2)
        w => q(:, 3)
    end if

    if (imode_ibm == 1) then
        call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
    end if

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLab_Write_ASCII(lfile, sRes)

        if (iread_scal) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        end if

        if (iread_flow) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, 0, q)
        end if

        if (imode_ibm == 1) then
            call IBM_BCS_FIELD_COMBINED(0, q)
            if (scal_on) call IBM_INITIALIZE_SCAL(0, s)
        end if

        call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

        ! -------------------------------------------------------------------
        ! Calculate intermittency
        ! -------------------------------------------------------------------
        if (opt_cond == 1) then ! External file
            write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname)); params_size = 2
            call IO_READ_FIELD_INT1(fname, 1, imax, jmax, kmax, itime, params_size, params, gate)
            igate_size = int(params(2))

            if (opt_main == 2) rtime = params(1)

        else if (opt_cond > 1) then
            opt_cond_scal = 1
            if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                opt_cond_scal = inb_scal_array
            end if

            call TLab_Write_ASCII(lfile, 'Calculating partition...')
            call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                         imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate)

            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE_INT1(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, gate, wrk1d)
            end if

        end if

        ! -------------------------------------------------------------------
        ! Type of averages
        ! -------------------------------------------------------------------
        select case (opt_main)

            ! ###################################################################
            ! Conventional statistics
            ! ###################################################################
        case (1)
            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 9), txc(1, 1), txc(1, 2), txc(1, 4), DCMP_TOTAL)
            end if

            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    txc(1:isize_field, 6) = txc(1:isize_field, 9) ! Pass the pressure in tmp6
                    call AVG_SCAL_XZ(is, q, s, s(1, is), &
                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)
                end do

                ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
                if (buoyancy%type == EQNS_BOD_QUADRATIC .or. buoyancy%type == EQNS_BOD_BILINEAR .or. &
                    imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                    if (buoyancy%type == EQNS_EXPLICIT) then
                        call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, txc(1, 7))
                    else
                        wrk1d(1:jmax, 1) = 0.0_wp
                        call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, txc(1, 7), wrk1d)
                    end if
                    dummy = 1.0_wp/froude
                    txc(1:isize_field, 7) = txc(1:isize_field, 7)*dummy

                    txc(1:isize_field, 6) = txc(1:isize_field, 9) ! Pass the pressure in tmp6
                    call AVG_SCAL_XZ(is, q, s, txc(1, 7), &
                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)

                end if

                if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                    if (imixture == MIXT_TYPE_AIRWATER) then
                        is = is + 1
                        call THERMO_ANELASTIC_THETA_L(imax, jmax, kmax, s, txc(1, 7))
                        !                 CALL THERMO_ANELASTIC_STATIC_CONSTANTCP(imax,jmax,kmax, s, txc(1,7))
                        txc(1:isize_field, 6) = txc(1:isize_field, 9) ! Pass the pressure in tmp6
                        call AVG_SCAL_XZ(is, q, s, txc(1, 7), &
                                         txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)
                    end if
                end if

            end if

            ! Lagrange Liquid and Liquid without diffusion
            if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then
                write (fname, *) itime; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
                call IO_READ_PARTICLE(fname, l_g, l_q)

                l_txc = 1.0_wp; ! We want density
                call PARTICLE_TO_FIELD(l_q, l_txc, txc(1, 7), wrk3d)

                txc(:, 7) = txc(:, 7) + small_wp
                idummy = inb_part - 3 ! # scalar properties solved in the lagrangian
                do is = inb_scal_array + 1 + 1, inb_scal_array + 1 + idummy
                    schmidt(is) = schmidt(1)
                    l_txc(:, 1) = l_q(:, 3 + is - inb_scal_array - 1) !!! DO WE WANT l_txc(:,is) ???
                    call PARTICLE_TO_FIELD(l_q, l_txc, txc(1, 8), wrk3d)
                    txc(:, 8) = txc(:, 8)/txc(:, 7)
                    txc(1:isize_field, 6) = txc(1:isize_field, 9) ! Pass the pressure in tmp6
                    call AVG_SCAL_XZ(is, q, s, txc(1, 8), &
                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)
                end do
            end if

            if (flow_on) then
                txc(1:isize_field, 3) = txc(1:isize_field, 9) ! Pass the pressure in tmp3
                call AVG_FLOW_XZ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), &
                                 txc(1, 7), txc(1, 8), txc(1, 9), mean)
            end if

            ! ###################################################################
            ! Partition of field
            ! ###################################################################
        case (2)
            do is = 1, igate_size
                write (varname(is), *) is; varname(is) = 'Partition'//trim(adjustl(varname(is)))
            end do
            write (fname, *) itime; fname = 'int'//trim(adjustl(fname))
            call INTER_N_XZ(fname, itime, rtime, imax, jmax, kmax, igate_size, varname, gate, y, mean)

            if (opt_cond > 1) then ! write only if the gate information has not been read
                write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname))
                params(1) = rtime; params(2) = real(igate_size, wp); params_size = 2
                call IO_WRITE_FIELD_INT1(fname, i1, imax, jmax, kmax, itime, params_size, params, gate)

                do is = 1, igate_size
                    gate_level = int(is, KIND=1)
                    call BOUNDARY_LOWER_INT1(imax, jmax, kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1, 2))
                    do k = 1, kmax ! rearranging
                        ij = (k - 1)*imax + 1
                        surface(1:imax, is, k) = wrk2d(ij:ij + imax - 1, 1)
                    end do
                    call BOUNDARY_UPPER_INT1(imax, jmax, kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1, 2))
                    do k = 1, kmax ! rearranging
                        ij = (k - 1)*imax + 1
                        surface(1:imax, is + igate_size, k) = wrk2d(ij:ij + imax - 1, 1)
                    end do
                end do
                varname = ''
                write (fname, *) itime; fname = 'envelopesJ.'//trim(adjustl(fname))
                call IO_WRITE_SUBARRAY(io_aux(IO_SUBARRAY_ENVELOPES), fname, varname, surface, io_sizes)

            end if

            ! ###################################################################
            ! Momentum equation
            ! ###################################################################
        case (3)
            write (fname, *) itime; fname = 'avgMom'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'U'
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'W'

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Uy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Uyy'
            call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), q(1, 1), txc(1, 2), txc(1, 1))
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Wy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Wyy'
            call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), q(1, 3), txc(1, 4), txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'VU)y'
            txc(1:isize_field, 6) = q(1:isize_field, 2)*q(1:isize_field, 1)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 6), txc(1, 5))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'VUy'
            txc(1:isize_field, 6) = q(1:isize_field, 2)*txc(1:isize_field, 1)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'UUx'
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(1, 1), txc(1, 7))
            txc(1:isize_field, 7) = q(1:isize_field, 1)*txc(1:isize_field, 7)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'WUz'
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(1, 1), txc(1, 8))
            txc(1:isize_field, 8) = q(1:isize_field, 3)*txc(1:isize_field, 8)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'WV)y'
            txc(1:isize_field, 10) = q(1:isize_field, 2)*q(1:isize_field, 3)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 10), txc(1, 9))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 10); vars(ifield)%tag = 'VWy'
            txc(1:isize_field, 10) = q(1:isize_field, 2)*txc(1:isize_field, 3)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 11); vars(ifield)%tag = 'UWx'
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(1, 3), txc(1, 11))
            txc(1:isize_field, 11) = q(1:isize_field, 1)*txc(1:isize_field, 11)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 12); vars(ifield)%tag = 'WWz'
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(1, 3), txc(1, 12))
            txc(1:isize_field, 12) = q(1:isize_field, 3)*txc(1:isize_field, 12)

            ! ###################################################################
            ! Main variables
            ! ###################################################################
        case (4)
            write (fname, *) itime; fname = 'avgMain'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ifield = ifield + 1; vars(ifield)%field => u(:); vars(ifield)%tag = 'U'
            ifield = ifield + 1; vars(ifield)%field => v(:); vars(ifield)%tag = 'V'
            ifield = ifield + 1; vars(ifield)%field => w(:); vars(ifield)%tag = 'W'

            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), DCMP_TOTAL)
                ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'P'
            else
                ifield = ifield + 1; vars(ifield)%field => q(:, 5); vars(ifield)%tag = 'R'
                ifield = ifield + 1; vars(ifield)%field => q(:, 6); vars(ifield)%tag = 'P'
                ifield = ifield + 1; vars(ifield)%field => q(:, 7); vars(ifield)%tag = 'T'
            end if

            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    ifield = ifield + 1; vars(ifield)%field => s(:, is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'Scalar'//trim(adjustl(vars(ifield)%tag))
                end do
            end if

            ! ###################################################################
            ! Enstrophy equation
            ! ###################################################################
        case (5)
            write (fname, *) itime; fname = 'avgW2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ! result vector in txc4, txc5, txc6
            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                if (buoyancy%type == EQNS_NONE) then
                    txc(:, 4) = 0.0_wp; txc(:, 5) = 0.0_wp; txc(:, 6) = 0.0_wp
                else
                    if (buoyancy%type == EQNS_EXPLICIT) then
                        call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, wrk3d)
                    else
                        wrk1d(1:jmax, 1) = 0.0_wp
                        call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, wrk3d, wrk1d)
                    end if
                    do ij = 1, isize_field
                        s(ij, 1) = wrk3d(ij)*buoyancy%vector(2)
                    end do

                    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 4))
                    txc(:, 4) = -txc(:, 4)
                    txc(:, 5) = 0.0_wp
                    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 6))
                end if

            else
                call FI_VORTICITY_BAROCLINIC(imax, jmax, kmax, q(1, 5), q(1, 6), txc(1, 4), txc(1, 3), txc(1, 7))
            end if

            call FI_CURL(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 7))
            txc(1:isize_field, 8) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 5) + txc(1:isize_field, 3)*txc(1:isize_field, 6)

            call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                         txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, u, v, w, txc(1, 2), &
                                        txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = visc*txc(1:isize_field, 2)

            call FI_VORTICITY(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5))  ! Enstrophy
            call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 4), txc(1, 5))  ! Dilatation

            txc(1:isize_field, 6) = txc(1:isize_field, 4)*txc(1:isize_field, 3) ! -w^2 div(u)
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3) ! production rate
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))                  ! ln(w^2)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'EnstrophyW_iW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnEnstrophyW_iW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionW_iW_jS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuW_iLapW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'DilatationMsW_iW_iDivU'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Baroclinic'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'RateAN_iN_jS_ij'

            ! ###################################################################
            ! Strain equation
            ! ###################################################################
        case (6)
            write (fname, *) itime; fname = 'avgS2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), DCMP_TOTAL)
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            else
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, q(1, 6), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            end if
            txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_PRODUCTION(imax, jmax, kmax, u, v, w, &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_DIFFUSION(imax, jmax, kmax, u, v, w, &
                                     txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
            txc(1:isize_field, 3) = 2.0_wp*visc*txc(1:isize_field, 3)

            call FI_STRAIN(imax, jmax, kmax, u, v, w, txc(1, 4), txc(1, 5), txc(1, 6))
            txc(1:isize_field, 4) = 2.0_wp*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = log(txc(1:isize_field, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Strain2S_ijS_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'LnStrain2S_ijS_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'ProductionMs2S_ijS_jkS_ki'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'DiffusionNuS_ijLapS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Pressure2S_ijP_ij'

            ! ###################################################################
            ! Scalar gradient equation
            ! ###################################################################
        case (7)
            write (fname, *) itime; fname = 'avgG2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s, u, v, w, &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s, &   ! array u used as auxiliar
                                       txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), u)
            txc(1:isize_field, 2) = txc(1:isize_field, 2)*visc/schmidt(inb_scal)

            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 3), txc(1, 4))
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3)
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GradientG_iG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnGradientG_iG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionMsG_iG_jS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuG_iLapG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'StrainAMsN_iN_jS_ij'

            ! ###################################################################
            ! Velocity gradient invariants
            ! ###################################################################
        case (8)
            write (fname, *) itime; fname = 'avgInv'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_INVARIANT_R(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_INVARIANT_Q(imax, jmax, kmax, u, v, w, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'InvariantP'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'InvariantQ'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'InvariantR'

            ! ###################################################################
            ! Scalar gradient components
            ! ###################################################################
        case (9)
            write (fname, *) itime; fname = 'avgGi'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 1))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s, txc(1, 2))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 3))
            do ij = 1, isize_field                       ! Angles; s array is overwritten to save space
                dummy = txc(ij, 2)/sqrt(txc(ij, 1)*txc(ij, 1) + txc(ij, 2)*txc(ij, 2) + txc(ij, 3)*txc(ij, 3))
                txc(ij, 4) = asin(dummy)                  ! with Oy
                s(ij, 1) = atan2(txc(ij, 3), txc(ij, 1))    ! with Ox in plane xOz
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'GradientX'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'GradientY'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GradientZ'
            ifield = ifield + 1; vars(ifield)%field => s(:, 1); vars(ifield)%tag = 'Theta'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Phi'

            ! ###################################################################
            ! eigenvalues of rate-of-strain tensor
            ! ###################################################################
        case (10)
            write (fname, *) itime; fname = 'avgEig'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_STRAIN_TENSOR(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'Lambda1'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Lambda2'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'Lambda3'

            ! ###################################################################
            ! eigenframe of rate-of-strain tensor
            ! ###################################################################
        case (11)
            write (fname, *) itime; fname = 'avgCos'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_STRAIN_TENSOR(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7))  ! txc7-txc9
            call TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1, 1), txc(1, 7))   ! txc1-txc6

            call FI_CURL(imax, jmax, kmax, u, v, w, txc(1, 7), txc(1, 8), txc(1, 9), txc(1, 10))
            do ij = 1, isize_field                                             ! local direction cosines of vorticity vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                u(ij) = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                v(ij) = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                w(ij) = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
            end do

            ifield = ifield + 1; vars(ifield)%field => u; vars(ifield)%tag = 'cos(w,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => v; vars(ifield)%tag = 'cos(w,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => w; vars(ifield)%tag = 'cos(w,lambda3)'

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 7))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s, txc(1, 8))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 9))
            do ij = 1, isize_field                                             ! local direction cosines of scalar gradient vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                cos1 = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                cos2 = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                cos3 = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
                txc(ij, 7) = cos1; txc(ij, 8) = cos2; txc(ij, 9) = cos3
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'cos(G,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'cos(G,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'cos(G,lambda3)'

            ! ###################################################################
            ! longitudinal velocity derivatives
            ! ###################################################################
        case (12)
            write (fname, *) itime; fname = 'avgDer'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, txc(1, 1))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, txc(1, 2))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'dudx'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'dvdy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'dwdz'

            ! ###################################################################
            ! Vertical fluxes
            ! ###################################################################
        case (13)
            ifield = 0
            write (fname, *) itime; fname = 'avgFluxY'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, txc(:, 1))
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, txc(:, 2))
            txc(:, 1) = (txc(:, 1) + txc(:, 2))*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'tauyx'

            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, txc(:, 2))
            txc(:, 2) = txc(:, 2)*2.0_wp*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'tauyy'

            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, txc(:, 3))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, txc(:, 4))
            txc(:, 3) = (txc(:, 3) + txc(:, 4))*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'tauyz'

            do is = 1, inb_scal_array
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), txc(:, 3 + is))
                txc(:, 3 + is) = txc(:, 3 + is)*visc/schmidt(inb_scal)
                ifield = ifield + 1; vars(ifield)%field => txc(:, 3 + is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'tauy'//trim(adjustl(vars(ifield)%tag))
            end do

            u = u*v
            ifield = ifield + 1; vars(ifield)%field => u; vars(ifield)%tag = 'vu'
            ! I need v below
            ifield = ifield + 1; vars(ifield)%field => v; vars(ifield)%tag = 'vv'
            w = w*v
            ifield = ifield + 1; vars(ifield)%field => w; vars(ifield)%tag = 'vw'
            do is = 1, inb_scal_array
                s(:, is) = s(:, is)*v
                ifield = ifield + 1; vars(ifield)%field => s(:, is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'v'//trim(adjustl(vars(ifield)%tag))
            end do
            v = v*v ! I need v above for the scalar fluxes

            ! ###################################################################
            ! Hydrostatic and dynamic pressure
            ! ###################################################################
        case (14)
            write (fname, *) itime; fname = 'avgP'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), DCMP_TOTAL)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'P'

            q = 0.0_wp
            call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), DCMP_TOTAL)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Psta'

            txc(:, 3) = txc(:, 1) - txc(:, 2)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Pdyn'

            ! ###################################################################
            ! Dissipation
            ! ###################################################################
        case (15)
            write (fname, *) itime; fname = 'avgEps'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_DISSIPATION(i1, imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Eps'

            ! ###################################################################
            ! Covariances among scalars
            ! ###################################################################
        case (16)
            write (fname, *) itime; fname = 'avgSiCov'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, s(1, 1))
            call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, s(1, 2))

            txc(1:isize_field, 1) = s(1:isize_field, 1)*s(1:isize_field, 2)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)*s(1:isize_field, 1)
            txc(1:isize_field, 3) = txc(1:isize_field, 1)*s(1:isize_field, 2)

            ifield = 0
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 's1s2'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 's1s2s1'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 's1s2s2'

            ! ###################################################################
            ! Potential vorticity
            ! ###################################################################
        case (17)
            write (fname, *) itime; fname = 'avgPV'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 6) = txc(1:isize_field, 1)*txc(1:isize_field, 1) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 2) &
                                    + txc(1:isize_field, 3)*txc(1:isize_field, 3) ! Enstrophy
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 3)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b

            txc(1:isize_field, 5) = sqrt(txc(1:isize_field, 5) + small_wp)
            txc(1:isize_field, 6) = sqrt(txc(1:isize_field, 6) + small_wp)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)/(txc(1:isize_field, 5)*txc(1:isize_field, 6)) ! Cosine of angle between 2 vectors

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'PV'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Cos'

        case (18)
            call AvgPhaseSpace(wrk2d, inb_flow, it, 0, 0, 1)
            call IO_Write_AvgPhase(1, inb_flow, IO_FLOW, 0, PhAvg%stride, avgu_name, 1, avg_flow, itime_vec(it))

            call AvgPhaseSpace(wrk2d, inb_scal, it, 0, 0, 2)
            call IO_Write_AvgPhase(1, inb_scal, IO_SCAL, 0, PhAvg%stride, avgp_name, 2, avg_scal, itime_vec(it))

            p => txc(:, 9) !makes sure to only pass the address, not the entire array
            call AvgPhaseSpace(wrk2d, 1, it, 0, 0, p)
            call IO_Write_AvgPhase(1, 1, IO_SCAL, 0, PhAvg%stride, avgs_name, 4, avg_p, itime_vec(it))

            call AvgPhaseStress(q, it, 0, 0)

            call AvgPhaseResetVariable()

        end select

        if (opt_main > 2) then
            if (nfield < ifield) then ! Check
                call TLab_Write_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
                call TLab_Stop(DNS_ERROR_WRKSIZE)
            end if

            if (jmax_aux*opt_block /= g(2)%size) then
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
            end if

            call AVG_N_XZ(fname, itime, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_order, vars, gate_level, gate, y_aux, mean)

        end if

    end do

    call TLab_Stop(0)
end program AVERAGES

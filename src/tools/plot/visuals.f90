#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "VISUALS"

!########################################################################
!#
!# Creating data blocks for visualizations. Derived from ensight.f
!# Partition to be incorporated via a MASK routine before VISUALS_WRITE
!#
!########################################################################
program VISUALS

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_pro
    use TLAB_MPI_PROCS
#endif
    use THERMO_VARS, only: imixture
    use THERMO_VARS, only: NSP, THERMO_SPNAME
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS
    use IBM_VARS
    use IO_FIELDS
    use OPR_FOURIER
    use OPR_PARTIAL

    implicit none

#include "integers.h"

    ! Parameter definitions
    TINTEGER, parameter :: itime_size_max = 3000
    TINTEGER, parameter :: iopt_size_max = 20
    TINTEGER, parameter :: igate_size_max = 8
    TINTEGER, parameter :: params_size_max = 2

    ! Additional local arrays
    integer(1), dimension(:), allocatable, save :: gate

    ! -------------------------------------------------------------------
    ! Local variables
    ! -------------------------------------------------------------------
    character*512 sRes
    character*32 fname, bakfile
    character*32 flow_file, scal_file, part_file, plot_file, time_str
    character*64 str

    TINTEGER opt_format
    TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
    TINTEGER ij, is, bcs(2, 2)
    TINTEGER iscal_offset, iread_flow, iread_scal, iread_part, idummy, MaskSize
    TREAL diff, dummy
    TINTEGER subdomain(6)

    ! Gates for the definition of the intermittency function (partition of the fields)
    TINTEGER igate_size
    TREAL gate_threshold(igate_size_max)

    TINTEGER itime_size, it
    TINTEGER itime_vec(itime_size_max)

    TINTEGER iopt_size, iv
    TINTEGER opt_vec(iopt_size_max)
    TREAL opt_vec2(iopt_size_max)

    TINTEGER params_size
    TREAL params(params_size_max)

    !########################################################################
    !########################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call PARTICLE_READ_GLOBAL(ifile)

    ! -------------------------------------------------------------------
    ! IBM status (before TLAB_MPI_INITIALIZE!)
    ! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, ifile, 'IBMParameter', 'Status', 'off', sRes)
    if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
    else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
    else
        call TLAB_WRITE_ASCII(efile, 'VISUALS. Wrong IBM Status option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    ! -------------------------------------------------------------------
    ! Initialize MPI
    ! -------------------------------------------------------------------
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    ! -------------------------------------------------------------------
    ! File names
    ! -------------------------------------------------------------------
#include "dns_read_times.h"

    ! -------------------------------------------------------------------
    ! Read local options
    ! -------------------------------------------------------------------
    opt_format = 2 ! default values

    if (imixture == MIXT_TYPE_NONE) then; iscal_offset = 9    ! define iscal_offset
    else; iscal_offset = 9 + NSP
    end if

    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamVisuals', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_INTEGER(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, '(A)') 'Option?'
        write (*, '(A)') ' 0. Grid'
        write (*, '(A)') ' 1. VelocityX'
        write (*, '(A)') ' 2. VelocityY'
        write (*, '(A)') ' 3. VelocityZ'
        write (*, '(A)') ' 4. VelocityVector'
        write (*, '(A)') ' 5. Velocity V_iV_i'
        write (*, '(A)') ' 6. Density'
        write (*, '(A)') ' 7. Temperature'
        write (*, '(A)') ' 8. Pressure'
        write (*, '(A)') ' 9. Scalars'
        if (imixture /= MIXT_TYPE_NONE) then
            do is = 1, NSP
                write (*, '(I2,A)') 9 + is, '. '//THERMO_SPNAME(is)
            end do
        end if

        write (*, '(I2,A)') iscal_offset + 1, '. ScalarGradientVector'
        write (*, '(I2,A)') iscal_offset + 2, '. ScalarGradient G_iG_i (Ln)'
        write (*, '(I2,A)') iscal_offset + 3, '. ScalarGradientEquation'
        write (*, '(I2,A)') iscal_offset + 4, '. VorticityVector'
        write (*, '(I2,A)') iscal_offset + 5, '. Enstrophy W_iW_i (Ln)'
        write (*, '(I2,A)') iscal_offset + 6, '. EnstrophyEquation'
        write (*, '(I2,A)') iscal_offset + 7, '. StrainTensor'
        write (*, '(I2,A)') iscal_offset + 8, '. Strain 2S_ijS_ij (Ln)'
        write (*, '(I2,A)') iscal_offset + 9, '. StrainEquation'
        write (*, '(I2,A)') iscal_offset + 10, '. VelocityGradientInvariants'
        write (*, '(I2,A)') iscal_offset + 11, '. Space partition'
        write (*, '(I2,A)') iscal_offset + 12, '. Buoyancy'
        write (*, '(I2,A)') iscal_offset + 14, '. HorizontalDivergence'
        write (*, '(I2,A)') iscal_offset + 15, '. Turbulent quantities'
        write (*, '(I2,A)') iscal_offset + 16, '. Radiative forcing'
        write (*, '(I2,A)') iscal_offset + 17, '. Relative humidity'
        write (*, '(I2,A)') iscal_offset + 18, '. Particle Density'
        write (*, '(I2,A)') iscal_offset + 19, '. Thermodynamic quantities'
        write (*, '(I2,A)') iscal_offset + 20, '. Analysis of B and V'
        read (*, '(A512)') sRes
#endif
    end if
    iopt_size = iopt_size_max
    call LIST_INTEGER(sRes, iopt_size, opt_vec)

    if (opt_vec(1) < 0) then ! Check
        call TLAB_WRITE_ASCII(efile, 'VISUALS. Missing input [PostProcessing.ParamVisuals] in dns.ini.')
        call TLAB_STOP(DNS_ERROR_INVALOPT)
    end if

    ! -------------------------------------------------------------------
    iread_flow = 0
    iread_scal = 0
    iread_part = 0
    inb_txc = 0

    do iv = 1, iopt_size
        if (opt_vec(iv) == 1) then; iread_flow = 1; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 2) then; iread_flow = 1; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 3) then; iread_flow = 1; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 4) then; iread_flow = 1; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == 5) then; iread_flow = 1; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 6) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == 7) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == 8) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 7); end if
        if (opt_vec(iv) == 9) then; iread_scal = 1; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) > 9 .and. opt_vec(iv) <= iscal_offset) then
            iread_scal = 1; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 1) then; iread_scal = 1; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 2) then; iread_scal = 1; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 3) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 4) then; iread_flow = 1; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 5) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 7); end if
        if (opt_vec(iv) == iscal_offset + 6) then; iread_flow = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 7) then; iread_flow = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 8) then; iread_flow = 1; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 9) then; iread_flow = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 10) then; iread_flow = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 12) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 14) then; iread_flow = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 15) then; iread_flow = 1; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 16) then; iread_scal = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 17) then; iread_scal = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 18) then; iread_part = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 19) then; iread_scal = 1; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 20) then; iread_flow = 1; iread_scal = 1; inb_txc = max(inb_txc, 7); end if
    end do

    ! check if enough memory is provided for the IBM
    if (imode_ibm == 1) then
#ifdef IBM_DEBUG
        inb_txc = max(inb_txc, 6)
#else
        inb_txc = max(inb_txc, 3)
#endif
    end if

    ! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Subdomain', '-1', sRes)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Subdomain limits ?'
        read (*, '(A64)') sRes
#endif
    end if
    idummy = 6
    call LIST_INTEGER(sRes, idummy, subdomain)

    if (idummy < 6) then ! default
        subdomain(1) = 1; subdomain(2) = g(1)%size
        subdomain(3) = 1; subdomain(4) = g(2)%size
        subdomain(5) = 1; subdomain(6) = g(3)%size
    end if

    ! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Format', '-1', sRes)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'File Format ?'
        write (*, *) ' 0. General restart format'
        write (*, *) ' 1. Ensight'
        write (*, *) ' 2. Raw, single precision, no header'
        read (*, '(A64)') sRes
#endif
    end if
    if (len_trim(adjustl(sRes)) > 0) then
        if (trim(adjustl(sRes)) == 'general') then; opt_format = 0
        else if (trim(adjustl(sRes)) == 'ensight') then; opt_format = 1
        else if (trim(adjustl(sRes)) == 'single') then; opt_format = 2
        else
            read (sRes, *) opt_format
        end if
    end if

    if (opt_format < 0) opt_format = 2 ! default is single precission, no header

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

    do iv = 1, iopt_size
        if (opt_vec(iv) == iscal_offset + 11) then
#include "dns_read_partition.h"
            if (opt_cond > 1) inb_txc = max(inb_txc, 5)
            exit
        end if
    end do

    ! -------------------------------------------------------------------
    ! Allocating memory space
    ! -------------------------------------------------------------------
    allocate (gate(isize_field))

    isize_wrk3d = isize_txc_field
#ifdef USE_MPI
    isize_wrk3d = isize_wrk3d + isize_field ! more space in wrk3d array needed in IO_WRITE_VISUALS
#endif
    if (part%type /= PART_TYPE_NONE) then
        isize_wrk3d = max(isize_wrk3d, (imax + 1)*jmax*(kmax + 1))
    end if

    call TLAB_ALLOCATE(C_FILE_LOC)

    if (iread_part == 1) then ! Particle variables
        inb_part_txc = max(inb_part_txc, 1)
        call PARTICLE_ALLOCATE(C_FILE_LOC)
    end if

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(C_FILE_LOC)
    end if

    ! -------------------------------------------------------------------
    ! Initialize
    ! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE(wrk1d) ! Initialize thermodynamic quantities

    if (ifourier == 1 .and. inb_txc >= 1) then ! For Poisson solver
        call OPR_FOURIER_INITIALIZE(txc, wrk1d, wrk2d, wrk3d)
    end if

    if (iread_flow == 1 .and. inb_txc >= 3) then ! We need array space
        call OPR_CHECK(imax, jmax, kmax, q, txc, wrk2d, wrk3d)
    end if

    if (imode_ibm == 1) then
        call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
    end if

#ifdef USE_MPI
    call VISUALS_MPIO_AUX(opt_format, subdomain)
#else
    io_aux(:)%offset = 0
#endif

    MaskSize = 6

    ! ###################################################################
    ! Grid
    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        do iv = 1, iopt_size
            if (opt_vec(iv) == 0) then
                call ENSIGHT_GRID('grid.ensight', g(1)%size, g(2)%size, g(3)%size, subdomain, g(1)%nodes, g(2)%nodes, g(3)%nodes)
            end if
        end do
#ifdef USE_MPI
    end if
#endif

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))//'.'
        call TLAB_WRITE_ASCII(lfile, sRes)

        if (icalc_scal == 1 .and. iread_scal == 1) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_READ_FIELDS(scal_file, IO_SCAL, imax, jmax, kmax, inb_scal, i0, s, wrk3d)
        elseif (icalc_scal == 0) then
            s = C_0_R
        end if

        if (iread_flow == 1) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_READ_FIELDS(flow_file, IO_FLOW, imax, jmax, kmax, inb_flow, i0, q, wrk3d)
        end if

        if (imode_ibm == 1) then
            call IBM_BCS_FIELD_COMBINED(i0, q)
            if (icalc_scal == 1) call IBM_INITIALIZE_SCAL(i0, s)
        end if

        call FI_DIAGNOSTIC(imax, jmax, kmax, q, s, wrk3d)

        if (iread_part == 1) then ! Particle variables
            write (part_file, *) itime; part_file = trim(adjustl(tag_part))//trim(adjustl(part_file))
        end if

        write (sRes, 100) rtime; sRes = 'Physical time '//trim(adjustl(sRes))
        call TLAB_WRITE_ASCII(lfile, sRes)

        ! -------------------------------------------------------------------
        ! Calculate intermittency
        ! -------------------------------------------------------------------
        if (opt_cond == 1) then ! Read external file
            write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname)); params_size = 2
            call IO_READ_FIELD_INT1(fname, i1, imax, jmax, kmax, itime, params_size, params, gate)
            igate_size = int(params(2))

        else if (opt_cond > 1) then
            opt_cond_scal = 1 ! Scalar field to use for the conditioning
            if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                opt_cond_scal = inb_scal_array
            end if

            call TLAB_WRITE_ASCII(lfile, 'Calculating partition...')
            call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                         imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate, wrk2d, wrk3d)
        end if

        ! -------------------------------------------------------------------
        ! define time string
        ! -------------------------------------------------------------------
        do ij = MaskSize, 1, -1
            time_str(ij:ij) = '0'
        end do
        write (plot_file, '(I10)') itime
        time_str(MaskSize - len_trim(adjustl(plot_file)) + 1:Masksize) = trim(adjustl(plot_file))

        ! -------------------------------------------------------------------
        ! Loop over options
        ! -------------------------------------------------------------------
        do iv = 1, iopt_size

            ! ###################################################################
            ! Velocities
            ! ###################################################################
            if (opt_vec(iv) == 1) then
                txc(1:isize_field, 1) = q(1:isize_field, 1)
                plot_file = 'VelocityX'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            else if (opt_vec(iv) == 2) then
                txc(1:isize_field, 1) = q(1:isize_field, 2)
                plot_file = 'VelocityY'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            else if (opt_vec(iv) == 3) then
                txc(1:isize_field, 1) = q(1:isize_field, 3)
                plot_file = 'VelocityZ'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            else if (opt_vec(iv) == 4) then
                txc(1:isize_field, 1:3) = q(1:isize_field, 1:3)
                plot_file = 'VelocityVector'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i3, subdomain, txc(1, 1), wrk3d)

            else if (opt_vec(iv) == 5) then
                txc(1:isize_field, 1) = sqrt(q(1:isize_field, 1)**2 + q(1:isize_field, 2)**2 + q(1:isize_field, 3)**2)
                plot_file = 'VelocityMagnitude'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            end if

            ! ###################################################################
            ! Thermodynamic state
            ! ###################################################################
            if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
                if (opt_vec(iv) == 6) then ! density
                    plot_file = 'Density'//time_str(1:MaskSize)
                    if (buoyancy%type == EQNS_EXPLICIT) then
                        call THERMO_ANELASTIC_DENSITY(imax, jmax, kmax, s, epbackground, pbackground, txc(1, 1))
                    else
                        wrk1d(1:jmax, 1) = C_0_R
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, txc(1, 1), wrk1d)
                        dummy = C_1_R/froude
                        txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy + C_1_R
                    end if
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                else if (opt_vec(iv) == 7 .and. imixture == MIXT_TYPE_AIRWATER) then ! temperature
                    plot_file = 'Temperature'//time_str(1:MaskSize)
                    call THERMO_ANELASTIC_TEMPERATURE(imax, jmax, kmax, s, epbackground, txc(1, 1))
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                    if (damkohler(1) > C_0_R) then ! Supersaturated liquid; this is wrong
                        plot_file = 'Supsat'//time_str(1:MaskSize)
                        txc(1:isize_field, 1:2) = s(1:isize_field, 1:2)
                        call THERMO_AIRWATER_PH(imax, jmax, kmax, txc(1, 2), txc(1, 1), epbackground, pbackground)
                        txc(1:isize_field, 3) = (s(1:isize_field, 3) - txc(1:isize_field, 3))/s(1, 3)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 3), wrk3d)
                    end if

                else if (opt_vec(iv) == 8) then ! pressure
                    plot_file = 'Pressure'//time_str(1:MaskSize)
                    call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk1d, wrk2d, wrk3d)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                    plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
                    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc(1, 1), txc(1, 2), wrk3d, wrk2d, wrk3d)
                    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 1), txc(1, 3), wrk3d, wrk2d, wrk3d)
                    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc(1, 1), txc(1, 4), wrk3d, wrk2d, wrk3d)
                    txc(1:isize_field, 2) = -(txc(1:isize_field, 2)*q(1:isize_field, 1) &
                                              + txc(1:isize_field, 3)*q(1:isize_field, 2) &
                                              + txc(1:isize_field, 4)*q(1:isize_field, 3))
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                    call TLAB_WRITE_ASCII(lfile, 'Computing pressure-strain correlation...')
                    txc(1:isize_field, 2) = txc(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 2))

                    plot_file = 'PressureStrainX'//time_str(1:MaskSize)
                    txc(1:isize_field, 3) = q(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc(1, 3), txc(1, 4), wrk3d, wrk2d, wrk3d)
                    txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 3), wrk3d)

                    plot_file = 'PressureStrainY'//time_str(1:MaskSize)
                    txc(1:isize_field, 3) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 3), txc(1, 4), wrk3d, wrk2d, wrk3d)
                    txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 3), wrk3d)

                    plot_file = 'PressureStrainZ'//time_str(1:MaskSize)
                    txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc(1, 3), txc(1, 4), wrk3d, wrk2d, wrk3d)
                    txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 3), wrk3d)

                    plot_file = 'PressureHydrostatic'//time_str(1:MaskSize)
                    q = C_0_R
                    call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk1d, wrk2d, wrk3d)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                    plot_file = 'PressureHydrodynamic'//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = txc(1:isize_field, 1) - txc(1:isize_field, 2)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                end if

            else ! compressible
                if (opt_vec(iv) == 6) then ! density
                    plot_file = 'Density'//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = q(1:isize_field, 5)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                else if (opt_vec(iv) == 7) then ! temperature
                    plot_file = 'Temperature'//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = q(1:isize_field, 7)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                else if (opt_vec(iv) == 8) then ! pressure
                    plot_file = 'Pressure'//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = q(1:isize_field, 6)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                end if

            end if

            ! ###################################################################
            ! Scalars
            ! ###################################################################
            if (opt_vec(iv) == 9) then ! All prognostic scalars
                do is = 1, inb_scal_array
                    write (str, *) is; plot_file = 'Scalar'//trim(adjustl(str))//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = s(1:isize_field, is)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
                end do

            else if (opt_vec(iv) > 9 .and. opt_vec(iv) <= iscal_offset) then ! Individual and diagnostic scalars
                if (imixture == MIXT_TYPE_AIRWATER) then ! s(1,inb_scal+1) contains liquid mass fraction
                    if (opt_vec(iv) == 10) then ! vapor water mass fraction
                        plot_file = trim(adjustl(THERMO_SPNAME(1)))//time_str(1:MaskSize)
                        s(:, 1) = s(:, inb_scal) - s(:, inb_scal + 1)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s, wrk3d)

                    else if (opt_vec(iv) == 11) then ! air mass fraction
                        plot_file = trim(adjustl(THERMO_SPNAME(2)))//time_str(1:MaskSize)
                        s(:, 1) = C_1_R - s(:, inb_scal)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s, wrk3d)

                    else if (opt_vec(iv) == 12) then ! liquid mass fraction
                        plot_file = trim(adjustl(THERMO_SPNAME(3)))//time_str(1:MaskSize)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s(1, inb_scal + 1), wrk3d)

                    end if

                else ! Plot the chosen species
                    is = opt_vec(iv) - 9; plot_file = trim(adjustl(THERMO_SPNAME(is)))//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = s(1:isize_field, is)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                end if

            end if

            ! ###################################################################
            ! Scalar Derivatives
            ! ###################################################################
            if (opt_vec(iv) >= iscal_offset + 1 .and. opt_vec(iv) <= iscal_offset + 3) then
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    if (opt_vec(iv) == iscal_offset + 1) then
                        plot_file = trim(adjustl(str))//'GradientVector'//time_str(1:MaskSize)
                        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, is), txc(1, 1), wrk3d, wrk2d, wrk3d)
                        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, is), txc(1, 2), wrk3d, wrk2d, wrk3d)
                        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, is), txc(1, 3), wrk3d, wrk2d, wrk3d)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i3, subdomain, txc(1, 1), wrk3d)
                    end if

                    if (opt_vec(iv) == iscal_offset + 2 .or. opt_vec(iv) == iscal_offset + 3) then ! Gradient magnitude
                        plot_file = trim(adjustl(str))//'Gradient'//time_str(1:MaskSize)
                        call FI_GRADIENT(imax, jmax, kmax, s(1, is), txc(1, 1), txc(1, 2), wrk2d, wrk3d)
                        if (opt_vec(iv) == iscal_offset + 2) then
                            plot_file = 'Ln'//trim(adjustl(plot_file))
                            txc(1:isize_field, 1) = log(txc(1:isize_field, 1) + C_SMALL_R)
                        end if
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
                    end if

                    if (opt_vec(iv) == iscal_offset + 3 .and. is <= inb_scal) then ! Scalar gradient equation
                        if (idiffusion == EQNS_NONE) then; diff = C_0_R
                        else; diff = visc/schmidt(is)
                        end if

                        call TLAB_WRITE_ASCII(lfile, 'Computing scalar gradient production...')
                        plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
                        call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s(1, is), q(1, 1), q(1, 2), q(1, 3), &
                                                    txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                        call TLAB_WRITE_ASCII(lfile, 'Computing scalar gradient diffusion...')
                        plot_file = trim(adjustl(str))//'GradientDiffusion'//time_str(1:MaskSize)
                        call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s(1, is), &
                                                   txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                        txc(1:isize_field, 1) = diff*txc(1:isize_field, 1)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                    end if

                end do

            end if

            ! ###################################################################
            ! Velocity Derivatives
            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 4) then ! VorticityVector
                plot_file = 'VorticityVector'//time_str(1:MaskSize)
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i3, subdomain, txc(1, 1), wrk3d)
            end if

            if (opt_vec(iv) == iscal_offset + 5 .or. opt_vec(iv) == iscal_offset + 6) then ! Enstrophy
                plot_file = 'Enstrophy'//time_str(1:MaskSize)
                call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, wrk3d)
                if (opt_vec(iv) == iscal_offset + 5) then ! Natural log
                    plot_file = 'Ln'//trim(adjustl(plot_file))
                    txc(1:isize_field, 1) = log(txc(1:isize_field, 1) + C_SMALL_R)
                end if
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'LnPotentialEnstrophy'//time_str(1:MaskSize)
                if (buoyancy%type == EQNS_EXPLICIT) then
                    call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, txc(1, 4))
                else
                    wrk1d(1:jmax, 1) = C_0_R
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, txc(1, 4), wrk1d)
                end if
                dummy = C_1_R/froude
                txc(1:isize_field, 4) = txc(1:isize_field, 4)*dummy
                call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc(1, 4), txc(1, 1), wrk3d, wrk2d, wrk3d)
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 4), txc(1, 2), wrk3d, wrk2d, wrk3d)
                call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc(1, 4), txc(1, 3), wrk3d, wrk2d, wrk3d)
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), wrk2d, wrk3d)
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                                        + txc(1:isize_field, 2)*txc(1:isize_field, 5) &
                                        + txc(1:isize_field, 3)*txc(1:isize_field, 6)
                txc(1:isize_field, 1) = log(txc(1:isize_field, 1)*txc(1:isize_field, 1) + C_SMALL_R)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
            end if

            if (opt_vec(iv) == iscal_offset + 6) then ! EnstrophyEquation
                call TLAB_WRITE_ASCII(lfile, 'Computing enstrophy production...')
                plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
                call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                             txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                call TLAB_WRITE_ASCII(lfile, 'Computing enstrophy diffusion...')
                plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
                call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                            txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                txc(1:isize_field, 1) = visc*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            end if

            ! -------------------------------------------------------------------
            if (opt_vec(iv) == iscal_offset + 7) then ! Strain Tensor
                plot_file = 'StrainTensor'//time_str(1:MaskSize)
     CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i6, subdomain, txc(1, 1), wrk3d)
            end if

            if (opt_vec(iv) == iscal_offset + 8 .or. opt_vec(iv) == iscal_offset + 9) then ! Strain
                plot_file = 'Strain'//time_str(1:MaskSize)
                call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, wrk3d)
                txc(1:isize_field, 1) = C_2_R*txc(1:isize_field, 1)
                if (opt_vec(iv) == iscal_offset + 8) then ! Natural log
                    plot_file = 'Ln'//trim(adjustl(plot_file))
                    txc(1:isize_field, 1) = log(txc(1:isize_field, 1) + C_SMALL_R)
                end if
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
            end if

            if (opt_vec(iv) == iscal_offset + 9) then ! StrainEquation (I need the pressure)
                call TLAB_WRITE_ASCII(lfile, 'Computing strain pressure...')
                plot_file = 'StrainPressure'//time_str(1:MaskSize)
                if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
                    call TLAB_WRITE_ASCII(efile, 'VISUALS. Strain eqn for incompressible undeveloped.')
                    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
                else
                    txc(:, 6) = q(:, 6)
                end if
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 6), &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
                txc(1:isize_field, 1) = C_2_R*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                call TLAB_WRITE_ASCII(lfile, 'Computing strain production...')
                plot_file = 'StrainProduction'//time_str(1:MaskSize)
                call FI_STRAIN_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                          txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                txc(1:isize_field, 1) = C_2_R*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                call TLAB_WRITE_ASCII(lfile, 'Computing strain diffusion...')
                plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
                call FI_STRAIN_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                         txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                txc(1:isize_field, 1) = C_2_R*visc*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            end if

            ! -------------------------------------------------------------------
            if (opt_vec(iv) == iscal_offset + 10) then ! Invariants
                plot_file = 'InvariantP'//time_str(1:MaskSize)
                call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'InvariantQ'//time_str(1:MaskSize)
          call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'InvariantR'//time_str(1:MaskSize)
                call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                    txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            end if

            ! ###################################################################
            ! Partition
            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 11) then
                call TLAB_WRITE_ASCII(efile, 'VISUALS. Partition undevelop.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            ! ###################################################################
            ! Buoyancy
            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 12) then
                plot_file = 'Buoyancy'//time_str(1:MaskSize)
                if (buoyancy%type == EQNS_EXPLICIT) then
                    call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, txc(1, 1))
                else
                    wrk1d(1:jmax, 1) = C_0_R
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, txc(1, 1), wrk1d)
                end if
                dummy = C_1_R/froude
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'Fvb'//time_str(1:MaskSize)     ! buoyancy flux along Oy
                txc(1:isize_field, 2) = txc(1:isize_field, 1)*q(1:isize_field, 2)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                plot_file = 'bPrime'//time_str(1:MaskSize)  ! buoyancy fluctuation
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 1))
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'Cvb'//time_str(1:MaskSize)     ! Covariance between b and v
                txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 2) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                plot_file = 'LnBuoyancySource'//time_str(1:MaskSize)
                if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                    call THERMO_AIRWATER_LINEAR_SOURCE(imax, jmax, kmax, s, txc(1, 1), txc(1, 2), txc(1, 3))
                    call FI_GRADIENT(imax, jmax, kmax, txc(1, 1), txc(1, 2), txc(1, 4), wrk2d, wrk3d)
                    dummy = buoyancy%parameters(inb_scal_array)
                    txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 3)*dummy
                else
                    call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 1), txc(1, 2), wrk2d, wrk3d)
                    call FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, txc(1, 1), txc(1, 2))
                end if
                dummy = visc/schmidt(1)/froude
                txc(1:isize_field, 1) = txc(1:isize_field, 2)*dummy
                txc(1:isize_field, 1) = log(abs(txc(1:isize_field, 1)) + C_SMALL_R)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 14) then
                plot_file = 'HorizontalDivergence'//time_str(1:MaskSize)
                call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(1, 1), txc(1, 2), wrk3d, wrk2d, wrk3d)
                call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(1, 3), txc(1, 1), wrk3d, wrk2d, wrk3d)
                txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 15) then ! Turbulent quantities
                plot_file = 'LnDissipation'//time_str(1:MaskSize)
                call FI_DISSIPATION(i1, imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                    txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk1d, wrk2d, wrk3d)
                txc(1:isize_field, 1) = log(txc(1:isize_field, 1) + C_SMALL_R)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'Tke'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = q(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 1))
                txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                txc(1:isize_field, 4) = C_05_R*(txc(1:isize_field, 1)**2 + txc(1:isize_field, 2)**2 + txc(1:isize_field, 3)**2)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 4), wrk3d)

                plot_file = 'ReynoldsTensor'//time_str(1:MaskSize)
                txc(1:isize_field, 4) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                txc(1:isize_field, 5) = txc(1:isize_field, 1)*txc(1:isize_field, 3)
                txc(1:isize_field, 6) = txc(1:isize_field, 2)*txc(1:isize_field, 3)
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 1)
                txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 2)
                txc(1:isize_field, 3) = txc(1:isize_field, 3)*txc(1:isize_field, 3)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i6, subdomain, txc(1, 1), wrk3d)

            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 16) then
                do is = 1, inb_scal

                    if (radiation%active(is)) then
                        write (str, *) is; plot_file = 'Radiation'//trim(adjustl(str))//time_str(1:MaskSize)
                        if (imode_eqns == DNS_EQNS_ANELASTIC) then
                         call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, s(1, radiation%scalar(is)), txc(1, 2))
                            call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), txc(1, 2), txc(1, 1), wrk1d, wrk3d)
                            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, ribackground, txc(1, 1))
                        else
                           call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(1, radiation%scalar(1)), txc(1, 1), wrk1d, wrk3d)
                        end if
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
                    end if

                end do
            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 17) then
                plot_file = 'RelativeHumidity'//time_str(1:MaskSize)
                call THERMO_ANELASTIC_RELATIVEHUMIDITY(imax, jmax, kmax, s, epbackground, pbackground, wrk3d, txc(1, 1))
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 18) then ! Particle information
                plot_file = 'ParticleDensity'//time_str(1:MaskSize)
                call IO_READ_PARTICLE(part_file, l_g, l_q)
                l_txc = C_1_R; ! We want density
                call PARTICLE_TO_FIELD(l_q, l_txc, txc(1, 1), wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then
                    txc(:, 1) = txc(:, 1) + 0.00000001
                    do is = 1, 2
                        plot_file = trim(adjustl(part_spname(is)))//time_str(1:MaskSize)
                        call PARTICLE_TO_FIELD(l_q, l_q(1, 3 + is), txc(1, 2), wrk3d)
                        txc(:, 2) = txc(:, 2)/txc(:, 1)
                        call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)
                    end do
                end if

                if (part%type == PART_TYPE_BIL_CLOUD_4) then
                    !inb_part_array is the last component -> residence times in bil_cloud_4
                    plot_file = trim(adjustl(part_spname(3)))//time_str(1:MaskSize)
                    call PARTICLE_TO_FIELD(l_q, l_q(1, inb_part_array - 1), txc(1, 2), wrk3d)
                    call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)
                end if

            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 19) then ! Thermodynamic quantities
                call VISUALS_FUNCTION1(imax, jmax, kmax, s, txc)

                plot_file = 'Enthalpy'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s(1, 1), wrk3d)
                plot_file = 'TotalWater'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s(1, 2), wrk3d)
                plot_file = 'LiquidWater'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, s(1, 3), wrk3d)
                plot_file = 'Temperature'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)
                plot_file = 'Density'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 20) then
                plot_file = 'LaplacianV'//time_str(1:MaskSize)
                call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), q(1, 2), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
                call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), q(1, 2), txc(1, 3), txc(1, 5), wrk2d, wrk3d)
                call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), q(1, 2), txc(1, 2), txc(1, 5), wrk2d, wrk3d)
                txc(1:isize_field, 2) = txc(1:isize_field, 2) + txc(1:isize_field, 3) + txc(1:isize_field, 4)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                plot_file = 'Buoyancy'//time_str(1:MaskSize)
                if (buoyancy%type == EQNS_EXPLICIT) then
                    call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, txc(1, 1))
                else
                    wrk1d(1:jmax, 1) = C_0_R
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, txc(1, 1), wrk1d)
                end if
                dummy = C_1_R/froude
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'LaplacianB'//time_str(1:MaskSize)
                call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), txc(1, 1), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
                call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), txc(1, 1), txc(1, 3), txc(1, 5), wrk2d, wrk3d)
                call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), txc(1, 1), txc(1, 2), txc(1, 5), wrk2d, wrk3d)
                txc(1:isize_field, 2) = txc(1:isize_field, 2) + txc(1:isize_field, 3) + txc(1:isize_field, 4)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)

                plot_file = 'Pressure'//time_str(1:MaskSize)
                bbackground = C_0_R
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk1d, wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 1), wrk3d)

                plot_file = 'PressureGradientY'//time_str(1:MaskSize)
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 1), txc(1, 2), wrk3d, wrk2d, wrk3d)
                call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, 2), wrk3d)
            end if

        end do

    end do

100 format(G_FORMAT_R)
    call TLAB_STOP(0)
end program VISUALS

!########################################################################
!# DESCRIPTION
!#
!# Calculate thermodynamic information for THERMO_AIRWATER_LINEAR
!#
!########################################################################
subroutine VISUALS_FUNCTION1(nx, ny, nz, s, txc)

    use TLAB_VARS, only: isize_txc_field
    use TLAB_VARS, only: epbackground, pbackground
    use THERMO_VARS, only: imixture

    implicit none

    TINTEGER nx, ny, nz
    TREAL, dimension(nx*ny*nz, *) :: s
    TREAL, dimension(isize_txc_field, *) :: txc

    ! -----------------------------------------------------------------------
    TREAL qt_0, qt_1, h_0, h_1, p, T_0, C_0, PsiRef

    ! #######################################################################
    imixture = MIXT_TYPE_AIRWATER

    qt_0 = 9.0d-3; qt_1 = 1.5d-3
    h_0 = 0.955376d0; h_1 = 0.981965d0
    p = 0.940d0
    T_0 = 0.952181d0 ! 283.75 / TREF
    C_0 = 1.0089
    PsiRef = 6.57d-4

    s(:, 3) = h_0 + s(:, 1)*(h_1 - h_0) + s(:, 2)*C_0*T_0*PsiRef ! enthalpy
    s(:, 2) = qt_0 + s(:, 1)*(qt_1 - qt_0)                            ! total water, space for q_l
    s(:, 1) = s(:, 3)

    epbackground = C_0_R                                    ! potential energy
    pbackground = p                                        ! pressure

    call THERMO_AIRWATER_PH(nx, ny, nz, s(1, 2), s(1, 1), epbackground, pbackground)
    call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s(1, 1), epbackground, txc(1, 1))
    call THERMO_ANELASTIC_DENSITY(nx, ny, nz, s(1, 1), epbackground, pbackground, txc(1, 2))

    return
end subroutine VISUALS_FUNCTION1

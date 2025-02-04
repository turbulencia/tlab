#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# Implementation of relaxation terms in buffer regions of the computational domain
!#
!# Data is read in dns_read_local
!#
!# Note that with the energy formulation the forcing terms have been
!# made proportional to differences on the conserved variables, i.e. rho, rho*V, rho*(e+v^2/2).
!#
!# The buffer files need to written without header information, so that
!# the header global variables (like itime) are not updated within this routine
!#
!########################################################################
module BOUNDARY_BUFFER
    use TLab_Constants, only: tag_flow, tag_scal, wfile, efile, lfile, MAX_VARS, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
#endif
    use NavierStokes, only: nse_eqns
    use TLab_WorkFlow, only: imode_sim
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal, isize_field
    use FDM, only: g
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermodynamics, only: CRATIO_INV
    use IO_FIELDS
    use OPR_FILTERS
    use Averages, only: COV2V1D, COV2V2D

#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS
    use TLabMPI_PROCS, only: TLabMPI_Panic
    use TLabMPI_Transpose, only: TLabMPI_Trp_TypeK_Create, ims_trp_plan_k
#endif

    implicit none
    save
    private

    integer(wi), parameter :: FORM_POWER_MIN = 1
    integer(wi), parameter :: FORM_POWER_MAX = 2

    type buffer_dt
        sequence
        integer type                                  ! relaxation, filter...
        integer(wi) size, total_size                      ! # points in buffer layer on this task and in total
        integer(wi) offset                                ! position in absolute grid
        integer(wi) nfields                               ! number of fields to apply tosition in absolute grid
        logical active(MAX_VARS), hard
        real(wp) strength(MAX_VARS), sigma(MAX_VARS)      ! parameters in relaxation term
        integer(wi) form                                  ! form of function of relaxation term
        real(wp) hardvalues(MAX_VARS)                     ! Fixed reference values
        real(wp), allocatable, dimension(:, :) :: tau  ! relaxation timescale for each field
        real(wp), allocatable, dimension(:, :, :, :) :: ref  ! reference fields
        integer(wi) id_subarray
    end type buffer_dt

    integer(wi), public :: BuffType
    logical, public :: BuffLoad
    type(buffer_dt), public :: BuffFlowImin, BuffFlowImax, BuffFlowJmin, BuffFlowJmax
    type(buffer_dt), public :: BuffScalImin, BuffScalImax, BuffScalJmin, BuffScalJmax
    !  TYPE(filter_dt), DIMENSION(3) :: FilterBuffer
    ! BufferFilter should then be a block in tlab.ini as [Filter], which is read in TLab_Initialize_Parameters.

    public :: BOUNDARY_BUFFER_READBLOCK
    public :: BOUNDARY_BUFFER_INITIALIZE
    public :: BOUNDARY_BUFFER_RELAX_FLOW
    public :: BOUNDARY_BUFFER_RELAX_SCAL
    public :: BOUNDARY_BUFFER_RELAX_SCAL_I
    public :: BOUNDARY_BUFFER_FILTER

    integer(wi) j, jloc, i, iloc, iq, is, idummy
    real(wp) dummy

contains

    ! ###################################################################
    ! ###################################################################
    subroutine BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, tag, variable, nfields)
        character(len=*) bakfile, inifile, tag
        type(buffer_dt) variable
        integer(wi) nfields

        real(wp) dummies(nfields + 1)
        character(len=512) sRes
        character(len=4) str

        ! ###################################################################
        variable%active(:) = .false.; variable%hard = .false.
        if (variable%size > 0) then
            call ScanFile_Char(bakfile, inifile, 'BufferZone', 'Parameters'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'void') then
                str = trim(adjustl(tag))
                call ScanFile_Char(bakfile, inifile, 'BufferZone', 'Parameters'//str(1:1), '1.0,2.0', sRes)
            end if
            is = nfields + 1; call LIST_REAL(sRes, is, dummies)
            if (is == 1) then
                variable%strength(:) = dummies(1)
                variable%sigma(:) = 2.0_wp
            else if (is == 2) then
                variable%strength(:) = dummies(1)
                variable%sigma(:) = dummies(2)
            else if (is == nfields + 1) then
                variable%strength(1:nfields) = dummies(1:nfields)
                variable%sigma(:) = dummies(nfields + 1)
            else
                call TLab_Write_ASCII(wfile, 'DNS_READ_LOCAL. Wrong number of values in BufferZone.ParametersUImin.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
            do is = 1, nfields
                if (variable%strength(is) /= 0.0_wp) variable%active(is) = .true.
            end do

            call ScanFile_Char(bakfile, inifile, 'BufferZone', 'HardValues'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                is = nfields; call LIST_REAL(sRes, is, variable%hardvalues)
                if (is == nfields) then
                    variable%hard = .true.
                else
                    call TLab_Write_ASCII(wfile, 'DNS_READ_LOCAL. Wrong number of values in BufferZone.HardValues.'//trim(adjustl(tag))//'.')
                    call TLab_Stop(DNS_ERROR_OPTION)
                end if
            end if

        end if

    end subroutine BOUNDARY_BUFFER_READBLOCK

    ! ###################################################################
    ! ###################################################################
    subroutine BOUNDARY_BUFFER_INITIALIZE(q, s, txc)
        ! Careful, txc is here contiguous
        real(wp), dimension(isize_field, *), intent(IN) :: q, s
        real(wp), dimension(isize_field, 6), intent(INOUT) :: txc

        ! ###################################################################
        BuffFlowImin%offset = 0; BuffFlowImax%offset = g(1)%size - BuffFlowImax%size
        BuffFlowJmin%offset = 0; BuffFlowJmax%offset = g(2)%size - BuffFlowJmax%size
        BuffScalImin%offset = 0; BuffScalImax%offset = g(1)%size - BuffScalImax%size
        BuffScalJmin%offset = 0; BuffScalJmax%offset = g(2)%size - BuffScalJmax%size

        BuffFlowImin%total_size = BuffFlowImin%size; BuffFlowImax%total_size = BuffFlowImax%size
        BuffFlowJmin%total_size = BuffFlowJmin%size; BuffFlowJmax%total_size = BuffFlowJmax%size
        BuffScalImin%total_size = BuffScalImin%size; BuffScalImax%total_size = BuffScalImax%size
        BuffScalJmin%total_size = BuffScalJmin%size; BuffScalJmax%total_size = BuffSCalJmax%size
#ifdef USE_MPI
        BuffFlowImin%size = min(imax, max(0, BuffFlowImin%total_size - ims_offset_i))
        BuffFlowImax%size = min(imax, max(0, ims_offset_i + imax - BuffFlowImax%offset))
        BuffScalImin%size = min(imax, max(0, BuffScalImin%total_size - ims_offset_i))
        BuffScalImax%size = min(imax, max(0, ims_offset_i + imax - BuffScalImax%offset))
#endif
        BuffFlowImin%nfields = inb_flow; BuffFlowImax%nfields = inb_flow
        BuffFlowJmin%nfields = inb_flow; BuffFlowJmax%nfields = inb_flow
        BuffScalImin%nfields = inb_scal; BuffScalImax%nfields = inb_scal
        BuffScalJmin%nfields = inb_scal; BuffScalJmax%nfields = inb_scal

        BuffFlowImin%form = FORM_POWER_MIN; BuffFlowImax%form = FORM_POWER_MAX
        BuffFlowJmin%form = FORM_POWER_MIN; BuffFlowJmax%form = FORM_POWER_MAX
        BuffScalImin%form = FORM_POWER_MIN; BuffScalImax%form = FORM_POWER_MAX
        BuffScalJmin%form = FORM_POWER_MIN; BuffScalJmax%form = FORM_POWER_MAX

        ! ###################################################################
        txc(:, 2:inb_flow + 1) = q(:, 1:inb_flow)
        if (nse_eqns == DNS_EQNS_TOTAL) then
            txc(:, 1) = q(:, 5) ! Density
            dummy = 0.5_wp*CRATIO_INV
            txc(:, 5) = q(:, 4) + dummy*(q(:, 1)*q(:, 1) + q(:, 2)*q(:, 2) + q(:, 3)*q(:, 3))
            txc(:, 6) = 1.0_wp
        else if (nse_eqns == DNS_EQNS_INTERNAL) then
            txc(:, 1) = q(:, 5) ! Density
            txc(:, 6) = 1.0_wp
        else
            txc(:, 1) = 1.0_wp ! Density
        end if

        if (BuffFlowImin%total_size > 0) call INI_BLOCK(1, trim(adjustl(tag_flow))//'bcs.imin', BuffFlowImin, txc(1, 1), txc(1, 2))
        if (BuffFlowImax%total_size > 0) call INI_BLOCK(1, trim(adjustl(tag_flow))//'bcs.imax', BuffFlowImax, txc(1, 1), txc(1, 2))
        if (BuffFlowJmin%total_size > 0) call INI_BLOCK(2, trim(adjustl(tag_flow))//'bcs.jmin', BuffFlowJmin, txc(1, 1), txc(1, 2))
        if (BuffFlowJmax%total_size > 0) call INI_BLOCK(2, trim(adjustl(tag_flow))//'bcs.jmax', BuffFlowJmax, txc(1, 1), txc(1, 2))

        if (BuffScalImin%total_size > 0) call INI_BLOCK(1, trim(adjustl(tag_scal))//'bcs.imin', BuffScalImin, txc(1, 1), s)
        if (BuffScalImax%total_size > 0) call INI_BLOCK(1, trim(adjustl(tag_scal))//'bcs.imax', BuffScalImax, txc(1, 1), s)
        if (BuffScalJmin%total_size > 0) call INI_BLOCK(2, trim(adjustl(tag_scal))//'bcs.jmin', BuffScalJmin, txc(1, 1), s)
        if (BuffScalJmax%total_size > 0) call INI_BLOCK(2, trim(adjustl(tag_scal))//'bcs.jmax', BuffScalJmax, txc(1, 1), s)

        return
    end subroutine BOUNDARY_BUFFER_INITIALIZE

    ! ###################################################################
    ! ###################################################################
    subroutine INI_BLOCK(idir, tag, item, rho, a)
        integer(wi), intent(IN) :: idir             ! Wall-normal direction of buffer zone
        character(len=*), intent(IN) :: tag         ! File name information
        type(buffer_dt), intent(INOUT) :: item
        real(wp), intent(IN) :: rho(isize_field), a(isize_field, item%nfields)

        ! -------------------------------------------------------------------
        integer(wi) io_sizes(5), id
        real(wp), dimension(2) :: var_minmax

        character*32 str, varname(item%nfields)
        character*128 line
#ifdef USE_MPI
        integer sa_comm_color
        integer, parameter :: sa_ndims = 3
        integer, dimension(sa_ndims) :: sa_size, sa_locsize, sa_offset
        real(wp), dimension(2) :: dummy2
#endif

        ! ###################################################################
        ! Reference fields

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING INI_BLOCK (boundary_buffer.f90)')
#endif

        if (item%size > 0) then
            if (idir == 1) allocate (item%ref(item%size, jmax, kmax, item%nfields))
            if (idir == 2) allocate (item%ref(imax, item%size, kmax, item%nfields))
        end if

        do iq = 1, item%nfields; write (varname(iq), *) iq; end do

        select case (idir)
        case (1)
#ifdef USE_MPI
            id = IO_SUBARRAY_BUFFER_ZOY
            io_aux(id)%offset = 0
            if (ims_npro_i > 1 .and. item%total_size == imax) then
                ! Buffer lives on exactly one PE in x; we can use ims_comm_z
                if (item%size > 0) then
                    io_aux(id)%active = .true.
                    io_aux(id)%subarray = IO_CREATE_SUBARRAY_ZOY(jmax*item%total_size, kmax, MPI_REAL8)
                    io_aux(id)%communicator = ims_comm_z
                else
                    io_aux(id)%active = .false.
                end if
            else ! Buffer occupies more or less than one PE --> need new subarray and communicator
                sa_comm_color = MPI_UNDEFINED
                if (item%size > 0) then
                    io_aux(id)%active = .true.
                    sa_size = [item%total_size, jmax, kmax*ims_npro_k]
                    sa_locsize = [item%size, jmax, kmax]
                    if (item%offset == 0) then
                        sa_offset = [ims_offset_i, 0, kmax*ims_pro_k]
                    else
                        sa_offset = [max(0, ims_offset_i - item%offset), 0, kmax*ims_pro_k]
                    end if

                    call MPI_Type_create_subarray(sa_ndims, sa_size, sa_locsize, sa_offset, MPI_ORDER_FORTRAN, MPI_REAL8, io_aux(id)%subarray, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                    call MPI_TYPE_COMMIT(io_aux(id)%subarray, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)

                    sa_comm_color = 1
                else
                    io_aux(id)%active = .false.
                end if
                io_aux(id)%communicator = MPI_UNDEFINED
                call MPI_Comm_Split(MPI_COMM_WORLD, sa_comm_color, ims_pro, io_aux(id)%communicator, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)

            end if
#endif
            idummy = item%size*jmax*kmax; 
            io_sizes = [idummy, 1, idummy, 1, item%nfields]
        case (2)
            id = IO_SUBARRAY_BUFFER_XOZ
            io_aux(id)%offset = 0
            io_aux(id)%precision = IO_TYPE_DOUBLE
#ifdef USE_MPI
            io_aux(id)%active = .true.
            io_aux(id)%communicator = MPI_COMM_WORLD
            io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ(imax, item%size, kmax, MPI_REAL8)
#endif
            idummy = imax*item%size*kmax; io_sizes = (/idummy, 1, idummy, 1, item%nfields/)
        end select

        if (BuffLoad) then
            call IO_READ_SUBARRAY(io_aux(id), tag, varname, item%ref, io_sizes)

        else
            select case (idir)
            case (1)
                do iq = 1, item%nfields
                    do j = 1, jmax
                        do iloc = 1, item%size
                            i = item%offset + iloc
                            if (.not. item%hard) item%hardvalues(iq) = COV2V1D(imax, jmax, kmax, i, j, rho, a(1, iq))
                            item%ref(iloc, j, :, iq) = item%hardvalues(iq)
                        end do
                    end do
                end do

            case (2)
                if (imode_sim == DNS_MODE_TEMPORAL) then
                    do iq = 1, item%nfields
                        do jloc = 1, item%size
                            j = item%offset + jloc
                            if (.not. item%hard) item%hardvalues(iq) = COV2V2D(imax, jmax, kmax, j, rho, a(1, iq))
                            item%ref(:, jloc, :, iq) = item%hardvalues(iq)
                        end do
                    end do
                else if (imode_sim == DNS_MODE_SPATIAL) then
                    do iq = 1, item%nfields
                        do jloc = 1, item%size
                            j = item%offset + jloc
                            do i = 1, imax
                                if (.not. item%hard) item%hardvalues(iq) = COV2V1D(imax, jmax, kmax, i, j, rho, a(1, iq))
                                item%ref(i, jloc, :, iq) = item%hardvalues(iq)
                            end do
                        end do
                    end do
                end if

            end select

            write (str, *) itime; str = trim(adjustl(tag))//'.'//trim(adjustl(str))
            call IO_WRITE_SUBARRAY(io_aux(id), str, varname, item%ref, io_sizes)

        end if

        do iq = 1, item%nfields ! Control
#ifdef USE_MPI
            if (io_aux(id)%active) then
                var_minmax = [minval(item%ref(:, :, :, iq)), -maxval(item%ref(:, :, :, iq))]
            else
                var_minmax = [huge(dummy), huge(dummy)]
            end if

            call MPI_ALLREDUCE(var_minmax, dummy2, 2, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
            var_minmax = dummy2; var_minmax(2) = -var_minmax(2)
#else
            var_minmax = [minval(item%ref(:, :, :, iq)), maxval(item%ref(:, :, :, iq))]
#endif
            write (line, 10) var_minmax(2)
            write (str, 10) var_minmax(1)

            line = trim(adjustl(str))//' and '//trim(adjustl(line))
            line = 'Checking bounds of field '//trim(adjustl(tag))//'.'//trim(adjustl(varname(iq)))//': '//trim(adjustl(line))
            call TLab_Write_ASCII(lfile, line)
        end do

        ! -----------------------------------------------------------------------
        ! Strength of the relaxation terms
        allocate (item%tau(item%size, item%nfields))
        if (item%size > 1) then
            do iq = 1, item%nfields
                dummy = 1.0_wp/(g(idir)%nodes(item%offset + item%size) - g(idir)%nodes(item%offset + 1)) ! Inverse of segment length
                do jloc = 1, item%size
                    j = item%offset + jloc
                    if (item%form == FORM_POWER_MAX) &
                        item%tau(jloc, iq) = item%strength(iq)*((g(idir)%nodes(j) - g(idir)%nodes(item%offset + 1))*dummy)**item%sigma(iq)
                    if (item%form == FORM_POWER_MIN) &
                        item%tau(jloc, iq) = item%strength(iq)*((g(idir)%nodes(item%offset + item%size) - g(idir)%nodes(j))*dummy)**item%sigma(iq)
                end do
            end do
        end if

        ! -----------------------------------------------------------------------
        ! Filters at boundaries; nseeds to be checked
#ifdef USE_MPI
        if (item%type == DNS_BUFFER_FILTER) then
            select case (idir)
            case (1)
                ! call TLab_Write_ASCII(lfile, 'Initialize MPI types for Ox BCs explicit filter.')
                ! id = TLAB_MPI_TRP_K_OUTBCS
                idummy = item%size*jmax
                ! call TLabMPI_TypeK_Create(ims_npro_k, kmax, idummy, 1, 1, 1, 1, id)
                ims_trp_plan_k(TLAB_MPI_TRP_K_OUTBCS) = TLabMPI_Trp_TypeK_Create(kmax, idummy, 1, 1, 1, 1, 'Ox BCs explicit filter.')

            case (2)
                ! call TLab_Write_ASCII(lfile, 'Initialize MPI types for Oy BCs explicit filter.')
                ! id = TLAB_MPI_TRP_K_TOPBCS
                idummy = imax*item%size
                ! call TLabMPI_TypeK_Create(ims_npro_k, kmax, idummy, 1, 1, 1, 1, id)
                ims_trp_plan_k(TLAB_MPI_TRP_K_TOPBCS) = TLabMPI_Trp_TypeK_Create(kmax, idummy, 1, 1, 1, 1, 'Oy BCs explicit filter.')

            end select
        end if
#endif

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'LEAVING INI_BLOCK (boundary_buffer.f90)')
#endif
        return
10      format(G_FORMAT_R)
    end subroutine INI_BLOCK

    ! ###################################################################
    ! ###################################################################
    subroutine BOUNDARY_BUFFER_RELAX_SCAL()
        use TLab_Arrays, only: q, s
        use DNS_ARRAYS, only: hs

        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL)
            if (BuffScalImin%size > 0) call RELAX_BLOCK_RHO(1, BuffScalImin, s, hs, q(:, 5))
            if (BuffScalImax%size > 0) call RELAX_BLOCK_RHO(1, BuffScalImax, s, hs, q(:, 5))
            if (BuffScalJmin%size > 0) call RELAX_BLOCK_RHO(2, BuffScalJmin, s, hs, q(:, 5))
            if (BuffScalJmax%size > 0) call RELAX_BLOCK_RHO(2, BuffScalJmax, s, hs, q(:, 5))
        case DEFAULT
            if (BuffScalImin%size > 0) call RELAX_BLOCK(1, BuffScalImin, s, hs)
            if (BuffScalImax%size > 0) call RELAX_BLOCK(1, BuffScalImax, s, hs)
            if (BuffScalJmin%size > 0) call RELAX_BLOCK(2, BuffScalJmin, s, hs)
            if (BuffScalJmax%size > 0) call RELAX_BLOCK(2, BuffScalJmax, s, hs)
        end select

        return
    end subroutine BOUNDARY_BUFFER_RELAX_SCAL

    ! ###################################################################
    ! ###################################################################
    subroutine BOUNDARY_BUFFER_RELAX_SCAL_I(is, s, hs)
        implicit none

        integer(wi), intent(IN) :: is ! field to which relaxation is applied
        real(wp), intent(IN) :: s(isize_field)
        real(wp), intent(OUT) :: hs(isize_field)

        ! ###################################################################
        if (BuffScalImin%size > 0) call RELAX_BLOCK_I(1, BuffScalImin, s, hs, is)
        if (BuffScalImax%size > 0) call RELAX_BLOCK_I(1, BuffScalImax, s, hs, is)
        if (BuffScalJmin%size > 0) call RELAX_BLOCK_I(2, BuffScalJmin, s, hs, is)
        if (BuffScalJmax%size > 0) call RELAX_BLOCK_I(2, BuffScalJmax, s, hs, is)

        return
    end subroutine BOUNDARY_BUFFER_RELAX_SCAL_I

    ! ###################################################################
    ! ###################################################################
    subroutine BOUNDARY_BUFFER_RELAX_FLOW()
        use TLab_Arrays, only: q
        use DNS_ARRAYS, only: hq

        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL)
            if (BuffFlowImin%size > 0) call RELAX_BLOCK_CF(1, BuffFlowImin, q, hq)
            if (BuffFlowImax%size > 0) call RELAX_BLOCK_CF(1, BuffFlowImax, q, hq)
            if (BuffFlowJmin%size > 0) call RELAX_BLOCK_CF(2, BuffFlowJmin, q, hq)
            if (BuffFlowJmax%size > 0) call RELAX_BLOCK_CF(2, BuffFlowJmax, q, hq)
        case DEFAULT
            if (BuffFlowImin%size > 0) call RELAX_BLOCK(1, BuffFlowImin, q, hq)
            if (BuffFlowImax%size > 0) call RELAX_BLOCK(1, BuffFlowImax, q, hq)
            if (BuffFlowJmin%size > 0) call RELAX_BLOCK(2, BuffFlowJmin, q, hq)
            if (BuffFlowJmax%size > 0) call RELAX_BLOCK(2, BuffFlowJmax, q, hq)
        end select

        return
    end subroutine BOUNDARY_BUFFER_RELAX_FLOW

    ! ###################################################################
    ! ###################################################################
    subroutine RELAX_BLOCK(idir, item, s, hs)
        integer(wi), intent(IN) :: idir
        type(buffer_dt), intent(IN) :: item
        real(wp), intent(IN) :: s(imax, jmax, kmax, item%nfields)
        real(wp), intent(OUT) :: hs(imax, jmax, kmax, item%nfields)

        ! ###################################################################
        select case (idir)
        case (1)
            do iq = 1, item%nfields
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hs(j, :, :, iq) = hs(j, :, :, iq) - item%Tau(jloc, iq)*(s(j, :, :, iq) - item%Ref(jloc, :, :, iq))
                end do
            end do

        case (2)
            do iq = 1, item%nfields
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hs(:, j, :, iq) = hs(:, j, :, iq) - item%Tau(jloc, iq)*(s(:, j, :, iq) - item%Ref(:, jloc, :, iq))
                end do
            end do

        end select

        return
    end subroutine RELAX_BLOCK

    ! ###################################################################
    ! ###################################################################
    subroutine RELAX_BLOCK_RHO(idir, item, s, hs, rho)
        integer(wi), intent(IN) :: idir
        type(buffer_dt), intent(IN) :: item
        real(wp), intent(IN) :: s(imax, jmax, kmax, item%nfields)
        real(wp), intent(OUT) :: hs(imax, jmax, kmax, item%nfields)
        real(wp), intent(IN) :: rho(imax, jmax, kmax)

        ! ###################################################################
        select case (idir)
        case (1)
            do iq = 1, item%nfields
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hs(j, :, :, iq) = hs(j, :, :, iq) - item%Tau(jloc, iq)*(rho(j, :, :)*s(j, :, :, iq) - item%Ref(jloc, :, :, iq))
                end do
            end do

        case (2)
            do iq = 1, item%nfields
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hs(:, j, :, iq) = hs(:, j, :, iq) - item%Tau(jloc, iq)*(rho(:, j, :)*s(:, j, :, iq) - item%Ref(:, jloc, :, iq))
                end do
            end do

        end select

        return
    end subroutine RELAX_BLOCK_RHO

    ! ###################################################################
    ! ###################################################################
    subroutine RELAX_BLOCK_I(idir, item, s, hs, iq)
        integer(wi), intent(IN) :: idir, iq
        type(buffer_dt), intent(IN) :: item
        real(wp), intent(IN) :: s(imax, jmax, kmax)
        real(wp), intent(OUT) :: hs(imax, jmax, kmax)

        ! ###################################################################
        select case (idir)
        case (1)
            do jloc = 1, item%size
                j = item%offset + jloc
                hs(j, :, :) = hs(j, :, :) - item%Tau(jloc, iq)*(s(j, :, :) - item%Ref(jloc, :, :, iq))
            end do

        case (2)
            do jloc = 1, item%size
                j = item%offset + jloc
                hs(:, j, :) = hs(:, j, :) - item%Tau(jloc, iq)*(s(:, j, :) - item%Ref(:, jloc, :, iq))
            end do

        end select

        return
    end subroutine RELAX_BLOCK_I

    ! ###################################################################
    ! ###################################################################
    subroutine RELAX_BLOCK_CF(idir, item, q, hq)  ! Compressible flow case
        implicit none

        integer(wi), intent(IN) :: idir
        type(buffer_dt), intent(IN) :: item
        real(wp), intent(IN) :: q(imax, jmax, kmax, item%nfields)
        real(wp), intent(OUT) :: hq(imax, jmax, kmax, item%nfields)

        ! ###################################################################
        if (nse_eqns == DNS_EQNS_TOTAL) then
            dummy = 0.5_wp*CRATIO_INV
        end if

        select case (idir)
        case (1)
            do iq = 1, item%nfields - 2
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(j, :, :, iq) = hq(j, :, :, iq) - item%Tau(jloc, iq)*(q(j, :, :, 5)*q(j, :, :, iq) - item%Ref(jloc, :, :, iq))
                end do
            end do
            iq = item%nfields - 1       ! Energy variable
            if (nse_eqns == DNS_EQNS_TOTAL) then
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(j, :, :, iq) = hq(j, :, :, iq) - item%Tau(jloc, iq)*(q(j, :, :, 5)*(q(j, :, :, iq) &
                                + dummy*(q(j, :, :, 1)*q(j, :, :, 1) + q(j, :, :, 2)*q(j, :, :, 2) + q(j, :, :, 3)*q(j, :, :, 3))) &
                                                                            - item%Ref(jloc, :, :, iq))
                end do
            else
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(j, :, :, iq) = hq(j, :, :, iq) - item%Tau(jloc, iq)*(q(j, :, :, 5)*q(j, :, :, iq) - item%Ref(jloc, :, :, iq))
                end do
            end if
            iq = item%nfields         ! Density, should be iq = 5
            do jloc = 1, item%size
                j = item%offset + jloc
                hq(j, :, :, iq) = hq(j, :, :, iq) - item%Tau(jloc, iq)*(q(j, :, :, iq) - item%Ref(jloc, :, :, iq))
            end do

        case (2)
            do iq = 1, item%nfields - 2
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(:, j, :, iq) = hq(:, j, :, iq) - item%Tau(jloc, iq)*(q(:, j, :, 5)*q(:, j, :, iq) - item%Ref(:, jloc, :, iq))
                end do
            end do
            iq = item%nfields - 1       ! Energy variable
            if (nse_eqns == DNS_EQNS_TOTAL) then
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(:, j, :, iq) = hq(:, j, :, iq) - item%Tau(jloc, iq)*(q(:, j, :, 5)*(q(:, j, :, iq) &
                                + dummy*(q(:, j, :, 1)*q(:, j, :, 1) + q(:, j, :, 2)*q(:, j, :, 2) + q(:, j, :, 3)*q(:, j, :, 3))) &
                                                                            - item%Ref(:, jloc, :, iq))
                end do
            else
                do jloc = 1, item%size
                    j = item%offset + jloc
                    hq(:, j, :, iq) = hq(:, j, :, iq) - item%Tau(jloc, iq)*(q(:, j, :, 5)*q(:, j, :, iq) - item%Ref(:, jloc, :, iq))
                end do
            end if
            iq = item%nfields         ! Density, should be iq = 5
            do jloc = 1, item%size
                j = item%offset + jloc
                hq(:, j, :, iq) = hq(:, j, :, iq) - item%Tau(jloc, iq)*(q(:, j, :, iq) - item%Ref(:, jloc, :, iq))
            end do

        end select

        return
    end subroutine RELAX_BLOCK_CF

    !########################################################################
    !########################################################################
    subroutine BOUNDARY_BUFFER_FILTER(rho, u, v, w, e, z1, txc1, txc2, txc3, txc4, txc5)
        real(wp), dimension(imax, jmax, kmax) :: rho, u, v, w, e
        real(wp), dimension(imax, jmax, kmax, *) :: z1
        real(wp), dimension(BuffFlowImax%size, jmax, kmax) :: txc1, txc2, txc3, txc4, txc5

        ! -------------------------------------------------------------------
        integer(wi) id, k, buff_imax!, ibc_x(4), ibc_y(4), ibc_z(4)
        real(wp) eta, delta, amp, ampr, rho_ratio

        ! ###################################################################
!!! Routines OPR_FILTER have been changed. This routine needs to be updates !!!
        call TLab_Write_ASCII(efile, 'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)

        ! BCs for the filters (see routine FILTER)
        ! ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
        ! ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0
        ! ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0

        ! ###################################################################
        ! Bottom boundary
        ! ###################################################################
        if (BuffFlowJmin%size > 1) then
            call TLab_Write_ASCII(efile, 'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! ###################################################################
        ! Top boundary
        ! ###################################################################
        if (BuffFlowJmax%size > 1) then
            call TLab_Write_ASCII(efile, 'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! ###################################################################
        ! Outflow boundary
        ! ###################################################################
        if (BuffFlowImax%size > 1) then
            id = TLAB_MPI_TRP_K_OUTBCS
            buff_imax = imax - BuffFlowImax%size + iloc
            ! -------------------------------------------------------------------
            ! Flow
            ! -------------------------------------------------------------------
            do k = 1, kmax
                do j = 1, jmax
                    do iloc = 1, BuffFlowImax%size ! Outflow boundary
                        i = imax - BuffFlowImax%size + iloc
                        txc1(iloc, j, k) = rho(i, j, k)
                        txc2(iloc, j, k) = rho(i, j, k)*u(i, j, k)
                        txc3(iloc, j, k) = rho(i, j, k)*v(i, j, k)
                        txc4(iloc, j, k) = rho(i, j, k)*w(i, j, k)
                        txc5(iloc, j, k) = rho(i, j, k)*e(i, j, k)
                        !              txc2(iloc,j,k) = u(i,j,k)
                        !              txc3(iloc,j,k) = v(i,j,k)
                        !              txc4(iloc,j,k) = w(i,j,k)
                        !              txc5(iloc,j,k) = e(i,j,k)
                    end do
                end do
            end do

            ! To be checked for the new formulation of opr_filter
            !   CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc1, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
            !   CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
            !   CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc3, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
            !   CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc4, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
            !   CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc5, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

            ! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
            delta = (g(1)%nodes(imax) - g(1)%nodes(buff_imax))/16.0_wp
            do i = buff_imax, imax
                iloc = i - buff_imax + 1

                eta = g(1)%nodes(i) - 0.5_wp*(g(1)%nodes(imax) + g(1)%nodes(buff_imax))
                amp = 0.5_wp*(1.0_wp + tanh(0.5_wp*eta/delta))
                ampr = 1.0_wp - amp

                do k = 1, kmax
                    do j = 1, jmax
                        rho_ratio = rho(i, j, k)
                        rho(i, j, k) = ampr*rho(i, j, k) + amp*txc1(iloc, j, k)
                        rho_ratio = rho_ratio/rho(i, j, k)

                        u(i, j, k) = ampr*rho_ratio*u(i, j, k) + amp*txc2(iloc, j, k)/rho(i, j, k)
                        v(i, j, k) = ampr*rho_ratio*v(i, j, k) + amp*txc3(iloc, j, k)/rho(i, j, k)
                        w(i, j, k) = ampr*rho_ratio*w(i, j, k) + amp*txc4(iloc, j, k)/rho(i, j, k)
                        e(i, j, k) = ampr*rho_ratio*e(i, j, k) + amp*txc5(iloc, j, k)/rho(i, j, k)

                        ! store rho_ratio to be used in scalar section
                        txc1(iloc, j, k) = rho_ratio

                        !              rho(i,j,k) = ampr*rho(i,j,k) + amp*txc1(iloc,j,k)
                        !              u(i,j,k)   = ampr*u(i,j,k)   + amp*txc2(iloc,j,k)
                        !              v(i,j,k)   = ampr*v(i,j,k)   + amp*txc3(iloc,j,k)
                        !              w(i,j,k)   = ampr*w(i,j,k)   + amp*txc4(iloc,j,k)
                        !              e(i,j,k)   = ampr*e(i,j,k)   + amp*txc5(iloc,j,k)

                    end do
                end do
            end do

            ! -------------------------------------------------------------------
            ! Scalar
            ! -------------------------------------------------------------------
            !     IF ( scal_on ) THEN
            do is = 1, inb_scal

                do k = 1, kmax
                    do j = 1, jmax
                        do iloc = 1, BuffFlowImax%size ! Outflow boundary
                            i = imax - BuffFlowImax%size + iloc
                            txc2(iloc, j, k) = rho(i, j, k)*z1(i, j, k, is)
                            !                    txc2(iloc,j,k) = z1(i,j,k,is)
                        end do
                    end do
                end do

                ! To be checked for the new formulation of opr_filter
                ! CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id,  txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

                ! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
                delta = (g(1)%nodes(imax) - g(1)%nodes(buff_imax))/16.0_wp
                do i = buff_imax, imax
                    iloc = i - buff_imax + 1

                    eta = g(1)%nodes(i) - 0.5_wp*(g(1)%nodes(imax) + g(1)%nodes(buff_imax))
                    amp = 0.5_wp*(1.0_wp + tanh(0.5_wp*eta/delta))
                    ampr = 1.0_wp - amp

                    do k = 1, kmax
                        do j = 1, jmax
                            z1(i, j, k, is) = ampr*txc1(iloc, j, k)*z1(i, j, k, is) + &
                                              amp*txc2(iloc, j, k)/rho(i, j, k)
                            !                    z1(i,j,k,is) = ampr*z1(i,j,k,is) + amp*txc2(iloc,j,k)
                        end do
                    end do

                end do

            end do
            !     ENDIF

        end if

        return
    end subroutine BOUNDARY_BUFFER_FILTER
end module BOUNDARY_BUFFER

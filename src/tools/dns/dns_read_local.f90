#include "dns_error.h"
#include "dns_const.h"
#include "avgij_map.h"

subroutine DNS_READ_LOCAL(inifile)

    use TLAB_TYPES, only: MAX_MODES
    use TLAB_CONSTANTS, only: wp, wi, big_wp, efile, lfile, wfile
    use TLAB_VARS
    use TLAB_PROCS
    use PARTICLE_VARS
    use DNS_LOCAL
    use TIME, only: rkm_mode, dtime, cfla, cfld, cflr
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use BOUNDARY_INFLOW
    use STATISTICS
    use PLANES

    implicit none

    character*(*) inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*64 lstr
    character*32 bakfile
    integer is, inb_scal_local1
    integer(wi) idummy
    real(wp) dummy(inb_flow + inb_scal + 1)

    integer :: bcs_visc_imin, bcs_visc_imax
    integer :: bcs_visc_jmin, bcs_visc_jmax
    integer :: bcs_visc_kmin, bcs_visc_kmax

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'

    call TLAB_WRITE_ASCII(lfile, 'Reading local input data.')

! ###################################################################
! Main Section
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Main]')
    call TLAB_WRITE_ASCII(bakfile, '#TimeOrder=<RungeKuttaExplicit3/RungeKuttaExplicit4/RungeKuttaDiffusion3>')
    call TLAB_WRITE_ASCII(bakfile, '#TimeStep=<value (used if CFL is negative)>')
    call TLAB_WRITE_ASCII(bakfile, '#TimeCFL=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#TimeDiffusiveCFL=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#TimeReactiveCFL=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#TermDivergence=<none/remove>')
    call TLAB_WRITE_ASCII(bakfile, '#RhsMode=<split/combined/nonblocking>')

    call SCANINICHAR(bakfile, inifile, 'Main', 'TimeOrder', 'dummy', sRes)
    if (trim(adjustl(sRes)) == 'rungekuttaexplicit3') then; rkm_mode = RKM_EXP3; lstr = '0.6'; 
    elseif (trim(adjustl(sRes)) == 'rungekuttaexplicit4') then; rkm_mode = RKM_EXP4; lstr = '1.2'; 
    elseif (trim(adjustl(sRes)) == 'rungekuttadiffusion3') then; rkm_mode = RKM_IMP3_DIFFUSION; lstr = '0.6'; 
!  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttasource3'    ) THEN; rkm_mode = RKM_IMP3_SOURCE;
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong TimeOrder option.')
        call TLAB_STOP(DNS_ERROR_RKORDER)
    end if

    ! Default cfla value set in lstr while reading TimeOrder
    call SCANINIREAL(bakfile, inifile, 'Main', 'TimeCFL', trim(adjustl(lstr)), cfla)
    write (lstr, *) 0.25_wp*cfla ! Default value for diffusive CFL
    call SCANINIREAL(bakfile, inifile, 'Main', 'TimeDiffusiveCFL', trim(adjustl(lstr)), cfld)
    write (lstr, *) 0.5_wp*cfla ! Default value for reactive CFL
    call SCANINIREAL(bakfile, inifile, 'Main', 'TimeReactiveCFL', trim(adjustl(lstr)), cflr)
    call SCANINIREAL(bakfile, inifile, 'Main', 'TimeStep', '0.05', dtime)

! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Main', 'TermDivergence', 'remove', sRes)
    if (trim(adjustl(sRes)) == 'none') then; remove_divergence = .false.
    else if (trim(adjustl(sRes)) == 'remove') then; remove_divergence = .true.
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong TermDivergence option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'RhsMode', 'combined', sRes)
    if (trim(adjustl(sRes)) == 'split') then; imode_rhs = EQNS_RHS_SPLIT
    else if (trim(adjustl(sRes)) == 'combined') then; imode_rhs = EQNS_RHS_COMBINED
#ifdef USE_PSFFT
    else if (trim(adjustl(sRes)) == 'nonblocking') then; imode_rhs = EQNS_RHS_NONBLOCKING
#endif
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong RhsMode option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! ###################################################################
! Iteration Section
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Iteration]')
    call TLAB_WRITE_ASCII(bakfile, '#Start=<integral start time>')
    call TLAB_WRITE_ASCII(bakfile, '#End=<integral stop time>')
    call TLAB_WRITE_ASCII(bakfile, '#Restart=<restart time step>')
    call TLAB_WRITE_ASCII(bakfile, '#Statistics=<statistics time step>')
    call TLAB_WRITE_ASCII(bakfile, '#IteraLog=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Saveplanes=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#RunAvera=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#Runtime=<seconds>')
    call TLAB_WRITE_ASCII(bakfile, '#ObsLog=<None/Ekman>')

    call SCANINIINT(bakfile, inifile, 'Iteration', 'Start', '0', nitera_first)
    call SCANINIINT(bakfile, inifile, 'Iteration', 'End', '0', nitera_last)
    call SCANINIINT(bakfile, inifile, 'Iteration', 'Restart', '50', nitera_save)
    call SCANINIINT(bakfile, inifile, 'Iteration', 'Statistics', '50', nitera_stats)
    call SCANINIINT(bakfile, inifile, 'Iteration', 'IteraLog', '10', nitera_log)
    call SCANINIINT(bakfile, inifile, 'Iteration', 'Saveplanes', '-1', nitera_pln)
    call SCANINIREAL(bakfile,inifile, 'Iteration', 'Runtime', '10000000', nruntime_sec)

! Accumulate statistics in spatial mode
    call SCANINIINT(bakfile, inifile, 'Iteration', 'SaveStats', '-1', nitera_stats_spa)

! Domain Filter (Should we move it to Iteration?)
    call SCANINIINT(bakfile, inifile, 'Filter', 'Step', '0', nitera_filter)
    if (nitera_filter == 0) FilterDomain(:)%type = DNS_FILTER_NONE


    call SCANINICHAR(bakfile, inifile, 'Iteration', 'ObsLog', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; dns_obs_log = OBS_TYPE_NONE
    else if (trim(adjustl(sRes)) == 'ekman') then; dns_obs_log = OBS_TYPE_EKMAN
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. ObsLog.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if


! ###################################################################
! Control Limits
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Control]')
    call TLAB_WRITE_ASCII(bakfile, '#FlowLimit=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#MinPressure=<pressure>')
    call TLAB_WRITE_ASCII(bakfile, '#MaxPressure=<pressure>')
    call TLAB_WRITE_ASCII(bakfile, '#MinDensity=<density>')
    call TLAB_WRITE_ASCII(bakfile, '#MaxDensity=<density>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalLimit=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#MinScalar=<scalar>')
    call TLAB_WRITE_ASCII(bakfile, '#MaxScalar=<scalar>')

    bound_r%active = .false.
    bound_p%active = .false.
    call SCANINICHAR(bakfile, inifile, 'Control', 'FlowLimit', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then
        bound_r%active = .true.
        bound_p%active = .true.
    end if
! Final check in last section of the routine
    call SCANINIREAL(bakfile, inifile, 'Control', 'MinPressure', '-1.0', bound_p%min)
    call SCANINIREAL(bakfile, inifile, 'Control', 'MaxPressure', '-1.0', bound_p%max)
    call SCANINIREAL(bakfile, inifile, 'Control', 'MinDensity', '-1.0', bound_r%min)
    call SCANINIREAL(bakfile, inifile, 'Control', 'MaxDensity', '-1.0', bound_r%max)

    bound_d%active = .false.
    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
        bound_d%active = .true.
        bound_d%max = big_wp ! default
        call SCANINICHAR(bakfile, inifile, 'Control', 'MaxDilatation', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            idummy = 1
            call LIST_REAL(sRes, idummy, dummy)
            bound_d%max = dummy(1)
        end if
    end if

    bound_s(:)%active = .false.
    call SCANINICHAR(bakfile, inifile, 'Control', 'ScalLimit', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then
        bound_s(:)%active = .true.
    end if

    bound_s(:)%min = 0.0_wp; inb_scal_local1 = MAX_NSP
    if (any(bound_s(:)%active)) then
        call SCANINICHAR(bakfile, inifile, 'Control', 'MinScalar', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call LIST_REAL(sRes, inb_scal_local1, dummy)
            if (inb_scal_local1 /= inb_scal) then ! Consistency check
                call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. MinScalar size does not match inb_scal.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            else
                bound_s(1:inb_scal)%min = dummy(1:inb_scal)
            end if
        end if
    end if

    bound_s%max = 1.0_wp; inb_scal_local1 = MAX_NSP
    if (any(bound_s(:)%active)) then
        call SCANINICHAR(bakfile, inifile, 'Control', 'MaxScalar', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call LIST_REAL(sRes, inb_scal_local1, dummy)
            if (inb_scal_local1 /= inb_scal) then ! Consistency check
                call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. MaxScalar size does not match inb_scal.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            else
                bound_s(1:inb_scal)%max = dummy(1:inb_scal)
            end if
        end if
    end if

! ###################################################################
! Boundary Conditions
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[BoundaryConditions]')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarImin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarJmin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarSfcTypeJmin=<static/linear>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarSfcTypeJmax=<static/linear>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarCouplingJmin=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarCouplingJmax=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ScalarKmin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#VelocityImin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#VelocityJmin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#VelocityKmin=<none/dirichlet/neumman>')
    call TLAB_WRITE_ASCII(bakfile, '#ViscousI=<none/inflow/outflow>')
    call TLAB_WRITE_ASCII(bakfile, '#ViscousJ=<none/inflow/outflow>')
    call TLAB_WRITE_ASCII(bakfile, '#ViscousK=<none/inflow/outflow>')
    call TLAB_WRITE_ASCII(bakfile, '#SigmaOut=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#SigmaInf=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#BetaTransverse=<value>')

! -------------------------------------------------------------------
! Scalar terms (including surface model at vertical boundaries)
! -------------------------------------------------------------------
    BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
    if (.not. g(1)%periodic) then
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Imin', BcsScalImin)
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Imax', BcsScalImax)
    end if

    BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
    if (.not. g(2)%periodic) then
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Jmin', BcsScalJmin)
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Jmax', BcsScalJmax)
    end if

    BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
    if (.not. g(3)%periodic) then
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Kmin', BcsScalKmin)
        call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Kmax', BcsScalKmax)
    end if

! -------------------------------------------------------------------
! Velocity terms / Euler part in compressible mode
! -------------------------------------------------------------------
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Imin', BcsFlowImin)
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Imax', BcsFlowImax)
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Jmin', BcsFlowJmin)
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Jmax', BcsFlowJmax)
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Kmin', BcsFlowKmin)
    call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Kmax', BcsFlowKmax)

! -------------------------------------------------------------------
! Viscous terms, used only t odefine bcs_inf and bcs_out
! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousI', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; bcs_visc_imin = 0; bcs_visc_imax = 0
    else if (trim(adjustl(sRes)) == 'inflow') then; bcs_visc_imin = 1; bcs_visc_imax = 2
    else if (trim(adjustl(sRes)) == 'outflow') then; bcs_visc_imin = 2; bcs_visc_imax = 2
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousI.')
        call TLAB_STOP(DNS_ERROR_IVSICBC)
    end if

    call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousJ', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; bcs_visc_jmin = 0; bcs_visc_jmax = 0
    else if (trim(adjustl(sRes)) == 'inflow') then; bcs_visc_jmin = 1; bcs_visc_jmax = 2
    else if (trim(adjustl(sRes)) == 'outflow') then; bcs_visc_jmin = 2; bcs_visc_jmax = 2
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousJ.')
        call TLAB_STOP(DNS_ERROR_JVSICBC)
    end if

    call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousK', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; bcs_visc_kmin = 0; bcs_visc_kmax = 0
    else if (trim(adjustl(sRes)) == 'inflow') then; bcs_visc_kmin = 1; bcs_visc_kmax = 2
    else if (trim(adjustl(sRes)) == 'outflow') then; bcs_visc_kmin = 2; bcs_visc_kmax = 2
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousK.')
        call TLAB_STOP(DNS_ERROR_KVSICBC)
    end if

! -------------------------------------------------------------------
! Relaxation coefficients towards reference values in characteristic formulation
! -------------------------------------------------------------------
    BcsDrift = .false.

! Inflow terms
    call SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaInf', '-1.0', dummy(1))
    if (dummy(1) >= 0.0_wp) BcsDrift = .true.
    BcsFlowImin%cinf = dummy(1); BcsFlowImax%cinf = dummy(1) ! so far, all of them the same
    BcsFlowJmin%cinf = dummy(1); BcsFlowJmax%cinf = dummy(1)
    BcsFlowKmin%cinf = dummy(1); BcsFlowKmax%cinf = dummy(1)
    BcsScalImin%cinf = dummy(1); BcsScalImax%cinf = dummy(1)
    BcsScalJmin%cinf = dummy(1); BcsScalJmax%cinf = dummy(1)
    BcsScalKmin%cinf = dummy(1); BcsScalKmax%cinf = dummy(1)

! Outflow terms
    call SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaOut', '-1.0', dummy(1))
    if (dummy(1) >= 0.0_wp) BcsDrift = .true.
    BcsFlowImin%cout = dummy(1); BcsFlowImax%cout = dummy(1) ! so far, all of them the same
    BcsFlowJmin%cout = dummy(1); BcsFlowJmax%cout = dummy(1)
    BcsFlowKmin%cout = dummy(1); BcsFlowKmax%cout = dummy(1)
    BcsScalImin%cout = dummy(1); BcsScalImax%cout = dummy(1)
    BcsScalJmin%cout = dummy(1); BcsScalJmax%cout = dummy(1)
    BcsScalKmin%cout = dummy(1); BcsScalKmax%cout = dummy(1)

! Transverse terms
    call SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'BetaTransverse', '-1.0', dummy(1))
    if (dummy(1) >= 0.0_wp) BcsDrift = .true.
    BcsFlowImin%ctan = dummy(1); BcsFlowImax%ctan = dummy(1) ! so far, all of them the same
    BcsFlowJmin%ctan = dummy(1); BcsFlowJmax%ctan = dummy(1)
    BcsFlowKmin%ctan = dummy(1); BcsFlowKmax%ctan = dummy(1)
    BcsScalImin%ctan = dummy(1); BcsScalImax%ctan = dummy(1)
    BcsScalJmin%ctan = dummy(1); BcsScalJmax%ctan = dummy(1)
    BcsScalKmin%ctan = dummy(1); BcsScalKmax%ctan = dummy(1)

! ###################################################################
! IBM Status Parameter
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[IBMParameter]')
    call TLAB_WRITE_ASCII(bakfile, '#Status=<on/off>')
    call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'Status', 'off', sRes)
    if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
    else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong IBM Status option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! ###################################################################
! Buffer Zone Parameters
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[BufferZone]')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<none/relaxation/filter/both>')
    call TLAB_WRITE_ASCII(bakfile, '#LoadBuffer=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#PointsJmin=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#PointsJmax=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#PointsImin=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#PointsImax=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ParametersJmin=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#ParametersJmax=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#ParametersImin=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#ParametersImax=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#HardValuesJmin=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#HardValuesJmax=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#HardValuesImin=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#HardValuesImax=<values>')

    call SCANINICHAR(bakfile, inifile, 'BufferZone', 'Type', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; BuffType = DNS_BUFFER_NONE
    else if (trim(adjustl(sRes)) == 'relaxation') then; BuffType = DNS_BUFFER_RELAX
    else if (trim(adjustl(sRes)) == 'filter') then; BuffType = DNS_BUFFER_FILTER
    else if (trim(adjustl(sRes)) == 'both') then; BuffType = DNS_BUFFER_BOTH
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong BufferType option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! Load buffer if used also by BCs
    call SCANINICHAR(bakfile, inifile, 'BufferZone', 'LoadBuffer', 'no', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; BuffLoad = .true.
    else; BuffLoad = .false.
    end if

! Sizes; read always because allocation checks if # points is zero
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUImin', '0', BuffFlowImin%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUImax', '0', BuffFlowImax%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUJmin', '0', BuffFlowJmin%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUJmax', '0', BuffFlowJmax%size)

    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSImin', '0', BuffScalImin%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSImax', '0', BuffScalImax%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSJmin', '0', BuffScalJmin%size)
    call SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSJmax', '0', BuffScalJmax%size)

    BuffScalImin%type = BuffType; BuffFlowImin%type = BuffType ! So far, all the same
    BuffScalImax%type = BuffType; BuffFlowImax%type = BuffType
    BuffScalJmin%type = BuffType; BuffFlowJmin%type = BuffType
    BuffScalJmax%type = BuffType; BuffFlowJmax%type = BuffType

    if (BuffScalImin%size /= BuffFlowImin%size .or. &
        BuffScalImax%size /= BuffFlowImax%size .or. &
        BuffScalJmin%size /= BuffFlowJmin%size .or. &
        BuffScalJmax%size /= BuffFlowJmax%size) then ! Because of io_subarray
        call TLAB_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Buffer zone sizes must be equal in flow and scal.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! Parameters
    if (BuffType /= DNS_BUFFER_NONE) then
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UImin', BuffFlowImin, inb_flow)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UImax', BuffFlowImax, inb_flow)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UJmin', BuffFlowJmin, inb_flow)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UJmax', BuffFlowJmax, inb_flow)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SImin', BuffScalImin, inb_scal)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SImax', BuffScalImax, inb_scal)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SJmin', BuffScalJmin, inb_scal)
        call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SJmax', BuffScalJmax, inb_scal)
    end if

! ###################################################################
! Viscosity Control
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[ViscChange]')
    call TLAB_WRITE_ASCII(bakfile, '#Time=<time>')

    call SCANINIREAL(bakfile, inifile, 'ViscChange', 'Time', '0.0', visc_time)

! ###################################################################
! Save planes to disk
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[SavePlanes]')

    call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesI', iplanes)
    call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesJ', jplanes)
    call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesK', kplanes)

! ###################################################################
! Save lines to disk
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[SaveTowers]')
    call TLAB_WRITE_ASCII(bakfile, 'Stride=<value_i,value_j,value_k>')

    call SCANINICHAR(bakfile, inifile, 'SaveTowers', 'Stride', '0,0,0', sRes)
    idummy = 3; call LIST_INTEGER(sRes, idummy, tower_stride)
    if (idummy /= 3) then
        tower_stride(:) = 0
        call TLAB_WRITE_ASCII(bakfile, 'Stride=0,0,0')
        call TLAB_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Cannot read stride for towers; set to 0,0,0.')
    end if

! ###################################################################
! Statistics Control
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Statsitics]')
    call TLAB_WRITE_ASCII(bakfile, '#Averages=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#Pdfs=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ConditionalAverages=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#Intermittency=<yes/no>')

    call SCANINICHAR(bakfile, inifile, 'Statistics', 'Averages', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stats_averages = .true.
    else; stats_averages = .false.; end if

    call SCANINICHAR(bakfile, inifile, 'Statistics', 'Pdfs', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stats_pdfs = .true.
    else; stats_pdfs = .false.; end if

    call SCANINICHAR(bakfile, inifile, 'Statistics', 'Intermittency', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stats_intermittency = .true.
    else; stats_intermittency = .false.; end if

! ###################################################################
! Inflow forcing conditions
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Inflow]')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<None/Discrete/Broadband/Both>')
    call TLAB_WRITE_ASCII(bakfile, '#Adapt=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Imax=<imax>')
    call TLAB_WRITE_ASCII(bakfile, '#Jmax=<jmax>')
    call TLAB_WRITE_ASCII(bakfile, '#Kmax=<kmax>')

    call SCANINICHAR(bakfile, inifile, 'Inflow', 'Type', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; inflow_mode = 0
    elseif (trim(adjustl(sRes)) == 'discrete') then; inflow_mode = 1
    elseif (trim(adjustl(sRes)) == 'broadbandperiodic') then; inflow_mode = 2
    elseif (trim(adjustl(sRes)) == 'broadbandsequential') then; inflow_mode = 3
    elseif (trim(adjustl(sRes)) == 'both') then; inflow_mode = 4
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Inflow.Type.')
        call TLAB_STOP(DNS_ERROR_INFTYPE)
    end if

    call SCANINIREAL(bakfile, inifile, 'Inflow', 'Adapt', '0.0', inflow_adapt)

! Broadband forcing: Grid size of the inflow domain
    g_inf(:)%size = 1       ! default
    g_inf(:)%periodic = g(:)%periodic
    g_inf(:)%uniform = g(:)%uniform
    if (inflow_mode == 2 .or. inflow_mode == 3 .or. inflow_mode == 4) then
        call SCANINIINT(bakfile, inifile, 'Inflow', 'Imax', '0', idummy)
        if (idummy > 0) then; g_inf(1)%size = idummy
        else
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Inflow.Imax.')
            call TLAB_STOP(DNS_ERROR_INFTYPE)
        end if

        call SCANINIINT(bakfile, inifile, 'Inflow', 'Jmax', '0', idummy)
        if (idummy > 0) then; g_inf(2)%size = idummy
        else
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Inflow.Jmax.')
            call TLAB_STOP(DNS_ERROR_INFTYPE)
        end if

        call SCANINIINT(bakfile, inifile, 'Inflow', 'Kmax', '0', idummy)
        if (idummy == g(3)%size) then; g_inf(3)%size = idummy
        else
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Inflow.Kmax.')
            call TLAB_STOP(DNS_ERROR_INFTYPE)
        end if

    end if
    g_inf(1)%inb_grid = g(1)%inb_grid
    g_inf(2)%inb_grid = 1
    g_inf(3)%inb_grid = 1

    if (inflow_mode == 2) then
        g_inf(1)%periodic = .true.
        g_inf(1)%uniform = .true.
    end if

    idummy = max(g_inf(1)%size, max(g_inf(2)%size, g_inf(3)%size))
    isize_wrk1d = max(isize_wrk1d, idummy)

    idummy = max(g_inf(1)%size*g_inf(2)%size, max(g_inf(1)%size*g_inf(3)%size, g_inf(2)%size*g_inf(3)%size))
    isize_wrk2d = max(isize_wrk2d, idummy)

    ! -------------------------------------------------------------------
    ! Discrete Forcing
    ! -------------------------------------------------------------------
    call DISCRETE_READBLOCK(bakfile, inifile, 'Discrete', fp)
    ! Modulation type in fp%type; varicose or sinuous

    ! Parameter 1 is the y-thickness of the perturbation
    ! Parameter 2 is the x-length of the inflow domain

! ###################################################################
! Final initialization and control statements
! ###################################################################
    if (nitera_first > nitera_last) then
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Not started because nitera_first > nitera_last.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! Avoid dividing by zero in time_integration routine
    if (nitera_save <= 0) nitera_save = nitera_last - nitera_first + 1
    if (nitera_stats <= 0) nitera_stats = nitera_last - nitera_first + 1
    if (nitera_log <= 0) nitera_log = nitera_last - nitera_first + 1
    if (nitera_pln <= 0) nitera_pln = nitera_last - nitera_first + 1
    if (nitera_filter <= 0) nitera_filter = nitera_last - nitera_first + 1

    if (imode_sim == DNS_MODE_TEMPORAL) nitera_stats_spa = -1 ! Never call avg_spatial routines
    if (nitera_stats_spa <= 0) nitera_stats_spa = nitera_last - nitera_first + 1

! -------------------------------------------------------------------
! Boundary conditions
! -------------------------------------------------------------------
! Make sure periodic BCs are not modified
    if (g(1)%periodic) then; 
        BcsFlowImin%type(:) = DNS_BCS_NONE; BcsFlowImax%type(:) = DNS_BCS_NONE
        BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
        bcs_visc_imin = DNS_BCS_NONE; bcs_visc_imax = DNS_BCS_NONE
    end if
    if (g(2)%periodic) then; 
        BcsFlowJmin%type(:) = DNS_BCS_NONE; BcsFlowJmax%type(:) = DNS_BCS_NONE
        BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
        bcs_visc_jmin = DNS_BCS_NONE; bcs_visc_jmax = DNS_BCS_NONE
    end if
    if (g(3)%periodic) then; 
        BcsFlowKmin%type(:) = DNS_BCS_NONE; BcsFlowKmax%type(:) = DNS_BCS_NONE
        BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
        bcs_visc_kmin = DNS_BCS_NONE; bcs_visc_kmax = DNS_BCS_NONE
    end if

! BCs for OPR_PARTIAL at xmin (1,*) and xmax (2,*)
    bcs_inf = 0                                    ! default is biased non-zero; if 1, set to zero
    if (bcs_visc_imin == 1) bcs_inf(1, 2, 1) = 1 ! Inflow conditions
    if (bcs_visc_imax == 1) bcs_inf(2, 2, 1) = 1
    if (bcs_visc_jmin == 1) bcs_inf(1, 2, 2) = 1
    if (bcs_visc_jmax == 1) bcs_inf(2, 2, 2) = 1
    if (bcs_visc_kmin == 1) bcs_inf(1, 2, 3) = 1
    if (bcs_visc_kmax == 1) bcs_inf(2, 2, 3) = 1

    bcs_out = 0
    if (bcs_visc_imin == 2) bcs_out(1, 2, 1) = 1 ! Outflow conditions
    if (bcs_visc_imax == 2) bcs_out(2, 2, 1) = 1
    if (bcs_visc_jmin == 2) bcs_out(1, 2, 2) = 1
    if (bcs_visc_jmax == 2) bcs_out(2, 2, 2) = 1
    if (bcs_visc_kmin == 2) bcs_out(1, 2, 3) = 1
    if (bcs_visc_kmax == 2) bcs_out(2, 2, 3) = 1

! Make sure there is array space for reference mean drift
    if (BcsDrift) then
        BuffFlowJmin%size = max(BuffFlowJmin%size, 1); BuffFlowJmax%size = max(BuffFlowJmax%size, 1)
        BuffScalJmin%size = BuffFlowJmin%size; BuffScalJmax%size = BuffFlowJmax%size
        if (imode_sim == DNS_MODE_SPATIAL) then
            BuffFlowImin%size = max(BuffFlowImin%size, 1); BuffFlowImax%size = max(BuffFlowImax%size, 1)
        end if
        BuffScalImin%size = BuffFlowImin%size; BuffScalImax%size = BuffFlowImax%size
    end if

! -------------------------------------------------------------------
! Interactive Boundary conditions
! -------------------------------------------------------------------
    do is = 1, inb_scal
        if (BcsScalJmin%type(is) /= DNS_BCS_DIRICHLET .and. &
            BcsScalJmin%SfcType(is) /= DNS_SFC_STATIC) then
            call TLAB_WRITE_ASCII(efile, &
                                  'DNS_READ_LOCAL. Interactive BC at jmin not implemented for non-Dirichlet BC')
            call TLAB_STOP(DNS_ERROR_JBC)
        end if
        if (BcsScalJmax%type(is) /= DNS_BCS_DIRICHLET .and. &
            BcsScalJmax%SfcType(is) /= DNS_SFC_STATIC) then
            write (*, *) BcsScalJmax%type(is), BcsScalJmax%SfcType(is), BcsScalJmax%cpl(is)
            call TLAB_WRITE_ASCII(efile, &
                                  'DNS_READ_LOCAL. Interactive BC at jmax not implemented for non-Dirichlet BC')
            call TLAB_STOP(DNS_ERROR_JBC)
        end if
    end do

! -------------------------------------------------------------------
! Implicit RKM part
! -------------------------------------------------------------------
    if (rkm_mode == RKM_IMP3_DIFFUSION) then

! Check if Neumann BCs for scalar are present and warn if so
        do is = 1, inb_scal
            if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN .or. &
                BcsScalJmax%type(is) == DNS_BCS_NEUMANN) then
                write (sRes, *) is; sRes = 'DNS_REAL_LOCAL. Scalar'//trim(adjustl(sRes))// &
 ': Finite flux BC not implemented for SEMI-IMPLICITE DIFFUSION'
                call TLAB_WRITE_ASCII(wfile, trim(adjustl(sRes)))
                write (sRes, *) is; sRes = 'DNS_REAL_LOCAL. Scalar'//trim(adjustl(sRes))// &
 ': Setting fluxes at boundary to zero'
                call TLAB_WRITE_ASCII(wfile, trim(adjustl(sRes)))
            end if
        end do

! Check if grid is non-uniform
        if (.not. g(1)%uniform .and. g(1)%mode_fdm1 /= FDM_COM6_DIRECT) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Non-uniform grid requires a direct FDM formulation.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        if (imode_eqns /= DNS_EQNS_INCOMPRESSIBLE) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Implicit formulation only available for incompressible case.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

    end if

! -------------------------------------------------------------------
! Towers information
! So far, the use or not use of tower information is
! set by the stride information.
! -------------------------------------------------------------------
    use_tower = .false.
    if (minval(tower_stride) > 0) then; use_tower = .true.; end if

! We need space
    if (use_tower) then
        idummy = tower_stride(1)*tower_stride(2)*tower_stride(3)
        if (idummy < 5) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Not enough space in wrk3d array to handle tower information. Increase strides.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    end if

! -------------------------------------------------------------------
! Pressure staggering and filtering
! -------------------------------------------------------------------
    if (stagger_on .or. any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        if (.not. (imode_rhs == EQNS_RHS_COMBINED)) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Horizontal pressure staggering or Pressure filter not implemented for this RHS type.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    end if

! -------------------------------------------------------------------
! Nonblocking formulation only valid for 2 scalars or less
! -------------------------------------------------------------------
    if (imode_rhs == EQNS_RHS_NONBLOCKING) then
        if (inb_scal > 2) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Nonblocking Communication not implemented >2 scalars')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        else if (inb_scal < 1) then
            call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Nonblocking Communication require at least 1 scalar')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    end if

    ! -------------------------------------------------------------------
    ! Array sizes
    ! -------------------------------------------------------------------
    inb_txc = 9

    if (imode_sim == DNS_MODE_SPATIAL) then ! because of the statistics
        inb_txc = max(inb_txc, 7)

        if (stats_averages) then
            idummy = max(MA_MOMENTUM_SIZE, MS_SCALAR_SIZE)
            if (mod(nstatavg*jmax*idummy, isize_txc_field) > 0) then
                idummy = nstatavg*jmax*idummy/isize_txc_field + 1
            else
                idummy = nstatavg*jmax*idummy/isize_txc_field
            end if
            inb_txc = max(inb_txc, idummy)
        end if

    end if

#ifdef USE_PSFFT
    if (imode_rhs == EQNS_RHS_NONBLOCKING) inb_txc = max(inb_txc, 15)
#endif

    isize_wrk3d = max(imax, g_inf(1)%size)*max(jmax, g_inf(2)%size)*kmax
    isize_wrk3d = max(isize_wrk3d, isize_txc_field)
    if (part%type /= PART_TYPE_NONE) then
        isize_wrk3d = max(isize_wrk3d, (imax + 1)*jmax*(kmax + 1))
        isize_wrk3d = max(isize_wrk3d, (jmax*(kmax + 1)*inb_part_interp*2))
        isize_wrk3d = max(isize_wrk3d, (jmax*(imax + 1)*inb_part_interp*2))
    end if
    if (use_tower) then
        isize_wrk3d = max(isize_wrk3d, nitera_save*(g(2)%size + 2))
    end if

    return
end subroutine DNS_READ_LOCAL

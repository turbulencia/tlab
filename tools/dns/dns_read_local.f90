#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2009/11/24 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Extracted from DNS_READ_GLOBAL
!#
!########################################################################
SUBROUTINE DNS_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile, wfile
  USE DNS_GLOBAL,    ONLY : p_init, mean_rho
  USE DNS_GLOBAL,    ONLY : imode_sim, icalc_scal, inb_flow, inb_scal
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, iunify, imode_fdm, imode_eqns
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 bakfile, lstr
  TINTEGER is,idummy,inb_scal_local1

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL IO_WRITE_ASCII(lfile, 'Reading local input data.')

! ###################################################################
! Main Section
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Main]')
  CALL IO_WRITE_ASCII(bakfile, '#TimeOrder=<RungeKuttaExplicit3/RungeKuttaExplicit4/RungeKuttaDiffusion3>')
  CALL IO_WRITE_ASCII(bakfile, '#TimeStep=<value (used if CFL is negative)>')
  CALL IO_WRITE_ASCII(bakfile, '#TimeCFL=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#TermDivergence=<none/remove>')
  CALL IO_WRITE_ASCII(bakfile, '#RhsMode=<split/combined/nonblocking>') 

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TimeOrder', 'dummy', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttaexplicit3'  ) THEN; rkm_mode = RKM_EXP3;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttaexplicit4'  ) THEN; rkm_mode = RKM_EXP4;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttadiffusion3' ) THEN; rkm_mode = RKM_IMP3_DIFFUSION;
!  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttasource3'    ) THEN; rkm_mode = RKM_IMP3_SOURCE;
  ELSE                                                         ! Old format
     CALL IO_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. TimeOrder obsolete.')

     CALL SCANINIINT(bakfile,  inifile, 'Main', 'TimeOrder', '4', rkm_mode)
     IF ( rkm_mode .LT. 3 .OR. rkm_mode .GT. 4 ) THEN
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Runge-Kutta Order equal to 3 or 4')
        CALL DNS_STOP(DNS_ERROR_RKORDER)
     ENDIF
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Main', 'TimeCFL',  '0.75', cfl  )
  CALL SCANINIREAL(bakfile, inifile, 'Main', 'TimeStep', '0.05', dtime)

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermDivergence', 'remove', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'   ) THEN; idivergence = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'remove' ) THEN; idivergence = EQNS_DIVERGENCE
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong TermDivergence option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'RhsMode', 'combined', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'split'       ) THEN; imode_rhs = EQNS_RHS_SPLIT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'combined'    ) THEN; imode_rhs = EQNS_RHS_COMBINED
#ifdef USE_PSFFT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'nonblocking' ) THEN; imode_rhs = EQNS_RHS_NONBLOCKING
#endif
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong RhsMode option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! ###################################################################
! Iteration Section
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Iteration]')
  CALL IO_WRITE_ASCII(bakfile, '#Start=<integral start time>')
  CALL IO_WRITE_ASCII(bakfile, '#End=<integral stop time>')
  CALL IO_WRITE_ASCII(bakfile, '#Restart=<restart time step>')
  CALL IO_WRITE_ASCII(bakfile, '#Statistics=<statistics time step>')
  CALL IO_WRITE_ASCII(bakfile, '#IteraLog=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Saveplanes=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#RunAvera=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#RunLines=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#RunPlane=<yes/no>')

  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start',      '0',  nitera_first)
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'End',        '0',  nitera_last )
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Restart',    '50', nitera_save )
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Statistics', '50', nitera_stats)
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'IteraLog',   '10', nitera_log  )
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Saveplanes', '-1', nitera_pln  )

  CALL SCANINICHAR(bakfile, inifile, 'Iteration', 'RunAvera', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; frunstat=1
  ELSE;                                       frunstat=0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Iteration', 'RunLines', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; frunline=1
  ELSE;                                       frunline=0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Iteration', 'RunPlane', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; frunplane=1
  ELSE;                                       frunplane=0; ENDIF

! ###################################################################
! Control Limits
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Control]')
  CALL IO_WRITE_ASCII(bakfile, '#FlowLimit=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#MinPressure=<pressure>')
  CALL IO_WRITE_ASCII(bakfile, '#MaxPressure=<pressure>')
  CALL IO_WRITE_ASCII(bakfile, '#MinDensity=<density>')
  CALL IO_WRITE_ASCII(bakfile, '#MaxDensity=<density>')
  CALL IO_WRITE_ASCII(bakfile, '#ScalLimit=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#MinScalar=<scalar>')
  CALL IO_WRITE_ASCII(bakfile, '#MaxScalar=<scalar>')

  CALL SCANINICHAR(bakfile, inifile, 'Control', 'FlowLimit', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; ilimit_flow=1
  ELSE;                                       ilimit_flow=0; ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'Control', 'ScalLimit', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; ilimit_scal=1
  ELSE;                                       ilimit_scal=0; ENDIF

! Final check in last section of the routine
  CALL SCANINIREAL(bakfile, inifile, 'Control', 'MinPressure', '-1.0', p_bound_min)
  CALL SCANINIREAL(bakfile, inifile, 'Control', 'MaxPressure', '-1.0', p_bound_max)
  CALL SCANINIREAL(bakfile, inifile, 'Control', 'MinDensity',  '-1.0', r_bound_min)
  CALL SCANINIREAL(bakfile, inifile, 'Control', 'MaxDensity',  '-1.0', r_bound_max)

  d_bound_max = C_BIG_R ! default
  CALL SCANINICHAR(bakfile, inifile, 'Control', 'MaxDilatation', 'void', sRes)
  IF ( sRes .NE. 'void' ) THEN
     idummy = 1
     CALL LIST_REAL(sRes, idummy, d_bound_max)
  ENDIF

  s_bound_min(:) = C_0_R; inb_scal_local1 = MAX_NSP
  IF ( ilimit_scal .EQ. 1 ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'Control', 'MinScalar',  'void', sRes)
     IF ( sRes .NE. 'void' ) THEN
        CALL LIST_REAL(sRes, inb_scal_local1, s_bound_min)
        IF ( inb_scal_local1 .NE. inb_scal ) THEN ! Consistency check
           CALL IO_WRITE_ASCII(efile,'DNS_READ_LOCAL. MinScalar size does not match inb_scal.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
     ENDIF
  ENDIF
     
  s_bound_max(:) = C_1_R; inb_scal_local1 = MAX_NSP
  IF ( ilimit_scal .EQ. 1 ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'Control', 'MaxScalar',  'void', sRes)
     IF ( sRes .NE. 'void' ) THEN
        CALL LIST_REAL(sRes, inb_scal_local1, s_bound_max)
        IF ( inb_scal_local1 .NE. inb_scal ) THEN ! Consistency check
           CALL IO_WRITE_ASCII(efile,'DNS_READ_LOCAL. MaxScalar size does not match inb_scal.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
     ENDIF
  ENDIF

! ###################################################################
! Boundary Conditions
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[BoundaryConditions]')
  CALL IO_WRITE_ASCII(bakfile, '#EulerI=<none/nonreflective/inflow/freeslip/noslip>')
  CALL IO_WRITE_ASCII(bakfile, '#EulerJ=<none/nonreflective/inflow/freeslip/noslip>')
  CALL IO_WRITE_ASCII(bakfile, '#EulerK=<none/nonreflective/inflow/freeslip/noslip>')
  CALL IO_WRITE_ASCII(bakfile, '#ViscousI=<none/inflow/outflow>')
  CALL IO_WRITE_ASCII(bakfile, '#ViscousJ=<none/inflow/outflow>')
  CALL IO_WRITE_ASCII(bakfile, '#ViscousK=<none/inflow/outflow>')
  CALL IO_WRITE_ASCII(bakfile, '#SigmaOut=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#SigmaInfImin=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#SigmaInfImax=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#SigmaInfJ=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#BetaTransverse=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ScalarImin=<none/dirichlet/neumman>')
  CALL IO_WRITE_ASCII(bakfile, '#ScalarJmin=<none/dirichlet/neumman>')
  CALL IO_WRITE_ASCII(bakfile, '#ScalarKmin=<none/dirichlet/neumman>')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityImin=<none/dirichlet/neumman>')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityJmin=<none/dirichlet/neumman>')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityKmin=<none/dirichlet/neumman>')

! -------------------------------------------------------------------
! Euler terms
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'EulerI', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; bcs_euler_imin = DNS_BCS_NONE;      bcs_euler_imax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'nonreflective' ) THEN; bcs_euler_imin = DNS_BCS_NR;        bcs_euler_imax = DNS_BCS_NR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'        ) THEN; bcs_euler_imin = DNS_BCS_INFLOW;    bcs_euler_imax = DNS_BCS_NR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip'      ) THEN; bcs_euler_imin = DNS_BCS_NEUMANN;   bcs_euler_imax = DNS_BCS_NEUMANN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'        ) THEN; bcs_euler_imin = DNS_BCS_DIRICHLET; bcs_euler_imax = DNS_BCS_DIRICHLET
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.EulerI.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'EulerJ', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; bcs_euler_jmin = DNS_BCS_NONE;      bcs_euler_jmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'nonreflective' ) THEN; bcs_euler_jmin = DNS_BCS_NR;        bcs_euler_jmax = DNS_BCS_NR 
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'        ) THEN; bcs_euler_jmin = DNS_BCS_INFLOW;    bcs_euler_jmax = DNS_BCS_NR 
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip'      ) THEN; bcs_euler_jmin = DNS_BCS_NEUMANN;   bcs_euler_jmax = DNS_BCS_NEUMANN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'        ) THEN; bcs_euler_jmin = DNS_BCS_DIRICHLET; bcs_euler_jmax = DNS_BCS_DIRICHLET
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.EulerJ.')
     CALL DNS_STOP(DNS_ERROR_JBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'EulerK', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; bcs_euler_kmin = DNS_BCS_NONE;      bcs_euler_kmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'nonreflective' ) THEN; bcs_euler_kmin = DNS_BCS_NR;        bcs_euler_kmax = DNS_BCS_NR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'        ) THEN; bcs_euler_kmin = DNS_BCS_INFLOW;    bcs_euler_kmax = DNS_BCS_NR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip'      ) THEN; bcs_euler_kmin = DNS_BCS_NEUMANN;   bcs_euler_kmax = DNS_BCS_NEUMANN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'        ) THEN; bcs_euler_kmin = DNS_BCS_DIRICHLET; bcs_euler_kmax = DNS_BCS_DIRICHLET
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.EulerK.')
     CALL DNS_STOP(DNS_ERROR_KBC)
  ENDIF

! -------------------------------------------------------------------
! Viscous terms
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousI', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'    ) THEN; bcs_visc_imin = 0;  bcs_visc_imax = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'  ) THEN; bcs_visc_imin = 1;  bcs_visc_imax = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'outflow' ) THEN; bcs_visc_imin = 2;  bcs_visc_imax = 2
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousI.')
     CALL DNS_STOP(DNS_ERROR_IVSICBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousJ', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'    ) THEN; bcs_visc_jmin = 0;  bcs_visc_jmax = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'  ) THEN; bcs_visc_jmin = 1;  bcs_visc_jmax = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'outflow' ) THEN; bcs_visc_jmin = 2;  bcs_visc_jmax = 2
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousJ.')
     CALL DNS_STOP(DNS_ERROR_JVSICBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'ViscousK', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'    ) THEN; bcs_visc_kmin = 0;  bcs_visc_kmax = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'inflow'  ) THEN; bcs_visc_kmin = 1;  bcs_visc_kmax = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'outflow' ) THEN; bcs_visc_kmin = 2;  bcs_visc_kmax = 2
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.ViscousK.')
     CALL DNS_STOP(DNS_ERROR_KVSICBC)
  ENDIF

! -------------------------------------------------------------------
! Scalar terms
! -------------------------------------------------------------------
  bcs_scal_imin(:) = DNS_BCS_NONE; bcs_scal_imax(:) = DNS_BCS_NONE
  IF ( i1bc .NE. DNS_BCS_PERIODIC ) THEN
  DO is = 1,inb_scal
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Imin'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'none', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_imin(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_imin(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_imin(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Imax'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'none', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_imax(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_imax(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_imax(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
  ENDDO
  ENDIF

  bcs_scal_jmin(:) = DNS_BCS_NONE; bcs_scal_jmax(:) = DNS_BCS_NONE
  IF ( j1bc .NE. DNS_BCS_PERIODIC ) THEN
  DO is = 1,inb_scal
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Jmin'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'void', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_jmin(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_jmin(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_jmin(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Jmax'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'void', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_jmax(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_jmax(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_jmax(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
  ENDDO
  ENDIF

  bcs_scal_kmin(:) = DNS_BCS_NONE; bcs_scal_kmax(:) = DNS_BCS_NONE
  IF ( k1bc .NE. DNS_BCS_PERIODIC ) THEN
  DO is = 1,inb_scal
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Kmin'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'none', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_kmin(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_kmin(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_kmin(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
     WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Kmax'
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'none', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; bcs_scal_kmax(is) = DNS_BCS_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; bcs_scal_kmax(is) = DNS_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; bcs_scal_kmax(is) = DNS_BCS_NEUMANN
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
        CALL DNS_STOP(DNS_ERROR_IBC)
     ENDIF
  ENDDO
  ENDIF

! -------------------------------------------------------------------
! Velocity terms
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityImin', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_imin = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_imin = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_imin = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityImin.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityImax', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_imax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_imax = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_imax = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityImax.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_jmin = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_jmin = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_jmin = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityJmin.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmax', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_jmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_jmax = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_jmax = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityJmax.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityKmin', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_kmin = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_kmin = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_kmin = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityKmin.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityKmax', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_kmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_kmax = DNS_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_kmax = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. BoundaryConditions.VelocityKmax.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

! -------------------------------------------------------------------
! Reference values in characteristic formulation
! -------------------------------------------------------------------
  bcs_euler_drift = 0
  CALL SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaOut', '-1.0', bcs_sigma_out)
  IF ( bcs_sigma_out .LE. C_0_R    ) THEN; bcs_sigma_out = C_0_R
  ELSE;                                    bcs_euler_drift = MAX(bcs_euler_drift,i1); ENDIF 
! we need it for the bcs array (see end of routine)

  CALL SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaInfJ', '-1.0', bcs_sigma_inf_j)
  IF ( bcs_sigma_inf_j .LE. C_0_R   ) THEN; bcs_sigma_inf_j = C_0_R
  ELSE;                                     bcs_euler_drift = MAX(bcs_euler_drift,i1); ENDIF 
! we need it for the bcs array (see end of routine)

! implemented s.t. SigmaInfImax cannot be used w/o SigmaInfImin
  CALL SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaInfImin', '-1.0', bcs_sigma_inf_imin)
  IF ( bcs_sigma_inf_imin .LE. C_0_R  ) THEN; bcs_sigma_inf_imin = C_0_R
  ELSE;                                       bcs_euler_drift = MAX(bcs_euler_drift,i1); ENDIF 
! we need it for the bcs array (see end of routine)

  CALL SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'SigmaInfImax', '-1.0', bcs_sigma_inf_imax)
  IF ( bcs_sigma_inf_imax .LE. C_0_R ) THEN; bcs_sigma_inf_imax = C_0_R
  ELSE;                                      bcs_euler_drift = MAX(bcs_euler_drift,i1); ENDIF 
! we need it for the bcs array (see end of routine)

! Relaxation coefficient for the transverse terms
  CALL SCANINIREAL(bakfile, inifile, 'BoundaryConditions', 'BetaTransverse', '-1.0', bcs_sigma_trans)
  IF ( bcs_sigma_trans .LE. C_0_R ) bcs_sigma_trans = C_0_R

! ###################################################################
! Buffer Zone Parameters
! inb_scal is needed
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[BufferZone]')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<none/relaxation/filter/both>')
  CALL IO_WRITE_ASCII(bakfile, '#LoadBuffer=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsUJmin=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsUJmax=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsEJmin=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsEJmax=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsSJmin=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsSJmax=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsImin=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PointsImax=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ParametersU=<values>')
  CALL IO_WRITE_ASCII(bakfile, '#ParametersS=<values>')
  CALL IO_WRITE_ASCII(bakfile, '#VFree=<value>')

  CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'Type', 'none', sRes)

  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'       ) THEN; buff_type = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'relaxation' ) THEN; buff_type = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'filter'     ) THEN; buff_type = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'both'       ) THEN; buff_type = 3
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong BufferType option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! load buffer if used also by BCs
  CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'LoadBuffer', 'no', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_load = 1
  ELSE;                                       buff_load = 0; ENDIF

  IF ( buff_type .GT. 0 ) THEN
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsImin',  '0',  buff_nps_imin)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsImax',  '0',  buff_nps_imax)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUJmin', '0',  buff_nps_u_jmin)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsUJmax', '0',  buff_nps_u_jmax)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsEJmin', '0',  buff_nps_e_jmin)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsEJmax', '0',  buff_nps_e_jmax)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSJmin','-1',  buff_nps_s_jmin)
     CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'PointsSJmax','-1',  buff_nps_s_jmax)
     ! control 
     IF ( buff_nps_s_jmin .LT. 0 ) THEN
        buff_nps_s_jmin = buff_nps_u_jmin
        CALL IO_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Field PointsSJmin default to PointsUJmin.') 
     ENDIF
     IF ( buff_nps_s_jmax .LT. 0 ) THEN
        buff_nps_s_jmax = buff_nps_u_jmax
        CALL IO_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Field PointsSJmax default to PointsUJmax.') 
     ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'ParametersU', '1.0,2.0',sRes)
     is = 2; CALL LIST_REAL(sRes, is, buff_param_u)
     IF ( is .EQ. 1 ) buff_param_u(2) = C_2_R ! force default if not given
     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'ParametersS', 'dummy',  sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'dummy' ) THEN
        buff_param_s(:) = buff_param_u(:)
        CALL IO_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Field ParametersS default to ParametersU.') 
     ELSE
        is = 2; CALL LIST_REAL(sRes, is, buff_param_u)
     ENDIF
     ! control
     CALL SCANINICHAR(bakfile,inifile, 'BufferZone', 'Sigma','void',sRes)
     IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
        CALL IO_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Field BufferZone/Sigma obsolete. Check manual.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'VFree', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_v_free = 1
     ELSE;                                       buff_v_free = 0; ENDIF

! -------------------------------------------------------------------
! Fixed values
! -------------------------------------------------------------------
     CALL IO_WRITE_ASCII(bakfile, '#HardU=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#U(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#U(T)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#HardV=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#V(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#V(T)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#HardW=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#W(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#W(T)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#HardE=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#E(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#E(T)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#HardRho=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#Rho(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#Rho(T)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#HardZ=<yes/no>')
     CALL IO_WRITE_ASCII(bakfile, '#Z(B)=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#Z(T)=<value>')

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardU', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(1) = 1
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'U(T)', '0.0', buff_hard(1,1))
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'U(B)', '0.0', buff_hard(1,2))
     ELSE;                                       buff_hard_on(1) = 0; ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardV', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(2) = 1
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'V(T)', '0.0', buff_hard(2,1))
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'V(B)', '0.0', buff_hard(2,2))
     ELSE;                                       buff_hard_on(2) = 0; ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardW', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(3) = 1
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'W(T)', '0.0', buff_hard(3,1))
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'W(B)', '0.0', buff_hard(3,2))
     ELSE;                                       buff_hard_on(3) = 0; ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardE', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(4) = 1
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'E(T)', '0.0', buff_hard(4,1))
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'E(B)', '0.0', buff_hard(4,2))
     ELSE;                                       buff_hard_on(4) = 0; ENDIF

     CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardRho', 'no', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(5) = 1
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'Rho(T)', '1.0', buff_hard(5,1))
        CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'Rho(B)', '1.0', buff_hard(5,2))
     ELSE;                                       buff_hard_on(5) = 0; ENDIF

     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal
           IF ( inb_scal .EQ. 1 ) THEN; lstr = 'HardZ'
           ELSE;                        WRITE(lstr,'(A5,I1)') 'HardZ', is; ENDIF
           CALL SCANINICHAR(bakfile, inifile, 'BufferZone', TRIM(ADJUSTL(lstr)), 'no', sRes)
           IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; buff_hard_on(inb_flow+is) = 1
              CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'Z(T)', '0.0', buff_hard(inb_flow+is,1))
              CALL SCANINIREAL(bakfile, inifile, 'BufferZone', 'Z(B)', '0.0', buff_hard(inb_flow+is,2))
           ELSE;                                       buff_hard_on(inb_flow+is) = 0; ENDIF
        ENDDO
     ENDIF

  ELSE
     buff_nps_imin   = 0; buff_nps_imax   = 0
     buff_nps_u_jmin = 0; buff_nps_u_jmax = 0
     buff_nps_e_jmin = 0; buff_nps_e_jmax = 0

  ENDIF

! ###################################################################
! Viscosity Control
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[ViscChange]')
  CALL IO_WRITE_ASCII(bakfile, '#Time=<time>')

  CALL SCANINIREAL(bakfile, inifile, 'ViscChange', 'Time', '0.0', visctime)

! ###################################################################
! Filter
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Filter]')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<none/compact/explicit6/adm>')
  CALL IO_WRITE_ASCII(bakfile, '#Alpha=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Step=<filter step>')
  CALL IO_WRITE_ASCII(bakfile, '#Scalar=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#Ifilt=<set filter in i direction>')
  CALL IO_WRITE_ASCII(bakfile, '#Jfilt=<set filter in j direction>')
  CALL IO_WRITE_ASCII(bakfile, '#Kfilt=<set filter in k direction>')

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'Type', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; ifilt_domain = DNS_FILTER_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'compact'   ) THEN; ifilt_domain = DNS_FILTER_COMPACT
     CALL SCANINIREAL(bakfile, inifile, 'Filter', 'Alpha', '0.49', ifilt_alpha)
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit6' ) THEN; ifilt_domain = DNS_FILTER_6E  
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit4' ) THEN; ifilt_domain = DNS_FILTER_4E  
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'adm'       ) THEN; ifilt_domain = DNS_FILTER_ADM;     ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Filter', 'Step', '0', ifilt_step)
  IF ( ifilt_step   .EQ. 0               ) ifilt_domain = DNS_FILTER_NONE
  IF ( ifilt_domain .EQ. DNS_FILTER_NONE ) ifilt_step   = 0 
  
  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'Scalar', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; ifilt_scalar = 1
  ELSE;                                       ifilt_scalar = 0; ENDIF

! active/no active
  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveX', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'no' ) THEN; ifilt_x = 0
  ELSE;                                      ifilt_x = 1; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveY', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'no' ) THEN; ifilt_y = 0
  ELSE;                                      ifilt_y = 1; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveZ', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'no' ) THEN; ifilt_z = 0
  ELSE;                                      ifilt_z = 1; ENDIF

! -------------------------------------------------------------------
! Inflow Filter
! -------------------------------------------------------------------
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[InflowFilter]')
  CALL IO_WRITE_ASCII(bakfile, '#Active=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#IWidth=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#JWidth=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Step=<value>')

  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
  CALL SCANINICHAR(bakfile, inifile, 'InflowFilter', 'Active', 'no', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; ifilt_inflow = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; ifilt_inflow = 1
  ELSE
     CALL SCANINIINT(bakfile, inifile, 'InflowFilter', 'Active', '0', ifilt_inflow)
  ENDIF

  IF ( ifilt_inflow .EQ. 1 ) THEN
     CALL SCANINIINT(bakfile, inifile, 'InflowFilter', 'IWidth', '1', ifilt_inflow_iwidth)
     
     IF ( ifilt_inflow_iwidth .GT. imax ) THEN
        CALL IO_WRITE_ASCII(efile, 'Error: Inflow filter i width larger than imax')
        CALL DNS_STOP(DNS_ERROR_INFFLTDOM)
     ENDIF
        
     CALL SCANINIINT(bakfile, inifile, 'InflowFilter', 'JWidth', '1', ifilt_inflow_jwidth)
     
     IF ( ifilt_inflow_jwidth .GT. jmax ) THEN
        CALL IO_WRITE_ASCII(efile, 'Error: Inflow filter j width larger than jmax')
        CALL DNS_STOP(DNS_ERROR_INFFLTDOM)
     ENDIF
     
     CALL SCANINIINT(bakfile, inifile, 'InflowFilter', 'Step', '1', ifilt_inflow_step)
  ENDIF

  ELSE
     ifilt_inflow = 0

  ENDIF

! ###################################################################
! Save planes to disk
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[SavePlanes]')
  CALL IO_WRITE_ASCII(bakfile, '#PlanesI=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PlanesJ=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#PlanesK=<value>')

  CALL SCANINICHAR(bakfile, inifile, 'SavePlanes', 'PlanesI', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void'  ) THEN
     nplanes_i = 0; planes_i = 0
  ELSE 
     nplanes_i = MAX_SAVEPLANES; CALL LIST_INTEGER(sRes, nplanes_i, planes_i)
  ENDIF
  
  CALL SCANINICHAR(bakfile, inifile, 'SavePlanes', 'PlanesJ', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void'  ) THEN
     nplanes_j = 0; planes_j = 0
  ELSE 
     nplanes_j = MAX_SAVEPLANES; CALL LIST_INTEGER(sRes, nplanes_j, planes_j)
  ENDIF
  
  CALL SCANINICHAR(bakfile, inifile, 'SavePlanes', 'PlanesK', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void'  ) THEN
     nplanes_k = 0; planes_k = 0
  ELSE 
     nplanes_k = MAX_SAVEPLANES; CALL LIST_INTEGER(sRes, nplanes_k, planes_k)
  ENDIF

! ###################################################################
! Save lines to disk
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[SaveTowers]')
  CALL IO_WRITE_ASCII(bakfile, 'Stride=<value_i,value_j,value_k>')
  
  CALL SCANINICHAR(bakfile, inifile, 'SaveTowers', 'Stride', '0,0,0', sRes)
  idummy = 3; CALL LIST_INTEGER(sRes,idummy,tower_stride)  
  IF ( idummy .NE. 3 ) THEN 
     tower_stride(:) = 0  
     CALL IO_WRITE_ASCII(bakfile, 'Stride=0,0,0')
     CALL IO_WRITE_ASCII(wfile,   'DNS_READ_LOCAL. Cannot read stride for towers; set to 0,0,0.')
  ENDIF

! ###################################################################
! Statistics Control   
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Statsitics]')
  CALL IO_WRITE_ASCII(bakfile, '#Averages=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#Pdfs=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#ConditionalAverages=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#Intermittency=<yes/no>')

  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'Averages', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; fstavg = 1
  ELSE;                                       fstavg = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'Pdfs', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; fstpdf = 1
  ELSE;                                       fstpdf = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'Intermittency', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; fstinter = 1
  ELSE;                                       fstinter = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'FilterEnergy', 'no', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; ffltdmp = 1
  ELSE;                                       ffltdmp = 0; ENDIF

! ###################################################################
! Inflow forcing conditions
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Inflow]')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<None/Discrete/Broadband/Both>')
  CALL IO_WRITE_ASCII(bakfile, '#Adapt=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Imax=<imax>')
  CALL IO_WRITE_ASCII(bakfile, '#Jmax=<jmax>')
  CALL IO_WRITE_ASCII(bakfile, '#Kmax=<kmax>')

  CALL SCANINICHAR(bakfile, inifile, 'Inflow', 'Type', 'None', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'none'                ) THEN; ifrc_mode = 0
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'discrete'            ) THEN; ifrc_mode = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'broadbandperiodic'   ) THEN; ifrc_mode = 2
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'broadbandsequential' ) THEN; ifrc_mode = 3
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'both'                ) THEN; ifrc_mode = 4
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Inflow.Type.')
     CALL DNS_STOP(DNS_ERROR_INFTYPE)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Inflow', 'Adapt', '0.0', frc_adapt)

! Broadband forcing: Grid size of the inflow domain
  IF ( ifrc_mode .EQ. 2 .OR. ifrc_mode .EQ. 3 .OR. ifrc_mode .EQ. 4 ) THEN                  
     CALL SCANINIINT(bakfile, inifile, 'Inflow', 'Imax', '0', imax_inf)
     CALL SCANINIINT(bakfile, inifile, 'Inflow', 'Jmax', '0', jmax_inf)
     CALL SCANINIINT(bakfile, inifile, 'Inflow', 'Kmax', '0', kmax_inf)

     IF ( kmax_inf .NE. kmax_total ) THEN
        CALL IO_WRITE_ASCII(efile,  'DNS_READ_LOCAL. Inflow KMAX missmatch.')
        CALL DNS_STOP(DNS_ERROR_KMAXMISS)
     ENDIF

     kmax_inf = kmax

  ELSE
     imax_inf = 1
     jmax_inf = 1
     kmax_inf = 1
  ENDIF

! ###################################################################
! Discrete Forcing
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Discrete]')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian>')
  CALL IO_WRITE_ASCII(bakfile, '#2DAmpl=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DAmpl=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#2DPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DXPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DZPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Broadening=<value>')

  CALL SCANINICHAR(bakfile,inifile,'Discrete','2DPhi', '0.0',sRes)
  Phix2D(:)=C_0_R; nx2d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nx2d, Phix2D)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','2DAmpl','0.0',sRes)
  A2D(:)=C_0_R; nx2d = MAX_FRC_FREC ! The amplitude sets the value of nx2d
  CALL LIST_REAL(sRes, nx2d, A2D)
  
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DXPhi','0.0',sRes)
  Phix3D(:)=C_0_R; nx3d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nx3d, Phix3d)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DZPhi','0.0',sRes)
  Phiz3D(:)=C_0_R; nz3d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nz3d, Phiz3D)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DAmpl','0.0',sRes)
  A3D(:)=C_0_R; nx3d = MAX_FRC_FREC ! The amplitude sets the value of nx3d
  CALL LIST_REAL(sRes, nx3d, A3D)


  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Type', 'Varicose', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'varicose' ) THEN; ifrcdsc_mode = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'sinuous'  ) THEN; ifrcdsc_mode = 2
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; ifrcdsc_mode = 3
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Error in Discrete.Type.')
     CALL DNS_STOP(DNS_ERROR_INFDISCR)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'XLength',     '1.0', frc_length)
  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', frc_delta)

! ###################################################################
! Final initialization and control statements
! ###################################################################
  IF ( nitera_first .GE. nitera_last ) THEN
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Not started because nitera_first >= nitera_last.' )
     CALL DNS_STOP(DNS_ERROR_OPTION)
  END IF

! Avoid dividing by zero in time_integration routine
  IF ( ifilt_step   .LE. 0 ) ifilt_step   = nitera_last - nitera_first + 1
  IF ( nitera_save  .LE. 0 ) nitera_save  = nitera_last - nitera_first + 1
  IF ( nitera_stats .LE. 0 ) nitera_stats = nitera_last - nitera_first + 1
  IF ( nitera_log   .LE. 0 ) nitera_log   = nitera_last - nitera_first + 1
  IF ( nitera_pln   .LE. 0 ) nitera_pln   = nitera_last - nitera_first + 1

! -------------------------------------------------------------------
! Control limits
! I need mean_rho
! -------------------------------------------------------------------
  IF ( p_bound_min .LT. C_0_R ) p_bound_min = p_init*C_1EM6_R
  IF ( p_bound_max .LT. C_0_R ) p_bound_max = p_init/C_1EM6_R
  IF ( r_bound_min .LT. C_0_R ) r_bound_min = mean_rho*C_1EM6_R
  IF ( r_bound_max .LT. C_0_R ) r_bound_max = mean_rho/C_1EM6_R

! -------------------------------------------------------------------
! Buffer zone
! -------------------------------------------------------------------
! Maximum size
  buff_nps_jmin = MAX(buff_nps_u_jmin,buff_nps_e_jmin)
  buff_nps_jmin = MAX(buff_nps_jmin,  buff_nps_s_jmin)
  buff_nps_jmax = MAX(buff_nps_u_jmax,buff_nps_e_jmax)
  buff_nps_jmax = MAX(buff_nps_jmax,  buff_nps_s_jmax)

! Pointers
  buff_imin   = buff_nps_imin;   buff_imax   = imax - buff_nps_imax    + 1
  buff_u_jmin = buff_nps_u_jmin; buff_u_jmax = jmax - buff_nps_u_jmax  + 1
  buff_e_jmin = buff_nps_e_jmin; buff_e_jmax = jmax - buff_nps_e_jmax  + 1
  buff_s_jmin = buff_nps_s_jmin; buff_s_jmax = jmax - buff_nps_s_jmax  + 1

! -------------------------------------------------------------------
! Boundary conditions
! -------------------------------------------------------------------
! Make sure periodic BCs are not modified
  IF ( i1bc .EQ. DNS_BCS_PERIODIC ) THEN;
     bcs_euler_imin    = DNS_BCS_NONE; bcs_euler_imax    = DNS_BCS_NONE
     bcs_visc_imin     = DNS_BCS_NONE; bcs_visc_imax     = DNS_BCS_NONE
     bcs_flow_imin     = DNS_BCS_NONE; bcs_flow_imax     = DNS_BCS_NONE
     bcs_scal_imin(:)  = DNS_BCS_NONE; bcs_scal_imax(:)  = DNS_BCS_NONE
  ENDIF
  IF ( j1bc .EQ. DNS_BCS_PERIODIC ) THEN;
     bcs_euler_jmin    = DNS_BCS_NONE; bcs_euler_jmax    = DNS_BCS_NONE
     bcs_visc_jmin     = DNS_BCS_NONE; bcs_visc_jmax     = DNS_BCS_NONE
     bcs_flow_jmin     = DNS_BCS_NONE; bcs_flow_jmax     = DNS_BCS_NONE
     bcs_scal_jmin(:)  = DNS_BCS_NONE; bcs_scal_jmax(:)  = DNS_BCS_NONE
  ENDIF
  IF ( k1bc .EQ. DNS_BCS_PERIODIC ) THEN;
     bcs_euler_kmin    = DNS_BCS_NONE; bcs_euler_kmax    = DNS_BCS_NONE
     bcs_visc_kmin     = DNS_BCS_NONE; bcs_visc_kmax     = DNS_BCS_NONE
     bcs_flow_kmin     = DNS_BCS_NONE; bcs_flow_kmax     = DNS_BCS_NONE
     bcs_scal_kmin(:)  = DNS_BCS_NONE; bcs_scal_kmax(:)  = DNS_BCS_NONE
  ENDIF

! Make sure there is array space for reference mean drift
  IF ( bcs_euler_drift .EQ. 1 ) THEN
     IF      ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
        buff_nps_jmin = MAX(buff_nps_jmin,i1); buff_nps_jmax = MAX(buff_nps_jmax,i1)
     ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL  ) THEN
        buff_nps_jmin = MAX(buff_nps_jmin,i1); buff_nps_jmax = MAX(buff_nps_jmax,i1)
        buff_nps_imin = MAX(buff_nps_imin,i1); buff_nps_imax = MAX(buff_nps_imax,i1)
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Implicit RKM part
! -------------------------------------------------------------------
  IF ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) THEN 

! Check if Neumann BCs for scalar are present and warn if so 
     DO is=1,inb_scal
        IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN .OR. &
             bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) THEN 
           WRITE(sRes, *) is; sRes='DNS_REAL_LOCAL. Scalar'//TRIM(ADJUSTL(sRes))//&
                ': Finite flux BC not implemented for SEMI-IMPLICITE DIFFUSION' 
           CALL IO_WRITE_ASCII(wfile, TRIM(ADJUSTL(sRes))) 
           WRITE(sRes, *) is; sRes='DNS_REAL_LOCAL. Scalar'//TRIM(ADJUSTL(sRes))//&
                ': Setting fluxes at boundary to zero' 
           CALL IO_WRITE_ASCII(wfile, TRIM(ADJUSTL(sRes)))
        ENDIF
     ENDDO

! Check if grid is non-uniform
     IF ( iunify .EQ. 1 .AND. imode_fdm .NE. FDM_COM6_DIRECT ) THEN
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Non-uniform grid requires a direct FDM formulation.') 
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
     ENDIF

  ENDIF

! -------------------------------------------------------------------
! Towers information
! So far, the use or not use of tower information (tower_mode) is 
! set by the stride information.
! -------------------------------------------------------------------
  IF ( MINVAL(tower_stride).GT.0 ) THEN; tower_mode=1; ELSE;  tower_mode=0; ENDIF

! We need space
  IF ( tower_mode .EQ. 1 ) THEN
     idummy = tower_stride(1) *tower_stride(2) *tower_stride(3)
     IF ( idummy .LT. 5 ) THEN
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Not enough space in wrk3d array to handle tower information. Increase strides.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Nonblocking formulation only valid for 2 scalars or less
! -------------------------------------------------------------------
  IF ( imode_rhs .EQ. EQNS_RHS_NONBLOCKING .AND. inb_scal .NE. 2 ) THEN
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Nonblocking formulation only valid for 2 scalars.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
  ENDIF

! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN 
     IF ( imixture .EQ. MIXT_TYPE_SUPSAT .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        IF ( imode_rhs .NE. EQNS_RHS_SPLIT ) THEN
           CALL IO_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Airwater formulation only for RhsSplit mode.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
        ENDIF
     ENDIF
  ENDIF
  
  RETURN
END SUBROUTINE DNS_READ_LOCAL

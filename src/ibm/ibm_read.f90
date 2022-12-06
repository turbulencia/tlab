#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/07/21 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 before time integration starts,
!#   reads all relevant ibm parameters and with consistency check 
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_READ_INI(inifile)

  use TLAB_CONSTANTS, only : efile, lfile
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
  use IBM_VARS
  
  implicit none

  character*(*), intent(in) :: inifile

  character                 :: bakfile*32, sRes*512

  ! ================================================================== !
  ! initialization 
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'
  call TLAB_WRITE_ASCII(lfile, 'Reading IBM input data from dns.ini.')

  ! read IBM parameters
  call TLAB_WRITE_ASCII(bakfile, '#')
  call TLAB_WRITE_ASCII(bakfile, '#[IBMParameter]')
  call TLAB_WRITE_ASCII(bakfile, '#IBMScalar=<on/off>')
  call TLAB_WRITE_ASCII(bakfile, '#RestartGeometry=<yes/no>')
  call TLAB_WRITE_ASCII(bakfile, '#DataTypeGeometry=<real/int/bit>')
  call TLAB_WRITE_ASCII(bakfile, '#MaxNumberObj=<value>')
  call TLAB_WRITE_ASCII(bakfile, '#FluidPoints=<value>')

  call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'IBMScalar', 'off', sRes)
  if      (TRIM(ADJUSTL(sRes)) == 'off') then; imode_ibm_scal = 0
  else if (TRIM(ADJUSTL(sRes)) == 'on' ) then; imode_ibm_scal = 1
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBMScalar option.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if

  call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'RestartGeometry', 'no', sRes)
  if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; ibm_restart = .true.
  else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; ibm_restart = .false.
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBM Restart option.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if

  call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'DataTypeGeometry', 'int', sRes)
  if      ( TRIM(ADJUSTL(sRes)) == 'real' ) then; ibm_io = IBM_IO_REAL
  else if ( TRIM(ADJUSTL(sRes)) == 'int'  ) then; ibm_io = IBM_IO_INT
  else if ( TRIM(ADJUSTL(sRes)) == 'bit'  ) then; ibm_io = IBM_IO_BIT
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBM Data type option.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if

  call SCANINIINT(bakfile, inifile, 'IBMParameter', 'MaxNumberObj', '0', nob_max)

  call SCANINIINT(bakfile, inifile, 'IBMParameter', 'FluidPoints', '3', nflu)

  ! read geometry parameters
  call TLAB_WRITE_ASCII(bakfile, '#')
  call TLAB_WRITE_ASCII(bakfile, '#[IBMGeometry]')
  call TLAB_WRITE_ASCII(bakfile, '#Type=<none,XBars>')       
  call TLAB_WRITE_ASCII(bakfile, '#Mirrored=<yes/no>')       
  call TLAB_WRITE_ASCII(bakfile, '#Number=<value>') ! max number of elements in one spatial direction
  call TLAB_WRITE_ASCII(bakfile, '#Height=<value>')
  call TLAB_WRITE_ASCII(bakfile, '#Width=<value>')

  call SCANINICHAR(bakfile, inifile, 'IBMGeometry', 'Type', 'none', sRes)
  if      ( TRIM(ADJUSTL(sRes)) == 'none' )  then
    if ( .not. ibm_restart ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. No IBM geometry available.')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if
    continue
  else if ( TRIM(ADJUSTL(sRes)) == 'xbars' ) then; xbars_geo%name   = 'xbars'
    call SCANINICHAR(bakfile, inifile, 'IBMGeometry', 'Mirrored', 'no', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; xbars_geo%mirrored = .true.
    else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; xbars_geo%mirrored = .false.
    end if
    call SCANINIINT(bakfile, inifile, 'IBMGeometry',  'Number', '0', xbars_geo%number)
    call SCANINIINT(bakfile, inifile, 'IBMGeometry',  'Height', '0', xbars_geo%height)
    call SCANINIINT(bakfile, inifile, 'IBMGeometry',  'Width',  '0', xbars_geo%width)
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBMGeometryType option.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if

  return
end subroutine IBM_READ_INI

!########################################################################

subroutine IBM_READ_CONSISTENCY_CHECK(imode_rhs,                              &           
                                      BcsFlowJmin_type,                       &
                                      BcsScalJmin_type,    BcsScalJmax_type,  &
                                      BcsScalJmin_SfcType, BcsScalJmax_SfcType)

  use TLAB_CONSTANTS, only : efile, MAX_VARS, wi, wp
  use TLAB_VARS,      only : imax, g
  use TLAB_VARS,      only : imode_eqns, iadvection, inb_scal
  use TLAB_VARS,      only : radiation, transport, chemistry, subsidence
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
  use IBM_VARS
  
  implicit none

  integer(wi), intent(in) :: imode_rhs
  integer(wi), intent(in) :: BcsFlowJmin_type(MAX_VARS)
  integer(wi), intent(in) :: BcsScalJmin_type(MAX_VARS),    BcsScalJmax_type(MAX_VARS)
  integer(wi), intent(in) :: BcsScalJmin_SfcType(MAX_VARS), BcsScalJmax_SfcType(MAX_VARS)

  integer(wi)             :: is

  ! ================================================================== !
  ! consistency check of IBM input data
  if ( g(3)%size == 1 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Not implemented for 2D-Flow.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if
  if ( ibm_io == IBM_IO_BIT ) then
    if ( wp == 8 ) then
      if ( mod( imax, 8 ) /= 0 ) then
        call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not possible, restriction: mod(imax/8)=0.')
        call TLAB_STOP(DNS_ERROR_OPTION)
      end if
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not tested/changed yet for "working precision = single precision".')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if
  end if
  if ( nob_max <= 0 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Too no objects in flow, set number of objects correctly.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if 
  if ( nflu < 2 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Too less FluidPoints (nflu>2).')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if 
  if ( xbars_geo%name == 'xbars' ) then
    if ( ( mod(g(3)%size,2*xbars_geo%number) == 0 ) .and. ( mod(xbars_geo%width,2) /= 0 ) ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(jmax_total,(2*nbars))==0 & mod(wbar,2)==0.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    else if ( ( mod(g(3)%size,2*xbars_geo%number) /= 0 ) .and. & 
              ( mod(real(g(3)%size/(2*xbars_geo%number), wp),0.5) == 0 ) .and. &
              ( mod(xbars_geo%width,2) /= 1) ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(jmax_total/(2*nbars),0.5)==0 & mod(wbar,2)==1.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)    
    end if
    if ( xbars_geo%mirrored .and. imode_ibm_scal == 1 ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. No IBM for scalars possible with objects on upper domain.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)   
    end if
  else if ( ( xbars_geo%name == 'none' ) .and. .not. ibm_restart ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. No IBM geometry defined.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)    
  end if
  if ( .not. ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE ) ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM only implemented for incompressible mode.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if
  if ( .not. ( (iadvection == EQNS_CONVECTIVE) .or. (iadvection == EQNS_SKEWSYMMETRIC) ) ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM only implemented for convective advection scheme.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if
  if ( .not. ( imode_rhs == EQNS_RHS_COMBINED ) ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM only implemented for combined rhs mode.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if
  if ( ( radiation%type  /= EQNS_NONE) .or. &
       ( transport%type  /= EQNS_NONE) .or. &
       ( radiation%type  /= EQNS_NONE) .or. &
       ( chemistry%type  /= EQNS_NONE) .or. &
       ( subsidence%type /= EQNS_NONE)        ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM not implemented for radiation, transport, chemistry, subsidence.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if
  do is = 1,inb_scal
    if ( ( BcsScalJmin_type(is)    /= DNS_BCS_DIRICHLET ) .or. ( BcsScalJmax_type(is)    /= DNS_BCS_DIRICHLET ) .or. &
         ( BcsScalJmin_SfcType(is) /= DNS_SFC_STATIC    ) .or. ( BcsScalJmax_SfcType(is) /= DNS_SFC_STATIC    ) ) then 
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Wrong scalar BCs.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if
  end do
  do is = 1,3
    if ( BcsFlowJmin_type(is) /= DNS_BCS_DIRICHLET ) then 
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Wrong Flow BCs.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if
  end do

  return
end subroutine IBM_READ_CONSISTENCY_CHECK

!########################################################################
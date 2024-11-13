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

  use TLab_Constants, only : efile, lfile
  use TLab_WorkFlow,     only : TLab_Stop, TLab_Write_ASCII
  use IBM_VARS
  
  implicit none

  character*(*), intent(in) :: inifile

  character                 :: bakfile*32, sRes*512

  ! ================================================================== !
  ! initialization 
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'
  call TLab_Write_ASCII(lfile, 'Reading IBM input data from tlab.ini.')

! IBM Status Parameter
  call TLab_Write_ASCII(bakfile, '#')
  call TLab_Write_ASCII(bakfile, '#[IBMParameter]')
  call TLab_Write_ASCII(bakfile, '#Status=<on/off>')
  call ScanFile_Char(bakfile, inifile, 'IBMParameter', 'Status', 'off', sRes)
  if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
  else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
  else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. Wrong IBM Status option.')
      call TLab_Stop(DNS_ERROR_OPTION)
  end if
  
  ! proceed if the IBM is turned on
  if (imode_ibm == 1) then 
    ! read IBM parameters
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[IBMParameter]')
    call TLab_Write_ASCII(bakfile, '#IBMScalar=<on/off>')
    call TLab_Write_ASCII(bakfile, '#RestartGeometry=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#DataTypeGeometry=<real/int/bit>')
    call TLab_Write_ASCII(bakfile, '#MaxNumberObj=<value>')
    call TLab_Write_ASCII(bakfile, '#FluidPoints=<value>')

    call ScanFile_Char(bakfile, inifile, 'IBMParameter', 'IBMScalar', 'off', sRes)
    if      (TRIM(ADJUSTL(sRes)) == 'off') then; imode_ibm_scal = 0
    else if (TRIM(ADJUSTL(sRes)) == 'on' ) then; imode_ibm_scal = 1
    else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. Wrong IBMScalar option.')
      call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'IBMParameter', 'RestartGeometry', 'no', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; ibm_restart = .true.
    else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; ibm_restart = .false.
    else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. Wrong IBM Restart option.')
      call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'IBMParameter', 'DataTypeGeometry', 'int', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'real' ) then; ibm_io = IBM_IO_REAL
    else if ( TRIM(ADJUSTL(sRes)) == 'int'  ) then; ibm_io = IBM_IO_INT
    else if ( TRIM(ADJUSTL(sRes)) == 'bit'  ) then; ibm_io = IBM_IO_BIT
    else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. Wrong IBM Data type option.')
      call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Int(bakfile, inifile, 'IBMParameter', 'MaxNumberObj', '0', nob_max)

    call ScanFile_Int(bakfile, inifile, 'IBMParameter', 'FluidPoints', '3', nflu)

    ! read geometry parameters
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[IBMGeometry]')
    call TLab_Write_ASCII(bakfile, '#Type=<none,XBars,Hill,Valley,Box>')       
    call TLab_Write_ASCII(bakfile, '#Mirrored=<yes/no>')       
    call TLab_Write_ASCII(bakfile, '#Number=<value>') ! max number of elements in one spatial direction
    call TLab_Write_ASCII(bakfile, '#Height=<value>')
    call TLab_Write_ASCII(bakfile, '#Width=<value>')
    call TLab_Write_ASCII(bakfile, '#Alpha=<value>')
    
    call ScanFile_Char(bakfile, inifile, 'IBMGeometry', 'Type', 'none', sRes)
    if ( TRIM(ADJUSTL(sRes)) == 'none' )  then
      if ( .not. ibm_restart ) then
        call TLab_Write_ASCII(efile, 'IBM_READ_INI. No IBM geometry available.')
        call TLab_Stop(DNS_ERROR_OPTION)
      end if
      continue
    else if ( TRIM(ADJUSTL(sRes)) == 'xbars' ) then; ibm_geo%name = 'xbars'
    else if ( TRIM(ADJUSTL(sRes)) == 'hill'  ) then; ibm_geo%name = 'hill' 
    else if ( TRIM(ADJUSTL(sRes)) == 'valley') then; ibm_geo%name = 'valley'
    else if ( TRIM(ADJUSTL(sRes)) == 'box'   ) then; ibm_geo%name = 'box'
    else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. Wrong IBMGeometryType option.')
      call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'IBMGeometry', 'Mirrored', 'no', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; ibm_geo%mirrored = .true.
    else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; ibm_geo%mirrored = .false.; end if
    call ScanFile_Int(bakfile, inifile, 'IBMGeometry', 'Number', '0', ibm_geo%number)
    call ScanFile_Int(bakfile, inifile, 'IBMGeometry', 'Height', '0', ibm_geo%height)
    call ScanFile_Int(bakfile, inifile, 'IBMGeometry', 'Width',  '0', ibm_geo%width)
    call ScanFile_Int(bakfile, inifile, 'IBMGeometry', 'Alpha',  '0', ibm_geo%hill_slope)
  end if

  return
end subroutine IBM_READ_INI

!########################################################################

subroutine IBM_READ_CONSISTENCY_CHECK()

  use TLab_Constants, only : efile, MAX_VARS, wi, wp
  use TLAB_VARS,      only : imax, g
  use TLab_WorkFlow,     only : TLab_Stop, TLab_Write_ASCII
  use IBM_VARS
  
  implicit none

  ! ================================================================== !
  ! consistency check of IBM input data
  if ( g(3)%size == 1 ) then
    call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Not implemented for 2D-Flow.')
    call TLab_Stop(DNS_ERROR_OPTION)
  end if
  if ( ibm_io == IBM_IO_BIT ) then
    if ( wp == 8 ) then
      if ( mod( imax, 8 ) /= 0 ) then
        call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not possible, restriction: mod(imax/8)=0.')
        call TLab_Stop(DNS_ERROR_OPTION)
      end if
    else
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not tested/changed yet for "working precision = single precision".')
      call TLab_Stop(DNS_ERROR_OPTION)
    end if
  end if
  if ( nob_max <= 0 ) then
    call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Too no objects in flow, set number of objects correctly.')
    call TLab_Stop(DNS_ERROR_OPTION)
  end if 
  if ( nflu < 2 ) then
    call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Too less FluidPoints (nflu>=2).')
    call TLab_Stop(DNS_ERROR_OPTION)
  end if
  if ( ( ibm_geo%name == 'xbars' ) .or. &
       ( ibm_geo%name == 'hill'  ) .or. &
       ( ibm_geo%name == 'valey' ) .or. &
       ( ibm_geo%name == 'box' ) ) then
    if ( ( mod(g(3)%size,2*ibm_geo%number) == 0 ) .and. ( mod(ibm_geo%width,2) /= 0 ) ) then
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(kmax_total,(2*nbars))==0 & mod(wbar,2)==0.')
      call TLab_Stop(DNS_ERROR_UNDEVELOP)
    else if ( ( mod(g(3)%size,2*ibm_geo%number) /= 0 ) .and. & 
              ( mod(real(g(3)%size/(2*ibm_geo%number), wp),0.5) == 0 ) .and. &
              ( mod(ibm_geo%width,2) /= 1) ) then
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(kmax_total/(2*nbars),0.5)==0 & mod(wbar,2)==1.')
      call TLab_Stop(DNS_ERROR_UNDEVELOP)    
    end if
    if ( ibm_geo%mirrored .and. imode_ibm_scal == 1 ) then
      call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. No IBM for scalars possible with objects on upper domain.')
      call TLab_Stop(DNS_ERROR_UNDEVELOP)   
    end if
  else if ( ( ibm_geo%name == 'none' ) .and. .not. ibm_restart ) then
    call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. No IBM geometry defined.')
    call TLab_Stop(DNS_ERROR_UNDEVELOP)    
  end if

  return
end subroutine IBM_READ_CONSISTENCY_CHECK

!########################################################################
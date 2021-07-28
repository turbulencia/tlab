#include "types.h"

MODULE THERMO_GLOBAL
  USE TLAB_CONSTANTS, ONLY : MAX_NSP, MAX_PROF

  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_NCP =  7 ! Thermodynamic data
  TINTEGER, PARAMETER :: MAX_SAT = 10 ! Saturation pressure

! ###################################################################
! General options
! ###################################################################
  TINTEGER :: imixture

! ###################################################################
! Thermodynamics parameters
! ###################################################################
  TREAL    :: thermo_param(MAX_PROF)

! ###################################################################
! Nondimensional numbers
! ###################################################################
  TREAL :: gama0, GRATIO ! Specific heat ratio
  TREAL :: MRATIO        ! Compressible parameter


  TINTEGER                         :: NSP, NCP_CHEMKIN
  CHARACTER*16, DIMENSION(MAX_NSP) :: THERMO_SPNAME
  TREAL,        DIMENSION(MAX_NSP) :: WGHT, WGHT_INV
  TREAL                            :: WREF, TREF, RGAS
  TREAL                            :: THERMO_AI(MAX_NCP,2,MAX_NSP), THERMO_TLIM(3,MAX_NSP)

  TINTEGER                         :: iuse_chemkin    ! If using chemkin data
  CHARACTER*128                    :: chemkin_file

  TREAL                            :: WMEAN, dsmooth  ! Inifinitely fast 
  TREAL, DIMENSION(MAX_NSP)        :: YMASS

  TINTEGER                         :: NPSAT           ! Vapor saturation pressure
  TREAL                            :: THERMO_PSAT(MAX_SAT), NEWTONRAPHSON_ERROR 

END MODULE THERMO_GLOBAL

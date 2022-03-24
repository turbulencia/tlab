#include "types.h"

MODULE TLAB_CONSTANTS
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MajorVersion = 7
  TINTEGER, PARAMETER :: MinorVersion = 0

  TINTEGER, PARAMETER :: sp = KIND(1.0)
  TINTEGER, PARAMETER :: dp = KIND(1.0d0)
  ! !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
  ! integer, parameter :: sp = selected_real_kind(6, 37)
  ! !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
  ! integer, parameter :: dp = selected_real_kind(15, 307)

  TINTEGER, PARAMETER :: MAX_VARS = 20
  TINTEGER, PARAMETER :: MAX_PROF = 10
  TINTEGER, PARAMETER :: MAX_JETS =  5

  TINTEGER, PARAMETER :: MAX_NSP = 10 ! Species in the mixture

  TINTEGER, PARAMETER :: MAX_AVG_TEMPORAL  = 230
  TINTEGER, PARAMETER :: MAX_STATS_SPATIAL = 100 ! Running statistics

  CHARACTER*32, PARAMETER :: gfile = 'grid'
  CHARACTER*32, PARAMETER :: ifile = 'dns.ini'
  CHARACTER*32, PARAMETER :: ofile = 'dns.out'
  CHARACTER*32, PARAMETER :: lfile = 'dns.log'
  CHARACTER*32, PARAMETER :: efile = 'dns.err'
  CHARACTER*32, PARAMETER :: wfile = 'dns.war'
  CHARACTER*32, PARAMETER :: tfile = 'dns.trc'

  CHARACTER*32, PARAMETER :: tag_flow ='flow.'
  CHARACTER*32, PARAMETER :: tag_scal ='scal.'
  CHARACTER*32, PARAMETER :: tag_part ='part.'
  CHARACTER*32, PARAMETER :: tag_traj ='traj.'

END MODULE TLAB_CONSTANTS

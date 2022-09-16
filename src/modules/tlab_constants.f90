#include "types.h"

module TLAB_CONSTANTS
    use TLAB_TYPES, only: wp
    implicit none
    save

    TINTEGER, parameter :: MajorVersion = 7
    TINTEGER, parameter :: MinorVersion = 0

    TINTEGER, parameter :: MAX_VARS = 20
    TINTEGER, parameter :: MAX_PROF = 10
    TINTEGER, parameter :: MAX_JETS = 5

    TINTEGER, parameter :: MAX_NSP = 10 ! Species in the mixture

    TINTEGER, parameter :: MAX_AVG_TEMPORAL = 230
    TINTEGER, parameter :: MAX_STATS_SPATIAL = 100 ! Running statistics

    character*32, parameter :: gfile = 'grid'
    character*32, parameter :: ifile = 'dns.ini'
    character*32, parameter :: ofile = 'dns.out'
    character*32, parameter :: lfile = 'dns.log'
    character*32, parameter :: efile = 'dns.err'
    character*32, parameter :: wfile = 'dns.war'
    character*32, parameter :: tfile = 'dns.trc'

    character*32, parameter :: tag_flow = 'flow.'
    character*32, parameter :: tag_scal = 'scal.'
    character*32, parameter :: tag_part = 'part.'
    character*32, parameter :: tag_traj = 'traj.'

    real(wp), parameter :: pi_wp = 3.14159265358979323846_wp 

end module TLAB_CONSTANTS

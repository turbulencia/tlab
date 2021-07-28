#include "types.h"

MODULE DNS_LOCAL
  USE DNS_TYPES,  ONLY : filter_dt, grid_dt
  USE TLAB_VARS, ONLY : MAX_NSP
#ifdef USE_PSFFT
  USE NB3DFFT,    ONLY : NB3DFFT_SCHEDLTYPE
#endif
  IMPLICIT NONE
  SAVE

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: nitera_first, nitera_last, nitera_save, nitera_stats, nitera_log, nitera_pln
  TINTEGER :: nitera_stats_spa ! Accumulate statistics in spatial mode

  TINTEGER :: idivergence, imode_rhs

! ###################################################################
! Control
! ###################################################################
  TINTEGER :: ilimit_flow, ilimit_scal
  TREAL    :: p_bound_min, p_bound_max, r_bound_min, r_bound_max ! pressure and density
  TREAL    :: s_bound_min(MAX_NSP), s_bound_max(MAX_NSP)         ! scalars
  TREAL    :: d_bound_max                                        ! dilatation

! ###################################################################
! Variable viscosity
! ###################################################################
  LOGICAL :: flag_viscosity
  TREAL   :: visc_stop, visc_time, visc_rate

! ###########################################################
! Filters
! ###########################################################
  TINTEGER :: FilterDomainStep

! ###################################################################
! Output data
! ###################################################################
  TINTEGER, DIMENSION(3)              :: tower_stride           ! Towers
  TINTEGER                            :: tower_mode

! ###################################################################
  TREAL    :: logs_data(20)

! ###################################################################
! NB3DFFT library
! ###################################################################
#ifdef USE_PSFFT
  TYPE(NB3DFFT_SCHEDLTYPE), SAVE :: nbcsetup
#endif

END MODULE DNS_LOCAL

MODULE DNS_ARRAYS
  IMPLICIT NONE
  SAVE
  PRIVATE

  TREAL, ALLOCATABLE, PUBLIC :: hq(:,:)      ! Right-hand sides Eulerian fields
  TREAL, ALLOCATABLE, PUBLIC :: hs(:,:)      ! Right-hand sides Eulerian fields
  TREAL, ALLOCATABLE, PUBLIC :: l_hq(:,:)    ! Right-hand sides Lagrangian fields
  TREAL, ALLOCATABLE, PUBLIC :: l_comm(:)    ! Communication space for Lagrangian fields

END MODULE DNS_ARRAYS

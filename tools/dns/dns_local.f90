#include "types.h"

MODULE DNS_LOCAL
  USE DNS_TYPES,  ONLY : filter_dt, grid_dt
  USE DNS_GLOBAL, ONLY : MAX_NSP
#ifdef USE_PSFFT  
  USE NB3DFFT,    ONLY : NB3DFFT_SCHEDLTYPE
#endif
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_SAVEPLANES = 10

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: rkm_mode, rkm_substep, rkm_endstep
  TREAL    :: cfl, dtime
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
! Statistics
! ###################################################################
  TINTEGER :: fstavg, fstpdf, fstinter, ffltdmp

! ###################################################################
! Variable viscosity
! ###################################################################
  TREAL    :: visctime, viscstart, viscstop
  TINTEGER :: iviscchg

! ###################################################################
! Boundary conditions
! ###################################################################
! Compressible
  TINTEGER :: bcs_inf(2,2,3), bcs_out(2,2,3) ! 1. index: lower and upper values
                                             ! 2. index: derivative order
                                             ! 3. index: direction
  TINTEGER :: bcs_euler_drift
  TREAL    :: bcs_sigma_out
  TREAL    :: bcs_sigma_inf_imin, bcs_sigma_inf_imax, bcs_sigma_inf_j
  TREAL    :: bcs_sigma_trans

! ###########################################################
! Filters
! ###########################################################
  TINTEGER :: FilterDomainStep

! ###################################################################
! Output data
! ###################################################################
  TINTEGER                            :: nplanes_i, nplanes_j, nplanes_k, nplanes_j_aux ! Planes
  TINTEGER, DIMENSION(MAX_SAVEPLANES) :: planes_i,  planes_j,  planes_k
  
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

! ###################################################################
! vaux array
! ###################################################################
#ifdef LES
  TINTEGER, PARAMETER :: vindex_size=22
#else
  TINTEGER, PARAMETER :: vindex_size=1
#endif

  TINTEGER, DIMENSION(vindex_size) :: vindex, vsize

  INTEGER VA_MEAN_WRK
  PARAMETER (VA_MEAN_WRK=1)

! -------------------------------------------------------------------
#ifdef LES
  INTEGER VA_LES_FLT0X, VA_LES_FLT0Y, VA_LES_FLT0Z
  PARAMETER (VA_LES_FLT0X=2)
  PARAMETER (VA_LES_FLT0Y=3)
  PARAMETER (VA_LES_FLT0Z=4)
  INTEGER VA_LES_FLT1X, VA_LES_FLT1Y, VA_LES_FLT1Z
  PARAMETER (VA_LES_FLT1X=5)
  PARAMETER (VA_LES_FLT1Y=6)
  PARAMETER (VA_LES_FLT1Z=7)
  INTEGER VA_LES_FLT2X, VA_LES_FLT2Y, VA_LES_FLT2Z
  PARAMETER (VA_LES_FLT2X=8)
  PARAMETER (VA_LES_FLT2Y=9)
  PARAMETER (VA_LES_FLT2Z=10)
  INTEGER VA_LES_SGSCOEFU, VA_LES_SGSCOEFE, VA_LES_SGSCOEFZ
  PARAMETER (VA_LES_SGSCOEFU=11)
  PARAMETER (VA_LES_SGSCOEFE=12)
  PARAMETER (VA_LES_SGSCOEFZ=13)
  INTEGER VA_LES_ARM_UF, VA_LES_ARM_PF, VA_LES_ARM_ZF
  PARAMETER (VA_LES_ARM_UF=14)
  PARAMETER (VA_LES_ARM_PF=15)
  PARAMETER (VA_LES_ARM_ZF=16)
  INTEGER VA_LES_ARM_UH, VA_LES_ARM_PH, VA_LES_ARM_ZH
  PARAMETER (VA_LES_ARM_UH=17)
  PARAMETER (VA_LES_ARM_PH=18)
  PARAMETER (VA_LES_ARM_ZH=19)
  INTEGER VA_ARM_WRK, VA_ARM_C0, VA_LES_FDF_BS
  PARAMETER (VA_ARM_WRK   =20)
  PARAMETER (VA_ARM_C0    =21)
  PARAMETER (VA_LES_FDF_BS=22) 
#endif
  
END MODULE DNS_LOCAL

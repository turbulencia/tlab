#include "types.h"

MODULE DNS_LOCAL
  USE DNS_GLOBAL, ONLY : MAX_NSP, MAX_VARS 
#ifdef USE_PSFFT  
  USE NB3DFFT,    ONLY : NB3DFFT_SCHEDLTYPE
#endif
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_SAVEPLANES = 10
  TINTEGER, PARAMETER :: MAX_FRC_FREC   = 32

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: rkm_mode, rkm_substep, rkm_endstep
  TREAL    :: cfl, dtime
  TINTEGER :: nitera_first, nitera_last, nitera_save, nitera_stats, nitera_log, nitera_pln
  TINTEGER :: frunstat, frunline, frunplane

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
  TINTEGER :: bcs_euler_imin, bcs_euler_imax, bcs_visc_imin, bcs_visc_imax
  TINTEGER :: bcs_euler_jmin, bcs_euler_jmax, bcs_visc_jmin, bcs_visc_jmax
  TINTEGER :: bcs_euler_kmin, bcs_euler_kmax, bcs_visc_kmin, bcs_visc_kmax

  TINTEGER :: bcs_euler_drift
  TREAL    :: bcs_sigma_out
  TREAL    :: bcs_sigma_inf_imin, bcs_sigma_inf_imax, bcs_sigma_inf_j
  TREAL    :: bcs_sigma_trans
  TREAL    :: bcs_p_imin, bcs_p_imax, bcs_p_jmin, bcs_p_jmax, bcs_p_kmin, bcs_p_kmax

! Incompressible
  TINTEGER :: bcs_scal_imin(MAX_NSP), bcs_scal_imax(MAX_NSP), bcs_flow_imin, bcs_flow_imax
  TINTEGER :: bcs_scal_jmin(MAX_NSP), bcs_scal_jmax(MAX_NSP), bcs_flow_jmin, bcs_flow_jmax
  TINTEGER :: bcs_scal_kmin(MAX_NSP), bcs_scal_kmax(MAX_NSP), bcs_flow_kmin, bcs_flow_kmax

! ################################################################
! Buffer zone
! ################################################################
  TINTEGER :: buff_type

  TINTEGER :: buff_nps_imin, buff_nps_imax
  TINTEGER :: buff_nps_jmin, buff_nps_jmax
  TINTEGER :: buff_nps_u_jmin, buff_nps_u_jmax, buff_nps_e_jmin, buff_nps_e_jmax, buff_nps_s_jmin, buff_nps_s_jmax
  TINTEGER :: buff_imin, buff_imax
  TINTEGER :: buff_u_jmin,buff_u_jmax, buff_e_jmin,buff_e_jmax, buff_s_jmin,buff_s_jmax

! Relaxation type
  TREAL    :: buff_param_u(2), buff_param_s(2)
  TINTEGER :: buff_v_free, buff_load, buff_hard_on(MAX_VARS)
  TREAL    :: buff_hard(MAX_VARS,2)

! ###########################################################
! Filters
! ###########################################################
  TINTEGER :: ifilt_domain, ifilt_step, ifilt_scalar, ifilt_x, ifilt_y, ifilt_z
  TREAL    :: ifilt_alpha
  TINTEGER :: ifilt_inflow, ifilt_inflow_iwidth, ifilt_inflow_jwidth, ifilt_inflow_step

! ###################################################################
! Inflow field in spatial mode
! ###################################################################
  TINTEGER :: imax_inf,jmax_inf,kmax_inf
  TREAL    :: scalex_inf,scaley_inf,scalez_inf

  TINTEGER :: ifrc_mode, ifrc_ifield
  TREAL    :: frc_length, frc_adapt

! ###################################################################
! Discrete forcing
! ###################################################################
  TINTEGER :: ifrcdsc_mode
  TREAL    :: frc_delta
  
  TINTEGER :: nx2d, nx3d, nz3d
  TREAL    :: A2D(MAX_FRC_FREC), Phix2d(MAX_FRC_FREC)
  TREAL    :: A3D(MAX_FRC_FREC), Phix3d(MAX_FRC_FREC), Phiz3d(MAX_FRC_FREC)

! ###################################################################
! Output data
! ###################################################################
  TINTEGER                            :: nplanes_i, nplanes_j, nplanes_k ! Planes
  TINTEGER, DIMENSION(MAX_SAVEPLANES) :: planes_i,  planes_j,  planes_k
  
  TINTEGER, DIMENSION(3)              :: tower_stride           ! Towers
  TINTEGER                            :: tower_mode  

! ###################################################################
  TREAL    :: logs_data(20)

! ###################################################################
! vaux array
! ###################################################################
#ifdef LES
  TINTEGER, PARAMETER :: vindex_size=36
#else
  TINTEGER, PARAMETER :: vindex_size=15
#endif

  TINTEGER, DIMENSION(vindex_size) :: vindex, vsize

! -------------------------------------------------------------------
  INTEGER VA_FLT_CX, VA_FLT_CY, VA_FLT_CZ
  PARAMETER (VA_FLT_CX=1)
  PARAMETER (VA_FLT_CY=2)
  PARAMETER (VA_FLT_CZ=3)
  INTEGER VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI
  PARAMETER (VA_BUFF_HT=4)
  PARAMETER (VA_BUFF_HB=5)
  PARAMETER (VA_BUFF_VI=6)
  PARAMETER (VA_BUFF_VO=7)
  INTEGER VA_BCS_HT, VA_BCS_HB, VA_BCS_VO, VA_BCS_VI
  PARAMETER (VA_BCS_HT=8)
  PARAMETER (VA_BCS_HB=9)
  PARAMETER (VA_BCS_VI=10)
  PARAMETER (VA_BCS_VO=11)

  INTEGER VA_MEAN_WRK
  PARAMETER (VA_MEAN_WRK=12)
  INTEGER VA_LINE_SPA_WRK
  PARAMETER (VA_LINE_SPA_WRK=13)
  INTEGER VA_TIMES
  PARAMETER (VA_TIMES=14)
  INTEGER VA_PLANE_SPA_WRK
  PARAMETER (VA_PLANE_SPA_WRK=15)

! ###################################################################
! NB3DFFT library
! ###################################################################
#ifdef USE_PSFFT  
  TYPE(NB3DFFT_SCHEDLTYPE), SAVE :: nbcsetup
#endif 


! ###################################################################
! -------------------------------------------------------------------
#ifdef LES
  INTEGER VA_LES_FLT0X, VA_LES_FLT0Y, VA_LES_FLT0Z
  PARAMETER (VA_LES_FLT0X=16)
  PARAMETER (VA_LES_FLT0Y=17)
  PARAMETER (VA_LES_FLT0Z=18)
  INTEGER VA_LES_FLT1X, VA_LES_FLT1Y, VA_LES_FLT1Z
  PARAMETER (VA_LES_FLT1X=19)
  PARAMETER (VA_LES_FLT1Y=20)
  PARAMETER (VA_LES_FLT1Z=21)
  INTEGER VA_LES_FLT2X, VA_LES_FLT2Y, VA_LES_FLT2Z
  PARAMETER (VA_LES_FLT2X=22)
  PARAMETER (VA_LES_FLT2Y=23)
  PARAMETER (VA_LES_FLT2Z=24)
  INTEGER VA_LES_SGSCOEFU, VA_LES_SGSCOEFE, VA_LES_SGSCOEFZ
  PARAMETER (VA_LES_SGSCOEFU=25)
  PARAMETER (VA_LES_SGSCOEFE=26)
  PARAMETER (VA_LES_SGSCOEFZ=27)
  INTEGER VA_LES_ARM_UF, VA_LES_ARM_PF, VA_LES_ARM_ZF
  PARAMETER (VA_LES_ARM_UF=28)
  PARAMETER (VA_LES_ARM_PF=29)
  PARAMETER (VA_LES_ARM_ZF=30)
  INTEGER VA_LES_ARM_UH, VA_LES_ARM_PH, VA_LES_ARM_ZH
  PARAMETER (VA_LES_ARM_UH=31)
  PARAMETER (VA_LES_ARM_PH=32)
  PARAMETER (VA_LES_ARM_ZH=33)
  INTEGER VA_ARM_WRK, VA_ARM_C0, VA_LES_FDF_BS
  PARAMETER (VA_ARM_WRK   =34)
  PARAMETER (VA_ARM_C0    =35)
  PARAMETER (VA_LES_FDF_BS=36) 
#endif
  
END MODULE DNS_LOCAL

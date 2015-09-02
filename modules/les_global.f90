MODULE LES_GLOBAL
  IMPLICIT NONE
  SAVE

#include "types.h"

! ###################################################################
! File names
! ###################################################################
  CHARACTER*32 :: lesfile

! #################################################################
! # LES Variables
! #################################################################
TREAL sgs_csm, smg_prandtl, smg_schmidt, &
     smg_udiff, smg_ediff, smg_zdiff, smg_ucoef, smg_ecoef, smg_zcoef, &
     sgs_alpha, sgs_smagtrans, sgs_smagdelta, sgs_pdil

TINTEGER iles, iles_type, iles_type_chem, iles_inviscid,&
     sgs_devsmag, isgs_f0size, isgs_f1size, iles_mpitype,&
     iles_type_diss, iles_type_recchem,&
     iles_type_tran, iles_type_regu, iles_type_disZchem,&
     iles_type_dischem

TREAL arm_remin, arm_dre, arm_c0_inviscid,&
     arm_tact, arm_tflame, arm_zst, arm_smooth,&
     arm_ucoef, arm_pcoef, arm_zcoef

TINTEGER iarm_spc, iarm_nl, iarm_nre, iarmavg_pts

TINTEGER iavgles, ilesstat_maj_ver, ilesstat_min_ver

TREAL LES_CHW_Qp

TINTEGER les_fdf_bs_maxmean, les_fdf_bs_maxvar

CHARACTER*128 les_fdf_bs_file

END MODULE LES_GLOBAL

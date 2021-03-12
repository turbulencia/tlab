#ifndef LES_CONST_H_INCLUDED
#define LES_CONST_H_INCLUDED

! Transport terms
#define LES_TRAN_NONE     0
#define LES_TRAN_SSM      1
#define LES_TRAN_ARMRANS  2
#define LES_TRAN_ARMSPTL  3
#define LES_TRAN_ARMINVS  4

! Regularization
#define LES_REGU_NONE        0
#define LES_REGU_SMGSTA      1
#define LES_REGU_SMGDYN      2
#define LES_REGU_SMGSTARMS   3
#define LES_REGU_SMGDYNRMS   4

! SGS dissipation function models
#define LES_DISS_NONE       0
#define LES_DISS_DISSSG     1

! SGS chemistry models
#define LES_CHEM_NONE       0
#define LES_CHEM_STAT_BS    1
#define LES_CHEM_STAT_PS    2
#define LES_CHEM_QUASIBS    3

! Reconstruction level in chemical sources
#define LES_RECCHEM_SSM      1
#define LES_RECCHEM_ARMRANS  2
#define LES_RECCHEM_ARMSPTL  3
#define LES_RECCHEM_ARMINVS  4

! Conditional dissipation model in chemical sources
#define LES_DISZCHEM_MEAN     1
#define LES_DISZCHEM_ONEDIM   2

! Dissipation model in chemical sources
#define LES_DISCHEM_GRAD     1
#define LES_DISCHEM_DISSSG   2

#endif

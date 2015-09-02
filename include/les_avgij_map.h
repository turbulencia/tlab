! LES of momentum equation: 
! Scale-similarity part & isotropic (Smagorinsky) part

#define LA_SMG_CS_U(j)    mean1d_les(j,1,1)
#define LA_SMG_RNU(j)     mean1d_les(j,1,2)

#define LA_SMG_TAUxx(j)   mean1d_les(j,1,3)
#define LA_SMG_TAUyy(j)   mean1d_les(j,1,4)
#define LA_SMG_TAUzz(j)   mean1d_les(j,1,5)
#define LA_SMG_TAUxy(j)   mean1d_les(j,1,6)
#define LA_SMG_TAUxz(j)   mean1d_les(j,1,7) 
#define LA_SMG_TAUyz(j)   mean1d_les(j,1,8)
#define LA_SMG_EPS_U(j)   mean1d_les(j,1,9)

#define LA_SSM_TAUxx(j)   mean1d_les(j,1,10)
#define LA_SSM_TAUyy(j)   mean1d_les(j,1,11)
#define LA_SSM_TAUzz(j)   mean1d_les(j,1,12)
#define LA_SSM_TAUxy(j)   mean1d_les(j,1,13)
#define LA_SSM_TAUxz(j)   mean1d_les(j,1,14) 
#define LA_SSM_TAUyz(j)   mean1d_les(j,1,15)
#define LA_SSM_EPS_U(j)   mean1d_les(j,1,16)

#define LA_SMG_CS_P(j)    mean1d_les(j,1,17)
#define LA_SMG_RALPHA(j)  mean1d_les(j,1,18)

#define LA_SMG_QPx(j)     mean1d_les(j,1,19)
#define LA_SMG_QPy(j)     mean1d_les(j,1,20)
#define LA_SMG_QPz(j)     mean1d_les(j,1,21)

#define LA_SSM_QPx(j)     mean1d_les(j,1,22)
#define LA_SSM_QPy(j)     mean1d_les(j,1,23)
#define LA_SSM_QPz(j)     mean1d_les(j,1,24)

#define LA_SMG_FORCE_U(j) mean1d_les(j,1,25)
#define LA_SMG_FORCE_V(j) mean1d_les(j,1,26)
#define LA_SSM_FORCE_U(j) mean1d_les(j,1,27)
#define LA_SSM_FORCE_V(j) mean1d_les(j,1,28)

#define LA_SMG_TRA_U(j)   mean1d_les(j,1,29)
#define LA_SSM_TRA_U(j)   mean1d_les(j,1,30)

#define LA_SMG_CON_Ux(j)  mean1d_les(j,1,31)
#define LA_SMG_CON_Uy(j)  mean1d_les(j,1,32)
#define LA_SSM_CON_Ux(j)  mean1d_les(j,1,33)
#define LA_SSM_CON_Uy(j)  mean1d_les(j,1,34)

#define LA_SMG_EPS_P(j)   mean1d_les(j,1,35)
#define LA_SSM_EPS_P(j)   mean1d_les(j,1,36)
#define LA_SMG_FORCE_P(j) mean1d_les(j,1,37)
#define LA_SSM_FORCE_P(j) mean1d_les(j,1,38)
#define LA_SMG_TRA_P(j)   mean1d_les(j,1,39)
#define LA_SSM_TRA_P(j)   mean1d_les(j,1,40)

#define LA_MOMENTUM_SIZE 40

! LES of n scalar equations: 
! Scale-similarity part & isotropic (Smagorinsky) part

#define LS_SMG_CS_Z(j)    mean1d_les_sc(j,1,1)
#define LS_SMG_RD(j)      mean1d_les_sc(j,1,2)

#define LS_SMG_Qx(j)      mean1d_les_sc(j,1,3)
#define LS_SMG_Qy(j)      mean1d_les_sc(j,1,4) 
#define LS_SMG_Qz(j)      mean1d_les_sc(j,1,5)
#define LS_SMG_EPS_Z(j)   mean1d_les_sc(j,1,6)

#define LS_SSM_Qx(j)      mean1d_les_sc(j,1,7)
#define LS_SSM_Qy(j)      mean1d_les_sc(j,1,8) 
#define LS_SSM_Qz(j)      mean1d_les_sc(j,1,9)
#define LS_SSM_EPS_Z(j)   mean1d_les_sc(j,1,10)

#define LS_SMG_FORCE_Z(j) mean1d_les_sc(j,1,11)
#define LS_SSM_FORCE_Z(j) mean1d_les_sc(j,1,12)

#define LS_SMG_TRA_Z(j)   mean1d_les_sc(j,1,13)
#define LS_SSM_TRA_Z(j)   mean1d_les_sc(j,1,14)

#define LS_SCALAR_SIZE 14

! LES of reacting terms of n scalar equations: 
! ARM model & PDF model

#define LC_ARM_c0(j)       mean1d_les_ch(j,1,1)
#define LC_ARM_Z2_R(j)     mean1d_les_ch(j,1,2)
#define LC_ARM_Z2_SG(j)    mean1d_les_ch(j,1,3)
#define LC_ARM_RZ2_SG(j)   mean1d_les_ch(j,1,4)
#define LC_ARM_CX_Z2_SG(j) mean1d_les_ch(j,1,5)
#define LC_ARM_CY_Z2_SG(j) mean1d_les_ch(j,1,6)
#define LC_ARM_Z4_R(j)     mean1d_les_ch(j,1,7)
#define LC_ARM_Z4_SG(j)    mean1d_les_ch(j,1,8)
#define LC_ARM_Z8_R(j)     mean1d_les_ch(j,1,9)
#define LC_ARM_Z8_SG(j)    mean1d_les_ch(j,1,10)
#define LC_ARM_TZ_R(j)     mean1d_les_ch(j,1,11)
#define LC_ARM_TZ_SG(j)    mean1d_les_ch(j,1,12)
#define LC_ARM_RZ_R(j)     mean1d_les_ch(j,1,13)
#define LC_ARM_RZ_SG(j)    mean1d_les_ch(j,1,14)
#define LC_ARM_WZ_R(j)     mean1d_les_ch(j,1,15)
#define LC_ARM_WZ_SG(j)    mean1d_les_ch(j,1,16)

#define LC_CHEM_SIZE 16

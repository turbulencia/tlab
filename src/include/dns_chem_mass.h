! Mass Fractions
  CHW_ETA = (MACRO_ZINPUT-XIST)/dsmooth
  YMASS(1) = CST(1)*                    CHW_DHALF*LOG((EXP(C_2_R*CHW_ETA)+C_1_R)/CHW_LN0)
  YMASS(2) = CST(2)*(MACRO_ZINPUT-C_1_R-CHW_DHALF*LOG((EXP(C_2_R*CHW_ETA)+C_1_R)/CHW_LN1))

  YMASS(3) = CST(3)*(MACRO_ZINPUT-YMASS(1)/YFUEL)
  YMASS(4) = CST(4)*(MACRO_ZINPUT-YMASS(1)/YFUEL)

  CHW_HETA  = C_05_R*(C_1_R+TANH(CHW_ETA))
  YMASSP(1) = CST(1)* CHW_HETA
  YMASSP(2) = CST(2)*(C_1_R-CHW_HETA)

  YMASSP(3) = CST(3)*(C_1_R-YMASSP(1)/YFUEL)
  YMASSP(4) = CST(4)*(C_1_R-YMASSP(1)/YFUEL)

  YMASS(NSP)  = C_1_R
  YMASSP(NSP) = C_0_R

! ################################################################
#ifdef VECTOR
  YMASS(NSP) = YMASS(NSP)-YMASS(1)
  YMASSP(NSP) = YMASSP(NSP)-YMASSP(1)
  YMASS(NSP) = YMASS(NSP)-YMASS(2)
  YMASSP(NSP) = YMASSP(NSP)-YMASSP(2)
  YMASS(NSP) = YMASS(NSP)-YMASS(3)
  YMASSP(NSP) = YMASSP(NSP)-YMASSP(3)
  YMASS(NSP) = YMASS(NSP)-YMASS(4)
  YMASSP(NSP) = YMASSP(NSP)-YMASSP(4)
#else 
  DO is=1,NSP-1
     YMASS(NSP)  = YMASS(NSP)-YMASS(is)
     YMASSP(NSP) = YMASSP(NSP)-YMASSP(is)
  ENDDO
#endif
! ################################################################
#ifdef VECTOR
  YMASS(1) = MAX(YMASS(1),C_0_R)
  YMASS(2) = MAX(YMASS(2),C_0_R)
  YMASS(3) = MAX(YMASS(3),C_0_R)
  YMASS(4) = MAX(YMASS(4),C_0_R)
  YMASS(NSP) = MAX(YMASS(NSP),C_0_R)
#else 
  DO is=1, NSP
     YMASS(is) = MAX(YMASS(is),C_0_R)
  ENDDO
#endif

! ################################################################
! Mean Molecular Weight
  WMEAN = C_0_R

#ifdef VECTOR
  WMEAN = WMEAN + YMASS(1)/WGHT(1)
  WMEAN = WMEAN + YMASS(2)/WGHT(2)
  WMEAN = WMEAN + YMASS(3)/WGHT(3)
  WMEAN = WMEAN + YMASS(4)/WGHT(4)
  WMEAN = WMEAN + YMASS(NSP)/WGHT(NSP)
#else 
  DO is=1, NSP
     WMEAN = WMEAN + YMASS(is)/WGHT(is)
  ENDDO
#endif

  WMEAN = C_1_R/(WMEAN*WREF)

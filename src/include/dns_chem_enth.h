! Enthalpy Terms
  CHW_CPTOTAL = C_0_R
  CH_G = C_0_R
  CH_F = C_0_R
  CH_H = C_0_R
  CHW_T = MACRO_TINPUT
  ! ################################################################
  DO is=1, NSP
     IF ( CHW_T .LT. THERMO_TLIM(3,is) ) THEN
        im = 2
     ELSE
        im = 1
     ENDIF
     CHW_HI = C_0_R
     CHW_CP = C_0_R
     DO icp=NCP, 1, -1
        CHW_HI = CHW_HI*CHW_T + THERMO_AI(icp,im,is)/M_REAL(icp)
        CHW_CP = CHW_CP*CHW_T + THERMO_AI(icp,im,is)
     ENDDO
     CHW_HI = CHW_HI*CHW_T + THERMO_AI(6,im,is)

     CHW_CPTOTAL = CHW_CPTOTAL + CHW_CP*YMASS(is)
     CH_F = CH_F + CHW_HI*YMASSP(is)
     CH_H = CH_H + CHW_HI*YMASS(is)
     CH_G = CH_G + YMASSP(is)/WGHT(is)
  ENDDO
! ################################################################
  CH_G = CH_G * WREF
  CPW = CHW_CPTOTAL*WMEAN

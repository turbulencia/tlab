! top-hat filter at level 0 (F)
IF ( iunifx .EQ. 0 ) THEN
   IF ( i1bc .EQ. 0 ) THEN
      vsize(VA_LES_FLT0X) = 1
      vsize(VA_LES_FLT1X) = 1
   ELSE
      vsize(VA_LES_FLT0X) = isgs_f0size
      vsize(VA_LES_FLT1X) = isgs_f1size
   ENDIF
ELSE
   vsize(VA_LES_FLT0X) = imax*(isgs_f0size + 1)
   vsize(VA_LES_FLT1X) = imax*(isgs_f1size + 1)
ENDIF
IF ( iunify .EQ. 0 ) THEN
   IF ( j1bc .EQ. 0 ) THEN
      vsize(VA_LES_FLT0Y) = 1
      vsize(VA_LES_FLT1Y) = 1
   ELSE
      vsize(VA_LES_FLT0Y) = isgs_f0size
      vsize(VA_LES_FLT1Y) = isgs_f1size
   ENDIF
ELSE
   vsize(VA_LES_FLT0Y) = jmax*(isgs_f0size + 1)
   vsize(VA_LES_FLT1Y) = jmax*(isgs_f1size + 1)
ENDIF
IF ( iunifz .EQ. 0 ) THEN
   IF ( k1bc .EQ. 0 ) THEN
      vsize(VA_LES_FLT0Z) = 1
      vsize(VA_LES_FLT1Z) = 1
   ELSE
      vsize(VA_LES_FLT0Z) = isgs_f0size
      vsize(VA_LES_FLT1Z) = isgs_f1size
   ENDIF
ELSE
   vsize(VA_LES_FLT0Z) = kmax_total*(isgs_f0size + 1)
   vsize(VA_LES_FLT1Z) = kmax_total*(isgs_f1size + 1)
ENDIF

!     top-hat filter at level 1 for dynamic mixed case (G, and H=F*G)
IF ( iles_type_regu .EQ. LES_REGU_SMGDYN      .OR. &
     iles_type_regu .EQ. LES_REGU_SMGDYNRMS ) THEN
   IF ( iunifx .EQ. 0 ) THEN
      IF ( i1bc .EQ. 0 ) THEN
         vsize(VA_LES_FLT2X) = isgs_f0size + isgs_f1size + 1
      ELSE
         vsize(VA_LES_FLT2X) = (isgs_f0size+isgs_f1size+1)*&
              (isgs_f0size/2+isgs_f1size/2+1)
      ENDIF
   ELSE
      vsize(VA_LES_FLT2X) = imax*(isgs_f0size+isgs_f1size+1)
   ENDIF
   IF ( iunify .EQ. 0 ) THEN
      IF ( j1bc .EQ. 0 ) THEN
         vsize(VA_LES_FLT2Y) = isgs_f0size+isgs_f1size+1
      ELSE
         vsize(VA_LES_FLT2Y) = (isgs_f0size+isgs_f1size+1)*&
                       (isgs_f0size/2+isgs_f1size/2+1)
      ENDIF
   ELSE
      vsize(VA_LES_FLT2Y) = jmax*(isgs_f0size+isgs_f1size+1)
   ENDIF
   IF ( iunifz .EQ. 0 ) THEN
      IF ( k1bc .EQ. 0 ) THEN
         vsize(VA_LES_FLT2Z) = isgs_f0size+isgs_f1size+1
      ELSE
         vsize(VA_LES_FLT2Z) = (isgs_f0size+isgs_f1size+1)*&
              (isgs_f0size/2+isgs_f1size/2+1)
      ENDIF
   ELSE
      vsize(VA_LES_FLT2Z) = kmax_total*&
           (isgs_f0size+isgs_f1size+1)
   ENDIF
ELSE
   vsize(VA_LES_FLT2X) = 1
   vsize(VA_LES_FLT2Y) = 1
   vsize(VA_LES_FLT2Z) = 1
ENDIF

! SMG coefficients
IF ( isim_mode .EQ. DNS_MODE_TEMPORAL ) THEN
   vsize(VA_LES_SGSCOEFU) = jmax
   vsize(VA_LES_SGSCOEFE) = jmax
   vsize(VA_LES_SGSCOEFZ) = jmax*inb_scalars
ELSE IF ( isim_mode .EQ. DNS_MODE_SPATIAL ) THEN
   vsize(VA_LES_SGSCOEFU) = imax*jmax
   vsize(VA_LES_SGSCOEFE) = imax*jmax
   vsize(VA_LES_SGSCOEFZ) = imax*jmax*inb_scalars
ENDIF

! ARM model for dynamic case
IF ( isim_mode .EQ. DNS_MODE_TEMPORAL ) THEN
   n_vaux_arm_size = jmax
ELSE IF ( isim_mode .EQ. DNS_MODE_SPATIAL ) THEN
   n_vaux_arm_size = imax*jmax
ENDIF

IF ( iles_type_tran .EQ. LES_TRAN_ARMRANS    .OR. &
     iles_type_tran .EQ. LES_TRAN_ARMSPTL    .OR. &
     iles_type_tran .EQ. LES_TRAN_ARMINVS  ) THEN
   vsize(VA_LES_ARM_UF) = n_vaux_arm_size
   vsize(VA_LES_ARM_PF) = n_vaux_arm_size
   vsize(VA_LES_ARM_ZF) = n_vaux_arm_size
   IF ( iles_type_regu .EQ. LES_REGU_SMGDYN      .OR. &
        iles_type_regu .EQ. LES_REGU_SMGDYNRMS ) THEN
      vsize(VA_LES_ARM_UH) = n_vaux_arm_size
      vsize(VA_LES_ARM_PH) = n_vaux_arm_size
      vsize(VA_LES_ARM_ZH) = n_vaux_arm_size
   ELSE
      vsize(VA_LES_ARM_UH) = 1
      vsize(VA_LES_ARM_PH) = 1
      vsize(VA_LES_ARM_ZH) = 1
   ENDIF
ELSE
   vsize(VA_LES_ARM_UF) = 1
   vsize(VA_LES_ARM_PF) = 1
   vsize(VA_LES_ARM_ZF) = 1
ENDIF

vsize(VA_ARM_WRK) = iarm_nl*iarm_nre + iarm_nre + iarm_nl*iarm_nre
vsize(VA_ARM_C0)  = n_vaux_arm_size


  IF ( bcs_visc_imin .EQ. 2 ) THEN
     i1vsout = 1
  ELSE
     i1vsout = 0
  ENDIF

  IF ( bcs_visc_imax .EQ. 2 ) THEN
     imxvsout = 1
  ELSE
     imxvsout = 0
  ENDIF

  IF ( bcs_visc_jmin .EQ. 2 ) THEN
     j1vsout = 1
  ELSE
     j1vsout = 0
  ENDIF

  IF ( bcs_visc_jmax .EQ. 2 ) THEN
     jmxvsout = 1
  ELSE
     jmxvsout = 0
  ENDIF

  IF ( bcs_visc_kmin .EQ. 2 ) THEN
     k1vsout = 1
  ELSE
     k1vsout = 0
  ENDIF

  IF ( bcs_visc_kmax .EQ. 2 ) THEN
     kmxvsout = 1
  ELSE
     kmxvsout = 0
  ENDIF


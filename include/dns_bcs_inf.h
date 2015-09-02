  IF ( bcs_visc_imin .EQ. 1 ) THEN
     i1vsin = 1
  ELSE
     i1vsin = 0
  ENDIF

  IF ( bcs_visc_imax .EQ. 1 ) THEN
     imxvsin = 1
  ELSE
     imxvsin = 0
  ENDIF

  IF ( bcs_visc_jmin .EQ. 1 ) THEN
     j1vsin = 1
  ELSE
     j1vsin = 0
  ENDIF

  IF ( bcs_visc_jmax .EQ. 1 ) THEN
     jmxvsin = 1
  ELSE
     jmxvsin = 0
  ENDIF

  IF ( bcs_visc_kmin .EQ. 1 ) THEN
     k1vsin = 1
  ELSE
     k1vsin = 0
  ENDIF

  IF ( bcs_visc_kmax .EQ. 1 ) THEN
     kmxvsin = 1
  ELSE
     kmxvsin = 0
  ENDIF


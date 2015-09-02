  CALL IO_READ_GRID(gfile, imax_total,jmax_total,kmax_total, scalex,scaley,scalez, x,y,z)
  CALL FDM_INITIALIZE(iunifx, imode_fdm, imax_total, i1bc, scalex, x, dx, wrk1d)
  CALL FDM_INITIALIZE(iunify, imode_fdm, jmax_total, j1bc, scaley, y, dy, wrk1d)
  CALL FDM_INITIALIZE(iunifz, imode_fdm, kmax_total, k1bc, scalez, z, dz, wrk1d)
  IF ( kmax_total .GT. 1 ) THEN; area = scalex*scalez             ! 3D case
  ELSE;                          area = scalex;       ENDIF       ! 2D case
  volume = area*scaley

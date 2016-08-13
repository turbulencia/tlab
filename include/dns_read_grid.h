  grid(1)%nodes => x(:,1); grid(1)%aux => x(:,2:)
  grid(2)%nodes => y(:,1); grid(2)%aux => y(:,2:)
  grid(3)%nodes => z(:,1); grid(3)%aux => z(:,2:)
  dx => x(:,2:) ! to be removed
  dy => y(:,2:)
  dz => z(:,2:)

!  CALL IO_READ_GRID(gfile, imax_total,jmax_total,kmax_total, scalex,scaley,scalez, x,y,z)
  CALL IO_READ_GRID(gfile, grid(1)%size,grid(2)%size,grid(3)%size, grid(1)%scale,grid(2)%scale,grid(3)%scale, &
		    grid(1)%nodes,grid(2)%nodes,grid(3)%nodes)
  scalex = grid(1)%scale ! to be removed 
  scaley = grid(2)%scale 
  scalez = grid(3)%scale
			     
!  CALL FDM_INITIALIZE(iunifx, imode_fdm, imax_total, i1bc, scalex, x, dx, wrk1d)
!  CALL FDM_INITIALIZE(iunify, imode_fdm, jmax_total, j1bc, scaley, y, dy, wrk1d)
!  CALL FDM_INITIALIZE(iunifz, imode_fdm, kmax_total, k1bc, scalez, z, dz, wrk1d)
  CALL FDM_INITIALIZE(iunifx, imode_fdm, imax_total, i1bc, scalex, grid(1)%nodes, grid(1)%aux, wrk1d)
  CALL FDM_INITIALIZE(iunify, imode_fdm, jmax_total, j1bc, scaley, grid(2)%nodes, grid(2)%aux, wrk1d)
  CALL FDM_INITIALIZE(iunifz, imode_fdm, kmax_total, k1bc, scalez, grid(3)%nodes, grid(3)%aux, wrk1d)
!  CALL FDM_INITIALIZE_NEW(imode_fdm, grid(1), wrk1d)
!  CALL FDM_INITIALIZE_NEW(imode_fdm, grid(2), wrk1d)
!  CALL FDM_INITIALIZE_NEW(imode_fdm, grid(3), wrk1d)
  IF ( kmax_total .GT. 1 ) THEN; area = scalex*scalez             ! 3D case
  ELSE;                          area = scalex;       ENDIF       ! 2D case
  volume = area*scaley

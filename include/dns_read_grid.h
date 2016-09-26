  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'

  g(1)%nodes => x(:,1); g(1)%aux => x(:,2:)
  g(2)%nodes => y(:,1); g(2)%aux => y(:,2:)
  g(3)%nodes => z(:,1); g(3)%aux => z(:,2:)

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes)
			     
  CALL FDM_INITIALIZE(imode_fdm, g(1)%size, g(1)%scale, g(1)%uniform, g(1)%periodic, g(1)%nodes, g(1)%aux, wrk1d, g(1))
  CALL FDM_INITIALIZE(imode_fdm, g(2)%size, g(2)%scale, g(2)%uniform, g(2)%periodic, g(2)%nodes, g(2)%aux, wrk1d, g(2))
  CALL FDM_INITIALIZE(imode_fdm, g(3)%size, g(3)%scale, g(3)%uniform, g(3)%periodic, g(3)%nodes, g(3)%aux, wrk1d, g(3))

  area = g(1)%scale
  IF ( g(3)%size .GT. 1 ) area = area *g(3)%scale ! 3D case
  volume = area *g(2)%scale

  dx => x(:,2:) ! to be removed
  dy => y(:,2:)
  dz => z(:,2:)
  scalex = g(1)%scale ! to be removed 
  scaley = g(2)%scale 
  scalez = g(3)%scale

  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z)
			     
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  area = g(1)%scale
  IF ( g(3)%size .GT. 1 ) area = area *g(3)%scale ! 3D case
  volume = area *g(2)%scale

  dx => x(:,2:) ! to be removed
  dy => y(:,2:)
  dz => z(:,2:)
  scalex = g(1)%scale ! to be removed 
  scaley = g(2)%scale 
  scalez = g(3)%scale

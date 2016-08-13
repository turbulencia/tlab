#include "types.h"

PROGRAM TRANSGRID

  USE DNS_TYPES, ONLY : grid_structure
  IMPLICIT NONE

  TYPE(grid_structure), DIMENSION(3) :: g, g_ref
  
  TINTEGER option, direction, n, isize_wrk1d
  CHARACTER*32 ifile,ffile,sfile, file_ref
  LOGICAL flag_exit
  
  TREAL, DIMENSION(:,:), ALLOCATABLE         :: wrk1d
  TREAL, DIMENSION(:),   ALLOCATABLE, TARGET :: x,y,z
  TREAL, DIMENSION(:),   ALLOCATABLE, TARGET :: x_ref,y_ref,z_ref
  TREAL offset, factor

! ###################################################################
! Initialize and read reference data
! ###################################################################
  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'
  
  WRITE(*,'("- Reference grid file ? ", $)')
  READ(*,*) ifile

  ffile = TRIM(ADJUSTL(ifile))//'.trn'
  sfile = TRIM(ADJUSTL(ffile))//'.sts'

  CALL RD_GRIDSIZE(ifile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale)

  isize_wrk1d = MAX(g(1)%size,MAX(g(2)%size,g(3)%size))

  ALLOCATE(x(2*g(1)%size)) ! Allocation of memory is doubled to allow introduction of planes
  ALLOCATE(y(2*g(2)%size))
  ALLOCATE(z(2*g(3)%size))
  g(1)%nodes => x
  g(2)%nodes => y
  g(3)%nodes => z
  ! ALLOCATE(g(1)%nodes(2*g(1)%size)) ! Allocation of memory is doubled to allow introduction of planes
  ! ALLOCATE(g(2)%nodes(2*g(2)%size))
  ! ALLOCATE(g(3)%nodes(2*g(3)%size))
  
  ALLOCATE(wrk1d(isize_wrk1d,3))

  CALL IO_READ_GRID(ifile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes)
  
! #######################################################################
! Main loop
! #######################################################################
  flag_exit = .FALSE.
  
  DO WHILE ( .NOT. flag_exit )
     
     WRITE(*,'(20(/),"TransGrid Main Menu",/)')
     WRITE(*,'(a)')   '0. Dump ASCII values to file'
     WRITE(*,'(a)')   '1. Offsets'
     WRITE(*,'(a)')   '2. Scaling'
     WRITE(*,'(a)')   '3. Drop planes'
     WRITE(*,'(a)')   '4. Introduce planes'
     WRITE(*,'(a)')   '5. Transfer grid between files'
     WRITE(*,'(a,/)') '9. Exit'
     READ(*,*) option 

     IF ( option .LT. 9 ) THEN
        WRITE(*,'(/,"Direction (1 for x, 2 for y, 3 for z) ? ", $)')
        READ(*,*) direction
     ENDIF
     
     SELECT CASE(option)
        
     CASE(0)
        OPEN(23,file=TRIM(ADJUSTL(g(direction)%name))//'.dat')
        DO n = 1,g(direction)%size; WRITE(23,*) g(direction)%nodes(n); ENDDO
        CLOSE(23)

     CASE(1)
        WRITE(*,'(/,"Offset ? ", $)')
        READ(*,*) offset
        g(direction)%nodes(:) = g(direction)%nodes(:) + offset

     CASE(2)
        WRITE(*,'(/,"Scaling factor ? ", $)')
        READ(*,*) factor
        g(direction)%nodes(:) = g(direction)%nodes(1) + ( g(direction)%nodes(:)-g(direction)%nodes(1)) *factor
        g(direction)%scale = g(direction)%scale *factor
        
     CASE(3) ! Dropping planes
        CALL TRANS_DROP_PLANES(g(direction)%size, g(direction)%nodes, g(direction)%scale)

     CASE(4) ! Introducing planes
        CALL TRANS_ADD_PLANES(g(direction)%size, g(direction)%nodes, g(direction)%scale)
        flag_exit = .TRUE. ! Exit directly to avoid memory allocation problems

     CASE(5) ! Transfering
        WRITE(*,'(/,"Reference file to copy the data from ? ", $)')
        READ(*,*) file_ref

        CALL RD_GRIDSIZE(file_ref, g_ref(1)%size,g_ref(2)%size,g_ref(3)%size, g_ref(1)%scale,g_ref(2)%scale,g_ref(3)%scale)
        DO n = 1,3
           IF ( g_ref(n)%size .GT. 2*g(n)%size ) THEN
              WRITE(*,*) 'Error. Reference grid too big.'
              STOP
           ENDIF
        ENDDO
        ALLOCATE(x_ref(g_ref(1)%size)); g_ref(1)%nodes => x_ref
        ALLOCATE(y_ref(g_ref(2)%size)); g_ref(2)%nodes => y_ref
        ALLOCATE(z_ref(g_ref(3)%size)); g_ref(3)%nodes => z_ref
           ! ELSE
           !    ALLOCATE(g_ref(n)%nodes(g_ref(n)%size))
        CALL IO_READ_GRID(ifile, g_ref(1)%size,g_ref(2)%size,g_ref(3)%size, g_ref(1)%scale,g_ref(2)%scale,g_ref(3)%scale, g_ref(1)%nodes,g_ref(2)%nodes,g_ref(3)%nodes)

        g(direction)%size  = g_ref(direction)%size
        g(direction)%scale = g_ref(direction)%scale
        g(direction)%nodes(1:g(direction)%size) = g_ref(direction)%nodes(1:g(direction)%size)

        flag_exit = .TRUE. ! Exit directly to avoid memory allocation problems

     CASE(9)
        flag_exit = .TRUE.
        
     END SELECT

  END DO
  
! #######################################################################
! Statistics and writing the data
! #######################################################################
  IF ( option .GT. 0 ) THEN
     DO direction = 1,3
        CALL TRANS_DATA(sfile, g(direction), wrk1d(:,1),wrk1d(:,2))
     ENDDO
     CALL IO_WRITE_GRID(ffile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes)
  ENDIF
  
  STOP
  
END PROGRAM TRANSGRID

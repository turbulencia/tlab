  ALLOCATE(pbackground(g(2)%size))
  ALLOCATE(rbackground(g(2)%size))
  ALLOCATE(bbackground(g(2)%size))
  ALLOCATE(tbackground(g(2)%size))
  ALLOCATE(epbackground(g(2)%size))

  WRITE(str,*) g(1)%inb_grid; line = 'Allocating array x of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) g(1)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(x(g(1)%size,g(1)%inb_grid),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for x.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  
  WRITE(str,*) g(2)%inb_grid; line = 'Allocating array y of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) g(2)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(y(g(2)%size,g(2)%inb_grid),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for y.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) g(3)%inb_grid; line = 'Allocating array z of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) g(3)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(z(g(3)%size,g(3)%inb_grid),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for z.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  IF ( iread_flow .EQ. 1 ) THEN
     WRITE(str,*) inb_flow_array; line = 'Allocating array flow  of size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(q(isize_field,inb_flow_array),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for flow.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

  IF ( iread_scal .EQ. 1 ) THEN
     WRITE(str,*) inb_scal_array; line = 'Allocating array scal  of size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(s(isize_field,inb_scal_array),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for scal.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

  IF ( inb_txc .GT. 0 ) THEN
     WRITE(str,*) inb_txc; line = 'Allocating array txc   of size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) isize_txc_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(txc(isize_txc_field,inb_txc),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for txc.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

  WRITE(str,*) isize_wrk3d; line = 'Allocating array wrk3d of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(wrk3d(isize_wrk3d),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for wrk3d.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

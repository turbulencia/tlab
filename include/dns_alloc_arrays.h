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

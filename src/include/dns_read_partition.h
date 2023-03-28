  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Partition', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec2)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Missing input [PostProcessing.Partition] in dns.ini.')
     CALL TLAB_STOP(DNS_ERROR_INVALOPT)
#else
     WRITE(*,*) 'Intermittency function for the conditioning ?'
     WRITE(*,*) ' 0. None'
     WRITE(*,*) ' 1. Based on external file'
     WRITE(*,*) ' 2. Based on scalar'
     WRITE(*,*) ' 3. Based on vorticity'
     WRITE(*,*) ' 4. Based on scalar gradient'
     WRITE(*,*) ' 5. Based on vertical velocity'
     WRITE(*,*) ' 6. Based on scalar fluctuation'
     WRITE(*,*) ' 7. Based on scalar vertical turbulent flux'
     WRITE(*,*) ' 8. Based on potential vorticity'
     READ(*,*) opt_cond

     IF ( opt_cond .EQ. 3 .OR. opt_cond .EQ. 4 ) THEN
        WRITE(*,*) 'Threshold based on absolute (0) or relative (1) values?'
        READ(*,*) opt_cond_relative
     ENDIF

     IF ( opt_cond .GE. 2 ) THEN
        WRITE(*,*) 'Intermittency thresholds ?'
        READ(*,'(A512)') sRes
        igate_size = igate_size_max
        CALL LIST_REAL(sRes, igate_size, gate_threshold)
     ENDIF

#endif
  ELSE
     idummy = 1
     opt_cond = INT(opt_vec2(idummy))
     IF ( opt_cond .EQ. 3 .OR. opt_cond .EQ. 4 ) THEN
        idummy = idummy + 1
        opt_cond_relative = INT(opt_vec2(idummy))
     ENDIF
     IF ( opt_cond .GE. 2 ) THEN
        idummy = idummy + 1
        igate_size                   = iopt_size - idummy + 1
        gate_threshold(1:igate_size) = opt_vec2(idummy:iopt_size)
     ENDIF

  ENDIF

  IF      ( opt_cond .EQ. 3 ) THEN; gate_threshold = gate_threshold**2.0_wp; ! using enstrophy
  ELSE IF ( opt_cond .EQ. 4 ) THEN; gate_threshold = gate_threshold**2.0_wp; ! using scalar dissipation
  ENDIF

  IF ( igate_size .GT. 0 ) THEN
     IF ( opt_cond .EQ. 7 ) THEN
        igate_size = INT(2.0**real(igate_size,wp)) ! double conditioning
     ELSE
        CALL SORT_REAL(igate_size,gate_threshold)
        igate_size = igate_size+1 ! # of gate levels is +1 number of thresholds between them
     ENDIF
  ENDIF

  IF ( igate_size .GT. igate_size_max ) THEN
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Not enough memory for gate_threshold.')
     CALL TLAB_STOP(DNS_ERROR_ALLOC)
  ENDIF

  IF      ( opt_cond .EQ. 2 ) THEN
     inb_txc    = MAX(inb_txc,1)
     iread_scal = .true.
  ELSE IF ( opt_cond .EQ. 3 ) THEN
     inb_txc    = MAX(inb_txc,3)
     iread_flow = .true.
  ELSE IF ( opt_cond .EQ. 4 ) THEN
     inb_txc    = MAX(inb_txc,3)
     iread_scal = .true.
  ELSE IF ( opt_cond .EQ. 5 ) THEN
     inb_txc    = MAX(inb_txc,1)
     iread_flow = .true.
  ELSE IF ( opt_cond .EQ. 6 ) THEN
     inb_txc    = MAX(inb_txc,1)
     iread_scal = .true.
  ELSE IF ( opt_cond .EQ. 7 ) THEN
     inb_txc    = MAX(inb_txc,1)
     iread_scal = .true.
     iread_flow = .true.
  ELSE IF ( opt_cond .EQ. 8 ) THEN
     inb_txc    = MAX(inb_txc,5)
     iread_scal = .true.
     iread_flow = .true.
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'Partition', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec2)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef PARALLEL
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Missing input [PostProcessing.Partition] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT) 
#else
     WRITE(*,*) 'Intermittency function for the conditioning ?'
     WRITE(*,*) ' 0. None'
     WRITE(*,*) ' 1. Based on external file'
     WRITE(*,*) ' 2. Based on scalar'
     WRITE(*,*) ' 3. Based on vorticity'
     WRITE(*,*) ' 4. Based on scalar gradient'
     WRITE(*,*) ' 5. Based on vertical velocity'
     READ(*,*) opt_cond

     IF ( opt_cond .EQ. 3 .OR. opt_cond .EQ. 4 ) THEN
        WRITE(*,*) 'Threshold based on relative (1) or absolute (2) values?'
        READ(*,*) opt_threshold
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
     opt_cond = DINT(opt_vec2(idummy))
     IF ( opt_cond .EQ. 3 .OR. opt_cond .EQ. 4 ) THEN
        idummy = idummy + 1
        opt_threshold = DINT(opt_vec2(idummy))
     ENDIF
     IF ( opt_cond .GE. 2 ) THEN
        idummy = idummy + 1
        igate_size                   = iopt_size - idummy + 1
        gate_threshold(1:igate_size) = opt_vec2(idummy:iopt_size)
     ENDIF

  ENDIF

  IF      ( opt_cond .EQ. 3 ) THEN; gate_threshold = gate_threshold**C_2_R; ! using enstrophy
  ELSE IF ( opt_cond .EQ. 4 ) THEN; gate_threshold = gate_threshold**C_2_R; ! using scalar dissipation
  ENDIF

  IF ( igate_size .GT. 0 ) THEN
     CALL SORT_REAL(igate_size,gate_threshold)
     igate_size = igate_size+1 ! # of gate levels is +1 number of thresholds between them
  ENDIF

  IF ( igate_size .GT. igate_size_max ) THEN
     CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Not enough memory for igate_vec.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF


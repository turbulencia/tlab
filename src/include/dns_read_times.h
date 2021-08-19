  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Files', '-1', sRes)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef PARALLEL
#else
     WRITE(*,*) 'Iteration numbers ?'
     READ(*,'(A512)') sRes
#endif
  ENDIF
  itime_size = itime_size_max
  CALL LIST_INTEGER(sRes, itime_size, itime_vec)

  IF ( itime_vec(1) .LT. 0 ) THEN ! Check
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Missing input [PostProcessing.Files] in dns.ini.')
     CALL TLAB_STOP(DNS_ERROR_INVALOPT)
  ENDIF

     IF      ( opt_cond .EQ. 2 ) THEN ! Based on scalar 
        txc(1:isize_field,1) = s(1:isize_field,opt_cond_scal)

     ELSE IF ( opt_cond .EQ. 3 ) THEN ! Based on vorticity
        CALL IO_WRITE_ASCII(lfile,'Calculating vorticity...')
        CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)

     ELSE IF ( opt_cond .EQ. 4 ) THEN ! Based on scalar gradient
        CALL IO_WRITE_ASCII(lfile,'Calculating scalar gradient...')
        CALL FI_GRADIENT(imax,jmax,kmax, s, txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)

     ELSE IF ( opt_cond .EQ. 5 ) THEN ! Based on vertical velocity
        txc(1:isize_field,1) = q(1:isize_field,2)
       
     ELSE IF ( opt_cond .EQ. 6 .OR. opt_cond .EQ. 7 ) THEN ! Based on scalar fluctuation
        CALL IO_WRITE_ASCII(lfile,'Calculating scalar fluctuation...')
        txc(1:isize_field,1) = s(1:isize_field,opt_cond_scal)
        CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(1,1))

     ELSE IF ( opt_cond .EQ. 8 ) THEN ! Based on potential vorticity
        CALL IO_WRITE_ASCII(lfile,'Calculating potential vorticity...')
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)

        txc(1:isize_field,4) = s(1:isize_field,opt_cond_scal)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) =                       txc(1:isize_field,1)*txc(1:isize_field,5) 
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,2)*txc(1:isize_field,5) 
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,3)*txc(1:isize_field,5) 

        txc(1:isize_field,1) = txc(1:isize_field,1)*txc(1:isize_field,1)
        txc(1:isize_field,1) = LOG(txc(1:isize_field,1)+C_SMALL_R)

     ENDIF

     IF      ( opt_cond .EQ. 1 ) THEN ! External file
        WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname)); params_size = 2
        CALL IO_READ_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
        igate_size = INT(params(2))
        IF ( igate_size .GT. igate_size_max ) THEN
           CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Not enough memory for igate_vec.')
           CALL DNS_STOP(DNS_ERROR_ALLOC)
        ENDIF

     ELSE IF ( opt_cond .EQ. 7 ) THEN ! double conditioning; flux
        DO ij = 1,isize_field
           IF      ( txc(ij,1) .GT. C_0_R .AND. q(ij,2) .GE. C_0_R ) THEN; gate(ij) = 1;
           ELSE IF ( txc(ij,1) .LE. C_0_R .AND. q(ij,2) .GT. C_0_R ) THEN; gate(ij) = 2;
           ELSE IF ( txc(ij,1) .LT. C_0_R .AND. q(ij,2) .LE. C_0_R ) THEN; gate(ij) = 3;
           ELSE;                                                           gate(ij) = 4; ENDIF
        ENDDO      

     ELSE                             ! Local file
        IF ( opt_threshold .EQ. 1 ) THEN ! case of threshold relative to maximum
           CALL MINMAX(imax,jmax,kmax, txc(1,1), umin,umax)
           WRITE(sRes,100) umin; sRes='Min/Max '//TRIM(ADJUSTL(sRes))
           WRITE(str,100) umax; sRes=TRIM(ADJUSTL(sRes))//'/'//TRIM(ADJUSTL(str))//'.'
           CALL IO_WRITE_ASCII(lfile,sRes)

           IF ( (umax-umin) .GT. C_0_R ) THEN
              txc(:,1) = (txc(:,1)-umin) / (umax-umin)
           ELSE
              CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Homogeneous field, nothing to partition.')
              CALL DNS_STOP(DNS_ERROR_INVALOPT)
           ENDIF
        ENDIF

! define gate field
        DO ij = 1,isize_field
           DO n = 1,igate_size-1
              IF ( txc(ij,1) .LT. gate_threshold(n) ) EXIT
           ENDDO
           gate(ij) = INT(n,KIND=1) ! note that gate can get -- correctly -- the value igate_size
        ENDDO
        
     ENDIF
     
     DO n = 1,igate_size
        igate_vec(n) = INT(n,KIND=1)
     ENDDO

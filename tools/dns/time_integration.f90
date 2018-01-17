#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"
#include "avgij_map.h"

!########################################################################
!# DESCRIPTION
!#
!# Performing the time integration over a given number of steps
!#
!########################################################################
SUBROUTINE TIME_INTEGRATION(q,hq, s,hs, q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
     l_q, l_hq, l_txc, l_tags, l_comm)
  
  USE DNS_CONSTANTS, ONLY : tag_flow, tag_scal, tag_part, tag_traj, lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_scal_array, inb_flow_array
  USE DNS_GLOBAL, ONLY : isize_particle, inb_particle, inb_particle_txc
  USE DNS_GLOBAL, ONLY : imode_sim, imode_eqns
  USE DNS_GLOBAL, ONLY : icalc_flow, icalc_scal, icalc_part
  USE DNS_GLOBAL, ONLY : rbackground, g
  USE DNS_GLOBAL, ONLY : itransport, visc
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : nstatavg, nstatavg_points
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL 
  USE DNS_TOWER
  USE LAGRANGE_GLOBAL, ONLY : itrajectory
  USE BOUNDARY_INFLOW
  USE PARTICLE_TRAJECTORIES
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
#include "integers.h"

  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(*)             :: txc, vaux
  TREAL, DIMENSION(*)             :: q_inf, s_inf
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  INTEGER(8), DIMENSION(isize_particle)                  :: l_tags
  TREAL,      DIMENSION(isize_particle,inb_particle    ) :: l_q, l_hq
  TREAL,      DIMENSION(isize_particle,inb_particle_txc) :: l_txc
  TREAL,      DIMENSION(*)                               :: l_comm

  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER is, iq, ip
  TINTEGER idummy, splanes_i(5), splanes_j(5), splanes_k(5)
  CHARACTER*32 fname, varname(1)
  CHARACTER*250 line1
  LOGICAL flag_save
  
! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w, e, rho, p, T, vis

! ###################################################################
! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)
     
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)

  ENDIF

! Sizes information for saving planes  
  idummy     = kmax*jmax*nplanes_i*(inb_flow_array+inb_scal_array)
  splanes_i  = (/idummy,1,idummy,1,1/)
  idummy     = imax*jmax*nplanes_k*(inb_flow_array+inb_scal_array)
  splanes_k  = (/idummy,1,idummy,1,1/)
  idummy     = imax*kmax*nplanes_j*(inb_flow_array+inb_scal_array)
  splanes_j  = (/idummy,1,idummy,1,1/)
  varname    = (/''/)
  
! ###################################################################
! Loop on iterations: itime counter
! ###################################################################
  itime = nitera_first
  
  WRITE(line1,*) itime; line1 = 'Starting time integration at It'//TRIM(ADJUSTL(line1))//'.'
  CALL IO_WRITE_ASCII(lfile,line1)

  DO WHILE ( itime .LT. nitera_last )

     CALL TIME_RUNGEKUTTA(q,hq, s,hs, q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
          l_q, l_hq, l_txc, l_tags, l_comm)

     itime = itime + 1
     rtime = rtime + dtime

! -----------------------------------------------------------------------
     IF ( MOD(itime-nitera_first,FilterDomainStep) .EQ. 0 ) THEN
        IF ( MOD(itime-nitera_first,nitera_stats)  .EQ. 0 ) THEN; flag_save = .TRUE.
        ELSE;                                                     flag_save = .FALSE.; ENDIF
        CALL DNS_FILTER(flag_save, q,s, txc, vaux, wrk1d,wrk2d,wrk3d)
     ENDIF

     ! This should be integrated into the inflow buffer, as the filter contribution
     IF ( MOD(itime-nitera_first,FilterInflowStep) .EQ. 0 ) THEN ! Inflow filter in spatial mode
        CALL BOUNDARY_INFLOW_FILTER(vaux(vindex(VA_BCS_VI)), q,s, txc, wrk1d,wrk2d,wrk3d)
     ENDIF

! -----------------------------------------------------------------------
     IF ( iviscchg .EQ. 1 ) THEN ! Change viscosity if necessary
        visc = visc - dtime*visctime
        IF ( ( (visc .LT. viscstop) .AND. (viscstart .GT. viscstop) )   .OR. &
             ( (visc .GT. viscstop) .AND. (viscstart .LT. viscstop) ) ) THEN
           iviscchg = 0; visc = viscstop
        ENDIF
     ENDIF

! -----------------------------------------------------------------------
     CALL TIME_COURANT(q,s, wrk3d)

! ###################################################################
! The rest: Logging, postprocessing and saving
! ###################################################################
     IF ( MOD(itime-nitera_first,nitera_log) .EQ. 0 .OR. INT(logs_data(1)) .NE. 0 ) THEN ! Log files
        CALL DNS_LOGS(i2)
#ifdef LES
        IF ( iles .EQ. 1 ) CALL LES_LOGS(i2)
#endif
     ENDIF

! -----------------------------------------------------------------------
! Accumulate statistics in spatially evolving cases
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. MOD(itime-nitera_first,nitera_stats_spa) .EQ. 0 ) THEN
        nstatavg_points = nstatavg_points + g(3)%size
        CALL DNS_SAVE_AVGIJ(rho,u,v,w,p,vis,T, hq,txc, vaux(vindex(VA_MEAN_WRK)), wrk2d,wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
           CALL DNS_SAVE_SCBDGIJ(rho,u,v,w,p,s,vis, hq,txc, vaux(vindex(VA_MEAN_WRK)+MA_MOMENTUM_SIZE*nstatavg*jmax), wrk2d,wrk3d)
        ENDIF
     ENDIF

! -----------------------------------------------------------------------
     IF ( tower_mode .EQ. 1 ) THEN 
        CALL DNS_TOWER_ACCUMULATE(q,1,wrk1d)
        CALL DNS_TOWER_ACCUMULATE(s,2,wrk1d)
     ENDIF

! -----------------------------------------------------------------------
     IF ( itrajectory .NE. LAG_TRAJECTORY_NONE ) THEN
        CALL PARTICLE_TRAJECTORIES_ACCUMULATE(q,s, txc, l_q,l_hq,l_txc,l_tags,l_comm, wrk2d,wrk3d)
     END IF

! -----------------------------------------------------------------------
     IF ( MOD(itime-nitera_first,nitera_stats) .EQ. 0 ) THEN ! Calculate statistics
        IF     ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
           CALL STATS_TEMPORAL_LAYER(q,s,hq, txc, vaux, wrk1d,wrk2d,wrk3d)
           IF ( icalc_part .EQ. 1 ) THEN
              CALL STATS_TEMPORAL_LAGRANGIAN(q,s,hq, l_q,l_hq,l_txc,l_tags,l_comm, txc, vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)
           ENDIF
        ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
           CALL STATS_SPATIAL_LAYER(vaux, txc, wrk1d,wrk2d)
        ENDIF
     ENDIF
     
! -----------------------------------------------------------------------
     IF ( MOD(itime-nitera_first,nitera_save) .EQ. 0 .OR. &      ! Save restart files
          itime .EQ. nitera_last .OR. INT(logs_data(1)) .NE. 0 ) THEN ! Secure that one restart file is saved
        
        IF ( icalc_flow .EQ. 1 ) THEN
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, isize_field, q, wrk3d)
        ENDIF
        IF ( icalc_scal .EQ. 1 ) THEN
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, isize_field, s, wrk3d)
        ENDIF
        
        IF ( tower_mode .EQ. 1 ) THEN
           CALL DNS_TOWER_WRITE(wrk3d) 
        ENDIF

        IF ( icalc_part .EQ. 1 ) THEN
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_PARTICLE(fname, l_tags, l_q)
           IF ( itrajectory .NE. LAG_TRAJECTORY_NONE ) THEN
              WRITE(fname,*) itime; fname =TRIM(ADJUSTL(tag_traj))//TRIM(ADJUSTL(fname))
              CALL PARTICLE_TRAJECTORIES_WRITE(fname, wrk3d)
           END IF
        END IF

        IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. nitera_stats_spa .GT. 0 ) THEN ! Spatial; running averages
           WRITE(fname,*) itime; fname='st'//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_AVGIJ(fname, vaux(vindex(VA_MEAN_WRK)))
        ENDIF
        
     ENDIF

! -----------------------------------------------------------------------
     IF ( MOD(itime-nitera_first,nitera_pln) .EQ. 0 ) THEN ! Save planes
        IF ( nplanes_k .GT. 0 ) THEN
           CALL REDUCE_Z_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, nplanes_k,planes_k, txc)
           WRITE(fname,*) itime; fname = 'planesK.'//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_SUBARRAY4(MPIO_SUBARRAY_PLANES_XOY, fname, varname, txc, splanes_k, hq)
        ENDIF

        IF ( nplanes_j .GT. 0 ) THEN
           IF ( nplanes_j_aux .GT. 0 ) THEN ! Calculate integrals
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
                 ip = 1
                 DO iq = 1,inb_flow_array
                    CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, q(1,iq), wrk3d(ip), wrk1d,hq(1,1)) !hq is aux variable
                    ip = ip + imax*kmax
                 ENDDO
                 DO is = 1,inb_scal_array
                    CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, s(1,is), wrk3d(ip), wrk1d,hq(1,1))
                    ip = ip + imax*kmax
                 ENDDO
              ELSE
                 wrk3d(1:(inb_flow_array+inb_scal_array)*imax*kmax) = C_0_R
              ENDIF
           ENDIF
           CALL REDUCE_Y_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, wrk3d, nplanes_j,nplanes_j_aux,planes_j, txc)
           WRITE(fname,*) itime; fname = 'planesJ.'//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_SUBARRAY4(MPIO_SUBARRAY_PLANES_XOZ, fname, varname, txc, splanes_j, hq)
        ENDIF

        IF ( nplanes_i .GT. 0 ) THEN
           CALL REDUCE_X_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, nplanes_i,planes_i, txc)
           WRITE(fname,*) itime; fname = 'planesI.'//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_SUBARRAY4(MPIO_SUBARRAY_PLANES_ZOY, fname, varname, txc, splanes_i, hq)
        ENDIF

     ENDIF
     
! -----------------------------------------------------------------------
     IF ( INT(logs_data(1)) .NE. 0 ) CALL DNS_STOP(INT(logs_data(1)))
     
  ENDDO

  RETURN
END SUBROUTINE TIME_INTEGRATION

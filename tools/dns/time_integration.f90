#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/06/12 - J.P. Mellado
!#              Derived from old version for the new energy formulation
!# 2007/08/02 - J.P. Mellado
!#              Filter has been studied. It has to be done in the
!#              pressure to get good dilatation evolution
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE TIME_INTEGRATION(x,y,z,dx,dy,dz, q,hq,s,hs, &
     x_inf,y_inf,z_inf,q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
     l_q, l_hq, l_txc, l_tags, l_comm, l_trajectories, l_trajectories_tags)
  
  USE DNS_CONSTANTS, ONLY : tag_flow, tag_scal, lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_vars, isize_particle, inb_particle
  USE DNS_GLOBAL, ONLY : i1bc,j1bc,k1bc, imode_fdm, imode_flow, imode_sim, imode_eqns
  USE DNS_GLOBAL, ONLY : icalc_flow, icalc_scal, icalc_particle
  USE DNS_GLOBAL, ONLY : itransport, visc
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : nspa_rest, nspa_step, nstatplnvars, iupdate_stat
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL 
  USE DNS_TOWER
  USE LAGRANGE_GLOBAL, ONLY : icalc_trajectories
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

  TREAL, DIMENSION(*)             :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(*)             :: txc, vaux
  TREAL, DIMENSION(*)             :: x_inf, y_inf, z_inf, q_inf, s_inf
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  TREAL,      DIMENSION(isize_particle,inb_particle) :: l_q, l_hq
  TREAL,      DIMENSION(*)                           :: l_comm, l_txc, l_trajectories
  INTEGER(8), DIMENSION(*)                           :: l_tags, l_trajectories_tags

  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER it_loc_data
  TINTEGER icount_stat, is
  TINTEGER idummy, sizes_l(4)
  CHARACTER*32 fname, varname(1)
  CHARACTER*250 line1
  LOGICAL flag_save
  
! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u, v, w, e, rho, p, T, vis

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

! Spatial mode
  it_loc_data  = nitera_first + nspa_rest*nspa_step ! Save statistical point information step control
  iupdate_stat = nitera_first + nspa_step           ! Update statistical information
  icount_stat = 1

! Information for saving planes  
  idummy = imax*kmax*npln_j*(inb_flow+inb_scal)
  sizes_l = (/idummy,idummy,0,1/)
  varname = (/'0'/)
  
! ###################################################################
! Loop on iterations: itime counter
! ###################################################################
  itime = nitera_first
  
  WRITE(line1,*) itime; line1 = 'Starting time integration at It'//TRIM(ADJUSTL(line1))//'.'
  CALL IO_WRITE_ASCII(lfile,line1)

  DO WHILE ( itime .LT. nitera_last )

     CALL TIME_RUNGEKUTTA(x,y,z,dx,dy,dz, q,hq,s,hs, &
          x_inf,y_inf,z_inf,q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
          l_q, l_hq, l_txc, l_tags, l_comm)

     itime = itime + 1
     rtime = rtime + dtime

     IF ( ifilt_step .GT. 0 .AND. MOD(itime-nitera_first,ifilt_step) .EQ. 0 ) THEN
        IF ( nitera_stats .GT. 0 .AND. MOD(itime-nitera_first,nitera_stats) .EQ. 0 ) THEN; flag_save = .TRUE.
        ELSE;                                                                              flag_save = .FALSE.; ENDIF
        CALL DNS_FILTER(flag_save, y, dx,dy,dz, q,s, txc, vaux, wrk1d,wrk2d,wrk3d)
     ENDIF

! -----------------------------------------------------------------------
! Specfics of spatially evolving cases
! -----------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN        
        IF ( ifilt_inflow .EQ. 1 .AND. MOD(itime,ifilt_inflow_step) .EQ. 0 ) THEN ! Inflow filter
           CALL BOUNDARY_INFLOW_FILTER(rho, txc, vaux(vindex(VA_BCS_VI))            , wrk1d,wrk2d,wrk3d)
           CALL BOUNDARY_INFLOW_FILTER(u,   txc, vaux(vindex(VA_BCS_VI)+jmax*kmax)  , wrk1d,wrk2d,wrk3d)
           CALL BOUNDARY_INFLOW_FILTER(v,   txc, vaux(vindex(VA_BCS_VI)+jmax*kmax*2), wrk1d,wrk2d,wrk3d)
           CALL BOUNDARY_INFLOW_FILTER(w,   txc, vaux(vindex(VA_BCS_VI)+jmax*kmax*3), wrk1d,wrk2d,wrk3d)
           CALL BOUNDARY_INFLOW_FILTER(p,   txc, vaux(vindex(VA_BCS_VI)+jmax*kmax*4), wrk1d,wrk2d,wrk3d)
           IF ( icalc_scal .EQ. 1 .AND. ifilt_scalar .EQ. 1 ) THEN
              DO is = 1,inb_scal
                 CALL BOUNDARY_INFLOW_FILTER(s(1,is),txc,vaux(vindex(VA_BCS_VI)+jmax*kmax*(4+is)), wrk1d,wrk2d,wrk3d)
              ENDDO
           ENDIF
           
! recalculation of p and T
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_AIRWATER_RP(imax,jmax,kmax, s, p, rho, T, wrk3d)
           ELSE
              CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, T)
           ENDIF
           CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, T, e)
           IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)
! This recalculation of T and p is made to make sure that the same numbers are
! obtained in statistics postprocessing as in the simulation; avg* files
! can then be compared with diff command.
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
              CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
           ENDIF
           
        ENDIF
     ENDIF
! -----------------------------------------------------------------------

     CALL TIME_COURANT(dx,dy,dz, q,s, wrk2d,wrk3d)

! ###################################################################
! The rest: Logging, postprocessing and saving
! ###################################################################
     IF ( nitera_log .GT. 0 .AND. MOD(itime-nitera_first,nitera_log) .EQ. 0 ) THEN ! Log files
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
             imode_eqns .EQ. DNS_EQNS_ANELASTIC     )THEN
           CALL FI_INVARIANT_P(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                dx,dy,dz, q(1,1),q(1,2),q(1,3), hq(1,1),hq(1,2), wrk1d,wrk2d,wrk3d)
           CALL MINMAX(imax,jmax,kmax, hq(1,1), logs_data(11),logs_data(10))
           logs_data(10)=-logs_data(10); logs_data(11)=-logs_data(11)
        ENDIF
        CALL DNS_LOGS(i2)
#ifdef LES
        IF ( iles .EQ. 1 ) CALL LES_LOGS(i2)
#endif
     ENDIF

! -----------------------------------------------------------------------
! Specfics of spatially evolving cases
! -----------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
        IF ( itime .EQ. iupdate_stat ) THEN ! Running statistics 
           CALL DNS_SPATIAL_STATS_RUN(icount_stat, dx,dy,dz, q,hq,s, txc, vaux, wrk1d,wrk2d,wrk3d)
           
           icount_stat  = icount_stat  + 1
           iupdate_stat = iupdate_stat + nspa_step
        ENDIF
        
! line and plane data may be written to disk at a different frequency from averages, 
! this latter is done along with flow/scalar restart files
        IF ( itime .EQ. it_loc_data ) THEN
           ! IF ( frunline .EQ. 1 ) THEN
           !    WRITE(fname,*) itime; fname='ln'//TRIM(ADJUSTL(fname))
           !    CALL DNS_WRITE_LINE(fname, isize_field, inb_vars, &
           !         vaux(vindex(VA_TIMES)), vaux(vindex(VA_LINE_SPA_WRK)), wrk3d)
           ! ENDIF
           ! IF ( frunplane .EQ. 1 ) THEN
           !    WRITE(fname,*) itime; fname='pl'//TRIM(ADJUSTL(fname))
           !    CALL DNS_WRITE_PLANE(fname, isize_field, nstatplnvars, &
           !         vaux(vindex(VA_TIMES)), vaux(vindex(VA_PLANE_SPA_WRK)), wrk3d)
           ! ENDIF

           icount_stat = 1
           it_loc_data = it_loc_data + nspa_rest*nspa_step
        ENDIF
      
     ENDIF
! -----------------------------------------------------------------------

     IF ( iviscchg .EQ. 1 ) THEN ! Change viscosity if necessary
        visc = visc - dtime*visctime
        IF ( ( (visc .LT. viscstop) .AND. (viscstart .GT. viscstop) )   .OR. &
             ( (visc .GT. viscstop) .AND. (viscstart .LT. viscstop) ) ) THEN
           iviscchg = 0; visc = viscstop
        ENDIF
     ENDIF

     IF ( tower_mode .EQ. 1 ) THEN 
        CALL DNS_TOWER_ACCUMULATE(q,1,dx,dy,dz,wrk1d)
        CALL DNS_TOWER_ACCUMULATE(s,2,dx,dy,dz,wrk1d)
     ENDIF

     IF ( icalc_particle .EQ. 1 ) THEN ! Lagrangian
        CALL PARTICLE_TIME_RESIDENCE(dtime, l_q, l_hq)
        CALL PARTICLE_TIME_LIQUID_CLIPPING(s,wrk1d,wrk2d,wrk3d,x ,y, z, l_txc, l_tags, l_hq, l_q)
     ENDIF

     IF ( icalc_trajectories .EQ. 1 ) THEN ! Lagrangian
        WRITE(fname,*) itime; fname = 'trajectories.'//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_TRAJECTORIES(fname,l_q,l_tags, l_trajectories, l_trajectories_tags, wrk3d,txc,itime, nitera_last, nitera_save, nitera_first)
     END IF

     IF ( nitera_stats .GT. 0 .AND. MOD(itime-nitera_first,nitera_stats) .EQ. 0 ) THEN ! Calculate statistics
        IF     ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
           IF ( imode_flow .EQ. DNS_FLOW_ISOTROPIC ) THEN ! TO BE DEVELOPED
           ELSE 
              CALL STATS_TEMPORAL_LAYER(x,y,z,dx,dy,dz, q,s,hq, txc, vaux, wrk1d,wrk2d,wrk3d)
              IF ( icalc_particle .EQ. 1 ) THEN ! Lagrangian
                 CALL STATS_TEMPORAL_LAGRANGIAN(x,y,z,dx,dy,dz, q,s, l_q,l_hq,l_txc,l_tags, txc, vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)
              ENDIF
           ENDIF
        ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
           CALL STATS_SPATIAL_LAYER(x,y,z,dx,dy,dz, vaux, txc, wrk1d,wrk2d,wrk3d)
        ENDIF
     ENDIF
     
     IF ( nitera_save .GT. 0 .AND. MOD(itime-nitera_first,nitera_save) .EQ. 0 ) THEN ! Save restart files
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

        IF ( icalc_particle .EQ. 1 ) THEN ! Lagrangian
           WRITE(fname,*) itime; fname = 'particle.'//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_PARTICLE(fname, l_q)
           WRITE(fname,*) itime; fname = 'particle_id.'//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_PARTICLE_TAGS(fname, l_tags)
        END IF

        IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. frunstat .EQ. 1 ) THEN ! Spatial; running averages
           WRITE(fname,*) itime; fname='st'//TRIM(ADJUSTL(fname))
           CALL DNS_WRITE_AVGIJ(fname, vaux(vindex(VA_MEAN_WRK)))
        ENDIF
        
     ENDIF

     IF ( nitera_pln .GT. 0 .AND. MOD(itime-nitera_first,nitera_pln) .EQ. 0 ) THEN ! Save restart files
        IF ( npln_j .GT. 0 ) THEN ! Saving horizontal planes
           CALL REDUCE_Y_ALL(imax,jmax,kmax, inb_flow,q, inb_scal,s, npln_j,pln_j, wrk3d)
           WRITE(fname,*) itime; fname = 'planes.'//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_SUBARRAY4(3,sizes_l,fname,varname,wrk3d,ims_subarray_pln_j,txc(1))
        ENDIF
     ENDIF
     
  ENDDO

  RETURN
END SUBROUTINE TIME_INTEGRATION

#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/03/25 - J.P. Mellado
!#              Arrays are expanded in direction perpendicular to the
!#              boundary.
!# 2007/06/14 - J.P. Mellado
!#              Moved to BOUNDARY_INIT, add reference pressure calculations
!#              and use energy formulation.
!# 2007/07/05 - J.P. Mellado
!#              Modification for AIRWATER case in ref. pressure calculations 
!# 2007/11/16 - J.P. Mellado
!#              Reference planes for BCs added.
!#
!########################################################################
!# DESCRIPTION
!#
!# The evolving buffer zone case should be revisited. Similarly, the
!# treatment of the VI zone should be checked and homogenized.
!#
!# The buffer regions are stored now in terms of rho, rho*v,
!# rho*(e+v^2/2) or rho*e, and rho*Z_i
!#
!# The buffer files need to written without header information, so that
!# the header global variables (like itime) are not updated within this
!# routine
!#
!########################################################################
SUBROUTINE BOUNDARY_INIT(buffer_ht, buffer_hb, buffer_vi, buffer_vo, &
     bcs_ht,bcs_hb,bcs_vi,bcs_vo, q,s, txc, wrk3d)

  USE DNS_CONSTANTS, ONLY : tag_flow,tag_scal, lfile, efile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal,inb_vars, isize_field
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : itime
  USE DNS_GLOBAL,    ONLY : mach, p_init
  USE DNS_GLOBAL,    ONLY : imode_eqns,imode_sim, icalc_scal
  USE DNS_GLOBAL,    ONLY : area, diam_u,ycoor_u
  USE THERMO_GLOBAL, ONLY : imixture, gama0
  USE DNS_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(imax,         buff_nps_jmin,kmax,*) :: buffer_hb
  TREAL, DIMENSION(imax,         buff_nps_jmax,kmax,*) :: buffer_ht
  TREAL, DIMENSION(buff_nps_imin,jmax,         kmax,*) :: buffer_vi
  TREAL, DIMENSION(buff_nps_imax,jmax,         kmax,*) :: buffer_vo

  TREAL, DIMENSION(imax,kmax,*)      :: bcs_hb, bcs_ht
  TREAL, DIMENSION(jmax,kmax,*)      :: bcs_vi, bcs_vo
  TREAL, DIMENSION(imax*jmax*kmax,*) :: q, s, txc
  TREAL, DIMENSION(*)                :: wrk3d

  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER i,j,k, iq,is, ip, idummy,io_sizes(5)
  TREAL AVG1V1D, AVG_IK, prefactor, dummy
  TREAL diam_loc, thick_loc, ycenter, r1, r05, FLOW_JET_TEMPORAL
  TREAL var_max, var_min

  CHARACTER*32 name, str, varname(inb_vars)
  CHARACTER*128 line 

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w, e, rho

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BONDARY_INIT' )
#endif

! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
  ENDIF

  prefactor = (gama0-C_1_R)*mach*mach
  r1        = C_1_R
  r05       = C_05_R

! Prepare data for subarrays
  DO iq = 1,inb_flow
     WRITE(varname(iq),*) iq;          varname(iq) = TRIM(ADJUSTL(varname(iq)))
  ENDDO
  DO is = 1,inb_scal
     WRITE(varname(is+inb_flow),*) is; varname(is) = TRIM(ADJUSTL(varname(is)))
  ENDDO

! ###################################################################
! Buffer zone treatment
! ###################################################################
  IF ( buff_type .GT. 0 .OR. bcs_euler_drift .EQ. 1 ) THEN

#ifdef USE_MPI
! I/O routines not yet developed for this particular case
  IF ( ims_npro_i .GT. 1 .AND. ( buff_nps_imin .GT. 0 .OR. buff_nps_imax .GT. 0 ) ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_INIT. I/O routines undeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)     
  ENDIF
#endif

! -------------------------------------------------------------------
! Read buffer zones
! -------------------------------------------------------------------
  IF ( buff_load .EQ. 1 ) THEN
        
     IF ( buff_nps_jmin .GT. 0 ) THEN
        name = 'flow.bcs.jmin'
        ! CALL DNS_READ_FIELDS(name, i0, imax,buff_nps_jmin,kmax, &
        !      inb_flow, i0, isize_field, buffer_hb,                   wrk3d)
        idummy = imax*buff_nps_jmin*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_READ_SUBARRAY8(i4, name, varname,             buffer_hb,                   io_sizes, wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        name = 'scal.bcs.jmin'
        ! CALL DNS_READ_FIELDS(name, i0, imax,buff_nps_jmin,kmax, &
        !      inb_scal, i0, isize_field, buffer_hb(1,1,1,inb_flow+1), wrk3d)
        io_sizes(5) = inb_scal
        CALL IO_READ_SUBARRAY8(i4, name, varname(inb_flow+1), buffer_hb(1,1,1,inb_flow+1), io_sizes, wrk3d)
        ENDIF
     ENDIF

     IF ( buff_nps_jmax .GT. 0 ) THEN 
        name = 'flow.bcs.jmax'
        ! CALL DNS_READ_FIELDS(name, i0, imax,buff_nps_jmax,kmax, &
        !      inb_flow, i0, isize_field, buffer_ht,                   wrk3d)
        idummy = imax*buff_nps_jmax*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_READ_SUBARRAY8(i5, name, varname,             buffer_ht,                   io_sizes, wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        name = 'scal.bcs.jmax'
        ! CALL DNS_READ_FIELDS(name, i0, imax,buff_nps_jmax,kmax, &
        !      inb_scal, i0, isize_field, buffer_ht(1,1,1,inb_flow+1), wrk3d)
        io_sizes(5) = inb_scal
        CALL IO_READ_SUBARRAY8(i5, name, varname(inb_flow+1), buffer_ht(1,1,1,inb_flow+1), io_sizes, wrk3d)
        ENDIF
     ENDIF

     IF ( buff_nps_imin .GT. 0 ) THEN 
        name = 'flow.bcs.imin'
        CALL DNS_READ_FIELDS(name, i0, buff_nps_imin,jmax,kmax, &
             inb_flow, i0, isize_field, buffer_vi,                   wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        name = 'scal.bcs.imin'
        CALL DNS_READ_FIELDS(name, i0, buff_nps_imin,jmax,kmax, &
             inb_scal, i0, isize_field, buffer_vi(1,1,1,inb_flow+1), wrk3d)
        ENDIF
     ENDIF

     IF ( buff_nps_imax .GT. 0 ) THEN 
        name = 'flow.bcs.imax'
        CALL DNS_READ_FIELDS(name, i0, buff_nps_imax,jmax,kmax, &
             inb_flow, i0, isize_field, buffer_vo,                   wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        name = 'scal.bcs.imax'
        CALL DNS_READ_FIELDS(name, i0, buff_nps_imax,jmax,kmax, &
             inb_scal, i0, isize_field, buffer_vo(1,1,1,inb_flow+1), wrk3d)
        ENDIF
     ENDIF

  ELSE
! -------------------------------------------------------------------
! Setting and saving buffer zones
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        DO j = 1,imax*jmax*kmax
           txc(j,2) = e(j) + C_05_R*prefactor*(u(j)*u(j)+v(j)*v(j)+w(j)*w(j))
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     ELSE
        txc(:,1) = C_1_R
     ENDIF

     IF ( buff_nps_jmin .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_HB(q,s, txc, buffer_hb)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.jmin.'//TRIM(ADJUSTL(name))
        ! CALL DNS_WRITE_FIELDS(name, i0, imax,buff_nps_jmin,kmax,&
        !      inb_flow, isize_field, buffer_hb,                   wrk3d)
        idummy = imax*buff_nps_jmin*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_WRITE_SUBARRAY8(i4, name, varname,             buffer_hb,                   io_sizes, wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.jmin.'//TRIM(ADJUSTL(name))
        ! CALL DNS_WRITE_FIELDS(name, i0, imax,buff_nps_jmin,kmax,&
        !      inb_scal, isize_field, buffer_hb(1,1,1,inb_flow+1), wrk3d)
        io_sizes(5) = inb_scal
        CALL IO_WRITE_SUBARRAY8(i4, name, varname(inb_flow+1), buffer_hb(1,1,1,inb_flow+1), io_sizes, wrk3d)
        ENDIF
     ENDIF

     IF ( buff_nps_jmax .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_HT(q,s, txc, buffer_ht)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.jmax.'//TRIM(ADJUSTL(name))
        ! CALL DNS_WRITE_FIELDS(name, i0, imax,buff_nps_jmax,kmax,&
        !      inb_flow, isize_field, buffer_ht,                   wrk3d)
        idummy = imax*buff_nps_jmax*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_WRITE_SUBARRAY8(i5, name, varname,             buffer_ht,                   io_sizes, wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.jmax.'//TRIM(ADJUSTL(name))
        ! CALL DNS_WRITE_FIELDS(name, i0, imax,buff_nps_jmax,kmax,&
        !      inb_scal, isize_field, buffer_ht(1,1,1,inb_flow+1), wrk3d)
        io_sizes(5) = inb_scal
        CALL IO_WRITE_SUBARRAY8(i5, name, varname(inb_flow+1), buffer_ht(1,1,1,inb_flow+1), io_sizes, wrk3d)
        ENDIF
     ENDIF
   
     IF ( buff_nps_imin .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_VI(q,s, txc, buffer_vi)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.imin.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, buff_nps_imin,jmax,kmax,&
             inb_flow, isize_field, buffer_vi,                   wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.imin.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, buff_nps_imin,jmax,kmax,&
             inb_scal, isize_field, buffer_vi(1,1,1,inb_flow+1), wrk3d)
        ENDIF
     ENDIF

     IF ( buff_nps_imax .GT. 0 ) THEN
        CALL BOUNDARY_INIT_VO(q,s, txc, buffer_vo)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.imax.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, buff_nps_imax,jmax,kmax,&
             inb_flow, isize_field, buffer_vo,                   wrk3d)
        IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.imax.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, buff_nps_imax,jmax,kmax,&
             inb_scal, isize_field, buffer_vo(1,1,1,inb_flow+1), wrk3d)
        ENDIF
     ENDIF

  ENDIF
     
! -------------------------------------------------------------------
! Control; so far only horizontal buffer layers...
! -------------------------------------------------------------------
  IF ( buff_nps_jmax .GT. 0 ) THEN
  DO is = 1, inb_flow + inb_scal
     var_max = MAXVAL(buffer_ht(:,:,:,is)); var_min = MINVAL(buffer_ht(:,:,:,is))
#ifdef USE_MPI
     CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     var_max = dummy
     CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
     var_min = dummy
#endif
     WRITE(line,10) var_max
     WRITE(str, 10) var_min
     line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
     WRITE(str,*) is
     line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at upper buffer are '//TRIM(ADJUSTL(line))
     CALL IO_WRITE_ASCII(lfile,line)
  ENDDO
  ENDIF

  IF ( buff_nps_jmin .GT. 0 ) THEN 
  DO is = 1, inb_flow + inb_scal
     var_max = MAXVAL(buffer_hb(:,:,:,is)); var_min = MINVAL(buffer_hb(:,:,:,is))
#ifdef USE_MPI
     CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     var_max = dummy
     CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
     var_min = dummy
#endif
     WRITE(line,10) var_max
     WRITE(str, 10) var_min
     line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
     WRITE(str,*) is
     line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at lower buffer are '//TRIM(ADJUSTL(line))
     CALL IO_WRITE_ASCII(lfile,line)
  ENDDO
  ENDIF

  ENDIF

! ###################################################################
! Compute reference pressure for Poinsot&Lele term in NR Bcs.
! Note that buffer_?(1,1,1,4) is now the energy, and we need p
! ###################################################################
  IF ( bcs_euler_drift .EQ. 1 ) THEN
! -------------------------------------------------------------------
! Using buffer fields; bottom
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        DO j = 1,imax*buff_nps_jmin*kmax
           dummy = C_05_R*prefactor*( buffer_hb(j,1,1,1)*buffer_hb(j,1,1,1) &
                                    + buffer_hb(j,1,1,2)*buffer_hb(j,1,1,2) &
                                    + buffer_hb(j,1,1,3)*buffer_hb(j,1,1,3) )/buffer_hb(j,1,1,5)
           txc(j,1) = ( buffer_hb(j,1,1,4) - dummy )/buffer_hb(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,imax*buff_nps_jmin*kmax
           txc(j,1) = buffer_hb(j,1,1,4)/buffer_hb(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,imax*buff_nps_jmin*kmax
           txc(j,4) = buffer_hb(j,1,1,inb_flow+1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,buff_nps_jmin,kmax, txc(1,4),&
             txc(1,1), buffer_hb(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,buff_nps_jmin,kmax, txc(1,4),&
             buffer_hb(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(imax,buff_nps_jmin,kmax, buffer_hb(1,1,1,inb_flow+1), &
             txc(1,1), buffer_hb(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,buff_nps_jmin,kmax, buffer_hb(1,1,1,inb_flow+1), &
             buffer_hb(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
     bcs_p_jmin = AVG_IK(imax,buff_nps_jmin,kmax, i1, txc(1,1), g(1)%jac,g(3)%jac, area)

! reference plane for BCS at bottom: rho, u_i, p, z_i
     DO k = 1,kmax; DO i = 1,imax
        bcs_hb(i,k,1) = buffer_hb(i,1,k,5)               ! density
        bcs_hb(i,k,2) = buffer_hb(i,1,k,2)/bcs_hb(i,k,1) ! normal velocity, in this case v
        bcs_hb(i,k,3) = buffer_hb(i,1,k,1)/bcs_hb(i,k,1)
        bcs_hb(i,k,4) = buffer_hb(i,1,k,3)/bcs_hb(i,k,1)
        ip = i + (k-1)*imax*buff_nps_jmin
        bcs_hb(i,k,5) = txc(ip,1)                        ! pressure
        DO is = inb_flow+1,inb_vars
           bcs_hb(i,k,is) = buffer_hb(i,1,k,is)/bcs_hb(i,k,1)
        ENDDO
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Using buffer fields; top
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,imax*buff_nps_jmax*kmax
           dummy = C_05_R*prefactor*( buffer_ht(j,1,1,1)*buffer_ht(j,1,1,1) &
                                    + buffer_ht(j,1,1,2)*buffer_ht(j,1,1,2) &
                                    + buffer_ht(j,1,1,3)*buffer_ht(j,1,1,3) )/buffer_ht(j,1,1,5)
           txc(j,1) = ( buffer_ht(j,1,1,4) - dummy )/buffer_ht(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,imax*buff_nps_jmax*kmax
           txc(j,1) = buffer_ht(j,1,1,4)/buffer_ht(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,imax*buff_nps_jmax*kmax
           txc(j,4) = buffer_ht(j,1,1,inb_flow+1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,buff_nps_jmax,kmax, txc(1,4),&
             txc(1,1), buffer_hb(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,buff_nps_jmax,kmax, txc(1,4),&
             buffer_ht(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(imax,buff_nps_jmax,kmax, buffer_ht(1,1,1,inb_flow+1), &
             txc(1,1), buffer_ht(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,buff_nps_jmax,kmax, buffer_ht(1,1,1,inb_flow+1), &
             buffer_ht(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
     bcs_p_jmax = AVG_IK(imax,buff_nps_jmax,kmax, buff_nps_jmax, txc(1,1), g(1)%jac,g(3)%jac, area)

! reference plane for BCS at top: rho, u_i, p, z_i
     DO k = 1,kmax; DO i = 1,imax
        bcs_ht(i,k,1) = buffer_ht(i,buff_nps_jmax,k,5)               ! density
        bcs_ht(i,k,2) = buffer_ht(i,buff_nps_jmax,k,2)/bcs_ht(i,k,1) ! normal velocity, in this case v
        bcs_ht(i,k,3) = buffer_ht(i,buff_nps_jmax,k,1)/bcs_ht(i,k,1)               
        bcs_ht(i,k,4) = buffer_ht(i,buff_nps_jmax,k,3)/bcs_ht(i,k,1)
        ip = i + (buff_nps_jmax-1)*imax + (k-1)*imax*buff_nps_jmax
        bcs_ht(i,k,5) = txc(ip,1)                                    ! pressure
        DO is = inb_flow+1,inb_vars
           bcs_ht(i,k,is) = buffer_ht(i,buff_nps_jmax,k,is)/bcs_ht(i,k,1)
        ENDDO
     ENDDO; ENDDO

! ###################################################################
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN

! -------------------------------------------------------------------
! Using buffer fields; left
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,buff_nps_imin*jmax*kmax
           dummy = C_05_R*prefactor*( buffer_vi(j,1,1,1)*buffer_vi(j,1,1,1) &
                                    + buffer_vi(j,1,1,2)*buffer_vi(j,1,1,2) &
                                    + buffer_vi(j,1,1,3)*buffer_vi(j,1,1,3) )/buffer_vi(j,1,1,5)
           txc(j,1) = ( buffer_vi(j,1,1,4) - dummy )/buffer_vi(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,buff_nps_imin*jmax*kmax
           txc(j,1) = buffer_vi(j,1,1,4)/buffer_vi(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,buff_nps_imin*jmax*kmax
           txc(j,4) = buffer_vi(j,1,1,inb_flow+1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(buff_nps_imin,jmax,kmax, txc(1,4),&
             txc(1,1), buffer_vi(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(buff_nps_imin,jmax,kmax, txc(1,4),&
             buffer_vi(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(buff_nps_imin,jmax,kmax, buffer_vi(1,1,1,inb_flow+1), &
             txc(1,1), buffer_vi(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(buff_nps_imin,jmax,kmax, buffer_vi(1,1,1,inb_flow+1), &
             buffer_vi(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
! reference value at the center of the vertical length
     j = jmax/2
     bcs_p_imin = AVG1V1D(buff_nps_imin,jmax,kmax, i1,j, i1, txc(1,1))

! reference plane for BCS at xmin: rho, u_i, p, z_i
     DO k = 1,kmax; DO j = 1,jmax
        bcs_vi(j,k,1) = buffer_vi(1,j,k,5)                    ! density
        bcs_vi(j,k,2) = buffer_vi(1,j,k,1)/bcs_vi(j,k,1)
        bcs_vi(j,k,3) = buffer_vi(1,j,k,2)/bcs_vi(j,k,1)
        bcs_vi(j,k,4) = buffer_vi(1,j,k,3)/bcs_vi(j,k,1)
        ip = 1 + (j-1)*buff_nps_imin + (k-1)*buff_nps_imin*jmax
        bcs_vi(j,k,5) = txc(ip,1)                             ! pressure
        DO is = inb_flow+1,inb_vars
           bcs_vi(j,k,is) = buffer_vi(1,j,k,is)/bcs_vi(j,k,1)
        ENDDO
     ENDDO; ENDDO

! shape factor
     diam_loc  = C_3_R*diam_u
     thick_loc = diam_u/C_8_R
     ycenter   = g(2)%nodes(1) + g(2)%scale *ycoor_u
     DO k = 1,kmax; DO j = 1,jmax
        bcs_vi(j,k,inb_vars+1) = FLOW_JET_TEMPORAL(i2, thick_loc, r1, r05, diam_loc, ycenter, dummy, g(2)%nodes(j))
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Using buffer fields; right
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,buff_nps_imax*jmax*kmax
           dummy = C_05_R*prefactor*( buffer_vo(j,1,1,1)*buffer_vo(j,1,1,1) &
                                    + buffer_vo(j,1,1,2)*buffer_vo(j,1,1,2) &
                                    + buffer_vo(j,1,1,3)*buffer_vo(j,1,1,3) )/buffer_vo(j,1,1,5)
           txc(j,1) = ( buffer_vo(j,1,1,4) - dummy )/buffer_vo(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,buff_nps_imax*jmax*kmax
           txc(j,1) = buffer_vo(j,1,1,4)/buffer_vo(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,buff_nps_imax*jmax*kmax
           txc(j,4) = buffer_vo(j,1,1,inb_flow+1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(buff_nps_imax,jmax,kmax, txc(1,4),&
             txc(1,1), buffer_vo(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(buff_nps_imax,jmax,kmax, txc(1,4),&
             buffer_vo(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(buff_nps_imax,jmax,kmax, buffer_vo(1,1,1,inb_flow+1), &
             txc(1,1), buffer_vo(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(buff_nps_imax,jmax,kmax, buffer_vo(1,1,1,inb_flow+1), &
             buffer_vo(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
! reference value at the center of the vertical length
     j = jmax/2
     bcs_p_imax = AVG1V1D(buff_nps_imax,jmax,kmax, buff_nps_imax,j, i1, txc(1,1))

! reference plane for BCS at xmax: rho, u_i, p, z_i
     DO k = 1,kmax; DO j = 1,jmax
!        bcs_vo(j,k,1) = buffer_vo(buff_nps_imax,j,k,1)
!        bcs_vo(j,k,2) = buffer_vo(buff_nps_imax,j,k,2)/buffer_vo(buff_nps_imax,j,k,1)
!        bcs_vo(j,k,3) = buffer_vo(buff_nps_imax,j,k,3)/buffer_vo(buff_nps_imax,j,k,1)
!        bcs_vo(j,k,4) = buffer_vo(buff_nps_imax,j,k,4)/buffer_vo(buff_nps_imax,j,k,1)
!        ip = buff_nps_imax + (j-1)*buff_nps_imin + (k-1)*buff_nps_imin*jmax
!        bcs_vo(j,k,5) = txc(ip,1)
!        DO is = 1,inb_scal
!           bcs_vo(j,k,5+is) = buffer_vo(buff_nps_imax,j,k,5+is)/buffer_vo(buff_nps_imax,j,k,1)
!        ENDDO
! try to use only the coflow values
        bcs_vo(j,k,1) = buffer_vo(buff_nps_imax,1,k,5)
        bcs_vo(j,k,2) = buffer_vo(buff_nps_imax,1,k,1)/bcs_vo(j,k,1)
        bcs_vo(j,k,3) = C_0_R
        bcs_vo(j,k,4) = buffer_vo(buff_nps_imax,1,k,3)/bcs_vo(j,k,1)
        ip = buff_nps_imax + (1-1)*buff_nps_imin + (k-1)*buff_nps_imin*jmax
        bcs_vo(j,k,5) = txc(ip,1)
        DO is = inb_flow+1,inb_vars
           bcs_vo(j,k,is) = buffer_vo(buff_nps_imax,1,k,is)/bcs_vo(j,k,1)
        ENDDO
     ENDDO; ENDDO
     
  ELSE
     bcs_p_imin = p_init; bcs_p_imax = p_init

  ENDIF

! rest
  bcs_p_kmin = p_init; bcs_p_kmax = p_init

! -------------------------------------------------------------------
! No buffer field available
! -------------------------------------------------------------------
  ELSE
     bcs_p_imin = p_init; bcs_p_imax = p_init
     bcs_p_jmin = p_init; bcs_p_jmax = p_init
     bcs_p_kmin = p_init; bcs_p_kmax = p_init

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BONDARY_INIT' )
#endif

  RETURN

10 FORMAT(G_FORMAT_R)

END SUBROUTINE BOUNDARY_INIT


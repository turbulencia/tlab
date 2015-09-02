!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2007/10/24 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Initializing inflow fields for broadband forcing case.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"

SUBROUTINE BOUNDARY_INFLOW_INIT(etime, y, x_inf,y_inf,z_inf, q_inf,z1_inf, txc, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : ifrc_mode, ifrc_ifield, frc_length, frc_adapt
  USE DNS_LOCAL, ONLY : imax_inf,jmax_inf,kmax_inf
  USE DNS_LOCAL, ONLY : scalex_inf,scaley_inf,scalez_inf
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_npro_i
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL etime
  TREAL y(jmax)
  TREAL x_inf(imax_inf)
  TREAL y_inf(jmax_inf)
  TREAL z_inf(kmax_total)
  TREAL q_inf(imax_inf*jmax_inf*kmax_inf,*), z1_inf(imax_inf*jmax_inf*kmax_inf,*)
  TREAL txc(imax_inf*jmax_inf*kmax_inf)
  TREAL wrk1d(imax_inf,*)
  TREAL wrk2d(*)
  TREAL wrk3d(*)

  TARGET :: q_inf

! -------------------------------------------------------------------
  TINTEGER ij, is
  TINTEGER joffset, jglobal, j, ibc_local, iwrk_size
  TREAL tolerance, dy
  TREAL visctmp, rtimetmp, itimetmp
  CHARACTER*32 fname, sname, str
  CHARACTER*128 line

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u_inf, v_inf, w_inf, p_inf, rho_inf

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_INIT')
#endif

#ifdef USE_MPI
! I/O routines not yet developed for this particular case
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_INIT. I/O routines undeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)     
  ENDIF
#endif

! Define pointers
  u_inf   => q_inf(:,1)
  v_inf   => q_inf(:,2)
  w_inf   => q_inf(:,3)
  p_inf   => q_inf(:,4)
  rho_inf => q_inf(:,5)

  iwrk_size = imax_inf*jmax_inf*kmax_inf
!  IF ( imax*jmax*kmax .LT. imax_inf*jmax_inf*kmax_inf ) THEN
!     CALL IO_WRITE_ASCII(efile,'BOUNDARY_INFLOW_INIT. Not enough space in array txc.')
!     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
!  ENDIF

! ###################################################################
  IF ( ifrc_mode .EQ. 2 .OR. ifrc_mode .EQ. 3 .OR. ifrc_mode .EQ. 4 ) THEN
     CALL IO_READ_GRID('inflow_grid', imax_inf,jmax_inf,kmax_total,&
          scalex_inf,scaley_inf,scalez_inf, x_inf,y_inf,z_inf)

     fname = 'inflow_flow '
     sname = 'inflow_scal '
     ifrc_ifield = INT(mean_u*etime/scalex_inf) + 1
     IF ( ifrc_mode .EQ. 3 ) THEN
        WRITE(str,*) ifrc_ifield
        fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(str))
        sname = TRIM(ADJUSTL(sname))//TRIM(ADJUSTL(str))
        line='Reading InflowFile'//TRIM(ADJUSTL(str))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDIF

     rtimetmp = rtime
     itimetmp = itime
     visctmp  = visc
     CALL DNS_READ_FIELDS(fname, i2, imax_inf,jmax_inf,kmax_inf, &
          inb_flow, i0, iwrk_size, q_inf, wrk3d)
     rtime = rtimetmp
     itime = itimetmp
     visc  = visctmp

     rtimetmp = rtime
     itimetmp = itime
     CALL DNS_READ_FIELDS(sname, i1, imax_inf,jmax_inf,kmax_inf, &
          inb_scal, i0, iwrk_size, z1_inf, wrk3d)
     rtime = rtimetmp
     itime = itimetmp

! array p contains the internal energy. Now we put in the pressure 
     CALL THERMO_CALORIC_TEMPERATURE&
          (imax_inf, jmax_inf, kmax_inf, z1_inf, p_inf, rho_inf, txc, wrk3d)
     CALL THERMO_THERMAL_PRESSURE&
          (imax_inf, jmax_inf, kmax_inf, z1_inf, rho_inf, txc, p_inf)

! ###################################################################
! Checking the matching 
! ###################################################################
     tolerance = C_1EM10_R
     joffset = (jmax-jmax_inf)/2
     DO j = 1,jmax_inf
        jglobal = joffset + j
        dy = ABS(y(jglobal)-y_inf(j))
        IF (dy.gt.tolerance) THEN
           CALL IO_WRITE_ASCII(efile, 'BOUNDARY_INFLOW. Inflow domain does not match')
           CALL DNS_STOP(DNS_ERROR_INFLOWDOMAIN)
        ENDIF
     ENDDO

! ###################################################################
! Performing the derivatives
! 
! Note that a field f_0(x) is convected with a velocity U along OX, i.e.
! f(x,t) = f_0(x-Ut), and therefore \partial f/\partial t = -U df_0/dx.
! ###################################################################
     IF ( ifrc_mode .EQ. 3 ) THEN
        ibc_local = 1
     ELSE
        ibc_local = 0 ! periodic BCs
     ENDIF

! space in wrk1d needs to be reviewed     
!     CALL FDM_INITIALIZE(i0, imode_fdm, imax_inf, ibc_local, scalex_inf, x_inf, wrk1d(1,1), wrk1d(1,2))

     IF ( icalc_flow .EQ. 1 ) THEN
        DO is = 1, inb_flow
           CALL PARTIAL_X(imode_fdm, imax_inf, jmax_inf, kmax_inf, ibc_local,&
                wrk1d(1,1), q_inf(1,is), txc, i0, i0, wrk1d(1,2), wrk2d, wrk3d)
           DO ij = 1,imax_inf*jmax_inf*kmax_inf
              q_inf(ij,is) = -txc(ij)*mean_u
           ENDDO
        ENDDO
     ENDIF

     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal
           CALL PARTIAL_X(imode_fdm, imax_inf, jmax_inf, kmax_inf, ibc_local,&
                wrk1d(1,1), z1_inf(1,is), txc, i0, i0, wrk1d(1,2), wrk2d, wrk3d)
           DO ij = 1,imax_inf*jmax_inf*kmax_inf
              z1_inf(ij,is) = -txc(ij)*mean_u
           ENDDO

        ENDDO
     ENDIF

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_INIT')
#endif

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_INIT

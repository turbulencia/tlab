#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!# DESCRIPTION
!#
!# Initializing inflow fields for broadband forcing case.
!#
!########################################################################
SUBROUTINE BOUNDARY_INFLOW_INIT(etime, q_inf,s_inf, txc, wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow, inb_scal, icalc_flow,icalc_scal
  USE DNS_GLOBAL,    ONLY : g, qbg
  USE DNS_GLOBAL,    ONLY : rtime,itime,visc
  USE DNS_LOCAL,     ONLY : ifrc_mode, ifrc_ifield
  USE DNS_LOCAL,     ONLY : g_inf
#ifdef USE_MPI
  USE DNS_LOCAL,     ONLY : FilterInflow
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL etime
  TREAL, DIMENSION(g_inf(1)%size*&
                   g_inf(2)%size*&
                   g_inf(3)%size,inb_flow), INTENT(INOUT) :: q_inf
  TREAL, DIMENSION(g_inf(1)%size*&
                   g_inf(2)%size*&
                   g_inf(3)%size,inb_scal), INTENT(INOUT) :: s_inf
  TREAL, DIMENSION(g_inf(1)%size*&
                   g_inf(2)%size*&
                   g_inf(3)%size),          INTENT(INOUT) :: txc, wrk3d
  TREAL, DIMENSION(*),                      INTENT(INOUT) :: wrk2d

  TARGET :: q_inf

! -------------------------------------------------------------------
  TINTEGER is, itimetmp, bcs(2,1)
  TINTEGER joffset, jglobal, j, iwrk_size
  TREAL tolerance, dy
  TREAL visctmp, rtimetmp
  CHARACTER*32 fname, sname, str
  CHARACTER*128 line

#ifdef USE_MPI
  TINTEGER isize_loc,id
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: p_inf, rho_inf

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

! #######################################################################
! Definining types for parallel mode
! #######################################################################
#ifdef USE_MPI
  IF ( FilterInflow(1)%type .NE. DNS_FILTER_NONE ) THEN !  Required for inflow explicit filter
     CALL IO_WRITE_ASCII(lfile,'Initialize MPI types for inflow filter.')
     id    = DNS_MPI_K_INFLOW
     isize_loc = FilterInflow(1)%size *FilterInflow(2)%size
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

! Define pointers
  p_inf   => q_inf(:,4)
  rho_inf => q_inf(:,5)

  iwrk_size = g_inf(1)%size *g_inf(2)%size *kmax
  IF ( imax*jmax*kmax .LT. iwrk_size ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_INFLOW_INIT. Not enough space in array txc.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! ###################################################################
  IF ( ifrc_mode .EQ. 2 .OR. ifrc_mode .EQ. 3 .OR. ifrc_mode .EQ. 4 ) THEN

! Checking the matching; we could move this outside... 
     tolerance = C_1EM10_R
     joffset = ( jmax - g_inf(2)%size )/2
     DO j = 1,g_inf(2)%size
        jglobal = joffset + j
        dy = ABS( g(2)%nodes(jglobal) -g_inf(2)%nodes(j) )
        IF (dy.gt.tolerance) THEN
           CALL IO_WRITE_ASCII(efile, 'BOUNDARY_INFLOW. Inflow domain does not match.')
           CALL DNS_STOP(DNS_ERROR_INFLOWDOMAIN)
        ENDIF
     ENDDO

! Reading fields
     fname = 'flow.inf'
     sname = 'scal.inf'
     ifrc_ifield = INT( qbg(1)%mean *etime /g_inf(1)%scale ) + 1
     IF ( ifrc_mode .EQ. 3 ) THEN
        WRITE(str,*) ifrc_ifield
        fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(str))
        sname = TRIM(ADJUSTL(sname))//TRIM(ADJUSTL(str))
        line='Reading InflowFile '//TRIM(ADJUSTL(str))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDIF

     rtimetmp = rtime
     itimetmp = itime
     visctmp  = visc
     CALL DNS_READ_FIELDS(fname, i2, g_inf(1)%size,g_inf(2)%size,kmax, inb_flow, i0, iwrk_size, q_inf, wrk3d)
     CALL DNS_READ_FIELDS(sname, i1, g_inf(1)%size,g_inf(2)%size,kmax, inb_scal, i0, iwrk_size, s_inf, wrk3d)
     rtime = rtimetmp
     itime = itimetmp
     visc  = visctmp

! array p contains the internal energy. Now we put in the pressure 
     CALL THERMO_CALORIC_TEMPERATURE&
          (g_inf(1)%size, g_inf(2)%size, kmax, s_inf, p_inf, rho_inf, txc, wrk3d)
     CALL THERMO_THERMAL_PRESSURE&
          (g_inf(1)%size, g_inf(2)%size, kmax, s_inf, rho_inf, txc, p_inf)

! ###################################################################
! Performing the derivatives
! 
! Note that a field f_0(x) is convected with a velocity U along OX, i.e.
! f(x,t) = f_0(x-Ut), and therefore \partial f/\partial t = -U df_0/dx.
! ###################################################################
     bcs = 0
     
     IF ( icalc_flow .EQ. 1 ) THEN
        DO is = 1, inb_flow
           CALL OPR_PARTIAL_X(OPR_P1, g_inf(1)%size,g_inf(2)%size,kmax, bcs, g_inf(1), q_inf(1,is), txc, wrk3d, wrk2d, wrk3d)
           q_inf(:,is) = -txc(:) *qbg(1)%mean
        ENDDO
     ENDIF
     
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal
           CALL OPR_PARTIAL_X(OPR_P1, g_inf(1)%size,g_inf(2)%size,kmax, bcs, g_inf(1), s_inf(1,is), txc, wrk3d, wrk2d, wrk3d)
           s_inf(:,is) = -txc(:) *qbg(1)%mean
        ENDDO
     ENDIF

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_INIT')
#endif

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_INIT

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Calculating RHS forcings at the inflow plane in spatially evolving cases
!#
!########################################################################
MODULE BOUNDARY_INFLOW

  USE TLAB_TYPES,     ONLY : filter_dt, grid_dt, discrete_dt
  USE TLAB_CONSTANTS, ONLY : efile, lfile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, inb_flow, inb_scal, inb_flow_array,inb_scal_array, icalc_flow,icalc_scal
  USE TLAB_VARS,    ONLY : imode_eqns, itransport
  USE TLAB_VARS,    ONLY : g, qbg, epbackground, pbackground
  USE TLAB_VARS,    ONLY : rtime,itime
  USE TLAB_VARS,    ONLY : visc,damkohler
  USE TLAB_PROCS
  USE THERMO_VARS, ONLY : imixture
  USE IO_FIELDS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_k
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE
  SAVE
  PRIVATE

  TYPE(grid_dt),      PUBLIC :: g_inf(3)
  TREAL, ALLOCATABLE         :: x_inf(:,:), y_inf(:,:), z_inf(:,:)
  TREAL, ALLOCATABLE         :: q_inf(:,:,:,:), s_inf(:,:,:,:)

  TINTEGER,           PUBLIC :: inflow_mode, inflow_ifield
  TREAL,              PUBLIC :: inflow_adapt
  TYPE(filter_dt),    PUBLIC :: FilterInflow(3)
  ! TINTEGER :: FilterInflowStep
  TYPE(discrete_dt),  PUBLIC :: fp ! Discrete forcing

  TARGET :: x_inf, y_inf, z_inf

  PUBLIC :: BOUNDARY_INFLOW_INITIALIZE
  PUBLIC :: BOUNDARY_INFLOW_BROADBAND
  PUBLIC :: BOUNDARY_INFLOW_DISCRETE
  PUBLIC :: BOUNDARY_INFLOW_FILTER

CONTAINS
  !########################################################################
  !########################################################################
  !# Initializing inflow fields for broadband forcing case.
  SUBROUTINE BOUNDARY_INFLOW_INITIALIZE(etime, txc, wrk1d,wrk2d,wrk3d)
    IMPLICIT NONE

#include "integers.h"

    TREAL etime
    TREAL, INTENT(INOUT) :: txc(g_inf(1)%size,g_inf(2)%size,g_inf(3)%size)
    TREAL, INTENT(INOUT) :: wrk1d(*),wrk2d(*),wrk3d(*)

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

    ! ###################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_INIT')
#endif

#ifdef USE_MPI
    ! I/O routines not yet developed for this particular case
    IF ( ims_npro_i .GT. 1 ) THEN
      CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_INIT. I/O routines undeveloped.')
      CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
    ENDIF
#endif

    IF ( .NOT. ALLOCATED(x_inf) ) ALLOCATE( x_inf(g_inf(1)%size,g_inf(1)%inb_grid) )
    IF ( .NOT. ALLOCATED(y_inf) ) ALLOCATE( y_inf(g_inf(2)%size,g_inf(2)%inb_grid) )
    IF ( .NOT. ALLOCATED(z_inf) ) ALLOCATE( z_inf(g_inf(3)%size,g_inf(3)%inb_grid) )
    IF ( .NOT. ALLOCATED(q_inf) ) ALLOCATE( q_inf(g_inf(1)%size,g_inf(2)%size,g_inf(3)%size,inb_flow_array) )
    IF ( .NOT. ALLOCATED(s_inf) ) ALLOCATE( s_inf(g_inf(1)%size,g_inf(2)%size,g_inf(3)%size,inb_scal_array) )

    IF ( g_inf(1)%size > 1 ) THEN ! Inflow fields for spatial simulations
      IF ( .NOT. ASSOCIATED(g_inf(1)%nodes) ) &
          CALL IO_READ_GRID('grid.inf', g_inf(1)%size, g_inf(2)%size, g_inf(3)%size,  &
          g_inf(1)%scale,g_inf(2)%scale,g_inf(3)%scale, x_inf,y_inf,z_inf)
      CALL FDM_INITIALIZE(x_inf, g_inf(1), wrk1d)
      IF ( .NOT. ASSOCIATED(g_inf(2)%nodes) ) g_inf(2)%nodes => y_inf(:,1)
      IF ( .NOT. ASSOCIATED(g_inf(3)%nodes) ) g_inf(3)%nodes => z_inf(:,1)
    ENDIF

    ! #######################################################################
    ! Definining types for parallel mode
    ! #######################################################################
#ifdef USE_MPI
    IF ( FilterInflow(1)%TYPE .NE. DNS_FILTER_NONE ) THEN !  Required for inflow explicit filter
      CALL TLAB_WRITE_ASCII(lfile,'Initialize MPI types for inflow filter.')
      id    = TLAB_MPI_K_INFLOW
      isize_loc = FilterInflow(1)%size *FilterInflow(2)%size
      CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
      FilterInflow(3)%mpitype = id
    ENDIF
#endif

    iwrk_size = g_inf(1)%size *g_inf(2)%size *kmax
    IF ( imax*jmax*kmax .LT. iwrk_size ) THEN
      CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_INFLOW_INIT. Not enough space in array txc.')
      CALL TLAB_STOP(DNS_ERROR_WRKSIZE)
    ENDIF

    ! ###################################################################
    IF ( inflow_mode .EQ. 2 .OR. inflow_mode .EQ. 3 .OR. inflow_mode .EQ. 4 ) THEN

      ! Checking the matching; we could move this outside...
      tolerance = C_1EM10_R
      joffset = ( jmax - g_inf(2)%size )/2
      DO j = 1,g_inf(2)%size
        jglobal = joffset + j
        dy = ABS( g(2)%nodes(jglobal) -g_inf(2)%nodes(j) )
        IF (dy.GT.tolerance) THEN
          CALL TLAB_WRITE_ASCII(efile, 'BOUNDARY_INFLOW. Inflow domain does not match.')
          CALL TLAB_STOP(DNS_ERROR_INFLOWDOMAIN)
        ENDIF
      ENDDO

      ! Reading fields
      fname = 'flow.inf'
      sname = 'scal.inf'
      inflow_ifield = INT( qbg(1)%mean *etime /g_inf(1)%scale ) + 1
      IF ( inflow_mode .EQ. 3 ) THEN
        WRITE(str,*) inflow_ifield
        fname = TRIM(ADJUSTL(fname))//TRIM(ADJUSTL(str))
        sname = TRIM(ADJUSTL(sname))//TRIM(ADJUSTL(str))
        line='Reading InflowFile '//TRIM(ADJUSTL(str))
        CALL TLAB_WRITE_ASCII(lfile,line)
      ENDIF

      rtimetmp = rtime
      itimetmp = itime
      visctmp  = visc
      CALL IO_READ_FIELDS(fname, IO_FLOW, g_inf(1)%size,g_inf(2)%size,kmax, inb_flow, i0, q_inf, wrk3d)
      CALL IO_READ_FIELDS(sname, IO_SCAL, g_inf(1)%size,g_inf(2)%size,kmax, inb_scal, i0, s_inf, wrk3d)
      rtime = rtimetmp
      itime = itimetmp
      visc  = visctmp

      ! array p contains the internal energy. Now we put in the pressure
      CALL THERMO_CALORIC_TEMPERATURE&
          (g_inf(1)%size, g_inf(2)%size, kmax, s_inf, q_inf(1,1,1,4), q_inf(1,1,1,5), txc, wrk3d)
      CALL THERMO_THERMAL_PRESSURE&
          (g_inf(1)%size, g_inf(2)%size, kmax, s_inf, q_inf(1,1,1,5), txc, q_inf(1,1,1,4))

      ! ###################################################################
      ! Performing the derivatives
      !
      ! Note that a field f_0(x) is convected with a velocity U along OX, i.e.
      ! f(x,t) = f_0(x-Ut), and therefore \partial f/\partial t = -U df_0/dx.
      ! ###################################################################
      bcs = 0

      IF ( icalc_flow .EQ. 1 ) THEN
        DO is = 1, inb_flow
          CALL OPR_PARTIAL_X(OPR_P1, g_inf(1)%size,g_inf(2)%size,kmax, bcs, g_inf(1), q_inf(1,1,1,is), txc, wrk3d, wrk2d, wrk3d)
          q_inf(:,:,:,is) = -txc(:,:,:) *qbg(1)%mean
        ENDDO
      ENDIF

      IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal
          CALL OPR_PARTIAL_X(OPR_P1, g_inf(1)%size,g_inf(2)%size,kmax, bcs, g_inf(1), s_inf(1,1,1,is), txc, wrk3d, wrk2d, wrk3d)
          s_inf(:,:,:,is) = -txc(:,:,:) *qbg(1)%mean
        ENDDO
      ENDIF

    ENDIF

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_INIT')
#endif

    RETURN
  END SUBROUTINE BOUNDARY_INFLOW_INITIALIZE

  !########################################################################
  !########################################################################
  SUBROUTINE BOUNDARY_INFLOW_BROADBAND(etime, inf_rhs, txc, wrk1d,wrk2d,wrk3d)
    IMPLICIT NONE

    TREAL etime
    TREAL, INTENT(  OUT) :: inf_rhs(jmax,kmax,inb_flow+inb_scal)
    TREAL, INTENT(INOUT) :: wrk1d(*), wrk2d(*), wrk3d(*), txc(*)

    ! -------------------------------------------------------------------
    TREAL xaux, dx_loc, vmult
    TINTEGER joffset, jglobal, ileft, iright, j, k, is, ip
    TREAL BSPLINES3P, BSPLINES3

    ! ###################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_BROADBAND' )
#endif

    ! Transient factor
    IF ( inflow_adapt .GT. C_0_R .AND. etime .LE. inflow_adapt ) THEN
      vmult = etime / inflow_adapt
    ELSE
      vmult = C_1_R
    ENDIF

    ! check if we need to read again inflow data
    IF ( inflow_mode .EQ. 3 .AND. INT(qbg(1)%mean*etime/g_inf(1)%scale)+1 .NE. inflow_ifield ) THEN
      CALL BOUNDARY_INFLOW_INITIALIZE(etime, txc, wrk1d,wrk2d,wrk3d)
    ENDIF

    ! ###################################################################
    ! Getting the position
    ! ###################################################################
    joffset = (jmax-g_inf(2)%size)/2

    xaux = qbg(1)%mean*etime
    ! Remove integral length scales of box
    xaux = xaux - INT(xaux/g_inf(1)%scale)*g_inf(1)%scale
    ! Set distance from box initial length
    xaux = g_inf(1)%scale-xaux

    dx_loc = g_inf(1)%nodes(2) - g_inf(1)%nodes(1)
    ! Get left index
    ileft = INT(xaux/dx_loc) +1
    ! Check bounds
    IF ( ileft .GT. g_inf(1)%size ) THEN
      ileft = 1
    ENDIF
    ! Set right index
    iright = ileft + 1
    ! Check bounds
    IF ( iright .GT. g_inf(1)%size ) THEN
      iright = 1
    ENDIF
    ! Get relative distance from left point
    xaux = (xaux-(g_inf(1)%nodes(ileft)-g_inf(1)%nodes(1)))/dx_loc

    ! ###################################################################
    ! Sampling the information
    ! ###################################################################
    ! -------------------------------------------------------------------
    ! Periodic
    ! -------------------------------------------------------------------
    IF ( inflow_mode .EQ. 2 ) THEN
      DO k = 1,kmax
        DO j = 1,g_inf(2)%size
          jglobal = joffset + j
          DO is = 1,inb_scal
            inf_rhs(jglobal,k,is) = inf_rhs(jglobal,k,is) + vmult *BSPLINES3P(q_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
          ENDDO

          IF ( icalc_scal .EQ. 1 ) THEN
            DO is = 1,inb_scal
              ip = inb_flow +is
              inf_rhs(jglobal,k,ip) = inf_rhs(jglobal,k,ip) + vmult *BSPLINES3P(s_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
            ENDDO
          ENDIF

        ENDDO
      ENDDO

      ! -------------------------------------------------------------------
      ! Sequential
      ! -------------------------------------------------------------------
    ELSE
      DO k = 1,kmax
        DO j = 1,g_inf(2)%size
          jglobal = joffset + j
          DO is = 1,inb_flow
            inf_rhs(jglobal,k,is) = inf_rhs(jglobal,k,is) + vmult *BSPLINES3(q_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
          ENDDO

          IF ( icalc_scal .EQ. 1 ) THEN
            DO is = 1,inb_scal
              ip = inb_flow +is
              inf_rhs(jglobal,k,ip) = inf_rhs(jglobal,k,ip) + vmult *BSPLINES3(s_inf(1,j,k,is), g_inf(1)%size, ileft, xaux)
            ENDDO
          ENDIF

        ENDDO
      ENDDO

    ENDIF

    ! ###################################################################
    ! Filling the rest
    ! ###################################################################
    DO j = 1,joffset
      inf_rhs(j,:,:) = inf_rhs(j,:,:) + C_0_R
    ENDDO
    DO j = jmax-joffset+1,jmax
      inf_rhs(j,:,:) = inf_rhs(j,:,:) + C_0_R
    ENDDO

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_BROADBAND' )
#endif
    RETURN
  END SUBROUTINE BOUNDARY_INFLOW_BROADBAND

  !########################################################################
  !########################################################################
  SUBROUTINE BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs, wrk1d, wrk2d)
    IMPLICIT NONE

    TREAL etime
    TREAL, INTENT(  OUT) :: inf_rhs(jmax,kmax,inb_flow+inb_scal)
    TREAL, INTENT(INOUT) :: wrk1d(jmax,2)
    TREAL, INTENT(INOUT) :: wrk2d(kmax,3)

    ! -------------------------------------------------------------------
    TINTEGER j, k, im, kdsp
    TREAL wx, wz, wx_1, wz_1, xaux, vmult, factorx, factorz, dummy

    TREAL PROFILES, ycenter, yr
    EXTERNAL PROFILES

    TREAL, DIMENSION(:), POINTER :: y,z

    ! ###################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_DISCRETE' )
#endif

    ! Define pointers
    y => g(2)%nodes
    z => g(3)%nodes

#ifdef USE_MPI
    kdsp = ims_offset_k
#else
    kdsp = 0
#endif

    xaux =-qbg(1)%mean *etime

    ! ###################################################################
    ! Shape function
    ! ###################################################################
    SELECT CASE ( fp%TYPE )
    CASE ( PROFILE_GAUSSIAN )
      ycenter = y(1) +g(2)%scale *qbg(1)%ymean
      DO j = 1,jmax
        yr = y(j)-ycenter
        wrk1d(j,1) = PROFILES( PROFILE_GAUSSIAN, fp%parameters(1), C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,2) = yr /( fp%parameters(1) **2 ) *wrk1d(j,1) ! Derivative of f
      ENDDO

    CASE (PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM)
      ycenter = y(1) +g(2)%scale *qbg(1)%ymean -C_05_R *qbg(1)%diam
      DO j = 1,jmax
        yr = y(j) - ycenter
        wrk1d(j,1) = PROFILES( PROFILE_GAUSSIAN, fp%parameters(1), C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,2) =-yr /( fp%parameters(1) **2 ) *wrk1d(j,1)
      ENDDO

      ycenter = y(1) +g(2)%scale *qbg(1)%ymean +C_05_R *qbg(1)%diam
      IF     ( fp%TYPE .EQ. PROFILE_GAUSSIAN_ANTISYM ) THEN; factorx =-C_1_R ! varicose
      ELSEIF ( fp%TYPE .EQ. PROFILE_GAUSSIAN_SYM     ) THEN; factorx = C_1_R ! Sinuous
      ENDIF
      DO j = 1,jmax
        yr = y(j) - ycenter
        dummy = factorx *PROFILES( PROFILE_GAUSSIAN, fp%parameters(1), C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,1) = wrk1d(j,1) +dummy
        wrk1d(j,2) = wrk1d(j,2) +yr /( fp%parameters(1) **2 ) *dummy
      ENDDO

    END SELECT

    ! ###################################################################
    ! Fourier series
    ! ###################################################################
    wx_1 = C_2_R * C_PI_R / fp%parameters(2) ! Fundamental wavelengths
    wz_1 = C_2_R * C_PI_R / g(3)%scale

    wrk2d = C_0_R
    DO im = 1,fp%size
      wx = M_REAL( fp%modex(im) ) *wx_1
      wz = M_REAL( fp%modez(im) ) *wz_1

      ! Factor to impose solenoidal constraint
      IF     ( fp%modex(im) .EQ. 0 .AND. fp%modez(im) .EQ. 0 ) THEN; EXIT
      ELSEIF (                           fp%modez(im) .EQ. 0 ) THEN; factorx= C_1_R /wx; factorz= C_0_R
      ELSEIF ( fp%modex(im) .EQ. 0                           ) THEN; factorx= C_0_R;     factorz= C_1_R /wz
      ELSE;                                                          factorx= C_05_R/wx; factorz= C_05_R/wz
      ENDIF

      DO k = 1,kmax
        wrk2d(k,2) = wrk2d(k,2) + fp%amplitude(im) *COS( wx *xaux +fp%phasex(im) ) *COS( wz *z(kdsp+k) +fp%phasez(im) )
        wrk2d(k,1) = wrk2d(k,1) + fp%amplitude(im) *SIN( wx *xaux +fp%phasex(im) ) *COS( wz *z(kdsp+k) +fp%phasez(im) ) *factorx
        wrk2d(k,3) = wrk2d(k,3) + fp%amplitude(im) *COS( wx *xaux +fp%phasex(im) ) *SIN( wz *z(kdsp+k) +fp%phasez(im) ) *factorz
      ENDDO

    ENDDO

    ! ###################################################################
    ! Forcing
    ! ###################################################################
    ! Transient factor
    IF ( inflow_adapt .GT. C_0_R .AND. etime .LE. inflow_adapt ) THEN
      vmult = etime / inflow_adapt
    ELSE
      vmult = C_1_R
    ENDIF

    DO k = 1,kmax
      DO j = 1,jmax
        inf_rhs(j,k,2) = inf_rhs(j,k,2) - vmult *qbg(1)%mean *wrk2d(k,1) *wrk1d(j,2) ! u
        inf_rhs(j,k,3) = inf_rhs(j,k,3) - vmult *qbg(1)%mean *wrk2d(k,2) *wrk1d(j,1) ! v
        inf_rhs(j,k,4) = inf_rhs(j,k,4) - vmult *qbg(1)%mean *wrk2d(k,3) *wrk1d(j,2) ! w
      ENDDO
    ENDDO

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_DISCRETE' )
#endif

    RETURN
  END SUBROUTINE BOUNDARY_INFLOW_DISCRETE

  !########################################################################
  !########################################################################
  ! Filter
  ! This should be integrated into the inflow buffer, as the filter contribution
  ! BufferFilter should then be a block in dns.ini as [Filter], which is read in io_read_global.

  SUBROUTINE BOUNDARY_INFLOW_FILTER(bcs_vi, bcs_vi_scal, q,s, txc, wrk1d,wrk2d,wrk3d)
    IMPLICIT NONE

    TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(INOUT) :: q,s
    TREAL, DIMENSION(jmax,kmax,*),      INTENT(IN)    :: bcs_vi, bcs_vi_scal
    TREAL, DIMENSION(imax*jmax*kmax,2), INTENT(INOUT) :: txc
    TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

    TARGET q

    ! -----------------------------------------------------------------------
    TINTEGER i,j,k,ip, iq, iq_loc(inb_flow), is
    TINTEGER j1, imx, jmx, ifltmx, jfltmx

    ! Pointers to existing allocated space
    TREAL, DIMENSION(:,:,:), POINTER :: e, rho, p, T, vis

    ! ###################################################################
    ! #######################################################################
    CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
    ! FilterInflow needs to be initiliazed
    CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

    ! Define pointers
    IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
      e   => q(:,:,:,4)
      rho => q(:,:,:,5)
      p   => q(:,:,:,6)
      T   => q(:,:,:,7)

      IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)

    ENDIF

    ! Define counters
    imx = FilterInflow(1)%size
    j1  = ( jmax -FilterInflow(2)%size )/2 +1
    jmx = ( jmax +FilterInflow(2)%size )/2
    j1  = MIN(MAX(j1,1),jmax)
    jmx = MIN(MAX(jmx,1),jmax)

    ifltmx = imx-1+1
    jfltmx = jmx-j1+1

    IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
      iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
    ELSE
      iq_loc = (/ 1,2,3 /)
    ENDIF

    ! #######################################################################
    DO iq = 1,inb_flow

      ! -----------------------------------------------------------------------
      ! Remove mean field
      ! -----------------------------------------------------------------------
      ip = 1
      DO k = 1,kmax
        DO j = j1,jmx
          DO i = 1,imx
            wrk3d(ip) = q(i,j,k,iq_loc(iq)) - bcs_vi(j,k,iq_loc(iq))
            ip = ip + 1
          ENDDO
        ENDDO
      ENDDO

      ! -----------------------------------------------------------------------
      CALL OPR_FILTER(ifltmx,jfltmx,kmax, FilterInflow, wrk3d, wrk1d,wrk2d,txc)

      ! -----------------------------------------------------------------------
      ! Add mean field
      ! -----------------------------------------------------------------------
      ip = 1
      DO k = 1,kmax
        DO j = j1,jmx
          DO i = 1,imx
            q(i,j,k,iq_loc(iq))=  wrk3d(ip) + bcs_vi(j,k,iq_loc(iq))
            ip = ip + 1
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    ! #######################################################################
    DO is = 1,inb_scal

      ! -----------------------------------------------------------------------
      ! Remove mean field
      ! -----------------------------------------------------------------------
      ip = 1
      DO k = 1,kmax
        DO j = j1,jmx
          DO i = 1,imx
            wrk3d(ip) = s(i,j,k,is) - bcs_vi_scal(j,k,is)
            ip = ip + 1
          ENDDO
        ENDDO
      ENDDO

      ! -----------------------------------------------------------------------
      CALL OPR_FILTER(ifltmx,jfltmx,kmax, FilterInflow, wrk3d, wrk1d,wrk2d,txc)

      ! -----------------------------------------------------------------------
      ! Add mean field
      ! -----------------------------------------------------------------------
      ip = 1
      DO k = 1,kmax
        DO j = j1,jmx
          DO i = 1,imx
            s(i,j,k,is) = wrk3d(ip) + bcs_vi_scal(j,k,is)
            ip = ip + 1
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    ! #######################################################################
    ! recalculation of diagnostic variables
    IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
      IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN
        CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,1,1,2), s(1,1,1,1), epbackground,pbackground)

      ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,1,1,inb_scal_array))

      ENDIF

    ELSE
      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_AIRWATER_RP(imax,jmax,kmax, s, p, rho, T, wrk3d)
      ELSE
        CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, T)
      ENDIF
      CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, T, e)

      ! This recalculation of T and p is made to make sure that the same numbers are
      ! obtained in statistics postprocessing as in the simulation; avg* files
      ! can then be compared with diff command.
      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
      ENDIF

      IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)

    ENDIF

    RETURN
  END SUBROUTINE BOUNDARY_INFLOW_FILTER

END MODULE BOUNDARY_INFLOW

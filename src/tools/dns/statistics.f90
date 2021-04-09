#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif
#ifdef LES
#include "les_const.h"
#endif

MODULE STATISTICS

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  IMPLICIT NONE
  SAVE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:)     :: mean
  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:,:)   :: mean_flow
  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:,:,:) :: mean_scal

  LOGICAL, PUBLIC :: stats_averages, stats_pdfs, stats_intermittency, stats_filter

  PUBLIC :: STATISTICS_INITIALIZE
  PUBLIC :: STATISTICS_TEMPORAL_LAYER
  PUBLIC :: STATISTICS_SPATIAL_LAYER

  PRIVATE

CONTAINS

! ###################################################################
! ###################################################################
SUBROUTINE STATISTICS_INITIALIZE()

  USE DNS_GLOBAL, ONLY : imode_sim, jmax, inb_scal, nstatavg

  IF      ( imode_sim .EQ. DNS_MODE_TEMPORAL) THEN
     ALLOCATE(mean(jmax,MAX_AVG_TEMPORAL))

  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     ALLOCATE(mean_flow(nstatavg,jmax,MA_MOMENTUM_SIZE))
     ALLOCATE(mean_scal(nstatavg,jmax,MS_SCALAR_SIZE,inb_scal))

  ENDIF

  RETURN
END SUBROUTINE STATISTICS_INITIALIZE

!########################################################################
!# Statistics that depend on time and the crosswise coordinate:
!# temporally evolving configurations (mixing layer, jet, rti...)
!# Isotropic case should be put in a different subroutine.
!########################################################################
SUBROUTINE STATISTICS_TEMPORAL_LAYER(q,s,hq, txc, wrk1d,wrk2d,wrk3d)

#ifdef TRACE_ON
  USE DNS_CONSTANTS, ONLY : tfile
#endif
  USE DNS_TYPES, ONLY : pointers_dt
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field, isize_txc_field, inb_scal_array
  USE DNS_GLOBAL,    ONLY : buoyancy, imode_eqns, icalc_scal
  USE DNS_GLOBAL,    ONLY : itransport, froude
  USE DNS_GLOBAL,    ONLY : epbackground, pbackground, rbackground
  USE DNS_GLOBAL,    ONLY : itime, rtime
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL, ONLY : iles, iles_type_chem
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field,    *), INTENT(IN)    :: q,s
  TREAL, DIMENSION(isize_field,    *), INTENT(INOUT) :: hq ! auxiliary array
  TREAL, DIMENSION(isize_txc_field,6), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                 INTENT(INOUT) :: wrk1d,wrk2d,wrk3d
!  INTEGER(1), DIMENSION(*)                :: wrk3d                ! see comments in header
!# wrk3d    In     Aux array, with space for REAL(8) 3D field; declared
!#                 as INTEGER(1) to be used as gate array.
!#                 PROBLEMS IN JUQUEEN WITH THIS TRICK!!! Maybe use hs?

  TARGET :: q,s, txc

! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER nbins, is, flag_buoyancy
  TINTEGER ibc(16), nfield
  TREAL amin(16), amax(16)
  CHARACTER*32 fname
  CHARACTER*64 str
  INTEGER(1) igate !, gate_levels(16)

  TYPE(pointers_dt), DIMENSION(16) :: data

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u,v,w,e,rho, p,T,vis

! #ifdef USE_MPI
!   TINTEGER id
! #endif
#ifdef LES
  TINTEGER ichi
#endif

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING STATS_TEMPORAL_LAYER' )
#endif

  flag_buoyancy = 0  ! default

! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)
  ELSE
     e   => txc(:,6) ! not used, but argument in LES routines
     rho => txc(:,6) ! not used
     p   => txc(:,3) ! to be used in AVG_FLOW_TEMPORAL_LAYER
! in case we need the buoyancy statistics
     IF ( buoyancy%type .EQ. EQNS_BOD_QUADRATIC   .OR. &
          buoyancy%type .EQ. EQNS_BOD_BILINEAR    .OR. &
          imixture .EQ. MIXT_TYPE_AIRWATER        .OR. &
          imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        flag_buoyancy = 1
     ENDIF
  ENDIF

! Calculate pressure
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,3), txc(1,1),txc(1,2), hq, wrk1d,wrk2d,wrk3d)

  ELSE
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)

  ENDIF

! ###################################################################
! Intermittency
! ###################################################################
!   IF ( stats_intermittency ) THEN
!      CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)

! ! calculate vorticity gate based on 1% threshold
!      CALL MINMAX(imax,jmax,kmax, txc(1,1), umin,umax)
!      umin = C_1EM3_R*C_1EM3_R*umax
!      DO ij = 1,isize_field
!         IF ( txc(ij,1) .GT. umin ) THEN; wrk3d(ij) = 1  ! gate array
!         ELSE;                            wrk3d(ij) = 0; ENDIF
!      ENDDO
!      varname(1) = 'Vorticity'; gate_levels(1) = 1

!      WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
!      CALL INTER_N_XZ(fname, itime,rtime, imax,jmax,kmax, i1, varname, wrk3d, y, mean)

!   ENDIF

! ###################################################################
! Unconditional plane PDFs
! ###################################################################
  IF ( stats_pdfs ) THEN
     nfield = 0
     nfield = nfield+1; data(nfield)%field => u(:); data(nfield)%tag = 'u'
     nfield = nfield+1; data(nfield)%field => v(:); data(nfield)%tag = 'v'
     nfield = nfield+1; data(nfield)%field => w(:); data(nfield)%tag = 'w'
     nfield = nfield+1; data(nfield)%field => p(:); data(nfield)%tag = 'p'
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        nfield = nfield+1; data(nfield)%field => rho(:); data(nfield)%tag = 'r'
        nfield = nfield+1; data(nfield)%field => T(:);   data(nfield)%tag = 't'
     ENDIF

     DO is = 1,inb_scal_array
       nfield = nfield+1; data(nfield)%field => s(:,is); data(nfield)%tag = 's'
       WRITE(str,*) is; data(nfield)%tag=TRIM(ADJUSTL(data(nfield)%tag))//TRIM(ADJUSTL(str))
     ENDDO

     ibc(1:nfield) = 2 ! BCs in the calculation of the PDFs
     igate = 0         ! no intermittency partition

     nbins = 32
     WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
     CALL PDF1V_N(fname, rtime, imax,jmax,kmax, &
      nfield, nbins, ibc, amin,amax,data, igate,wrk3d, g(2)%nodes, txc, wrk1d)

  ENDIF

! ###################################################################
! Plane averages
! ###################################################################
  IF ( stats_averages ) THEN
     CALL AVG_FLOW_XZ(q,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),hq(1,1),hq(1,2),hq(1,3),  &
          mean, wrk1d,wrk2d,wrk3d)

     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
           CALL AVG_SCAL_XZ(is, q,s, s(1,is), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
        ENDDO

! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
        IF ( flag_buoyancy .EQ. 1 ) THEN
           dummy = C_1_R/froude
           IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
              CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, hq(1,1))
           ELSE
              wrk1d(1:jmax) = C_0_R
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, hq(1,1), wrk1d) ! note that wrk3d is defined as integer.
           ENDIF
           hq(1:isize_field,1) = hq(1:isize_field,1) *dummy

           CALL AVG_SCAL_XZ(is, q,s, hq(1,1), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)

        ENDIF

        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              is = is + 1
              CALL THERMO_ANELASTIC_THETA_L(imax,jmax,kmax, s, epbackground,pbackground, hq(1,1))
              CALL AVG_SCAL_XZ(is, q,s, hq(1,1), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
           ENDIF
        ENDIF

     ENDIF

#ifdef LES
! -------------------------------------------------------------------
! LES
! -------------------------------------------------------------------
     IF ( iles .EQ. 1 ) THEN
        CALL LES_AVG_FLOW_TEMPORAL_LAYER&
             (rho,u,v,w,e, hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
             txc, vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d,&
             vaux(vindex(VA_LES_FLT0X)),vaux(vindex(VA_LES_FLT0Y)),vaux(vindex(VA_LES_FLT0Z)), &
             vaux(vindex(VA_LES_FLT1X)),vaux(vindex(VA_LES_FLT1Y)),vaux(vindex(VA_LES_FLT1Z)), &
             vaux(vindex(VA_LES_FLT2X)),vaux(vindex(VA_LES_FLT2Y)),vaux(vindex(VA_LES_FLT2Z)), &
             vaux(vindex(VA_LES_ARM_UF)), vaux(vindex(VA_LES_ARM_PF)),&
             vaux(vindex(VA_LES_ARM_UH)), vaux(vindex(VA_LES_ARM_PH)))
        IF ( iles_type_chem .EQ. LES_CHEM_QUASIBS ) THEN
           ichi = inb_txc-1; txc(:,ichi) = C_0_R
        ELSE
! this memory address is not used
           ichi = 1
        ENDIF
        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal
              CALL LES_AVG_SCAL_TEMPORAL_LAYER&
                   (is, rho,u,v,w, s(1,is), hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
                   txc, vaux(vindex(VA_MEAN_WRK)), txc(1,ichi), wrk1d,wrk2d,wrk3d, &
                   vaux(vindex(VA_LES_FLT0X)),vaux(vindex(VA_LES_FLT0Y)),vaux(vindex(VA_LES_FLT0Z)), &
                   vaux(vindex(VA_LES_FLT1X)),vaux(vindex(VA_LES_FLT1Y)),vaux(vindex(VA_LES_FLT1Z)), &
                   vaux(vindex(VA_LES_FLT2X)),vaux(vindex(VA_LES_FLT2Y)),vaux(vindex(VA_LES_FLT2Z)), &
                   vaux(vindex(VA_LES_ARM_UF)),vaux(vindex(VA_LES_ARM_ZF)),&
                   vaux(vindex(VA_LES_ARM_UH)),vaux(vindex(VA_LES_ARM_ZH)))
           ENDDO
        ENDIF
!        IF ( iles_type_chem .EQ. LES_CHEM_QUASIBS ) THEN
!           CALL LES_AVG_CHEM_BS(x, y, z, dx, dy, dz, rho, &
!                u, s, gama, txc, vaux(vindex(VA_MEAN_WRK)), &
!                txc(1,ichi), wrk1d, wrk2d, wrk3d,&
!                vaux(vindex(VA_LES_FLT0X)), vaux(vindex(VA_LES_FLT0Y)), &
!                vaux(vindex(VA_LES_FLT0Z)), vaux(vindex(VA_LES_FDF_BS)))
!        ENDIF
     ENDIF
#endif

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING STATS_TEMPORAL_LAYER' )
#endif

  RETURN
END SUBROUTINE STATISTICS_TEMPORAL_LAYER

! ###################################################################
! ###################################################################
SUBROUTINE STATISTICS_SPATIAL_LAYER(txc, wrk1d,wrk2d)

#ifdef TRACE_ON
  USE DNS_CONSTANTS, ONLY : tfile
#endif
  USE DNS_GLOBAL
  USE DNS_LOCAL
  USE BOUNDARY_BUFFER
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(*) :: txc, wrk1d, wrk2d

! -----------------------------------------------------------------------
  TINTEGER is, buff_u_jmin, buff_u_jmax, isize_txc

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING STATS_SPATIAL_LAYER' )
#endif

! #######################################################################
! Averages
! #######################################################################
  IF ( stats_averages ) THEN
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif
        isize_txc = inb_txc*isize_txc_field

        buff_u_jmin = BuffFlowJmax%size
        buff_u_jmax = jmax -BuffFlowJmax%size +1
        CALL AVG_FLOW_SPATIAL_LAYER(isize_txc, buff_u_jmin,buff_u_jmax, &
             mean_flow, txc, wrk1d,wrk2d)

        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal
              CALL AVG_SCAL_SPATIAL_LAYER(is, isize_txc, buff_u_jmin,buff_u_jmax, &
                   mean_flow, mean_scal(1,1,1,is), txc, wrk1d)
           ENDDO
        ENDIF

#ifdef LES
        IF ( iles .EQ. 1 ) THEN
           CALL LES_AVG_SPATIAL_LAYER(isize_txc, x,y, vaux(vindex(VA_MEAN_WRK)), txc, wrk1d,wrk2d)
        ENDIF
#endif

#ifdef USE_MPI
     ENDIF
#endif
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING STATS_SPATIAL_LAYER' )
#endif

  RETURN
END SUBROUTINE STATISTICS_SPATIAL_LAYER

END MODULE STATISTICS

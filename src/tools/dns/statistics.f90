#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"

MODULE STATISTICS

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  IMPLICIT NONE
  SAVE

  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:)     :: mean
  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:,:)   :: mean_flow
  TREAL, ALLOCATABLE, PUBLIC, DIMENSION(:,:,:,:) :: mean_scal

  PRIVATE

  LOGICAL, PUBLIC :: stats_averages, stats_pdfs, stats_intermittency, stats_buoyancy

  PUBLIC :: STATISTICS_INITIALIZE, STATISTICS_TEMPORAL, STATISTICS_SPATIAL

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE STATISTICS_INITIALIZE()

    USE TLAB_VARS, ONLY : imode_sim, jmax, inb_scal, nstatavg

    IF      ( imode_sim .EQ. DNS_MODE_TEMPORAL) THEN
      ALLOCATE(mean(jmax,MAX_AVG_TEMPORAL))

    ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
      ALLOCATE(mean_flow(nstatavg,jmax,MA_MOMENTUM_SIZE))
      ALLOCATE(mean_scal(nstatavg,jmax,MS_SCALAR_SIZE,inb_scal))

    END IF

    RETURN
  END SUBROUTINE STATISTICS_INITIALIZE

  !########################################################################
  !########################################################################
  SUBROUTINE STATISTICS_TEMPORAL()

#ifdef TRACE_ON
    USE DNS_CONSTANTS, ONLY : tfile
#endif
    USE TLAB_TYPES,     ONLY : pointers_dt
    USE TLAB_VARS,    ONLY : g
    USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field, inb_scal_array
    USE TLAB_VARS,    ONLY : buoyancy, imode_eqns, icalc_scal
    USE TLAB_VARS,    ONLY : froude
    USE TLAB_VARS,    ONLY : epbackground, pbackground, rbackground
    USE TLAB_VARS,    ONLY : itime, rtime
    USE TLAB_VARS,      ONLY : inb_part, icalc_part
    USE TLAB_VARS,      ONLY : sbg, schmidt
    USE TLAB_ARRAYS
    USE DNS_ARRAYS
    USE THERMO_GLOBAL, ONLY : imixture
    USE LAGRANGE_GLOBAL, ONLY : l_g, ilagrange
    USE LAGRANGE_GLOBAL, ONLY : icalc_part_pdf
    USE LAGRANGE_ARRAYS

    IMPLICIT NONE

    ! -------------------------------------------------------------------
    TREAL dummy, amin(16), amax(16)
    TINTEGER ij, is, idummy, nbins, ibc(16), nfield
    TYPE(pointers_dt) vars(16)
    CHARACTER*32 fname, gatename(1)
    CHARACTER*64 str
    INTEGER(1) igate
    INTEGER(1), ALLOCATABLE, SAVE :: gate(:)

    ! ###################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING STATS_TEMPORAL_LAYER' )
#endif

    stats_buoyancy = .FALSE.  ! default

    ! in case we need the buoyancy statistics
    IF ( buoyancy%TYPE .EQ. EQNS_BOD_QUADRATIC   .OR. &
        buoyancy%TYPE .EQ. EQNS_BOD_BILINEAR    .OR. &
        imixture .EQ. MIXT_TYPE_AIRWATER        .OR. &
        imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
      stats_buoyancy = .TRUE.
    END IF

    ! Calculate pressure
    IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
      CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,3), txc(1,1),txc(1,2), txc(1,4), wrk1d,wrk2d,wrk3d)
    END IF

    ! ###################################################################
    ! Intermittency
    ! ###################################################################
    IF ( stats_intermittency ) THEN
      ALLOCATE(gate(isize_field))
      CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,4), wrk2d,wrk3d)

      ! calculate vorticity gate based on 1% threshold
      CALL MINMAX(imax,jmax,kmax, txc(1,1), amin(1),amax(1))
      amin(1) = C_1EM3_R*C_1EM3_R*amax(1)
      DO ij = 1,isize_field
        IF ( txc(ij,1) .GT. amin(1) ) THEN; gate(ij) = 1  ! gate array
        ELSE;                               gate(ij) = 0
        END IF
      END DO
      nfield = 1; gatename(1) = 'Vorticity'

      WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
      CALL INTER_N_XZ(fname, itime,rtime, imax,jmax,kmax, nfield, gatename, gate, g(2)%nodes, mean)

      DEALLOCATE(gate)
    END IF

    ! ###################################################################
    ! Unconditional plane PDFs
    ! ###################################################################
    IF ( stats_pdfs ) THEN
      nfield = 0
      nfield = nfield+1; vars(nfield)%field => q(:,1); vars(nfield)%tag = 'u'
      nfield = nfield+1; vars(nfield)%field => q(:,2); vars(nfield)%tag = 'v'
      nfield = nfield+1; vars(nfield)%field => q(:,3); vars(nfield)%tag = 'w'
      IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        nfield = nfield+1; vars(nfield)%field => txc(:,3); vars(nfield)%tag = 'p'
      ELSE
        nfield = nfield+1; vars(nfield)%field => q(:,6); vars(nfield)%tag = 'p'
        nfield = nfield+1; vars(nfield)%field => q(:,5); vars(nfield)%tag = 'r'
        nfield = nfield+1; vars(nfield)%field => q(:,7); vars(nfield)%tag = 't'
      END IF

      DO is = 1,inb_scal_array
        nfield = nfield+1; vars(nfield)%field => s(:,is); vars(nfield)%tag = 's'
        WRITE(str,*) is; vars(nfield)%tag=TRIM(ADJUSTL(vars(nfield)%tag))//TRIM(ADJUSTL(str))
      END DO

      ibc(1:nfield) = 2 ! BCs in the calculation of the PDFs
      igate = 0         ! no intermittency partition

      nbins = 32
      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
      CALL PDF1V_N(fname, rtime, imax,jmax,kmax, &
          nfield, nbins, ibc, amin,amax,vars, igate,wrk3d, g(2)%nodes, txc, wrk1d)

    END IF

    ! ###################################################################
    ! Plane averages
    ! ###################################################################
    IF ( stats_averages ) THEN
      IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
          hq(1:isize_field,3) = txc(1:isize_field,3) ! Pass the pressure
          CALL AVG_SCAL_XZ(is, q,s, s(1,is), &
              txc(1,1),txc(1,2),txc(1,4),txc(1,5),txc(1,6),hq(1,3), mean, wrk1d,wrk2d,wrk3d)
        END DO

        ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
        IF ( stats_buoyancy ) THEN
          IF ( buoyancy%TYPE .EQ. EQNS_EXPLICIT ) THEN
            CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, hq(1,1))
          ELSE
            wrk1d(1:jmax,1) = C_0_R
            CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, hq(1,1), wrk1d) ! note that wrk3d is defined as integer.
          END IF
          dummy = C_1_R/froude
          hq(1:isize_field,1) = hq(1:isize_field,1) *dummy

          hq(1:isize_field,3) = txc(1:isize_field,3) ! Pass the pressure
          CALL AVG_SCAL_XZ(is, q,s, hq(1,1), &
              txc(1,1),txc(1,2),txc(1,4),txc(1,5),txc(1,6),hq(1,3), mean, wrk1d,wrk2d,wrk3d)

        END IF

        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
          IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
            is = is + 1
            CALL THERMO_ANELASTIC_THETA_L(imax,jmax,kmax, s, epbackground,pbackground, hq(1,1))

            hq(1:isize_field,3) = txc(1:isize_field,3) ! Pass the pressure
            CALL AVG_SCAL_XZ(is, q,s, hq(1,1), &
                txc(1,1),txc(1,2),txc(1,4),txc(1,5),txc(1,6),hq(1,3), mean, wrk1d,wrk2d,wrk3d)
          END IF
        END IF

      END IF

      CALL AVG_FLOW_XZ(q,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),hq(1,1),hq(1,2),hq(1,3),  &
          mean, wrk1d,wrk2d,wrk3d)

      ! Lagrange Liquid and Liquid without diffusion
      IF ( icalc_part == 1 ) THEN
        IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
          l_txc(:,1) = C_1_R; ! We want density
          CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,5), wrk2d,wrk3d)

          hq(:,1) = hq(:,1) + C_SMALL_R
          idummy = inb_part - 3 ! # scalar properties solved in the lagrangian
          DO is = inb_scal_array +1 +1, inb_scal_array+1 +idummy
            sbg(is)%mean = C_1_R; sbg(is)%delta = C_0_R; sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)
            CALL PARTICLE_TO_FIELD(l_q, l_q(1,3+is-inb_scal_array-1), hq(1,2), wrk2d,wrk3d)
            hq(:,2) = hq(:,2) /hq(:,1)
            CALL AVG_SCAL_XZ(is, q,s, hq(1,2), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4), txc(1,5),txc(1,6),mean, wrk1d,wrk2d,wrk3d)
          END DO
        END IF

        IF ( icalc_part_pdf .EQ. 1) THEN                ! Save particle pathlines for particle_pdf
          WRITE(fname,*) itime; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
          CALL PARTICLE_PDF(fname,s, l_g,l_q,l_txc,l_comm, wrk3d)
        END IF

        IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN  ! Save particle residence times
          WRITE(fname,*) itime; fname = "residence_pdf."//TRIM(ADJUSTL(fname))
          CALL PARTICLE_RESIDENCE_PDF(fname, l_g%np, l_q)
        END IF

      END IF

    END IF

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING STATS_TEMPORAL_LAYER' )
#endif

    RETURN
  END SUBROUTINE STATISTICS_TEMPORAL

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE STATISTICS_SPATIAL()

#ifdef TRACE_ON
    USE DNS_CONSTANTS, ONLY : tfile
#endif
    USE TLAB_VARS
    USE TLAB_ARRAYS
    USE DNS_LOCAL
    USE BOUNDARY_BUFFER
#ifdef USE_MPI
    USE TLAB_MPI_VARS
#endif

    IMPLICIT NONE

    ! -----------------------------------------------------------------------
    TINTEGER is, buff_u_jmin, buff_u_jmax, isize_txc

    ! #######################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING STATS_SPATIAL_LAYER' )
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
          END DO
        END IF

#ifdef LES
        IF ( iles .EQ. 1 ) THEN
          CALL LES_AVG_SPATIAL_LAYER(isize_txc, x,y, vaux(vindex(VA_MEAN_WRK)), txc, wrk1d,wrk2d)
        END IF
#endif

#ifdef USE_MPI
      END IF
#endif
    END IF

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING STATS_SPATIAL_LAYER' )
#endif

    RETURN
  END SUBROUTINE STATISTICS_SPATIAL

END MODULE STATISTICS

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

MODULE BOUNDARY_BCS

  USE TLAB_CONSTANTS, ONLY : MAX_VARS
  USE TLAB_PROCS

  IMPLICIT NONE
  SAVE

  TYPE bcs_dt
     SEQUENCE
     TINTEGER type(MAX_VARS)                      ! dirichlet, neumann for incompressible
     TINTEGER SfcType(MAX_VARS)                   ! Type of Surface Model
     TREAL cpl(MAX_VARS)                          ! Coupling parameter for surface model
     TREAL cinf, cout, ctan                       ! characteristic formulation for compressible
     TREAL, ALLOCATABLE, DIMENSION(:,:,:) :: ref  ! reference fields
  END type bcs_dt

  TYPE(bcs_dt), PUBLIC :: BcsFlowImin,BcsFlowImax,BcsFlowJmin,BcsFlowJmax,BcsFlowKmin,BcsFlowKmax
  TYPE(bcs_dt), PUBLIC :: BcsScalImin,BcsScalImax,BcsScalJmin,BcsScalJmax,BcsScalKmin,BcsScalKmax

  LOGICAL, PUBLIC :: BcsDrift

! Compressible viscous
  TINTEGER, PUBLIC :: bcs_inf(2,2,3), bcs_out(2,2,3) ! 1. index: lower and upper values
                                             ! 2. index: derivative order
                                             ! 3. index: direction
  PUBLIC :: BOUNDARY_BCS_INITIALIZE

  PRIVATE

CONTAINS

! ###################################################################
! ###################################################################
SUBROUTINE BOUNDARY_BCS_INITIALIZE(wrk3d)

  USE TLAB_CONSTANTS, ONLY : tag_flow,tag_scal, lfile, efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_VARS,    ONLY : imode_eqns
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, inb_flow_array,inb_scal_array
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : mach, pbg, qbg
  USE THERMO_VARS, ONLY : gama0
  USE BOUNDARY_BUFFER
#ifdef USE_MPI
  USE MPI
  USE TLAB_VARS,    ONLY : inb_scal_array
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_VARS, ONLY : ims_bcs_imax, ims_bcs_jmax
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER j, is
  TREAL prefactor, param(5)
  TREAL diam_loc, thick_loc, ycenter, r1, r05, PROFILES

#ifdef USE_MPI
  CHARACTER*32 str
  TINTEGER isize_loc,id
#endif

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_INITIALLIZE' )
#endif

! ###################################################################
! Allocate memory space
! ###################################################################
! we use inb_flow_array and inb_scal_array to have space for aux fields
! we add space at the end for the shape factor
  ALLOCATE( BcsFlowImin%ref(jmax,kmax,inb_flow_array+1) )
  ALLOCATE( BcsFlowImax%ref(jmax,kmax,inb_flow_array+1) )
  ALLOCATE( BcsFlowJmin%ref(imax,kmax,inb_flow_array+1) )
  ALLOCATE( BcsFlowJmax%ref(imax,kmax,inb_flow_array+1) )
  ALLOCATE( BcsFlowKmin%ref(imax,jmax,inb_flow_array+1) ) ! not yet used
  ALLOCATE( BcsFlowKmax%ref(imax,jmax,inb_flow_array+1) )

  ALLOCATE( BcsScalImin%ref(jmax,kmax,inb_scal_array+1) )
  ALLOCATE( BcsScalImax%ref(jmax,kmax,inb_scal_array+1) )
  ALLOCATE( BcsScalJmin%ref(imax,kmax,inb_scal_array+1) )
  ALLOCATE( BcsScalJmax%ref(imax,kmax,inb_scal_array+1) )
  ALLOCATE( BcsScalKmin%ref(imax,jmax,inb_scal_array+1) ) ! not yet used
  ALLOCATE( BcsScalKmax%ref(imax,jmax,inb_scal_array+1) )

! #######################################################################
! Incompressible mode
! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
! So far, nothing to initialize

  ELSE
! #######################################################################
! Compressible mode
! #######################################################################
#ifdef USE_MPI
! -------------------------------------------------------------------
! Characteristic BCs
! -------------------------------------------------------------------
     IF ( .NOT. g(1)%periodic ) THEN ! Required for NRBCs in Ox
        id    = TLAB_MPI_K_NRBCX
        isize_loc = MOD(jmax,ims_npro_k)
        ims_bcs_imax = 2*(inb_flow+inb_scal_array)
        DO WHILE ( MOD(isize_loc*ims_bcs_imax,ims_npro_k) .GT. 0 )
           ims_bcs_imax = ims_bcs_imax + 1
        ENDDO
        WRITE(str,*) ims_bcs_imax
        str = 'Initialize MPI types for Ox BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
        CALL TLAB_WRITE_ASCII(lfile,str)
        isize_loc = ims_bcs_imax*jmax
        CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
             ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF

     IF ( .NOT. g(2)%periodic ) THEN ! Required for NRBCs in Oy
        id    = TLAB_MPI_K_NRBCY
        isize_loc = MOD(imax,ims_npro_k)
        ims_bcs_jmax = 2*(inb_flow+inb_scal_array)
        DO WHILE ( MOD(isize_loc*ims_bcs_jmax,ims_npro_k) .GT. 0 )
           ims_bcs_jmax = ims_bcs_jmax + 1
        ENDDO
        WRITE(str,*) ims_bcs_jmax
        str = 'Initialize MPI types for Oy BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
        CALL TLAB_WRITE_ASCII(lfile,str)
        isize_loc = imax*ims_bcs_jmax
        CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
             ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF
#endif

! ###################################################################
! Compute reference pressure for Poinsot&Lele term in NR Bcs.
! Note that buffer_?(1,1,1,4) is now the energy, and we need p
! ###################################################################
     BcsFlowImin%ref(:,:,5) = pbg%mean; BcsFlowImax%ref(:,:,5) = pbg%mean ! default values
     BcsFlowJmin%ref(:,:,5) = pbg%mean; BcsFlowJmax%ref(:,:,5) = pbg%mean
     BcsFlowKmin%ref(:,:,5) = pbg%mean; BcsFlowKmax%ref(:,:,5) = pbg%mean

     prefactor = C_05_R *(gama0-C_1_R) *mach*mach
     r1        = C_1_R
     r05       = C_05_R

! -------------------------------------------------------------------
! Using buffer fields; bottom
! -------------------------------------------------------------------
     IF ( BuffFlowJmin%size .GT. 0 ) THEN
        BcsFlowJmin%ref(:,:,1) = BuffFlowJmin%Ref(:,1,:,5)                          ! density
        BcsFlowJmin%ref(:,:,2) = BuffFlowJmin%Ref(:,1,:,2) / BcsFlowJmin%ref(:,:,1) ! normal velocity, in this case v
        BcsFlowJmin%ref(:,:,3) = BuffFlowJmin%Ref(:,1,:,1) / BcsFlowJmin%ref(:,:,1)
        BcsFlowJmin%ref(:,:,4) = BuffFlowJmin%Ref(:,1,:,3) / BcsFlowJmin%ref(:,:,1)
        BcsFlowJmin%ref(:,:,6) = BuffFlowJmin%Ref(:,1,:,4) / BcsFlowJmin%ref(:,:,1) ! energy, to get pressure into ref5 below
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           BcsFlowJmin%ref(:,:,6) = BcsFlowJmin%ref(:,:,6) &
                                  - prefactor *( BcsFlowJmin%ref(:,:,2) *BcsFlowJmin%ref(:,:,2) &
                                               + BcsFlowJmin%ref(:,:,3) *BcsFlowJmin%ref(:,:,3) &
                                               + BcsFlowJmin%ref(:,:,4) *BcsFlowJmin%ref(:,:,4) )
        ENDIF
        DO is = 1,inb_scal
           BcsScalJmin%ref(:,:,is) = BuffScalJmin%Ref(:,1,:,is) / BcsFlowJmin%ref(:,:,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,i1,kmax, BcsScalJmin%Ref, &
             BcsFlowJmin%ref(1,1,6), BcsFlowJmin%Ref(1,1,1), BcsFlowJmin%Ref(1,1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,i1,kmax, BcsScalJmin%Ref, &
             BcsFlowJmin%Ref(1,1,1), BcsFlowJmin%Ref(1,1,7), BcsFlowJmin%Ref(1,1,5))

! shape factor
        BcsFlowJmin%ref(:,:,inb_flow+1) = C_1_R
        BcsScalJmin%ref(:,:,inb_scal+1) = C_1_R

     ENDIF

! -------------------------------------------------------------------
! Using buffer fields; top
! -------------------------------------------------------------------
     IF ( BuffFlowJmax%size .GT. 0 ) THEN
        j = BuffFlowJmax%size
        BcsFlowJmax%ref(:,:,1) = BuffFlowJmax%Ref(:,j,:,5)                          ! density
        BcsFlowJmax%ref(:,:,2) = BuffFlowJmax%Ref(:,j,:,2) / BcsFlowJmax%ref(:,:,1) ! normal velocity, in this case v
        BcsFlowJmax%ref(:,:,3) = BuffFlowJmax%Ref(:,j,:,1) / BcsFlowJmax%ref(:,:,1)
        BcsFlowJmax%ref(:,:,4) = BuffFlowJmax%Ref(:,j,:,3) / BcsFlowJmax%ref(:,:,1)
        BcsFlowJmax%ref(:,:,6) = BuffFlowJmax%Ref(:,j,:,4) / BcsFlowJmax%ref(:,:,1) ! energy, to get pressure into ref5 below
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           BcsFlowJmax%ref(:,:,6) = BcsFlowJmax%ref(:,:,6) &
                                  - prefactor *( BcsFlowJmax%ref(:,:,2) *BcsFlowJmax%ref(:,:,2) &
                                               + BcsFlowJmax%ref(:,:,3) *BcsFlowJmax%ref(:,:,3) &
                                               + BcsFlowJmax%ref(:,:,4) *BcsFlowJmax%ref(:,:,4) )
        ENDIF
        DO is = 1,inb_scal
           BcsScalJmax%ref(:,:,is) = BuffScalJmax%Ref(:,j,:,is) / BcsFlowJmax%ref(:,:,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,i1,kmax, BcsScalJmax%Ref, &
             BcsFlowJmax%ref(1,1,6), BcsFlowJmax%Ref(1,1,1), BcsFlowJmax%Ref(1,1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,i1,kmax, BcsScalJmax%Ref, &
             BcsFlowJmax%Ref(1,1,1), BcsFlowJmax%Ref(1,1,7), BcsFlowJmax%Ref(1,1,5))

! shape factor
        BcsFlowJmax%ref(:,:,inb_flow+1) = C_1_R
        BcsScalJmax%ref(:,:,inb_scal+1) = C_1_R
     ENDIF

! -------------------------------------------------------------------
! Using buffer fields; left
! -------------------------------------------------------------------
     IF ( BuffFlowImin%size .GT. 0 ) THEN
        BcsFlowImin%ref(:,:,1) = BuffFlowImin%Ref(:,:,1,5)                          ! density
        BcsFlowImin%ref(:,:,2) = BuffFlowImin%Ref(:,:,1,1) / BcsFlowImin%ref(:,:,1) ! normal velocity, in this case u
        BcsFlowImin%ref(:,:,3) = BuffFlowImin%Ref(:,:,1,2) / BcsFlowImin%ref(:,:,1)
        BcsFlowImin%ref(:,:,4) = BuffFlowImin%Ref(:,:,1,3) / BcsFlowImin%ref(:,:,1)
        BcsFlowImin%ref(:,:,6) = BuffFlowImin%Ref(:,:,1,4) / BcsFlowImin%ref(:,:,1) ! energy, to get pressure into ref5 below
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           BcsFlowImin%ref(:,:,6) = BcsFlowImin%ref(:,:,6) &
                                  - prefactor *( BcsFlowImin%ref(:,:,2) *BcsFlowImin%ref(:,:,2) &
                                               + BcsFlowImin%ref(:,:,3) *BcsFlowImin%ref(:,:,3) &
                                               + BcsFlowImin%ref(:,:,4) *BcsFlowImin%ref(:,:,4) )
        ENDIF
        DO is = 1,inb_scal
           BcsScalImin%ref(:,:,is) = BuffScalImin%Ref(:,:,1,is) / BcsFlowImin%ref(:,:,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,i1,kmax, BcsScalImin%Ref, &
             BcsFlowImin%ref(1,1,6), BcsFlowImin%Ref(1,1,1), BcsFlowImin%Ref(1,1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,i1,kmax, BcsScalImin%Ref, &
             BcsFlowImin%Ref(1,1,1), BcsFlowImin%Ref(1,1,7), BcsFlowImin%Ref(1,1,5))

! shape factor
        diam_loc  = C_3_R*qbg(1)%diam
        thick_loc = qbg(1)%diam/C_8_R
        ycenter   = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean
        param=C_0_R; param(5) = diam_loc
        BcsFlowImin%ref(:,:,inb_flow+1) = PROFILES(i2, thick_loc, r1, r05, ycenter, param, g(2)%nodes(j))
        BcsScalImin%ref(:,:,inb_scal+1) = PROFILES(i2, thick_loc, r1, r05, ycenter, param, g(2)%nodes(j))

     ENDIF

! -------------------------------------------------------------------
! Using buffer fields; right
! -------------------------------------------------------------------
     IF ( BuffFlowImax%size .GT. 0 ) THEN
        BcsFlowImax%ref(:,:,1) = BuffFlowImax%Ref(:,:,1,5)                          ! density
        BcsFlowImax%ref(:,:,2) = BuffFlowImax%Ref(:,:,1,1) / BcsFlowImax%ref(:,:,1) ! normal velocity, in this case u
        BcsFlowImax%ref(:,:,3) = BuffFlowImax%Ref(:,:,1,2) / BcsFlowImax%ref(:,:,1)
        BcsFlowImax%ref(:,:,4) = BuffFlowImax%Ref(:,:,1,3) / BcsFlowImax%ref(:,:,1)
        BcsFlowImax%ref(:,:,6) = BuffFlowImax%Ref(:,:,1,4) / BcsFlowImax%ref(:,:,1) ! energy, to get pressure into ref5 below
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           BcsFlowImax%ref(:,:,6) = BcsFlowImax%ref(:,:,6) &
                                  - prefactor *( BcsFlowImax%ref(:,:,2) *BcsFlowImax%ref(:,:,2) &
                                               + BcsFlowImax%ref(:,:,3) *BcsFlowImax%ref(:,:,3) &
                                               + BcsFlowImax%ref(:,:,4) *BcsFlowImax%ref(:,:,4) )
        ENDIF
        DO is = 1,inb_scal
           BcsScalImax%ref(:,:,is) = BuffScalImax%Ref(:,:,1,is) / BcsFlowImax%ref(:,:,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,i1,kmax, BcsScalImax%Ref, &
             BcsFlowImax%ref(1,1,6), BcsFlowImax%Ref(1,1,1), BcsFlowImax%Ref(1,1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,i1,kmax, BcsScalImax%Ref, &
             BcsFlowImax%Ref(1,1,1), BcsFlowImax%Ref(1,1,7), BcsFlowImax%Ref(1,1,5))

! try to use only the coflow values
        BcsFlowImax%ref(:,:,3) = C_0_R

! shape factor
        BcsFlowImax%ref(:,:,inb_flow+1) = C_1_R
        BcsScalImax%ref(:,:,inb_scal+1) = C_1_R

     ENDIF

  ENDIF

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_INITIALLIZE' )
#endif

  RETURN

END SUBROUTINE BOUNDARY_BCS_INITIALIZE

END MODULE BOUNDARY_BCS

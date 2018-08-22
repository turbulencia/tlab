#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

MODULE BOUNDARY_BCS

  IMPLICIT NONE
  SAVE
  
  TYPE bcs_dt
     SEQUENCE
     TINTEGER type                                ! dirichlet, neumann, nonreflective, inflow, outflow
     TREAL, ALLOCATABLE, DIMENSION(:,:,:) :: ref  ! reference fields
  END type bcs_dt
  
  TYPE(bcs_dt), PUBLIC :: BcsFlowImin,BcsFlowImax,BcsFlowJmin,BcsFlowJmax
  TYPE(bcs_dt), PUBLIC :: BcsScalImin,BcsScalImax,BcsScalJmin,BcsScalJmax

CONTAINS
  
! ###################################################################
! ###################################################################
SUBROUTINE BOUNDARY_BCS_INITIALIZE(bcs_ht,bcs_hb,bcs_vi,bcs_vo, q,s, txc, wrk3d)

  USE DNS_CONSTANTS, ONLY : tag_flow,tag_scal, lfile, efile
  USE DNS_GLOBAL,    ONLY : imode_eqns, imode_sim
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, inb_flow_array,inb_scal_array
  USE DNS_GLOBAL,    ONLY : g,area
  USE DNS_GLOBAL,    ONLY : mach, pbg, qbg
  USE THERMO_GLOBAL, ONLY : imixture, gama0
  USE DNS_LOCAL,     ONLY : bcs_euler_drift, bcs_p_imin,bcs_p_imax, bcs_p_jmin,bcs_p_jmax, bcs_p_kmin,bcs_p_kmax
  USE BOUNDARY_BUFFER
#ifdef USE_MPI
  USE DNS_GLOBAL,    ONLY : inb_scal_array
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
  
  TREAL, DIMENSION(imax,kmax,*)      :: bcs_hb, bcs_ht
  TREAL, DIMENSION(jmax,kmax,*)      :: bcs_vi, bcs_vo
  TREAL, DIMENSION(imax*jmax*kmax,*) :: q, s, txc
  TREAL, DIMENSION(*)                :: wrk3d

! -------------------------------------------------------------------
  TINTEGER i,j,k, is, ip
  TREAL AVG1V1D, AVG_IK, prefactor, dummy
  TREAL diam_loc, thick_loc, ycenter, r1, r05, FLOW_JET_TEMPORAL

#ifdef USE_MPI
  CHARACTER*32 str
  TINTEGER isize_loc,id
#endif
  
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_INITIALLIZE' )
#endif

! ###################################################################
! Allocate memory space
! ###################################################################
! we use inb_flow_array and inb_scal_array to have space for aux fields
  ALLOCATE( BcsFlowImin%ref(jmax,kmax,inb_flow_array) )
  ALLOCATE( BcsFlowImax%ref(jmax,kmax,inb_flow_array) )
  ALLOCATE( BcsFlowJmin%ref(imax,kmax,inb_flow_array) )
  ALLOCATE( BcsFlowJmax%ref(imax,kmax,inb_flow_array) )
  
  ALLOCATE( BcsScalImin%ref(jmax,kmax,inb_scal_array) )
  ALLOCATE( BcsScalImax%ref(jmax,kmax,inb_scal_array) )
  ALLOCATE( BcsScalJmin%ref(imax,kmax,inb_scal_array) )
  ALLOCATE( BcsScalJmax%ref(imax,kmax,inb_scal_array) )

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
        id    = DNS_MPI_K_NRBCX
        isize_loc = MOD(jmax,ims_npro_k)
        ims_bcs_imax = 2*(inb_flow+inb_scal_array)
        DO WHILE ( MOD(isize_loc*ims_bcs_imax,ims_npro_k) .GT. 0 ) 
           ims_bcs_imax = ims_bcs_imax + 1
        ENDDO
        WRITE(str,*) ims_bcs_imax
        str = 'Initialize MPI types for Ox BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
        CALL IO_WRITE_ASCII(lfile,str)
        isize_loc = ims_bcs_imax*jmax
        CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
             ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF
     
     IF ( .NOT. g(2)%periodic ) THEN ! Required for NRBCs in Oy
        id    = DNS_MPI_K_NRBCY
        isize_loc = MOD(imax,ims_npro_k)
        ims_bcs_jmax = 2*(inb_flow+inb_scal_array)
        DO WHILE ( MOD(isize_loc*ims_bcs_jmax,ims_npro_k) .GT. 0 ) 
           ims_bcs_jmax = ims_bcs_jmax + 1
        ENDDO
        WRITE(str,*) ims_bcs_jmax
        str = 'Initialize MPI types for Oy BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
        CALL IO_WRITE_ASCII(lfile,str)
        isize_loc = imax*ims_bcs_jmax
        CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
             ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF
#endif

! ###################################################################
! Compute reference pressure for Poinsot&Lele term in NR Bcs.
! ###################################################################
     bcs_p_imin = pbg%mean; bcs_p_imax = pbg%mean ! default values
     bcs_p_jmin = pbg%mean; bcs_p_jmax = pbg%mean
     bcs_p_kmin = pbg%mean; bcs_p_kmax = pbg%mean

     IF ( bcs_euler_drift .EQ. 1 ) THEN
        prefactor = C_05_R *(gama0-C_1_R) *mach*mach
        r1        = C_1_R
        r05       = C_05_R
        
! -------------------------------------------------------------------
! Note that buffer_?(1,1,1,4) is now the energy, and we need p
! -------------------------------------------------------------------

! -------------------------------------------------------------------
! Using buffer fields; bottom
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        DO j = 1,imax*BuffFlowJmin%size*kmax
           dummy = prefactor *( BuffFlowJmin%Ref(j,1,1,1)*BuffFlowJmin%Ref(j,1,1,1) &
                              + BuffFlowJmin%Ref(j,1,1,2)*BuffFlowJmin%Ref(j,1,1,2) &
                              + BuffFlowJmin%Ref(j,1,1,3)*BuffFlowJmin%Ref(j,1,1,3) )/BuffFlowJmin%Ref(j,1,1,5)
           txc(j,1) = ( BuffFlowJmin%Ref(j,1,1,4) - dummy )/BuffFlowJmin%Ref(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,imax*BuffFlowJmin%size*kmax
           txc(j,1) = BuffFlowJmin%Ref(j,1,1,4)/BuffFlowJmin%Ref(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,imax*BuffFlowJmin%size*kmax
           txc(j,4) = BuffScalJmin%Ref(j,1,1,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,BuffFlowJmin%size,kmax, txc(1,4),&
             txc(1,1), BuffFlowJmin%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,BuffFlowJmin%size,kmax, txc(1,4),&
             BuffFlowJmin%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(imax,BuffFlowJmin%size,kmax, BuffScalJmin%Ref, &
             txc(1,1), BuffFlowJmin%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,BuffFlowJmin%size,kmax, BuffScalJmin%Ref, &
             BuffFlowJmin%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
     bcs_p_jmin = AVG_IK(imax,BuffFlowJmin%size,kmax, i1, txc(1,1), g(1)%jac,g(3)%jac, area)

! reference plane for BCS at bottom: rho, u_i, p, z_i
     DO k = 1,kmax; DO i = 1,imax
        bcs_hb(i,k,1) = BuffFlowJmin%Ref(i,1,k,5)                ! density
        bcs_hb(i,k,2) = BuffFlowJmin%Ref(i,1,k,2) /bcs_hb(i,k,1) ! normal velocity, in this case v
        bcs_hb(i,k,3) = BuffFlowJmin%Ref(i,1,k,1) /bcs_hb(i,k,1)
        bcs_hb(i,k,4) = BuffFlowJmin%Ref(i,1,k,3) /bcs_hb(i,k,1)
        ip = i + (k-1)*imax*BuffFlowJmin%size
        bcs_hb(i,k,5) = txc(ip,1)                        ! pressure
        DO is = 1,inb_scal
           ip = inb_flow +is
           bcs_hb(i,k,ip) = BuffScalJmin%Ref(i,1,k,is) /bcs_hb(i,k,1)
        ENDDO
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Using buffer fields; top
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,imax*BuffFlowJmax%size*kmax
           dummy = prefactor *( BuffFlowJmax%Ref(j,1,1,1)*BuffFlowJmax%Ref(j,1,1,1) &
                              + BuffFlowJmax%Ref(j,1,1,2)*BuffFlowJmax%Ref(j,1,1,2) &
                              + BuffFlowJmax%Ref(j,1,1,3)*BuffFlowJmax%Ref(j,1,1,3) )/BuffFlowJmax%Ref(j,1,1,5)
           txc(j,1) = ( BuffFlowJmax%Ref(j,1,1,4) - dummy )/BuffFlowJmax%Ref(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,imax*BuffFlowJmax%size*kmax
           txc(j,1) = BuffFlowJmax%Ref(j,1,1,4)/BuffFlowJmax%Ref(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,imax*BuffFlowJmax%size*kmax
           txc(j,4) = BuffScalJmax%Ref(j,1,1,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(imax,BuffFlowJmax%size,kmax, txc(1,4),&
             txc(1,1), BuffFlowJmax%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,BuffFlowJmax%size,kmax, txc(1,4),&
             BuffFlowJmax%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(imax,BuffFlowJmax%size,kmax, BuffScalJmax%Ref, &
             txc(1,1), BuffFlowJmax%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(imax,BuffFlowJmax%size,kmax, BuffScalJmax%Ref, &
             BuffFlowJmax%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
     bcs_p_jmax = AVG_IK(imax,BuffFlowJmax%size,kmax, BuffFlowJmax%size, txc(1,1), g(1)%jac,g(3)%jac, area)

! reference plane for BCS at top: rho, u_i, p, z_i
     DO k = 1,kmax; DO i = 1,imax
        bcs_ht(i,k,1) = BuffFlowJmax%Ref(i,BuffFlowJmax%size,k,5)               ! density
        bcs_ht(i,k,2) = BuffFlowJmax%Ref(i,BuffFlowJmax%size,k,2)/bcs_ht(i,k,1) ! normal velocity, in this case v
        bcs_ht(i,k,3) = BuffFlowJmax%Ref(i,BuffFlowJmax%size,k,1)/bcs_ht(i,k,1)               
        bcs_ht(i,k,4) = BuffFlowJmax%Ref(i,BuffFlowJmax%size,k,3)/bcs_ht(i,k,1)
        ip = i + (BuffFlowJmax%size-1)*imax + (k-1)*imax*BuffFlowJmax%size
        bcs_ht(i,k,5) = txc(ip,1)                                    ! pressure
        DO is = 1,inb_scal
           ip = inb_flow +is
           bcs_ht(i,k,ip) = BuffScalJmax%Ref(i,BuffFlowJmax%size,k,is) /bcs_ht(i,k,1)
        ENDDO
     ENDDO; ENDDO

! ###################################################################
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN

! -------------------------------------------------------------------
! Using buffer fields; left
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,BuffFlowImin%size*jmax*kmax
           dummy = prefactor *( BuffFlowImin%Ref(j,1,1,1)*BuffFlowImin%Ref(j,1,1,1) &
                              + BuffFlowImin%Ref(j,1,1,2)*BuffFlowImin%Ref(j,1,1,2) &
                              + BuffFlowImin%Ref(j,1,1,3)*BuffFlowImin%Ref(j,1,1,3) )/BuffFlowImin%Ref(j,1,1,5)
           txc(j,1) = ( BuffFlowImin%Ref(j,1,1,4) - dummy )/BuffFlowImin%Ref(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,BuffFlowImin%size*jmax*kmax
           txc(j,1) = BuffFlowImin%Ref(j,1,1,4)/BuffFlowImin%Ref(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,BuffFlowImin%size*jmax*kmax
           txc(j,4) = BuffScalImin%Ref(j,1,1,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(BuffFlowImin%size,jmax,kmax, txc(1,4),&
             txc(1,1), BuffFlowImin%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(BuffFlowImin%size,jmax,kmax, txc(1,4),&
             BuffFlowImin%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(BuffFlowImin%size,jmax,kmax, BuffScalImin%Ref, &
             txc(1,1), BuffFlowImin%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(BuffFlowImin%size,jmax,kmax, BuffScalImin%Ref, &
             BuffFlowImin%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
! reference value at the center of the vertical length
     j = jmax/2
     bcs_p_imin = AVG1V1D(BuffFlowImin%size,jmax,kmax, i1,j, i1, txc(1,1))

! reference plane for BCS at xmin: rho, u_i, p, z_i
     DO k = 1,kmax; DO j = 1,jmax
        bcs_vi(j,k,1) = BuffFlowImin%Ref(1,j,k,5)                    ! density
        bcs_vi(j,k,2) = BuffFlowImin%Ref(1,j,k,1)/bcs_vi(j,k,1)
        bcs_vi(j,k,3) = BuffFlowImin%Ref(1,j,k,2)/bcs_vi(j,k,1)
        bcs_vi(j,k,4) = BuffFlowImin%Ref(1,j,k,3)/bcs_vi(j,k,1)
        ip = 1 + (j-1)*BuffFlowImin%size + (k-1)*BuffFlowImin%size*jmax
        bcs_vi(j,k,5) = txc(ip,1)                             ! pressure
        DO is = 1,inb_scal
           ip = inb_flow +is
           bcs_vi(j,k,ip) = BuffScalImin%Ref(1,j,k,is)/bcs_vi(j,k,1)
        ENDDO
     ENDDO; ENDDO

! shape factor
     diam_loc  = C_3_R*qbg(1)%diam
     thick_loc = qbg(1)%diam/C_8_R
     ycenter   = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean
     DO k = 1,kmax; DO j = 1,jmax
        bcs_vi(j,k,inb_flow+inb_scal+1) = FLOW_JET_TEMPORAL(i2, thick_loc, r1, r05, diam_loc, ycenter, dummy, g(2)%nodes(j))
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Using buffer fields; right
! -------------------------------------------------------------------
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO j = 1,BuffFlowImax%size*jmax*kmax
           dummy = prefactor *( BuffFlowImax%Ref(j,1,1,1)*BuffFlowImax%Ref(j,1,1,1) &
                              + BuffFlowImax%Ref(j,1,1,2)*BuffFlowImax%Ref(j,1,1,2) &
                              + BuffFlowImax%Ref(j,1,1,3)*BuffFlowImax%Ref(j,1,1,3) )/BuffFlowImax%Ref(j,1,1,5)
           txc(j,1) = ( BuffFlowImax%Ref(j,1,1,4) - dummy )/BuffFlowImax%Ref(j,1,1,5)
        ENDDO
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1,BuffFlowImax%size*jmax*kmax
           txc(j,1) = BuffFlowImax%Ref(j,1,1,4)/BuffFlowImax%Ref(j,1,1,5)
        ENDDO
     ENDIF
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! in AIRWATER case the s array is moved to txc4 because q_l in computed
! just after s(1) !!
        DO j = 1,BuffFlowImax%size*jmax*kmax
           txc(j,4) = BuffScalImax%Ref(j,1,1,1)
        ENDDO
        CALL THERMO_CALORIC_TEMPERATURE(BuffFlowImax%size,jmax,kmax, txc(1,4),&
             txc(1,1), BuffFlowImax%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(BuffFlowImax%size,jmax,kmax, txc(1,4),&
             BuffFlowImax%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(BuffFlowImax%size,jmax,kmax, BuffScalImax%Ref, &
             txc(1,1), BuffFlowImax%Ref(1,1,1,5), txc(1,2), txc(1,3))
        CALL THERMO_THERMAL_PRESSURE(BuffFlowImax%size,jmax,kmax, BuffScalImax%Ref, &
             BuffFlowImax%Ref(1,1,1,5), txc(1,2), txc(1,1))
     ENDIF
! reference value at the center of the vertical length
     j = jmax/2
     bcs_p_imax = AVG1V1D(BuffFlowImax%size,jmax,kmax, BuffFlowImax%size,j, i1, txc(1,1))

! reference plane for BCS at xmax: rho, u_i, p, z_i
     DO k = 1,kmax; DO j = 1,jmax
!        bcs_vo(j,k,1) = BuffFlowImax%Ref(BuffFlowImax%size,j,k,1)
!        bcs_vo(j,k,2) = BuffFlowImax%Ref(BuffFlowImax%size,j,k,2)/BuffFlowImax%Ref(BuffFlowImax%size,j,k,1)
!        bcs_vo(j,k,3) = BuffFlowImax%Ref(BuffFlowImax%size,j,k,3)/BuffFlowImax%Ref(BuffFlowImax%size,j,k,1)
!        bcs_vo(j,k,4) = BuffFlowImax%Ref(BuffFlowImax%size,j,k,4)/BuffFlowImax%Ref(BuffFlowImax%size,j,k,1)
!        ip = BuffFlowImax%size + (j-1)*BuffFlowImin%size + (k-1)*BuffFlowImin%size*jmax
!        bcs_vo(j,k,5) = txc(ip,1)
!        DO is = 1,inb_scal
!           ip = inb_flow +is
!           bcs_vo(j,k,ip) = BuffScalImax%Ref(BuffFlowImax%size,j,k,is)/BuffFlowImax%Ref(BuffFlowImax%size,j,k,1)
!        ENDDO
! try to use only the coflow values
        bcs_vo(j,k,1) = BuffFlowImax%Ref(BuffFlowImax%size,1,k,5)
        bcs_vo(j,k,2) = BuffFlowImax%Ref(BuffFlowImax%size,1,k,1)/bcs_vo(j,k,1)
        bcs_vo(j,k,3) = C_0_R
        bcs_vo(j,k,4) = BuffFlowImax%Ref(BuffFlowImax%size,1,k,3)/bcs_vo(j,k,1)
        ip = BuffFlowImax%size + (1-1)*BuffFlowImin%size + (k-1)*BuffFlowImin%size*jmax
        bcs_vo(j,k,5) = txc(ip,1)
        DO is = 1,inb_scal
           ip = inb_flow +is
           bcs_vo(j,k,ip) = BuffScalImax%Ref(BuffFlowImax%size,1,k,is)/bcs_vo(j,k,1)
        ENDDO
     ENDDO; ENDDO
     
  ENDIF

  ENDIF

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_INITIALLIZE' )
#endif

  RETURN

END SUBROUTINE BOUNDARY_BCS_INITIALIZE

END MODULE BOUNDARY_BCS

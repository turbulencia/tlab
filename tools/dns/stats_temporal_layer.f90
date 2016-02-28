#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif
#ifdef LES
#include "les_const.h"
#endif

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Checked to include RTI case
!# 2007/05/10 - J.P. Mellado
!#              Multispecies cases added
!# 2009/06/01 - J.P. Mellado
!#              Incompressible case added
!#
!########################################################################
!# DESCRIPTION
!#
!# Statistics that depend on time and the crosswise coordinate:
!# temporally evolving configurations (mixing layer, jet, rti...)
!# Isotropic case should be put in a different subroutine.
!#
!########################################################################
SUBROUTINE STATS_TEMPORAL_LAYER(x,y,z,dx,dy,dz, q,s,hq, txc, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : pointers_structure
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL, ONLY : fstavg, fstpdf, fstinter
  USE DNS_LOCAL, ONLY : vindex, VA_MEAN_WRK
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL, ONLY : iles, iles_type_chem
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*),                 INTENT(IN)    :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field,    *), INTENT(IN)    :: q,s
  TREAL, DIMENSION(isize_field,    *), INTENT(INOUT) :: hq ! auxiliary array
  TREAL, DIMENSION(isize_txc_field,6), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                 INTENT(INOUT) :: wrk1d,wrk2d,wrk3d, vaux
!  INTEGER(1), DIMENSION(*)                :: wrk3d                ! see comments in header
!# wrk3d    In     Aux array, with space for REAL(8) 3D field; declared
!#                 as INTEGER(1) to be used as gate array.
!#                 PROBLEMS IN JUQUEEN WITH THIS TRICK!!! Maybe use hs?
  
  TARGET :: q,s, txc

! -------------------------------------------------------------------
  TREAL umin,umax, dummy, s_aux(MAX_NSP)
  TINTEGER nbins, is, flag_buoyancy!, ij
  TINTEGER ibc(16), nfield
  TREAL amin(16), amax(16)
  CHARACTER*32 fname
  CHARACTER*32 varname(16)
  INTEGER(1) igate !, gate_levels(16)

  TYPE(pointers_structure), DIMENSION(16) :: data

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u,v,w,e,rho, p,T,vis

#ifdef USE_MPI
  TINTEGER id
#endif
#ifdef LES
  TINTEGER ichi
#endif

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING STATS_TEMPORAL_LAYER' )
#endif

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
  ENDIF

! in case we need the buoyancy statistics
  IF ( ibodyforce .EQ. EQNS_BOD_QUADRATIC          .OR. &
       ibodyforce .EQ. EQNS_BOD_BILINEAR           .OR. &
       ibodyforce .EQ. EQNS_BOD_PIECEWISE_LINEAR   .OR. &
       ibodyforce .EQ. EQNS_BOD_PIECEWISE_BILINEAR  .OR. &
       imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
     flag_buoyancy = 1
  ELSE 
     flag_buoyancy = 0   
  ENDIF

! ###################################################################
! Intermittency
! ###################################################################
!   IF ( fstinter .EQ. 1 ) THEN
!      CALL FI_VORTICITY(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
!           dx,dy,dz, u,v,w, txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)

! ! calculate vorticity gate based on 1% threshold
!      CALL MINMAX(imax,jmax,kmax, txc(1,1), umin,umax)
!      umin = C_1EM3_R*C_1EM3_R*umax
!      DO ij = 1,isize_field
!         IF ( txc(ij,1) .GT. umin ) THEN; wrk3d(ij) = 1  ! gate array
!         ELSE;                            wrk3d(ij) = 0; ENDIF
!      ENDDO
!      varname(1) = 'Vorticity'; gate_levels(1) = 1

!      WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
!      CALL INTER2D_N(fname, varname, rtime, imax,jmax,kmax, i1, y, wrk3d, gate_levels)

!   ENDIF

! ###################################################################
! Calculate pressure
! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN 
     CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, u,v,w,s,txc(:,3), txc(:,1),txc(:,2),txc(:,4), wrk1d,wrk2d,wrk3d)
     
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     p(:) = C_1_R ! to be developed
     
  ELSE
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)
     
  ENDIF

! ###################################################################
! Unconditional plane PDFs
! ###################################################################
  IF ( fstpdf .EQ. 1 ) THEN
     nfield = 0
     nfield = nfield+1; data(nfield)%field => u(:); varname(nfield) = 'Pu'
     nfield = nfield+1; data(nfield)%field => v(:); varname(nfield) = 'Pv'
     nfield = nfield+1; data(nfield)%field => w(:); varname(nfield) = 'Pw'
     nfield = nfield+1; data(nfield)%field => p(:); varname(nfield) = 'Pp'
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        nfield = nfield+1; data(nfield)%field => rho(:); varname(nfield) = 'Pr'
        nfield = nfield+1; data(nfield)%field => T(:);   varname(nfield) = 'Pt'
     ENDIF

     IF ( icalc_scal .EQ. 1 ) THEN
        nfield = nfield+1; data(nfield)%field => s(:,1);varname(nfield) = 'PScalar1'
     ENDIF

     ibc(1:nfield) = 2 ! BCs in the calculation of the PDFs
     igate = 0         ! no intermittency partition

     nbins = 32
     WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
     CALL PDF2D_N(fname, varname, igate, rtime, &
          imax,jmax,kmax, nfield, ibc, amin, amax, y, wrk3d, &
          data, nbins, isize_field, txc(:,1), wrk1d)
     
  ENDIF

! ###################################################################
! Plane averages
! ###################################################################
  IF ( fstavg .EQ. 1 ) THEN
     CALL AVG_FLOW_TEMPORAL_LAYER(y,dx,dy,dz, q, s,&
          txc(:,1),txc(:,2),txc(:,3),txc(1:,4),txc(:,5),txc(:,6),hq(:,1),hq(:,2),hq(:,3),  &
          vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)

     IF ( icalc_scal .EQ. 1 ) THEN
        IF ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN ! But this should be already done...
           CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array)) ! Update the liquid function
        ENDIF
        DO is = inb_scal+1,inb_scal_array ! Add diagnostic fields, if any
           mean_i(is) = C_1_R; delta_i(is) = C_0_R; ycoor_i(is) = ycoor_i(1); schmidt(is) = schmidt(1)
        ENDDO
        DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
           CALL AVG_SCAL_TEMPORAL_LAYER(is, y,dx,dy,dz, q,s, s(1,is), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4), vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)
        ENDDO

! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
        IF ( flag_buoyancy .EQ. 1 ) THEN
           wrk1d(1:jmax) = C_0_R
           CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, txc(:,1), wrk1d) ! note that wrk3d is defined as integer.
           dummy = C_1_R/froude
           txc(1:isize_field,1) = txc(1:isize_field,1)*dummy
! mean values
           s_aux(1:inb_scal) = mean_i(1:inb_scal) - C_05_R*delta_i(1:inb_scal)
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
              CALL THERMO_AIRWATER_LINEAR(i1,i1,i1, s_aux, s_aux(inb_scal_array), dummy)
           ENDIF
           dummy = C_0_R
           CALL FI_BUOYANCY(ibodyforce, i1,i1,i1, body_param, s_aux, umin, dummy)
           s_aux(1:inb_scal) = mean_i(1:inb_scal) + C_05_R*delta_i(1:inb_scal)
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
              CALL THERMO_AIRWATER_LINEAR(i1,i1,i1, s_aux, s_aux(inb_scal_array), dummy)
           ENDIF
           dummy = C_0_R
           CALL FI_BUOYANCY(ibodyforce, i1,i1,i1, body_param, s_aux, umax, dummy)
           mean_i(is) = (umax+umin)/froude; delta_i(is) = ABS(umax-umin)/froude; ycoor_i(is) = ycoor_i(1); schmidt(is) = schmidt(1)
           CALL AVG_SCAL_TEMPORAL_LAYER(is, y,dx,dy,dz, q,s, txc(1,1), &
                txc(1,2),txc(1,3),txc(1,4),txc(1,5), vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)
           
        ENDIF

     ENDIF

#ifdef LES
! -------------------------------------------------------------------
! LES
! -------------------------------------------------------------------
     IF ( iles .EQ. 1 ) THEN
        CALL LES_AVG_FLOW_TEMPORAL_LAYER&
             (x,y,z,dx,dy,dz, rho,u,v,w,e, hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
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
                   (is, x,y,z,dx,dy,dz, rho,u,v,w, s(1,is), hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
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
END SUBROUTINE STATS_TEMPORAL_LAYER

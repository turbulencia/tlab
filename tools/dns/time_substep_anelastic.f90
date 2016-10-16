#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2012/12/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Extracted from TIME_SUBSTEP_INCOMPRESSIBLE
!#
!########################################################################
!# ARGUMENTS 
!#
!# txc     Aux   3D array size 6
!#
!########################################################################
SUBROUTINE TIME_SUBSTEP_ANELASTIC(dte,etime, x,y,z,dx,dy,dz, q,hq,s,hs, &
     txc, vaux, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE

  TREAL dte,etime
  TREAL, DIMENSION(*)                 :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field, *)    :: q,hq, s,hs
  TREAL, DIMENSION(isize_txc_field,*) :: txc
  TREAL, DIMENSION(*)                 :: wrk1d,wrk2d,wrk3d, vaux

  TARGET :: q

! -----------------------------------------------------------------------
  TINTEGER is

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w

#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_SUBSTEP_ANELASTIC')
#endif

#ifdef USE_BLAS
  ilen = isize_field
#endif

! Define pointers
  u => q(:,1)
  v => q(:,2)
  w => q(:,3)

! #######################################################################
! Evaluate standard RHS of anelastic equations
! #######################################################################
! -----------------------------------------------------------------------
  IF ( iadvection .EQ. EQNS_CONVECTIVE .AND. &
       iviscous   .EQ. EQNS_EXPLICIT   .AND. & 
       idiffusion .EQ. EQNS_EXPLICIT         ) THEN
     CALL RHS_FLOW_GLOBAL_ANELASTIC_1(dte, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3), s, &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
          vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), &
          wrk1d,wrk2d,wrk3d)
     
     DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_ANELASTIC_1(is, dx,dy,dz, u,v,w,s(1,is),hs(1,is), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
     ENDDO
     
! -----------------------------------------------------------------------
  ELSE
     CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_ANELASTIC. Undeveloped formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     
  ENDIF

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     ! Flow part needs to be taken into account in the pressure
     DO is = 1,inb_scal
     CALL BOUNDARY_BUFFER_RELAXATION_SCAL(is,&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), x,y, q,s(1,is),hs(1,is))
     ENDDO
  ENDIF

! #######################################################################
! Perform the time stepping for anelastic equations
! wrk1d contains reference density profile
! #######################################################################
! not yet develop; you need the reference density profile r_ref

!     DO is = 1,inb_flow
!        DO ij = 1,isize_field
!           j = MOD(ij-1,imax*jmax)+1; j = (j-1)/imax+1
!           q(ij,is) = q(ij,is) + dte*hq(ij,is)/r_ref(j)
!        ENDDO
!     ENDDO
!
!     DO is = 1,inb_scal
!        DO ij = 1,isize_field
!           j = MOD(ij-1,imax*jmax)+1; j = (j-1)/imax+1
!           s(ij,is) = s(ij,is) + dte*hs(ij,is)/r_ref(j)
!        ENDDO
!     ENDDO

! ###################################################################
! Calculate other intensive thermodynamic variables
! ###################################################################
  IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_SUPSAT ) THEN
     IF ( damkohler(1) .LE. C_0_R )  THEN
        CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(1,2), p_init, s(1,1))
     ENDIF
     
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
     CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))

  ENDIF

! ###################################################################
! Simulation control
! ###################################################################
  CALL DNS_CONTROL_SCAL(s)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_SUBSTEP_ANELASTIC')
#endif

  RETURN
END SUBROUTINE TIME_SUBSTEP_ANELASTIC

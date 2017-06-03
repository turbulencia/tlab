#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# The case of internal energy adds the term p div u here to avoid the 
!# computation of the dilatation in RHS_FLOW_EULER_?, so it is
!# not only the viscous part in that case...
!# Internal energy eqn formulation does 18 derivatives.
!# Total energy eqn formulation does 21 derivatives.
!#
!########################################################################
SUBROUTINE RHS_FLOW_VISCOUS_DIVERGENCE&
     (vis, u,v,w,p, h1,h2,h3,h4, tau_xx,tau_xy,tau_xz,tau_yy,tau_yz,tau_zz,&
     tmp1,tmp2,tmp3, wrk2d,wrk3d)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field, imode_eqns
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : visc, mach
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,     ONLY : bcs_out, bcs_inf

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field), INTENT(IN)    :: vis, u,v,w,p
  TREAL, DIMENSION(isize_field), INTENT(OUT)   :: h1,h2,h3,h4
  TREAL, DIMENSION(isize_field), INTENT(OUT)   :: tau_xx,tau_xy,tau_xz,tau_yy,tau_yz,tau_zz
  TREAL, DIMENSION(isize_field), INTENT(INOUT) :: tmp1,tmp2,tmp3
  TREAL, DIMENSION(*),           INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1)
  TINTEGER i
  TREAL prefactor, dil, c23

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

  bcs = 0
  prefactor = (gama0-C_1_R)*mach*mach

! ###################################################################
! Define viscous stress tensor
! Add corresponding terms in the internal energy equation if needed
! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  c23 = C_2_R/C_3_R*visc
  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xx(i) = c23*vis(i)*( C_2_R*tmp1(i) - (tmp2(i) + tmp3(i)) )
        tau_yy(i) = c23*vis(i)*( C_2_R*tmp2(i) - (tmp1(i) + tmp3(i)) )
        tau_zz(i) = c23*vis(i)*( C_2_R*tmp3(i) - (tmp1(i) + tmp2(i)) )
        dil = tmp1(i) + tmp2(i) + tmp3(i)
        h4(i) = h4(i) + prefactor*( tau_xx(i)*tmp1(i) + tau_yy(i)*tmp2(i) + tau_zz(i)*tmp3(i) - p(i)*dil )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN ! total energy equation
     DO i = 1,imax*jmax*kmax
        tau_xx(i) = c23*vis(i)*( C_2_R*tmp1(i) - (tmp2(i) + tmp3(i)) )
        tau_yy(i) = c23*vis(i)*( C_2_R*tmp2(i) - (tmp1(i) + tmp3(i)) )
        tau_zz(i) = c23*vis(i)*( C_2_R*tmp3(i) - (tmp1(i) + tmp2(i)) )
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL   ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xy(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_xy(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN ! total energy equation
     tau_xy = visc *vis *( tmp1 + tmp2 )

  ENDIF

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)

  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_xz(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN ! total energy equation
     tau_xz = visc*vis*( tmp1 + tmp2 )
     
  ENDIF

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)

  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_yz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_yz(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN ! total energy equation
     tau_yz = visc*vis*( tmp1 + tmp2 )

  ENDIF

! ###################################################################
! Momentum equation
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,1), g(1), tau_xx, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_out(1,2,2), g(2), tau_xy, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_out(1,2,3), g(3), tau_xz, tmp3, wrk3d, wrk2d,wrk3d)
  h1 = h1 + tmp1 + tmp2 + tmp3

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_out(1,2,1), g(1), tau_xy, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,2), g(2), tau_yy, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_out(1,2,3), g(3), tau_yz, tmp3, wrk3d, wrk2d,wrk3d)
  h2 = h2 + tmp1 + tmp2 + tmp3

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_out(1,2,1), g(1), tau_xz, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_out(1,2,2), g(2), tau_yz, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,3), g(3), tau_zz, tmp3, wrk3d, wrk2d,wrk3d)
  h3 = h3 + tmp1 + tmp2 + tmp3

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     DO i = 1,imax*jmax*kmax
        tau_xx(i) = tau_xx(i)*u(i) + tau_xy(i)*v(i) + tau_xz(i)*w(i)
        tau_yy(i) = tau_xy(i)*u(i) + tau_yy(i)*v(i) + tau_yz(i)*w(i)
        tau_zz(i) = tau_xz(i)*u(i) + tau_yz(i)*v(i) + tau_zz(i)*w(i)
     ENDDO
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tau_xx, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tau_yy, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tau_zz, tmp1, wrk3d, wrk2d,wrk3d)
     h4 = h4 + prefactor*( tmp1+tmp2+tmp3 )

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_VISCOUS_DIVERGENCE

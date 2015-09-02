!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/08 - J.P. Mellado
!#              Created
!# 2007/08/16 - J.P. Mellado
!#              Case of internal energy formulation added
!# 2007/08/30 - J.P. Mellado
!#              Name change
!#
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
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE RHS_FLOW_VISCOUS_DIVERGENCE&
     (dx,dy,dz, vis, u,v,w,p, h1,h2,h3,h4, tau_xx,tau_xy,tau_xz,tau_yy,tau_yz,tau_zz,&
     tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL vis(*), u(*), v(*), w(*), p(*)
  TREAL h1(*), h2(*), h3(*), h4(*)
  TREAL tmp1(*), tmp2(*), tmp3(*)
  TREAL tau_xx(*), tau_xy(*), tau_xz(*), tau_yy(*), tau_yz(*), tau_zz(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsin, i1vsout, imxvsin, imxvsout
  TINTEGER j1vsin, j1vsout, jmxvsin, jmxvsout
  TINTEGER k1vsin, k1vsout, kmxvsin, kmxvsout
  TINTEGER i
  TREAL prefactor, dil, c23

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

  prefactor = (gama0-C_1_R)*mach*mach

#include "dns_bcs_inf.h"
#include "dns_bcs_out.h"

! ###################################################################
! Define viscous stress tensor
! Add corresponding terms in the internal energy equation if needed
! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  c23 = C_2_R/C_3_R*visc
  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xx(i) = c23*vis(i)*( C_2_R*tmp1(i) - (tmp2(i) + tmp3(i)) )
        tau_yy(i) = c23*vis(i)*( C_2_R*tmp2(i) - (tmp1(i) + tmp3(i)) )
        tau_zz(i) = c23*vis(i)*( C_2_R*tmp3(i) - (tmp1(i) + tmp2(i)) )
        dil = tmp1(i) + tmp2(i) + tmp3(i)
        h4(i) = h4(i) + prefactor*( &
             tau_xx(i)*tmp1(i) + tau_yy(i)*tmp2(i) + tau_zz(i)*tmp3(i) - p(i)*dil )
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
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL   ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xy(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_xy(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN ! total energy equation
     DO i = 1,imax*jmax*kmax
        tau_xy(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
     ENDDO

  ENDIF

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)

  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_xz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_xz(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN ! total energy equation
     DO i = 1,imax*jmax*kmax
        tau_xz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
     ENDDO
  ENDIF

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)

  IF      ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! internal energy equation
     DO i = 1,imax*jmax*kmax
        tau_yz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
        h4(i) = h4(i) + prefactor*( tau_yz(i)*(tmp1(i)+tmp2(i)) )
     ENDDO

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN ! total energy equation
     DO i = 1,imax*jmax*kmax
        tau_yz(i) = visc*vis(i)*( tmp1(i) + tmp2(i) )
     ENDDO
  ENDIF

! ###################################################################
! Momentum equation
! ###################################################################
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tau_xx, tmp1, i1vsin, imxvsin, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tau_xy, tmp2, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tau_xz, tmp3, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + tmp1(i) + tmp2(i) + tmp3(i)
  ENDDO

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tau_xy, tmp1, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tau_yy, tmp2, j1vsin, jmxvsin, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tau_yz, tmp3, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + tmp1(i) + tmp2(i) + tmp3(i)
  ENDDO

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tau_xz, tmp1, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tau_yz, tmp2, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tau_zz, tmp3, k1vsin, kmxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) + tmp1(i) + tmp2(i) + tmp3(i)
  ENDDO

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
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tau_xx, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tau_yy, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tau_zz, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + prefactor*( tmp1(i)+tmp2(i)+tmp3(i) )
     ENDDO
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_VISCOUS_DIVERGENCE

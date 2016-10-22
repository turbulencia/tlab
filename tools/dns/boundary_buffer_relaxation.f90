#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/12/01 - J.P. Mellado
!#              Modified to expand buffer_vo to 3D field
!# 2007/03/25 - J.P. Mellado
!#              Arrays are expanded in direction perpendicular to the
!#              remaining boundary.
!# 2007/04/03 - J.P. Mellado
!#              Spatial and temporal modes are unified: brach base on
!#              number of points in corresponding buffer
!# 2007/06/14 - J.P. Mellado
!#              Energy formulation. Instead of pressure, energy e.
!#
!########################################################################
!# DESCRIPTION
!#
!# Note that with the energy formulation the forcing terms have been
!# made proportional to differences on the conserved variables 
!# themselves, i.e. rho, rho*V, rho*(e+v^2/2). 
!# The scalars not yes, they still have the formulation in s.
!# In this respect, inflow and outflow cases must be reviewed.
!#
!########################################################################
SUBROUTINE BOUNDARY_BUFFER_RELAXATION_FLOW(buffer_ht,buffer_hb,buffer_vi,buffer_vo, q,hq)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imode_eqns
  USE DNS_GLOBAL, ONLY : mach
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL

  IMPLICIT NONE

  TREAL, DIMENSION(imax,         buff_nps_jmin,kmax,*) :: buffer_hb
  TREAL, DIMENSION(imax,         buff_nps_jmax,kmax,*) :: buffer_ht
  TREAL, DIMENSION(buff_nps_imin,jmax,         kmax,*) :: buffer_vi
  TREAL, DIMENSION(buff_nps_imax,jmax,         kmax,*) :: buffer_vo

  TREAL, DIMENSION(imax,jmax,kmax,*) :: q,hq

  TARGET q

! -------------------------------------------------------------------
  TINTEGER i,j, iloc,jloc, k
  TREAL lb, sigma, prefactor, dummy

  TREAL, DIMENSION(:,:,:), POINTER :: u, v, w, e, rho

! ###################################################################
  prefactor = (gama0-C_1_R)*mach*mach

  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,:,:,4)
     rho => q(:,:,:,5)
  ENDIF

! ###################################################################
! Bottom boundary
! ###################################################################
  IF ( buff_nps_jmin .GT. 1 ) THEN
     lb = g(2)%nodes(buff_u_jmin) - g(2)%nodes(1)
     DO j = 1,buff_u_jmin-1
        sigma = buff_param_u(1)*((g(2)%nodes(buff_u_jmin)-g(2)%nodes(j))/lb)**buff_param_u(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(:,j,:,1) = hq(:,j,:,1) - sigma*(rho(:,j,:)*u(:,j,:)-buffer_hb(:,j,:,1))
           IF ( buff_v_free .NE. 1 ) &
           hq(:,j,:,2) = hq(:,j,:,2) - sigma*(rho(:,j,:)*v(:,j,:)-buffer_hb(:,j,:,2))
           hq(:,j,:,3) = hq(:,j,:,3) - sigma*(rho(:,j,:)*w(:,j,:)-buffer_hb(:,j,:,3))
                   
        ELSE ! do not use the density
           hq(:,j,:,1) = hq(:,j,:,1) - sigma*(u(:,j,:)-buffer_hb(:,j,:,1))
           IF ( buff_v_free .NE. 1 ) &
           hq(:,j,:,2) = hq(:,j,:,2) - sigma*(v(:,j,:)-buffer_hb(:,j,:,2))
           hq(:,j,:,3) = hq(:,j,:,3) - sigma*(w(:,j,:)-buffer_hb(:,j,:,3))

        ENDIF

     ENDDO

! if compressible
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     lb = g(2)%nodes(buff_e_jmin) - g(2)%nodes(1)
     DO j = 1,buff_e_jmin-1
        sigma = buff_param_u(1)*((g(2)%nodes(buff_e_jmin)-g(2)%nodes(j))/lb)**buff_param_u(2)

        IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
           DO k = 1,kmax; DO i = 1,imax
              dummy = C_05_R*prefactor*rho(i,j,k)*&
                   ( u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
              hq(i,j,k,4) = hq(i,j,k,4) - sigma*( rho(i,j,k)*e(i,j,k)+dummy-buffer_hb(i,j,k,4) )
              hq(i,j,k,5) = hq(i,j,k,5) - sigma*( rho(i,j,k)               -buffer_hb(i,j,k,5) )
           ENDDO; ENDDO

        ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(:,j,:,4) = hq(:,j,:,4) - sigma*( rho(:,j,:)*e(:,j,:)-buffer_hb(:,j,:,4) )
           hq(:,j,:,5) = hq(:,j,:,5) - sigma*( rho(:,j,:)         -buffer_hb(:,j,:,5) )

        ENDIF

     ENDDO
     ENDIF

  ENDIF

! ###################################################################
! Top boundary
! ###################################################################
  IF ( buff_nps_jmax .GT. 1 ) THEN
     lb = g(2)%nodes(jmax)-g(2)%nodes(buff_u_jmax)
     DO j = buff_u_jmax+1,jmax
        jloc = j - buff_u_jmax + 1
        sigma = buff_param_u(1)*((g(2)%nodes(j)-g(2)%nodes(buff_u_jmax))/lb)**buff_param_u(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(:,j,:,1) = hq(:,j,:,1) - sigma*(rho(:,j,:)*u(:,j,:)-buffer_ht(:,jloc,:,1))
           IF ( buff_v_free .NE. 1 ) &
           hq(:,j,:,2) = hq(:,j,:,2) - sigma*(rho(:,j,:)*v(:,j,:)-buffer_ht(:,jloc,:,2))
           hq(:,j,:,3) = hq(:,j,:,3) - sigma*(rho(:,j,:)*w(:,j,:)-buffer_ht(:,jloc,:,3))

        ELSE ! do not use the density
           hq(:,j,:,1) = hq(:,j,:,1) - sigma*(u(:,j,:)-buffer_ht(:,jloc,:,1))
           IF ( buff_v_free .NE. 1 ) &
           hq(:,j,:,2) = hq(:,j,:,2) - sigma*(v(:,j,:)-buffer_ht(:,jloc,:,2))
           hq(:,j,:,3) = hq(:,j,:,3) - sigma*(w(:,j,:)-buffer_ht(:,jloc,:,3))

        ENDIF

     ENDDO

! if compressible
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     lb = g(2)%nodes(jmax)-g(2)%nodes(buff_e_jmax)
     DO j = buff_e_jmax+1,jmax
        jloc = j - buff_e_jmax + 1
        sigma = buff_param_u(1)*((g(2)%nodes(j)-g(2)%nodes(buff_e_jmax))/lb)**buff_param_u(2)

        IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
           DO k = 1,kmax; DO i=1,imax
              dummy = C_05_R*prefactor*rho(i,j,k)*&
                   ( u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
              hq(i,j,k,4) = hq(i,j,k,4) - sigma*( rho(i,j,k)*e(i,j,k)+dummy-buffer_ht(i,jloc,k,4) )
              hq(i,j,k,5) = hq(i,j,k,5) - sigma*( rho(i,j,k)               -buffer_ht(i,jloc,k,5) )
           ENDDO; ENDDO

        ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(:,j,:,4) = hq(:,j,:,4) - sigma*( rho(:,j,:)*e(:,j,:)-buffer_ht(:,jloc,:,4) )
           hq(:,j,:,5) = hq(:,j,:,5) - sigma*( rho(:,j,:)         -buffer_ht(:,jloc,:,5) )

        ENDIF

     ENDDO
     ENDIF

  ENDIF

! ###################################################################
! Inflow boundary
! ###################################################################
  IF ( buff_nps_imin .GT. 1 ) THEN
     DO i = 1,buff_imin-1
        lb = g(1)%nodes(buff_imin)-g(1)%nodes(1)
        sigma = ((g(1)%nodes(buff_imin)-g(1)%nodes(i))/lb)**buff_param_u(2)

        hq(i,:,:,1) = hq(i,:,:,1) - sigma*(rho(i,:,:)*u(i,:,:)-buffer_vi(i,:,:,1))
        IF ( buff_v_free .NE. 1 ) &
        hq(i,:,:,2) = hq(i,:,:,2) - sigma*(rho(i,:,:)*v(i,:,:)-buffer_vi(i,:,:,2))
        hq(i,:,:,3) = hq(i,:,:,3) - sigma*(rho(i,:,:)*w(i,:,:)-buffer_vi(i,:,:,3))

! if compressible
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(i,:,:,4) = hq(i,:,:,4) - sigma*(rho(i,:,:)*e(i,:,:)-buffer_vi(i,:,:,4))
           hq(i,:,:,5) = hq(i,:,:,5) - sigma*(rho(i,:,:)         -buffer_vi(i,:,:,5))
        ENDIF

     ENDDO

  ENDIF

! ###################################################################
! Outflow boundary
! ###################################################################
  IF ( buff_nps_imax .GT. 1 ) THEN
     DO i = buff_imax+1,imax
        iloc = i-buff_imax+1
        lb = g(1)%nodes(imax)-g(1)%nodes(buff_imax)
        sigma = buff_param_u(1)*((g(1)%nodes(i)-g(1)%nodes(buff_imax))/lb)**buff_param_u(2)

        hq(i,:,:,1) = hq(i,:,:,1) - sigma*(rho(i,:,:)*u(i,:,:)-buffer_vo(iloc,:,:,1))
        IF ( buff_v_free .NE. 1 ) &
        hq(i,:,:,2) = hq(i,:,:,2) - sigma*(rho(i,:,:)*v(i,:,:)-buffer_vo(iloc,:,:,2))
        hq(i,:,:,3) = hq(i,:,:,3) - sigma*(rho(i,:,:)*w(i,:,:)-buffer_vo(iloc,:,:,3))

! if compressible
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           hq(i,:,:,4) = hq(i,:,:,4) - sigma*(rho(i,:,:)*e(i,:,:)-buffer_vo(iloc,:,:,4))
           hq(i,:,:,5) = hq(i,:,:,5) - sigma*(rho(i,:,:)         -buffer_vo(iloc,:,:,5))
        ENDIF

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BUFFER_RELAXATION_FLOW

! #######################################################################
! Scalar
! #######################################################################
SUBROUTINE BOUNDARY_BUFFER_RELAXATION_SCAL(is, buffer_ht,buffer_hb,buffer_vi,buffer_vo, q, s,hs)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imode_eqns, inb_flow
  USE DNS_LOCAL

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(imax,         buff_nps_jmin,kmax,*) :: buffer_hb
  TREAL, DIMENSION(imax,         buff_nps_jmax,kmax,*) :: buffer_ht
  TREAL, DIMENSION(buff_nps_imin,jmax,         kmax,*) :: buffer_vi
  TREAL, DIMENSION(buff_nps_imax,jmax,         kmax,*) :: buffer_vo

  TREAL, DIMENSION(imax,jmax,kmax)           :: s,hs
  TREAL, DIMENSION(imax,jmax,kmax,*), TARGET :: q ! for the density

! -------------------------------------------------------------------
  TINTEGER i,j, iloc,jloc
  TREAL lb, sigma

  TREAL, DIMENSION(:,:,:), POINTER :: rho

! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     rho => q(:,:,:,5) ! density
  ENDIF

! ###################################################################
! Bottom boundary
! ###################################################################
  IF ( buff_nps_jmin .GT. 1 ) THEN
     lb = g(2)%nodes(buff_s_jmin) - g(2)%nodes(1)
     DO j = 1,buff_s_jmin-1
        sigma = buff_param_s(1)*((g(2)%nodes(buff_s_jmin)-g(2)%nodes(j))/lb)**buff_param_s(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        hs(:,j,:) = hs(:,j,:) - sigma*(rho(:,j,:)*s(:,j,:)-buffer_hb(:,j,:,inb_flow+is))
        
        ELSE ! do not use the density
        hs(:,j,:) = hs(:,j,:) - sigma*(           s(:,j,:)-buffer_hb(:,j,:,inb_flow+is))

        ENDIF

     ENDDO
  ENDIF

! ###################################################################
! Top boundary
! ###################################################################
  IF ( buff_nps_jmax .GT. 1 ) THEN
     lb = g(2)%nodes(jmax)-g(2)%nodes(buff_s_jmax)
     DO j = buff_s_jmax+1,jmax
        jloc = j - buff_s_jmax + 1
        sigma = buff_param_s(1)*((g(2)%nodes(j)-g(2)%nodes(buff_s_jmax))/lb)**buff_param_s(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        hs(:,j,:) = hs(:,j,:) - sigma*(rho(:,j,:)*s(:,j,:)-buffer_ht(:,jloc,:,inb_flow+is))

        ELSE ! do not use the density
        hs(:,j,:) = hs(:,j,:) - sigma*(           s(:,j,:)-buffer_ht(:,jloc,:,inb_flow+is))

        ENDIF

     ENDDO
  ENDIF

! ###################################################################
! Inflow boundary
! ###################################################################
  IF ( buff_nps_imin .GT. 1 ) THEN
     DO i = 1,buff_imin-1
        lb = g(1)%nodes(buff_imin)-g(1)%nodes(1)
        sigma = buff_param_s(1)*((g(1)%nodes(buff_imin)-g(1)%nodes(i))/lb)**buff_param_s(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        hs(i,:,:) = hs(i,:,:) - sigma*(rho(i,:,:)*s(i,:,:)-buffer_vi(i,:,:,inb_flow+is))

        ELSE ! do not use the density
        hs(i,:,:) = hs(i,:,:) - sigma*(           s(i,:,:)-buffer_vi(i,:,:,inb_flow+is))

        ENDIF

     ENDDO
  ENDIF

! ###################################################################
! Outflow boundary
! ###################################################################
  IF ( buff_nps_imax .GT. 1 ) THEN
     DO i = buff_imax+1,imax
        iloc = i-buff_imax+1
        lb = g(1)%nodes(imax)-g(1)%nodes(buff_imax)
        sigma = buff_param_s(1)*((g(1)%nodes(i)-g(1)%nodes(buff_imax))/lb)**buff_param_s(2)

        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        hs(i,:,:) = hs(i,:,:) - sigma*(rho(i,:,:)*s(i,:,:)-buffer_vo(iloc,:,:,inb_flow+is))

        ELSE ! do not use the density
        hs(i,:,:) = hs(i,:,:) - sigma*(           s(i,:,:)-buffer_vo(iloc,:,:,inb_flow+is))

        ENDIF

     ENDDO
  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BUFFER_RELAXATION_SCAL

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE AVG_FLOW_TEMPORAL_ISOTROPIC(dx,dy,dz, rho,u,v,w,p, &
     vis, dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, wrk1d,wrk2d,wrk3d)

! I have commmented out terms with temperature field T because
! I have removed this array (JP Mellado, 2008/11/20

! ##############################################################
! # Created from avgj.f for isotropic case.
! # 
! # Juan Pedro Mellado, Nov. 2001
! ##############################################################

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
!  TINTEGER itime, imax, jmax, kmax, i1bc, j1bc, k1bc, imode_fdm
!  TREAL rtime, vol, visc, gama0, cond, mean_rho, pbg%mean

  TREAL, DIMENSION(*)                :: dx, dy, dz
  TREAL, DIMENSION(imax,jmax,kmax)   :: u, v, w, rho, p, vis
  TREAL, DIMENSION(imax,jmax,kmax)   :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, wrk3d

  TREAL wrk1d(jmax)
  TREAL wrk2d(imax,kmax,5)

! -------------------------------------------------------------------
  TREAL rU
  TREAL rV
  TREAL rW
  TREAL rP
  TREAL rR

  TREAL fU
  TREAL fV
  TREAL fW

  TREAL rUf
  TREAL rVf
  TREAL rWf

  TREAL rP2

  TREAL Rxx
  TREAL Ryy
  TREAL Rzz
  TREAL Rxy
  TREAL Rxz
  TREAL Ryz

  TREAL Ep2
  TREAL rT 
  TREAL Tp2
  TREAL Sp2
  TREAL fT 

  TREAL Tau_xx 
  TREAL Tau_yy 
  TREAL Tau_zz 
  TREAL Tau_xy 
  TREAL Tau_xz 
  TREAL Tau_yz 

  TREAL rR2    

  TINTEGER i, k
  TREAL AVG1V3D
  TREAL lambda_x, lambda_y, lambda_z
  TREAL lk
  TREAL Exy, PIxy, Rxy_t
  TREAL Kin, Kin_t, Eps, Pi
  TREAL DV2, Phi
  TREAL Pp2, Pip2, PPhi
  TREAL entropy, entropy2
  TREAL skew_x, skew_y, skew_z
  TREAL rRP, rRT, rT2, fT2
  TREAL pp, vs, dil, relambda
  TREAL tau11, tau22, tau33, tau12, tau23, tau13

  TREAL ccnd, ccnd2, c23
  TREAL c2, rho_prime, p_prime, T_prime

  TREAL rho_ac, rho_en, T_ac, T_en, M_t
  CHARACTER*32 fname
#ifdef USE_MPI
  INTEGER ims_err, ims_pro
#endif

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 1')
#endif

  c23 = C_2_R/C_3_R
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; ccnd = C_0_R
  ELSE;                                  ccnd = visc/(mach*mach*prandtl); ENDIF
  ccnd2 = ccnd*C_2_R

  WRITE(fname,*) itime; fname='avg'//TRIM(ADJUSTL(fname))

#ifdef USE_MPI
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)

  IF ( ims_pro .EQ. 0 ) THEN
#endif

#ifdef USE_RECLEN
     OPEN(UNIT=23, RECL=1480, FILE=fname, STATUS='unknown')
#else
     OPEN(UNIT=23, FILE=fname, STATUS='unknown')
#endif

     WRITE(23, '(A8,E14.7E3)') 'RTIME = ', rtime

     WRITE(23, 1010) 'GROUP = MEAN rR rU rV rW rP rT fU fV fW fT'
     WRITE(23, 1010) 'GROUP = ACOUSTIC rP2 rR2 rT2 fT2 DV2 '&
          //'Rho_ac Rho_en T_ac T_en M_t'
     WRITE(23, 1010) 'GROUP = REYNOLDS Rxx Ryy Rzz Rxy Rxz Ryz'
     WRITE(23, 1010) 'GROUP = CROSS Rxy Rxy_t Exy PIxy'
     WRITE(23, 1010) 'GROUP = TKE Kin Kin_t Eps Pi'
     WRITE(23, 1010) 'GROUP = MOMENTS ReLambda '&
          //'Lambda_x Lambda_y Lambda_z '&
          //'Skew_x Skew_y Skew_z'

     WRITE(23, 1020)

1010 FORMAT(A)
1020 FORMAT('I J Lk ReLambda Lambda_x Lambda_y Lambda_z ',&
          'rR rU rV rW rP rT fU fV fW fT rP2 rR2 rT2 fT2 ',&
          'Rxx Ryy Rzz Rxy Rxz Ryz ',&
          'Rxy_t Exy PIxy ',&
          'Kin Kin_t Eps Pi ',&
          'DV2 Tau_xx Tau_yy Tau_zz Tau_xy Tau_xz Tau_yz ',&
          'Phi rRP rRT Pp2 Ep2 Pip2 PPhi Sp2 ',&
          'Skew_x Skew_y Skew_z Entropy Entropy2 ',&
          'Rho_ac Rho_en T_ac T_en M_t')

#ifdef USE_MPI
  ENDIF
#endif

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 2')
#endif

! #############################################################
! # Temporary Array Storage
! #
! # dudx = unused
! # dudy = unused
! # dudz = unused
! # dvdx = unused
! # dvdy = unused
! # dvdz = unused
! # dwdx = unused
! # dwdy = unused
! # dwdz = unused
! #
! #############################################################

! ###########################################################
! #    Calculate Mean                                       #
! ###########################################################

  rR = AVG1V3D(imax, jmax, kmax, i1, rho)
  rU = AVG1V3D(imax, jmax, kmax, i1, u)
  rV = AVG1V3D(imax, jmax, kmax, i1, v)
  rW = AVG1V3D(imax, jmax, kmax, i1, w)
  rP = AVG1V3D(imax, jmax, kmax, i1, p)
!  rT = AVG1V3D(imax, jmax, kmax, i1, T)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) = rho(i,1,k)*u(i,1,k)
        dudy(i,1,k) = rho(i,1,k)*v(i,1,k)
        dudz(i,1,k) = rho(i,1,k)*w(i,1,k)
!        dvdx(i,1,k) = rho(i,1,k)*T(i,1,k)
     ENDDO
  ENDDO

  fU = AVG1V3D(imax, jmax, kmax, i1, dudx)/rR
  fV = AVG1V3D(imax, jmax, kmax, i1, dudy)/rR
  fW = AVG1V3D(imax, jmax, kmax, i1, dudz)/rR
  fT = AVG1V3D(imax, jmax, kmax, i1, dvdx)/rR

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 3')
#endif

! #############################################################
! # Temporary Array Storage
! #
! # dudx = unused
! # dudy = unused
! # dudz = unused
! # dvdx = unused
! # dvdy = unused
! # dvdz = unused
! # dwdx = U_prime
! # dwdy = V_prime
! # dwdz = W_prime
! #
! #############################################################

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dwdx(i,1,k) = u(i,1,k) - fU
        dwdy(i,1,k) = v(i,1,k) - fV
        dwdz(i,1,k) = w(i,1,k) - fW
     ENDDO
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 4')
#endif

  rUf = AVG1V3D(imax, jmax, kmax, i1, dwdx)
  rVf = AVG1V3D(imax, jmax, kmax, i1, dwdy)
  rWf = AVG1V3D(imax, jmax, kmax, i1, dwdz)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) = rho(i,1,k)*dwdx(i,1,k)**2
        dudy(i,1,k) = rho(i,1,k)*dwdy(i,1,k)**2
        dudz(i,1,k) = rho(i,1,k)*dwdz(i,1,k)**2
        dvdx(i,1,k) = (p(i,1,k) - rP)**2
     ENDDO
  ENDDO

  Rxx = AVG1V3D(imax, jmax, kmax, i1, dudx)/rR
  Ryy = AVG1V3D(imax, jmax, kmax, i1, dudy)/rR
  Rzz = AVG1V3D(imax, jmax, kmax, i1, dudz)/rR
  rP2 = AVG1V3D(imax, jmax, kmax, i1, dvdx)

  Kin = C_05_R*(Rxx+Ryy+Rzz)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) = rho(i,1,k)*dwdx(i,1,k)*dwdy(i,1,k)
        dudy(i,1,k) = rho(i,1,k)*dwdx(i,1,k)*dwdz(i,1,k)
        dudz(i,1,k) = rho(i,1,k)*dwdy(i,1,k)*dwdz(i,1,k)
     ENDDO
  ENDDO

  Rxy = AVG1V3D(imax, jmax, kmax, i1, dudx)/rR
  Rxz = AVG1V3D(imax, jmax, kmax, i1, dudy)/rR
  Ryz = AVG1V3D(imax, jmax, kmax, i1, dudz)/rR

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 6')
#endif

! ##################################################################
! #           Pressure Dissipation                                 #
! ##################################################################

! #############################################################
! # Temporary Array Storage
! #
! # dudx = unused
! # dudy = P_prime
! # dudz = Rho_prime
! # dvdx = d T / d x
! # dvdy = d T / d y
! # dvdz = d T / d z
! # dwdx = d P_prime / d x
! # dwdy = d P_prime / d y
! # dwdz = d P_prime / d z
! #
! #############################################################

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudy(i,1,k) = p(i,1,k)   - rP
        dudz(i,1,k) = rho(i,1,k) - rR
     ENDDO
  ENDDO

!  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
!       dx, T, dvdx, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
!       dy, T, dvdy, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
!       dz, T, dvdz, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, dudy, dwdx, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, dudy, dwdy, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, dudy, dwdz, i0, i0, wrk1d, wrk2d, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) = (dwdx(i,1,k)**2&
             + dwdy(i,1,k)**2&
             + dwdz(i,1,k)**2)/rho(i,1,k)
     ENDDO
  ENDDO
  Ep2 = AVG1V3D(imax, jmax, kmax, i1, dudx)*ccnd2

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) =  (v(i,1,k)-rV)*dudy(i,1,k)**2&
             - ccnd2*(dudy(i,1,k)*dvdy(i,1,k))
     ENDDO
  ENDDO
  Tp2 = AVG1V3D(imax, jmax, kmax, i1, dudx)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dudx(i,1,k) = dvdx(i,1,k)*dwdx(i,1,k) &
             + dvdy(i,1,k)*dwdy(i,1,k)&
             + dvdz(i,1,k)*dwdz(i,1,k)
     ENDDO
  ENDDO

  Sp2 = Ep2 - ccnd2*AVG1V3D(imax, jmax, kmax, i1, dudx)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 7')
#endif

! ##################################################################
! #  Entropy Variation                                             #
! ##################################################################

  DO k = 1,kmax
     DO i = 1,imax*jmax
        IF ( p(i,1,k) .GT. C_0_R .AND. &
             rho(i,1,k) .GT. C_0_R ) THEN
           wrk3d(i,1,k) = log(p(i,1,k)/pbg%mean)&
                - gama0 * log(rho(i,1,k)/rbg%mean)
        ELSE
           wrk3d(i,1,k) = C_BIG_R
        ENDIF
     ENDDO
  ENDDO
  entropy = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        IF ( p(i,1,k) .GT. C_0_R .AND. &
             rho(i,1,k) .GT. C_0_R ) THEN
           wrk3d(i,1,k) = (log(p(i,1,k)/pbg%mean)&
                - gama0 * log(rho(i,1,k)/rbg%mean) &
                - entropy)**2
        ELSE
           wrk3d(i,1,k) = C_BIG_R
        ENDIF
     ENDDO
  ENDDO
  entropy2 = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 13')
#endif

! ##################################################################
! #  Acoustic and Entropy Density and Temperature fluctuations     #
! ##################################################################

  c2 = gama0*rP/rR

  DO k = 1,kmax
     DO i = 1,imax*jmax
        rho_prime = rho(i,1,k) - rR
        p_prime = p(i,1,k) - rP
!        T_prime = T(i,1,k) - rT
        rho_ac = p_prime/c2
        rho_en = rho_prime - rho_ac
        T_ac = rT*(p_prime/rP-rho_ac/rR)
        T_en = T_prime - T_ac

        dudx(i,1,k) = rho_ac*rho_ac
        dudy(i,1,k) = rho_en*rho_en
        dudz(i,1,k) = T_ac*T_ac
        dvdx(i,1,k) = T_en*T_en
     ENDDO
  ENDDO

  rho_ac = AVG1V3D(imax, jmax, kmax, i1, dudx)
  rho_en = AVG1V3D(imax, jmax, kmax, i1, dudy)
  T_ac = AVG1V3D(imax, jmax, kmax, i1, dudz)
  T_en = AVG1V3D(imax, jmax, kmax, i1, dvdx)

  IF ( rP .GT. C_0_R .and. rR .GT. C_0_R ) THEN
     M_t = SQRT( c23*Kin*rR/rP/gama0 )
  ELSE
     M_t = C_BIG_R
  ENDIF

! ##################################################################
! #           Reynolds Averages Derivatives Fluctuations           #
! ##################################################################

! #############################################################
! # Temporary Array Storage
! #
! # dudx = d U / d x
! # dudy = d U / d y
! # dudz = d U / d z
! # dvdx = d V / d x
! # dvdy = d V / d y
! # dvdz = d V / d z
! # dwdx = d W / d x
! # dwdy = d W / d y
! # dwdz = d W / d z
! #
! #############################################################

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, dudx, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, dudy, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, dudz, i0, i0, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, dvdx, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, dvdy, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, dvdz, i0, i0, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, dwdx, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, dwdy, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, dwdz, i0, i0, wrk1d, wrk2d, wrk3d)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 8')
#endif

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dil = vis(i,1,k)*(dudx(i,1,k)&
             +dvdy(i,1,k)+dwdz(i,1,k))*c23
        wrk3d(i,1,k) = vis(i,1,k)*C_2_R*dudx(i,1,k)-dil
     ENDDO
  ENDDO
  Tau_xx = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc
  DO k=1,kmax
     DO i=1,imax*jmax
        dil = vis(i,1,k)*(dudx(i,1,k)&
             +dvdy(i,1,k)+dwdz(i,1,k))*c23
        wrk3d(i,1,k) = vis(i,1,k)*C_2_R*dvdy(i,1,k)-dil
     ENDDO
  ENDDO
  Tau_yy = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc
  DO k=1,kmax
     DO i=1,imax*jmax
        dil = vis(i,1,k)*(dudx(i,1,k)&
             +dvdy(i,1,k)+dwdz(i,1,k))*c23
        wrk3d(i,1,k) = vis(i,1,k)*C_2_R*dwdz(i,1,k)-dil
     ENDDO
  ENDDO
  Tau_zz = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc
  DO k=1,kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = vis(i,1,k)*(dudy(i,1,k)+dvdx(i,1,k))
     ENDDO
  ENDDO
  Tau_xy = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc
  DO k=1,kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = vis(i,1,k)*(dudz(i,1,k)+dwdx(i,1,k))
     ENDDO
  ENDDO
  Tau_xz = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc
  DO k=1,kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = vis(i,1,k)*(dvdz(i,1,k)+dwdy(i,1,k))
     ENDDO
  ENDDO
  Tau_yz = AVG1V3D(imax, jmax, kmax, i1, wrk3d) * visc

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 10')
#endif

! ####################################
! # Calculate Taylor MicroScales     #
! ####################################

  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dudx(i,1,k)**2
     ENDDO
  ENDDO
  lambda_x = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dudx(i,1,k)**3
     ENDDO
  ENDDO
  skew_x = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  IF ( lambda_x .EQ. C_0_R ) THEN
     lambda_x = C_BIG_R
     skew_x = C_BIG_R
  ELSE
     skew_x = skew_x / lambda_x**C_1_5_R
     lambda_x = SQRT(Rxx/lambda_x)
  ENDIF

  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dvdy(i,1,k)**2
     ENDDO
  ENDDO
  lambda_y = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dvdy(i,1,k)**3
     ENDDO
  ENDDO
  skew_y = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  IF ( lambda_y .EQ. C_0_R ) THEN
     lambda_y = C_BIG_R
     skew_y = C_BIG_R
  ELSE
     skew_y = skew_y / lambda_y**C_1_5_R
     lambda_y = SQRT(Ryy/lambda_y)
  ENDIF

  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dwdz(i,1,k)**2
     ENDDO
  ENDDO
  lambda_z = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = dwdz(i,1,k)**3
     ENDDO
  ENDDO
  skew_z = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  IF ( lambda_z .EQ. C_0_R ) THEN
     lambda_z = C_BIG_R
     skew_z = C_BIG_R
  ELSE
     skew_z = skew_z / lambda_z**C_1_5_R
     lambda_z = SQRT(Rzz/lambda_z)
  ENDIF

  DO k=1, kmax
     DO i=1,imax*jmax
        wrk3d(i,1,k) = (rho(i,1,k)-rR)*(p(i,1,k)-rP)
     ENDDO
  ENDDO
  rRP = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  DO k=1, kmax
     DO i=1,imax*jmax
!        wrk3d(i,1,k) = (rho(i,1,k)-rR)*(T(i,1,k)-rT)
     ENDDO
  ENDDO
  rRT = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  IF ( rR2 .GT. C_0_R .AND. &
       rP2 .GT. C_0_R ) THEN
     rRP = rRP/sqrt(rR2*rP2)
  ELSE
     rRP = C_2_R
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 14')
#endif

! ##################################################################
! #  Dissipation Terms                                             #
! ##################################################################

! Epsilon
  DO k = 1,kmax
     DO i = 1,imax*jmax
        vs = visc*vis(i,1,k)
        dil = (dudx(i,1,k)+dvdy(i,1,k)+dwdz(i,1,k))*c23
        tau11 = vs*(C_2_R*dudx(i,1,k)-dil)-Tau_xx
        tau22 = vs*(C_2_R*dvdy(i,1,k)-dil)-Tau_yy
        tau33 = vs*(C_2_R*dwdz(i,1,k)-dil)-Tau_zz
        tau12 = vs*(dudy(i,1,k)+dvdx(i,1,k))-Tau_xy
        tau13 = vs*(dudz(i,1,k)+dwdx(i,1,k))-Tau_xz
        tau23 = vs*(dvdz(i,1,k)+dwdy(i,1,k))-Tau_yz
        wrk3d(i,1,k) = tau11*dudx(i,1,k)+tau12*dudy(i,1,k)&
             +tau13*dudz(i,1,k) +&
             tau12*dvdx(i,1,k)+tau22*dvdy(i,1,k)&
             +tau23*dvdz(i,1,k) +&
             tau13*dwdx(i,1,k)+tau23*dwdy(i,1,k)&
             +tau33*dwdz(i,1,k)
     ENDDO
  ENDDO
  Eps = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  Eps = Eps / rR

! Exy
  DO k = 1,kmax
     DO i = 1,imax*jmax
        vs = visc*vis(i,1,k)
        dil = (dudx(i,1,k)+dvdy(i,1,k)+dwdz(i,1,k))*c23
        tau11 = vs*(C_2_R*dudx(i,1,k)-dil)-Tau_xx
        tau22 = vs*(C_2_R*dvdy(i,1,k)-dil)-Tau_yy
        tau12 = vs*(dudy(i,1,k)+dvdx(i,1,k))-Tau_xy
        tau13 = vs*(dudz(i,1,k)+dwdx(i,1,k))-Tau_xz
        tau23 = vs*(dvdz(i,1,k)+dwdy(i,1,k))-Tau_yz
        wrk3d(i,1,k) = tau11*dvdx(i,1,k)+tau12*dvdy(i,1,k)&
             +tau13*dvdz(i,1,k)+tau12*dudx(i,1,k)&
             +tau22*dudy(i,1,k)+tau23*dudz(i,1,k)
     ENDDO
  ENDDO
  Exy = AVG1V3D(imax, jmax, kmax, i1, wrk3d)
  Exy = Exy / rR

!
! Mean Viscous Dissipation
! 
  DO k = 1,kmax
     DO i = 1,imax*jmax
        wrk3d(i,1,k) = vis(i,1,k)*(dudx(i,1,k)**2&
             +dvdy(i,1,k)**2+dwdz(i,1,k)**2&
             +C_05_R*((dudy(i,1,k)+dvdx(i,1,k))**2&
             +(dudz(i,1,k)+dwdx(i,1,k))**2&
             +(dvdz(i,1,k)+dwdy(i,1,k))**2) &
             -((dudx(i,1,k)+dvdy(i,1,k)&
             +dwdz(i,1,k))**2)/C_3_R)
     ENDDO
  ENDDO
  Phi = AVG1V3D(imax, jmax, kmax, i1, wrk3d) &
       * C_2_R * visc
! 
! Viscous Dissipation * Pressure Fluctuations
!
  DO k = 1,kmax
     DO i = 1,imax*jmax
        wrk3d(i,1,k) = vis(i,1,k)*(dudx(i,1,k)**2&
             +dvdy(i,1,k)**2+dwdz(i,1,k)**2&
             +C_05_R*((dudy(i,1,k)+dvdx(i,1,k))**2&
             +(dudz(i,1,k)+dwdx(i,1,k))**2&
             +(dvdz(i,1,k)+dwdy(i,1,k))**2) &
             -((dudx(i,1,k)+dvdy(i,1,k)&
             +dwdz(i,1,k))**2)/C_3_R)
        wrk3d(i,1,k) = (p(i,1,k)-rP)*wrk3d(i,1,k)
     ENDDO
  ENDDO

  PPhi = AVG1V3D(imax, jmax, kmax, i1, wrk3d) &
       * C_2_R * visc

  IF ( Eps .GT. C_0_R ) THEN

     DO k=1, kmax
        DO i=1, imax*jmax
           wrk3d(i,1,k) = visc*vis(i,1,k)/rho(i,1,k)
        ENDDO
     ENDDO

     lk = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

     relambda = C_2_R*Kin*SQRT(C_5_R/(lk*Eps))

     lk = (lk**3/Eps)**C_025_R
  ELSE
     lk = C_BIG_R
     relambda = C_0_R
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 15')
#endif

! ##################################################################
! #  Pressure Strain Terms                                         #
! ##################################################################

  DO k = 1,kmax
     DO i = 1,imax*jmax
        dil = dudx(i,1,k) + dvdy(i,1,k) + dwdz(i,1,k)
        dil = dil*dil
     ENDDO
  ENDDO
  DV2 = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        pp = p(i,1,k)-rP
        dil = dudx(i,1,k) + dvdy(i,1,k) + dwdz(i,1,k)
        wrk3d(i,1,k) = pp*dil
     ENDDO
  ENDDO
  Pi = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        pp = p(i,1,k)-rP
        wrk3d(i,1,k) = pp*(dudy(i,1,k)+dvdx(i,1,k))
     ENDDO
  ENDDO
  PIxy = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
        pp = p(i,1,k)-rP
        dil = dudx(i,1,k) + dvdy(i,1,k) + dwdz(i,1,k)
        wrk3d(i,1,k) = pp*pp*dil
     ENDDO
  ENDDO
  Pip2 = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 16')
#endif

! ##################################################################
! #  Temperature Variation                                         #
! ##################################################################

  DO k = 1,kmax
     DO i = 1,imax*jmax
!        wrk3d(i,1,k) = (T(i,1,k)- rT)**2
     ENDDO
  ENDDO
  rT2 = AVG1V3D(imax, jmax, kmax, i1, wrk3d)

  DO k = 1,kmax
     DO i = 1,imax*jmax
!        wrk3d(i,1,k) = rR*(T(i,1,k) - fT)**2
     ENDDO
  ENDDO
  fT2 = AVG1V3D(imax, jmax, kmax, i1, wrk3d)/rR

  IF ( rR2 .GT. C_0_R .AND. &
       rT2 .GT. C_0_R ) THEN
     rRT = rRT/sqrt(rR2*rT2)
  ELSE
     rRT = C_2_R
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_ISOTROPIC: Section 17')
#endif

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

     Rxy_t = - Exy + PIxy / rR

     Kin_t = - Eps + PI / rR

     WRITE(23,1024) 1, 1, &
          lk, relambda, lambda_x, lambda_y, lambda_z,&
          rR, rU, rV,&
          rW, rP, rT, &
          fU, fV, &
          fW, fT,&
          rP2, rR2, rT2, fT2,&
          Rxx, Ryy, Rzz, &
          Rxy, Rxz, Ryz,&
          Rxy_t, Exy, PIxy, &
          Kin, Kin_t, Eps, Pi, &
          DV2, &
          Tau_xx, Tau_yy, Tau_zz, &
          Tau_xy, Tau_xz, Tau_yz, &
          Phi, &
          rRP, rRT, &
          Pp2, Ep2,&
          Pip2, PPhi, Sp2,&
          skew_x, skew_y, skew_z, &
          entropy, entropy2, rho_ac, rho_en, &
          T_ac, T_en, M_t

1024 FORMAT(I3,1X,I3,57(1X,E12.5E3))

#ifdef USE_MPI
  ENDIF
#endif

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CLOSE(23)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE AVG_FLOW_TEMPORAL_ISOTROPIC


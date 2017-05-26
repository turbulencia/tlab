#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# The calculation of the normal direction is reapeted for different
!# groupd of variables to save memory
!#
!########################################################################
SUBROUTINE SL_NORMAL_VORTICITY(isl, ith, iavg, nmax, istep, kstep, nfield, itxc_size, &
     threshold, ibuffer_npy, u,v,w,p,z1, a, sl, profiles, txc, mean, wrk1d,wrk2d,wrk3d)
  
  USE DNS_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

#define L_NFIELDS_MAX 13

  TINTEGER isl, ith, nmax, istep, kstep, nfield, itxc_size, iavg, ibuffer_npy
  TREAL threshold
  TREAL u(*), v(*), w(*), p(*), z1(*), a(*), sl(imax,kmax)
  TREAL profiles(L_NFIELDS_MAX,nmax,imax/istep,kmax/kstep)
  TREAL mean(L_NFIELDS_MAX,nmax,2)
  TREAL txc(imax*jmax*kmax,7)
  TREAL wrk1d(*), wrk2d(3,imax,kmax), wrk3d(*)

! -------------------------------------------------------------------
  TREAL vmin, vmax, vmean, AVG_IK, diff, normal_factor, dn_u, norm
  TINTEGER ij, i, k, n, ifield, nfield_loc, ipfield, jmin_loc, jmax_loc
  CHARACTER*32 fname
#ifdef USE_MPI
  TINTEGER ioffset, ip
  INTEGER mpio_ip, mpio_locsize
  INTEGER status(MPI_STATUS_SIZE)
#endif
! ###################################################################
  jmin_loc = MAX(1,2*ibuffer_npy)
  jmax_loc = MIN(jmax,jmax - 2*ibuffer_npy +1)

  IF ( nfield .LT. L_NFIELDS_MAX ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_NORMAL_VORTICITY. Profiles array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ELSE
     nfield = L_NFIELDS_MAX
  ENDIF
  IF ( itxc_size .LT. imax*jmax*kmax*7 ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_NORMAL_VORTICITY. Txc array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! Calculate vorticiy field w_iw_i
  CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, a, txc(1,1),txc(1,2), wrk2d,wrk3d)

! -------------------------------------------------------------------
! Calculate boundaries
! -------------------------------------------------------------------
! threshold w.r.t w_max, therefore threshold^2 w.r.t. w^2_max
  IF ( ith .EQ. 1 ) THEN
     CALL MINMAX(imax,jmax,kmax, a, vmin,vmax)
     vmin = threshold*threshold*vmax
! threshold w.r.t w_mean, therefore threshold^2 w.r.t. w^2_mean
  ELSE IF ( ith .EQ. 2 ) THEN
     ij = jmax/2
     vmean = AVG_IK(imax,jmax,kmax, ij, a, dx, dz, area)
     vmin = threshold*threshold*vmean
  ENDIF
! upper or lower depending on flag isl
  IF ( isl .EQ. 1 ) THEN
     CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, y, a, txc(1,1), sl, wrk2d)
  ELSE IF ( isl .EQ. 2 ) THEN
     CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, y, a, txc(1,1), sl, wrk2d)
  ENDIF

! grid size factor for the normal direction
  normal_factor = C_05_R

! ###################################################################
! Normal analysis:
! txc1 ....: vorticity w_i w_i
! txc2 ....: scalar gradient G_i G_i
! txc3 ....: rate-of-strain 2 s_ij s_ij
! ###################################################################
  ipfield    = 1
  nfield_loc = 3

  CALL FI_STRAIN(imax,jmax,kmax, u,v,w, txc(1,3),txc(1,1),txc(1,2), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,3) = C_2_R*txc(ij,3)
  ENDDO
  CALL FI_GRADIENT(imax,jmax,kmax, z1,txc(1,2), txc(1,1), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,1) = a(ij)
  ENDDO

! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, txc(1,4), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, txc(1,5), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, txc(1,6), wrk3d, wrk2d,wrk3d)

  CALL SL_NORMAL_SAMPLE&
       (imax,jmax,kmax, kmax_total, nmax, istep, kstep, nfield_loc, nfield, &
       g(1)%scale, g(3)%scale, normal_factor, &
       g(1)%nodes,g(2)%nodes,g(3)%nodes, sl, txc, profiles(ipfield,1,1,1), txc(1,4), txc(1,5), txc(1,6))

! ###################################################################
! Normal analysis:
! txc1 ....: first invariant P
! txc2 ....: second invariant Q
! txc3 ....: third invariant R
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 3

  CALL FI_INVARIANT_R(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4),txc(1,5),txc(1,6),txc(1,1),txc(1,2), wrk2d,wrk3d)
  CALL FI_INVARIANT_Q(imax,jmax,kmax, u,v,w, txc(1,2), txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
  CALL FI_INVARIANT_P(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,4), wrk2d,wrk3d)

! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, txc(1,4), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, txc(1,5), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, txc(1,6), wrk3d, wrk2d,wrk3d)

  CALL SL_NORMAL_SAMPLE&
       (imax,jmax,kmax, kmax_total, nmax, istep, kstep, nfield_loc, nfield, &
       g(1)%scale, g(3)%scale, normal_factor,&
       g(1)%nodes,g(2)%nodes,g(3)%nodes, sl, txc, profiles(ipfield,1,1,1), txc(1,4), txc(1,5), txc(1,6))

! ###################################################################
! Normal analysis:
! txc1 ....: vorticity production
! txc2 ....: vorticity diffusion
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 2

  CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, u,v,w, txc(1,1),&
       txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
  CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, u,v,w, txc(1,2),&
       txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,2) = txc(ij,2)*visc
  ENDDO

! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, txc(1,4), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, txc(1,5), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, txc(1,6), wrk3d, wrk2d,wrk3d)

  CALL SL_NORMAL_SAMPLE&
       (imax,jmax,kmax, kmax_total, nmax, istep, kstep, nfield_loc, nfield, &
       g(1)%scale, g(3)%scale, normal_factor,&
       g(1)%nodes,g(2)%nodes,g(3)%nodes, sl, txc, profiles(ipfield,1,1,1), txc(1,4), txc(1,5), txc(1,6))

! ###################################################################
! Normal analysis:
! txc1 ....: scalar gradient production
! txc2 ....: scalar gradient diffusion
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 2

  CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, z1, u,v,w, &
       txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
  CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, z1, &
       txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
  diff = visc/schmidt(inb_scal)
  DO ij = 1,imax*jmax*kmax
     txc(ij,2) = txc(ij,2)*diff
  ENDDO

! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, txc(1,4), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, txc(1,5), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, txc(1,6), wrk3d, wrk2d,wrk3d)

  CALL SL_NORMAL_SAMPLE&
       (imax,jmax,kmax, kmax_total, nmax, istep, kstep, nfield_loc, nfield, &
       g(1)%scale, g(3)%scale, normal_factor, &
       g(1)%nodes,g(2)%nodes,g(3)%nodes, sl, txc, profiles(ipfield,1,1,1), txc(1,4), txc(1,5), txc(1,6))

! ###################################################################
! Normal analysis:
! txc1 ....: strain production
! txc2 ....: strain diffusion
! txc3 ....: strain-pressure
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 3

  CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, u,v,w, &
       txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,1) = txc(ij,1)*C_2_R
  ENDDO
  CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, u,v,w, &
       txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,2) = txc(ij,2)*visc*C_2_R
  ENDDO
  CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, u,v,w,p, &
       txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,3) = txc(ij,3)*C_2_R
  ENDDO

! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, txc(1,4), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, txc(1,5), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, txc(1,6), wrk3d, wrk2d,wrk3d)

  CALL SL_NORMAL_SAMPLE&
       (imax,jmax,kmax, kmax_total, nmax, istep, kstep, nfield_loc, nfield, &
       g(1)%scale, g(3)%scale, normal_factor, &
       g(1)%nodes,g(2)%nodes,g(3)%nodes, sl, txc, profiles(ipfield,1,1,1), txc(1,4), txc(1,5), txc(1,6))

! ###################################################################
! Output averages
! This option is note developed in parallel
! ###################################################################
  IF ( iavg .EQ. 1 ) THEN

     DO n = 1,nfield*nmax
! mean
        mean(n,1,1) = C_0_R
        DO k = 1,kmax/kstep
           DO i = 1,imax/istep
              mean(n,1,1) = mean(n,1,1) + profiles(n,1,i,k)
           ENDDO
        ENDDO
        mean(n,1,1) = mean(n,1,1)/M_REAL(imax/istep*kmax/kstep)

! standard deviation
        mean(n,1,2) = C_0_R
        DO k = 1,kmax/kstep
           DO i = 1,imax/istep
              mean(n,1,2) = mean(n,1,2) + profiles(n,1,i,k)*profiles(n,1,i,k)
           ENDDO
        ENDDO
        mean(n,1,2) = mean(n,1,2)/M_REAL(imax/istep*kmax/kstep) - mean(n,1,1)*mean(n,1,1)
        mean(n,1,2) = SQRT(mean(n,1,2))

     ENDDO

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
! mean value of grid spacing
     dn_u = (x(2)-x(1) + y(2)-y(1) + z(2)-z(1))/C_3_R*normal_factor

     WRITE(fname,*) itime; fname = 'avgSl'//TRIM(ADJUSTL(fname))

     OPEN(unit=21,file=fname)
     IF ( ith .EQ. 1 ) THEN
        WRITE(21,'(A12,E14.7E3)') '# w^2_max = ', vmax
     ELSE IF ( ith .EQ. 2 ) THEN
        WRITE(21,'(A12,E14.7E3)') '# w^2_avg = ', vmean
     ENDIF
     WRITE(21,'(A14,E14.7E3)') '# Threshold = ', vmin
     IF ( isl .EQ. 1 ) THEN
        WRITE(21,'(A)') '# Upper envelope surface '
     ELSE IF ( isl .EQ. 2 ) THEN
        WRITE(21,'(A)') '# Lower envelope surface '
     ENDIF
     WRITE(21,'(A8,E14.7E3)') 'RTIME = ', rtime
     WRITE(21,'(A)') 'GROUP = Mean '&
          //'rW2 rG2 r2S2 rP rQ rR rP_W rD_W rP_G rD_G r2P_S r2D_S r2SijPij'
     WRITE(21,'(A)') 'GROUP = Sigma '&
          //'sW2 sG2 s2S2 sP sQ sR sP_W sD_W sP_G sD_G s2P_S s2D_S s2SijPij'
     WRITE(21,'(A)') 'I J N '&
          //'rW2 rG2 r2S2 rP rQ rR rP_W rD_W rP_G rD_G r2P_S r2D_S r2SijPij '&
          //'sW2 sG2 s2S2 sP sQ sR sP_W sD_W sP_G sD_G s2P_S s2D_S s2SijPij'

     DO n = 1,nmax
        WRITE(21,1040) i1, i1, &
             M_REAL(n-1-nmax/2)*dn_u, &
             (mean(ifield,n,1),ifield=1,nfield),&
             (mean(ifield,n,2),ifield=1,nfield)
     ENDDO

! ###################################################################
! Output profiles
! ###################################################################
  ELSE
! sample the normal vector into array wrk2d for output file
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i3, i3, y, sl, txc(1,4), wrk2d)

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
#ifdef USE_MPI
     CALL DNS_MPI_TAGUPDT

     IF ( ims_pro .EQ. 0 ) THEN
#endif

! mean value of grid spacing
        dn_u = (x(2)-x(1) + y(2)-y(1) + z(2)-z(1))/C_3_R*normal_factor

        WRITE(str,*) itime; fname = 'slw'//TRIM(ADJUSTL(str))

        OPEN(unit=21,file=fname)
        IF ( ith .EQ. 1 ) THEN
           WRITE(21,'(A12,E14.7E3)') '# w^2_max = ', vmax
        ELSE IF ( ith .EQ. 2 ) THEN
           WRITE(21,'(A12,E14.7E3)') '# w^2_avg = ', vmean
        ENDIF
        WRITE(21,'(A14,E14.7E3)') '# Threshold = ', vmin
        IF ( isl .EQ. 1 ) THEN
           WRITE(21,'(A)') '# Upper envelope surface '
        ELSE IF ( isl .EQ. 2 ) THEN
           WRITE(21,'(A)') '# Lower envelope surface '
        ENDIF
        WRITE(21,'(A8,E14.7E3)') 'RTIME = ', rtime
        WRITE(21,'(A)') 'I J N W2 G2 2S2 P Q R '&
             //'P_W D_W P_G D_G 2P_S 2D_S 2SijPij Px Py Pz Nx Ny Nz'

        DO k = 1,kmax/kstep
           DO i = 1,imax/istep
              norm = SQRT(wrk2d(1,i,k)**2 + wrk2d(2,i,k)**2 + wrk2d(3,i,k)**2)
              DO n = 1,nmax-1
                 WRITE(21,1020) istep*i, kstep*k, &
                      M_REAL(n-1-nmax/2)*dn_u, &
                      (profiles(ifield,n,i,k),ifield=1,nfield)
              ENDDO
              n = nmax
              WRITE(21,1030) istep*i, kstep*k, &
                   M_REAL(n-1-nmax/2)*dn_u, &
                   (profiles(ifield,n,i,k),ifield=1,nfield),&
                   x(istep*i), sl(i,k), z(kstep*k), &
                   -wrk2d(1,i,k)/norm, -wrk2d(2,i,k)/norm, -wrk2d(3,i,k)/norm
           ENDDO
        ENDDO

#ifdef USE_MPI
        ioffset = 0
        DO ip = 2,ims_npro
           mpio_ip = ip-1
           mpio_locsize = nfield*(imax/istep)*nmax*(kmax/kstep)
           CALL MPI_RECV(profiles, mpio_locsize, MPI_REAL8, mpio_ip, &
                ims_tag, MPI_COMM_WORLD, status, ims_err)

           ioffset = ioffset + kmax
           DO k = 1,kmax/kstep
              DO i = 1,imax/istep
                 norm = SQRT(wrk2d(1,i,k)**2 + wrk2d(2,i,k)**2 + wrk2d(3,i,k)**2)
                 DO n = 1,nmax
                    WRITE(21,1020) istep*i, kstep*k+ioffset, &
                         M_REAL(n-1-nmax/2)*dn_u, &
                         (profiles(ifield,n,i,k),ifield=1,nfield)
                 ENDDO
                 n = nmax
                 WRITE(21,1030) istep*i, kstep*k+ioffset, &
                      M_REAL(n-1-nmax/2)*dn_u, &
                      (profiles(ifield,n,i,k),ifield=1,nfield),&
                      x(istep*i), sl(i,k), z(kstep*k+ioffset), &
                      -wrk2d(1,i,k)/norm, -wrk2d(2,i,k)/norm, -wrk2d(3,i,k)/norm
              ENDDO
           ENDDO
        ENDDO
#endif

        CLOSE(21)

#ifdef USE_MPI
     ELSE
        mpio_locsize = nfield*(imax/istep)*nmax*(kmax/kstep)
        CALL MPI_SEND&
             (profiles, mpio_locsize, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
     ENDIF
#endif

  ENDIF

1020 FORMAT(I3,1X,I3,1X,E10.3E3,L_NFIELDS_MAX(1X,E10.3E3))
1030 FORMAT(I3,1X,I3,1X,E10.3E3,L_NFIELDS_MAX(1X,E10.3E3),6(1X,E10.3E3))
1040 FORMAT(I3,1X,I3,1X,E10.3E3,L_NFIELDS_MAX(1X,E10.3E3),L_NFIELDS_MAX(1X,E10.3E3))

  RETURN
END SUBROUTINE SL_NORMAL_VORTICITY

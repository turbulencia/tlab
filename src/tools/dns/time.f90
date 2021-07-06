#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!#
!# Runge-Kutta explicit 3th order from Williamson 1980
!# Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
!# Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
!#
!########################################################################
MODULE TIME

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_flow,inb_scal, inb_flow_array,inb_scal_array
  USE DNS_GLOBAL, ONLY : icalc_flow,icalc_scal,icalc_part, imode_eqns
  USE DNS_GLOBAL, ONLY : isize_particle, inb_part,inb_part_array
  USE DNS_GLOBAL, ONLY : rtime, itime
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : itransport, visc, prandtl, schmidt
  USE DNS_LOCAL,  ONLY : nitera_first, nitera_log, logs_data
  USE LAGRANGE_GLOBAL, ONLY : l_g, ilagrange
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
  SAVE
  PRIVATE

  TINTEGER, PUBLIC :: rkm_mode              ! Type of Runge-Kutta scheme
  TINTEGER, PUBLIC :: rkm_endstep           ! number of substeps
  TINTEGER, PUBLIC :: rkm_substep           ! substep counter

  TREAL,    PUBLIC :: dtime                 ! time step
  TREAL :: dte                   ! time step of each substep
  TREAL :: etime                 ! time at each substep

  TREAL,    PUBLIC ::  cfla, cfld, cflr     ! CFL numbers

  TREAL kdt(5), kco(4), ktime(5)            ! explicit scheme coefficients
  TREAL kex(3), kim(3)                      ! implicit scheme coefficients

  TREAL schmidtfactor, dx2i
  TINTEGER i,j,k, kdsp,idsp
  TREAL dummy

  PUBLIC :: TIME_INITIALIZE
  PUBLIC :: TIME_RUNGEKUTTA
  PUBLIC :: TIME_COURANT

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TIME_INITIALIZE()
    IMPLICIT NONE

    ! ###################################################################
    ! RK coefficients
    SELECT CASE ( rkm_mode )
    CASE( RKM_EXP3 )
      rkm_endstep = 3

      kdt(1:3)   = (/ C_1_R/C_3_R, C_15_R/C_16_R,  C_8_R/C_15_R /)
      ktime(1:3) = (/ C_0_R,       C_1_R/C_3_R,    C_3_R/C_4_R  /)
      kco(1:2)   = (/-C_5_R/C_9_R,-C_153_R/C_128_R /)

    CASE( RKM_EXP4 )
      rkm_endstep = 5

      kdt(1) = C_1432997174477_R/C_9575080441755_R
      kdt(2) = C_5161836677717_R/C_13612068292357_R
      kdt(3) = C_1720146321549_R/C_2090206949498_R
      kdt(4) = C_3134564353537_R/C_4481467310338_R
      kdt(5) = C_2277821191437_R/C_14882151754819_R

      ktime(1) = C_0_R
      ktime(2) = C_1432997174477_R/C_9575080441755_R
      ktime(3) = C_2526269341429_R/C_6820363962896_R
      ktime(4) = C_2006345519317_R/C_3224310063776_R
      ktime(5) = C_2802321613138_R/C_2924317926251_R

      kco(1) = -C_567301805773_R/C_1357537059087_R
      kco(2) = -C_2404267990393_R/C_2016746695238_R
      kco(3) = -C_3550918686646_R/C_2091501179385_R
      kco(4) = -C_1275806237668_R/C_842570457699_R

    CASE( RKM_IMP3_DIFFUSION)
      rkm_endstep = 3

      kdt(1:3)   = (/   C_8_R/C_15_R,   C_5_R/C_12_R,  C_3_R/C_4_R /)

      kim(1:3)   = (/ C_111_R/C_256_R,  C_1_R/C_2_R,   C_2_R/C_9_R /)
      kex(1:3)   = (/ C_145_R/C_256_R, -C_9_R/C_50_R , C_2_R/C_9_R/)
      kco(1:3)   = (/ C_0_R,          -C_17_R/C_25_R, -C_5_R/C_9_R /)
      ! TO DO - calculate ktime from coefficients  ktime
      ktime(1:3) = (/ C_0_R,  C_0_R, C_0_R /)

      ! Coefficients from Spalart, Moser, Rogers (1991)
      ! kim = beta/gamma
      ! kex = alpha/gamma
      ! kco = zeta/gamma
      !
      ! alpha = (/ 29./96.,   -3./40,    1./6. /)
      ! beta  = (/ 37./160.,   5./24.,   1./6. /)
      ! gamma = (/  8./15.,    5./12.,   3./4. /)
      ! zeta  = (/  0.,      -17./60.,  -5./12./)

    END SELECT

    ! ###################################################################
    ! maximum diffusivities for TIME_COURANT
    schmidtfactor = C_1_R
    dummy         = C_1_R/prandtl
    schmidtfactor = MAX(schmidtfactor, dummy)
    dummy         = C_1_R/MINVAL(schmidt(1:inb_scal))
    schmidtfactor = MAX(schmidtfactor, dummy)
    schmidtfactor = schmidtfactor *visc

    ! -------------------------------------------------------------------
    ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TIME_COURANT
#ifdef USE_MPI
    idsp = ims_offset_i
    kdsp = ims_offset_k
#else
    idsp = 0
    kdsp = 0
#endif

    dx2i = C_0_R
    IF ( g(3)%size > 1 ) THEN
      DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        dummy = g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4)
        dx2i = MAX(dx2i,dummy)
      END DO; END DO; END DO
    ELSE
      DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        dummy = g(1)%jac(i+idsp,4) + g(2)%jac(j,4)
        dx2i = MAX(dx2i,dummy)
      END DO; END DO; END DO
    END IF

    RETURN
  END SUBROUTINE TIME_INITIALIZE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TIME_RUNGEKUTTA(q,hq,s,hs, txc, wrk1d,wrk2d,wrk3d, l_q,l_hq,l_txc,l_comm)
    IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

    TREAL q(isize_field,inb_flow_array)
    TREAL s(isize_field,inb_scal_array)
    TREAL hq(isize_field,inb_flow)
    TREAL hs(isize_field,inb_scal)
    TREAL, DIMENSION(*)             :: txc, wrk1d,wrk2d,wrk3d

    TREAL l_q(isize_particle,inb_part_array)
    TREAL l_hq(isize_particle,inb_part)
    TREAL, DIMENSION(*)                :: l_comm, l_txc

    ! -------------------------------------------------------------------
    TINTEGER is, flag_control
    TREAL alpha

    TINTEGER srt,END,siz ! Variables for OpenMP Paritioning
#ifdef USE_PROFILE
    TINTEGER t_srt,t_end,t_dif,idummy,PROC_CYCLES,MAX_CYCLES
    CHARACTER*256 time_string
#endif

#ifdef USE_BLAS
    INTEGER ilen
#endif

    !########################################################################
#ifdef USE_BLAS
    ilen = isize_field
#endif

    ! -------------------------------------------------------------------
    ! Initialize arrays to zero for the explcit low-storage algorithm
    ! -------------------------------------------------------------------
    IF ( rkm_mode == RKM_EXP3 .OR. rkm_mode == RKM_EXP4 ) THEN
      IF ( icalc_flow == 1 ) hq = C_0_R
      IF ( icalc_scal == 1 ) hs = C_0_R
      IF ( icalc_part == 1 ) l_hq = C_0_R
    END IF

    !########################################################################
    ! Loop over the sub-stages
    !########################################################################
    DO rkm_substep = 1,rkm_endstep

      ! -------------------------------------------------------------------
      ! Update transported (or prognostic) variables q and s
      ! -------------------------------------------------------------------
      dte   =         dtime   *kdt(rkm_substep)
      etime = rtime + dtime *ktime(rkm_substep)

#ifdef USE_PROFILE
      CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)
#endif
      SELECT CASE ( imode_eqns )
      CASE( DNS_EQNS_INCOMPRESSIBLE,DNS_EQNS_ANELASTIC )
        IF    ( rkm_mode == RKM_EXP3 .OR. rkm_mode == RKM_EXP4 ) THEN
          CALL TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT(&
              dte,etime, q,hq,s,hs,txc, wrk1d,wrk2d,wrk3d, l_q, l_hq, l_txc, l_comm)

        ELSE
          CALL TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT(&
              dte,etime, kex(rkm_substep), kim(rkm_substep), kco(rkm_substep), q,hq,s,hs,txc, wrk1d,wrk2d,wrk3d)
        END IF

      CASE( DNS_EQNS_INTERNAL,DNS_EQNS_TOTAL )
        CALL TIME_SUBSTEP_COMPRESSIBLE(dte,etime, q,hq,s,hs, txc, wrk1d,wrk2d,wrk3d)

      END SELECT

      CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

      ! -------------------------------------------------------------------
      ! Control updated values
      ! -------------------------------------------------------------------
      flag_control = MOD(rkm_substep,rkm_endstep) + MOD(itime+1-nitera_first,nitera_log) ! Only when datalogs are written
      CALL DNS_CONTROL(flag_control, q,s, txc, wrk2d,wrk3d)
      IF ( INT(logs_data(1)) /= 0 ) RETURN ! Error detected

      IF ( icalc_part == 1 .AND. ilagrange == LAG_TYPE_BIL_CLOUD_4 ) THEN
        CALL PARTICLE_TIME_RESIDENCE(dtime, l_g%np, l_q)
        CALL PARTICLE_TIME_LIQUID_CLIPPING(s, l_q,l_txc, l_comm, wrk3d)
      END IF

      ! -------------------------------------------------------------------
      ! Update RHS hq and hs in the explicit low-storage algorithm
      ! -------------------------------------------------------------------
      IF ( ( rkm_mode == RKM_EXP3 .OR. rkm_mode == RKM_EXP4 ) .AND. &
          rkm_substep < rkm_endstep ) THEN

#ifdef USE_BLAS
        !$omp parallel default(shared) &
        !$omp private (ilen,srt,end,siz,alpha,is)
#else
        !$omp parallel default(shared) &
        !$omp private (i,   srt,end,siz,alpha,is)
#endif

        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)
#ifdef USE_BLAS
        ILEN = siz
#endif

        alpha = kco(rkm_substep)

        IF ( icalc_flow == 1 ) THEN
          DO is = 1,inb_flow
#ifdef USE_BLAS
            CALL DSCAL(ILEN, alpha, hq(srt,is), 1)
#else
            hq(srt:END,is) = alpha *hq(srt:END,is)
#endif
          END DO
        END IF

        IF ( icalc_scal == 1 ) THEN
          DO is = 1,inb_scal
#ifdef USE_BLAS
            CALL DSCAL(ILEN, alpha, hs(srt,is), 1)
#else
            hs(srt:END,is) = alpha *hs(srt:END,is)
#endif
          END DO
        END IF
        !$omp end parallel

        IF ( icalc_part == 1 ) THEN
          DO is = 1,inb_part
            l_hq(1:l_g%np,is) = alpha*l_hq(1:l_g%np,is)
          END DO
        END IF

      END IF

      ! -------------------------------------------------------------------
      ! Profiling data
      ! -------------------------------------------------------------------
#ifdef USE_PROFILE
      CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)
      idummy = t_end-t_srt

#ifdef USE_MPI
      CALL MPI_REDUCE(idummy,t_dif,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD)
      IF ( ims_pro == 0 ) THEN
        WRITE(time_string,999) ims_npro, ims_npro_i, ims_npro_k, rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro
999     FORMAT(I5.5,' (ims_npro_i X ims_npro_k:',I4.4,'x',I4.4,1x,') RK-Substep',I1,':', E13.5,'s')
        CALL IO_WRITE_ASCII(lfile, time_string)
      END IF
#else
      t_dif = idummy
      WRITE(time_string,999) rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro
999   FORMAT('RK-Substep',I1,':', E13.5,'s')
      CALL IO_WRITE_ASCII(lfile,time_string)
#endif

#endif

    END DO

    RETURN
  END SUBROUTINE TIME_RUNGEKUTTA

  !########################################################################
  !# DESCRIPTION
  !#
  !# Determine the variable time step is a positive cfl is given
  !# If negative cfl is prescribed, then constant dtime and this routine
  !# calculates CFL and diffustion numbers for log files
  !#
  !# The reacting case should be reviewed in the new formulation in terms
  !# of energy.
  !#
  !# The diffusion number is fixed in terms of the CFL, which is the input.
  !# This depends on the scheme used. From Lele (1992), page 32, we have that, if
  !# the sixth order tridiagonal scheme is used, then the maximum CFL number
  !# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994
  !# used here we have 3.36/1.989, about 1.69.
  !# This holds for periodic case. A safety margin leads to the common value of 1.2.
  !#
  !# If second order finite different operator is used, then the maximum
  !# diffusion number is 2.9/6.857, about 0.42.
  !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
  !# If the extension by Lamballais et al is used, then the maximum
  !# diffusion number is 2.9/pi^2, about 0.29.
  !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
  !#
  !# If twice the first order finite difference operator is used, then the
  !# maximum diffusion number is 2.9/1.989^2, about 0.73.
  !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
  !#
  !# In incompressible mode the arrays rho, p and vis are not used
  !#
  !########################################################################
  SUBROUTINE TIME_COURANT(q, wrk3d)

    USE THERMO_GLOBAL, ONLY : gama0
#ifdef CHEMISTRY
    USE CHEM_GLOBAL, ONLY : TGFM
#endif

    IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

    TREAL, INTENT(IN   ) :: q(imax,jmax,kmax,inb_flow_array)
    TREAL, INTENT(INOUT) :: wrk3d(imax,jmax,kmax)

    ! -------------------------------------------------------------------
    TINTEGER ipmax
    TREAL dt_loc
    TREAL pmax(3), dtc,dtd,dtr
#ifdef CHEMISTRY
    TREAL mtgfm, zmin, zmax, umin, umax
    TREAL vmin, vmax, wmin, wmax
#endif
#ifdef USE_MPI
    TREAL pmax_aux(3)
#endif

    ! ###################################################################
#ifdef USE_MPI
    idsp = ims_offset_i
    kdsp = ims_offset_k
#else
    idsp = 0
    kdsp = 0
#endif

    dtc = C_BIG_R   ! So that the minimum non-zero determines dt at the end
    dtd = C_BIG_R
    dtr = C_BIG_R

    ipmax = 0       ! Initialize counter of time restrictions

    ! ###################################################################
    ! CFL number condition
    ! ###################################################################
#define rho(i,j,k) q(i,j,k,5)
#define p(i,j,k)   q(i,j,k,6)

    ipmax = ipmax +1

    ! -------------------------------------------------------------------
    ! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
    ! -------------------------------------------------------------------
    SELECT CASE( imode_eqns )
    CASE( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
      IF ( g(3)%size > 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
          wrk3d(i,j,k) = ABS(q(i,j,k,1)) *g(1)%jac(i+idsp,3) &
              + ABS(q(i,j,k,2)) *g(2)%jac(j,3)      &
              + ABS(q(i,j,k,3)) *g(3)%jac(k+kdsp,3)
        END DO; END DO; END DO
      ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
          wrk3d(i,j,k) = ABS(q(i,j,k,1)) *g(1)%jac(i+idsp,3) &
              + ABS(q(i,j,k,2)) *g(2)%jac(j,3)
        END DO; END DO; END DO
      END IF

      ! -------------------------------------------------------------------
      ! Compressible: Calculate global maximum of (u+c)/dx + (v+c)/dy + (w+c)/dz
      ! -------------------------------------------------------------------
    CASE( DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL )
      wrk3d = SQRT(gama0*p(:,:,:)/rho(:,:,:)) ! sound speed; positiveness of p and rho checked in routine DNS_CONTROL
      IF ( g(3)%size > 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
          wrk3d(i,j,k) = (ABS(q(i,j,k,1))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
              + (ABS(q(i,j,k,2))+wrk3d(i,j,k)) *g(2)%jac(j,3)      &
              + (ABS(q(i,j,k,3))+wrk3d(i,j,k)) *g(3)%jac(k+kdsp,3)
        END DO; END DO; END DO
      ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
          wrk3d(i,j,k) = (ABS(q(i,j,k,1))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
              + (ABS(q(i,j,k,2))+wrk3d(i,j,k)) *g(2)%jac(j,3)
        END DO; END DO; END DO
      END IF

    END SELECT

    pmax(1) = MAXVAL(wrk3d)

    ! ###################################################################
    ! Diffusion number condition
    ! ###################################################################
#define vis(i,j,k) q(i,j,k,8)

    ipmax = ipmax +1

    ! -------------------------------------------------------------------
    ! Incompressible: Calculate global maximum of \nu*(1/dx^2 + 1/dy^2 + 1/dz^2)
    ! -------------------------------------------------------------------
    SELECT CASE( imode_eqns )
    CASE( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
      pmax(2) = schmidtfactor *dx2i

      ! -------------------------------------------------------------------
      ! Compressible: Calculate global maximum of \mu/rho*(1/dx^2 + 1/dy^2 + 1/dz^2)
      ! -------------------------------------------------------------------
    CASE( DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL )
      IF ( itransport == EQNS_TRANS_POWERLAW ) THEN
        IF ( g(3)%size > 1 ) THEN
          DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
            wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4) ) *vis(i,j,k) /rho(i,j,k)
          END DO; END DO; END DO
        ELSE
          DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
            wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) ) *vis(i,j,k) /rho(i,j,k)
          END DO; END DO; END DO
        END IF

      ELSE ! constant dynamic viscosity
        IF ( g(3)%size > 1 ) THEN
          DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
            wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4) ) /rho(i,j,k)
          END DO; END DO; END DO
        ELSE
          DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
            wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) ) /rho(i,j,k)
          END DO; END DO; END DO
        END IF

      END IF

      pmax(2) = schmidtfactor *MAXVAL(wrk3d)

    END SELECT


#ifdef CHEMISTRY
    ! ###################################################################
    ! Reacting time step control
    ! ###################################################################
    ipmax = ipmax +1
    pmax(ipmax) = C_0_R

    ! -------------------------------------------------------------------
    ! Infinitely fast
    ! -------------------------------------------------------------------
    IF ( ireactive == CHEM_INFINITE ) THEN
      CALL MINMAX(imax,jmax,kmax, q(1,1,1,1), umin,umax)
      CALL MINMAX(imax,jmax,kmax, q(1,1,1,2), vmin,vmax)
      CALL MINMAX(imax,jmax,kmax, q(1,1,1,3), wmin,wmax)
      CALL MINMAX(imax,jmax,kmax, s(1,1,1,inb_scal), zmin,zmax)

      umax = umax-umin
      vmax = vmax-vmin
      wmax = wmax-wmin

      umax = MAX(umax, vmax)
      umax = MAX(umax, wmax)

      mtgfm = TGFM*(zmax-zmin)/umax

      IF ( itransport == EQNS_TRANS_POWERLAW ) THEN
        DO k=1, kmax
          DO ij=1, imax*jmax
            index = ij-1
            wrk2d(ij,1) = C_1_R/(g(1)%jac(MOD(index,imax)+1+idsp,1)**3)
            wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(g(2)%jac((ij-1)/imax+1,1)**3)
            wrk2d(ij,1) = (mtgfm*vis(ij,1,k)/rho(ij,1,k))*wrk2d(ij,1)
          END DO

          IF ( g(3)%size > 1 ) THEN
            DO ij=1, imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm*vis(ij,1,k)/rho(ij,1,k))/(g(3)%jac(k+kdsp,1)**3)
            END DO
          END IF

          bmax = MAXVAL(wrk2d(1:nxy,1))
          IF (bmax >= pmax(ipmax)) pmax(ipmax) = bmax
        END DO

      ELSE
        DO k=1, kmax
          DO ij=1, imax*jmax
            index = ij-1
            wrk2d(ij,1) = C_1_R/(g(1)%jac(MOD(index,imax)+1+idsp,1)**3)
            wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(g(2)%jac((ij-1)/imax+1,1)**3)
            wrk2d(ij,1) = (mtgfm/rho(ij,1,k))*wrk2d(ij,1)
          END DO

          IF ( g(3)%size > 1 ) THEN
            DO ij=1, imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm/rho(ij,1,k))/(g(3)%jac(k+kdsp,1)**3)
            END DO
          END IF

          bmax = MAXVAL(wrk2d(1:nxy,1))
          IF (bmax >= pmax(ipmax)) pmax(ipmax) = bmax
        END DO

      END IF

      ! -------------------------------------------------------------------
      ! Finite rate
      ! -------------------------------------------------------------------
    ELSE IF ( ireactive == CHEM_FINITE ) THEN
      ! I obtain TGFM/(\rho*c^2) (not TGFM/(\rho*u^2)) as eigenvalue if I
      ! assume that omega depends only on \rho as \rho*S, S constant.
      IF ( g(3)%size > 1 ) THEN
        wrk3d = rho *( q(:,:,:,1) *q(:,:,:,1) +q(:,:,:,2) *q(:,:,:,2) +q(:,:,:,3) *q(:,:,:,3) )
      ELSE
        wrk3d = rho *( q(:,:,:,1) *q(:,:,:,1) +q(:,:,:,2) *q(:,:,:,2) )
      END IF

      bmax = MAXVAL(wrk3d)
      IF (bmax >= pmax(ipmax)) pmax(ipmax) = bmax
      pmax(ipmax) = TGFM/pmax(ipmax)

    END IF

#endif

    ! ###################################################################
    ! Final operations
    ! ###################################################################
#ifdef USE_MPI
    CALL MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
    pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

    IF ( pmax(1) > C_0_R ) dtc = cfla /pmax(1) ! Set time step for the given CFL number
    IF ( pmax(2) > C_0_R ) dtd = cfld /pmax(2) ! Set time step for the given diffusion number
#ifdef CHEMISTRY
    IF ( pmax(ipmax) > C_0_R ) dtr = cflr/pmax(ipmax)
#endif

    ! -------------------------------------------------------------------
    IF ( cfla > C_0_R ) THEN
      IF ( rkm_mode == RKM_EXP3 .OR. rkm_mode == RKM_EXP4 ) THEN
        dt_loc = MIN(dtc, dtd)
      ELSE
        dt_loc = dtc
      END IF
      dt_loc = MIN(dt_loc, dtr)

      dtime = dt_loc

    END IF

    ! Real CFL and diffusion numbers being used, for the logfile
    logs_data(2) = dtime * pmax(1)
    logs_data(3) = dtime * pmax(2)
#ifdef CHEMISTRY
    IF      ( ireactive == CHEM_INFINITE ) THEN; logs_data(4) = (dtime**2) *pmax(3)
    ELSE IF ( ireactive == CHEM_FINITE   ) THEN; logs_data(4) =  dtime     *pmax(3)
    END IF
#endif

    RETURN

  END SUBROUTINE TIME_COURANT

END MODULE TIME

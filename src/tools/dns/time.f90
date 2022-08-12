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
  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, inb_flow,inb_scal, inb_flow_array,inb_scal_array
  USE TLAB_VARS, ONLY : icalc_flow,icalc_scal,icalc_part, imode_eqns
  USE TLAB_VARS, ONLY : isize_particle, inb_part,inb_part_array
  USE TLAB_VARS, ONLY : rtime, itime
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : itransport, visc, prandtl, schmidt
  USE DNS_LOCAL,  ONLY : nitera_first, nitera_log, logs_data
  USE TLAB_PROCS
  USE LAGRANGE_VARS, ONLY : l_g, ilagrange
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE
  SAVE
  PRIVATE

  TINTEGER, PUBLIC :: rkm_mode              ! Type of Runge-Kutta scheme
  TINTEGER, PUBLIC :: rkm_endstep           ! number of substeps
  TINTEGER, PUBLIC :: rkm_substep           ! substep counter

  TREAL,    PUBLIC :: cfla, cfld, cflr      ! CFL numbers
  TREAL,    PUBLIC :: dtime                 ! time step
  TREAL dte                                 ! time step of each substep
  TREAL etime                               ! time at each substep

  TREAL kdt(5), kco(4), ktime(5)            ! explicit scheme coefficients
  TREAL kex(3), kim(3)                      ! implicit scheme coefficients

  TREAL schmidtfactor, dx2i
  TINTEGER i,j,k, kdsp,idsp, is
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
  SUBROUTINE TIME_RUNGEKUTTA()
    USE TLAB_ARRAYS
    USE LAGRANGE_ARRAYS
    USE DNS_ARRAYS
    IMPLICIT NONE

    ! -------------------------------------------------------------------
    TINTEGER flag_control
    TREAL alpha

    TINTEGER ij_srt,ij_end,ij_siz ! Variables for OpenMP Paritioning
#ifdef USE_PROFILE
    TINTEGER t_srt,t_end,t_dif,idummy,PROC_CYCLES,MAX_CYCLES
    CHARACTER*256 time_string
#endif

#ifdef USE_BLAS
    INTEGER ij_len
#endif

    !########################################################################
#ifdef USE_BLAS
    ij_len = isize_field
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
          CALL TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT()
        ELSE
          CALL TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT()
        END IF

      CASE( DNS_EQNS_INTERNAL,DNS_EQNS_TOTAL )
        CALL TIME_SUBSTEP_COMPRESSIBLE()

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
!$omp private (ij_len,ij_srt,ij_end,ij_siz,alpha,is)
#else
!$omp parallel default(shared) &
!$omp private (i,   ij_srt,ij_end,ij_siz,alpha,is)
#endif

        CALL DNS_OMP_PARTITION(isize_field,ij_srt,ij_end,ij_siz)
#ifdef USE_BLAS
        ij_len = ij_siz
#endif

        alpha = kco(rkm_substep)

        IF ( icalc_flow == 1 ) THEN
          DO is = 1,inb_flow
#ifdef USE_BLAS
            CALL DSCAL(ij_len, alpha, hq(ij_srt,is), 1)
#else
            hq(ij_srt:ij_end,is) = alpha *hq(ij_srt:ij_end,is)
#endif
          END DO
        END IF

        IF ( icalc_scal == 1 ) THEN
          DO is = 1,inb_scal
#ifdef USE_BLAS
            CALL DSCAL(ij_len, alpha, hs(ij_srt,is), 1)
#else
            hs(ij_srt:ij_end,is) = alpha *hs(ij_srt:ij_end,is)
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
        CALL TLAB_WRITE_ASCII(lfile, time_string)
      END IF
#else
      t_dif = idummy
      WRITE(time_string,999) rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro
999   FORMAT('RK-Substep',I1,':', E13.5,'s')
      CALL TLAB_WRITE_ASCII(lfile,time_string)
#endif

#endif

    END DO

    RETURN
  END SUBROUTINE TIME_RUNGEKUTTA

  !########################################################################
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

    USE THERMO_VARS, ONLY : gama0
#ifdef CHEMISTRY
    USE CHEM_GLOBAL, ONLY : TGFM
#endif

    IMPLICIT NONE

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

#undef rho
#undef p
#undef vis

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

  !########################################################################
  !#
  !# Branching among different formulations of the RHS.
  !#
  !# Be careful to define here the pointers and to enter the RHS routines
  !# with individual fields. Definition of pointer inside of RHS routines
  !# decreased performance considerably (at least in JUGENE)
  !#
  !########################################################################
  SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT()
    USE TLAB_VARS, ONLY : iadvection
    USE TLAB_ARRAYS
    USE DNS_LOCAL, ONLY : imode_rhs
    USE DNS_ARRAYS
    USE LAGRANGE_ARRAYS
    USE BOUNDARY_BUFFER
    IMPLICIT NONE

    ! -----------------------------------------------------------------------
    TINTEGER ij_srt,ij_end,ij_siz    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
    INTEGER ij_len
#endif

    ! #######################################################################
#ifdef USE_BLAS
    ij_len = isize_field
#endif

    SELECT CASE ( iadvection )
    CASE( EQNS_DIVERGENCE )
      CALL FI_SOURCES_FLOW(q,s, hq, txc(1,1), wrk1d,wrk2d,wrk3d)
      CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3(dte, q(1,1),q(1,2),q(1,3),hq(1,1),hq(1,2),hq(1,3), &
          q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)

      CALL FI_SOURCES_SCAL(s, hs, txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d)
      DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3(is, q(1,1),q(1,2),q(1,3), s(1,is),hs(1,is), &
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      END DO

    CASE( EQNS_SKEWSYMMETRIC )
      CALL FI_SOURCES_FLOW(q,s, hq, txc(1,1), wrk1d,wrk2d,wrk3d)
      CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2(dte, q(1,1),q(1,2),q(1,3),hq(1,1),hq(1,2),hq(1,3), &
          q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)

      CALL FI_SOURCES_SCAL(s, hs, txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d)
      DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2(is, q(1,1),q(1,2),q(1,3), s(1,is),hs(1,is), &
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
      END DO

    CASE( EQNS_CONVECTIVE )
      SELECT CASE ( imode_rhs )
      CASE( EQNS_RHS_SPLIT )
        CALL FI_SOURCES_FLOW(q,s, hq, txc(1,1), wrk1d,wrk2d,wrk3d)
        CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1(dte, q(1,1),q(1,2),q(1,3),hq(1,1),hq(1,2),hq(1,3), &
            q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)

        CALL FI_SOURCES_SCAL(s, hs, txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d)
        DO is = 1,inb_scal
          CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1(is, q(1,1),q(1,2),q(1,3), s(1,is),hs(1,is), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
        END DO

      CASE( EQNS_RHS_COMBINED )
        CALL FI_SOURCES_FLOW(q,s, hq, txc(1,1),          wrk1d,wrk2d,wrk3d)
        CALL FI_SOURCES_SCAL(  s, hs, txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d)
        CALL RHS_GLOBAL_INCOMPRESSIBLE_1(dte, q(1,1),q(1,2),q(1,3),hq(1,1),hq(1,2),hq(1,3), &
            q,hq, s,hs, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)

      CASE( EQNS_RHS_NONBLOCKING )
#ifdef USE_PSFFT
        CALL RHS_GLOBAL_INCOMPRESSIBLE_NBC(dte, &
            q(1,1),q(1,2),q(1,3),s(1,1),&
            txc(1,1), txc(1,2), &
            txc(1,3), txc(1,4), txc(1,5), txc(1,6), txc(1,7), txc(1,8),txc(1,9),txc(1,10), &
            txc(1,11),txc(1,12),txc(1,13),txc(1,14),&
            hq(1,1),hq(1,2),hq(1,3), hs(1,1), &
            wrk1d,wrk2d,wrk3d)
#else
        CALL TLAB_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT. Need compiling flag -DUSE_PSFFT.')
        CALL TLAB_STOP(DNS_ERROR_PSFFT)
#endif
      END SELECT
    END SELECT

    IF ( icalc_part == 1 ) THEN
      CALL RHS_PARTICLE_GLOBAL(q,s, txc, l_q,l_hq,l_txc,l_comm, wrk1d,wrk2d,wrk3d)
    END IF

    IF ( BuffType == DNS_BUFFER_RELAX .OR. BuffType == DNS_BUFFER_BOTH ) THEN
      CALL BOUNDARY_BUFFER_RELAX_SCAL(s,hs, q) ! Flow part needs to be taken into account in the pressure
    END IF

    ! #######################################################################
    ! Perform the time stepping for incompressible equations
    ! #######################################################################
#ifdef USE_OPENMP
#ifdef USE_BLAS
!$omp parallel default(shared) &
!$omp private (ij_len,is,ij_srt,ij_end,ij_siz)
#else
!$omp parallel default(shared) &
!$omp private (ij,  is,ij_srt,ij_end,ij_siz)
#endif
#endif

    CALL DNS_OMP_PARTITION(isize_field,ij_srt,ij_end,ij_siz)
#ifdef USE_BLAS
    ij_len = ij_siz
#endif
    DO is = 1,inb_flow

#ifdef USE_BLAS
      CALL DAXPY(ij_len, dte, hq(ij_srt,is), 1, q(ij_srt,is), 1)
#else
      q(ij_srt:ij_end,is) = q(ij_srt:ij_end,is) + dte *hq(ij_srt:ij_end,is)
#endif
    END DO

    DO is = 1,inb_scal
#ifdef BLAS
      CALL DAXPY(ij_len, dte, hs(ij_srt,is), 1, s(ij_srt,is), 1)
#else
      s(ij_srt:ij_end,is) = s(ij_srt:ij_end,is) + dte *hs(ij_srt:ij_end,is)
#endif
    END DO
#ifdef USE_OPENMP
!$omp end parallel
#endif

    ! ######################################################################
    ! Particle POSTION UPDATED and  SEND/RECV TO THE NEW PROCESSOR
    ! ######################################################################
    IF ( icalc_part == 1 ) THEN
      CALL PARTICLE_TIME_SUBSTEP(dte, l_q, l_hq, l_comm )
    END IF

    RETURN
  END SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT

  !########################################################################
  !########################################################################
  SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT()
    USE TLAB_ARRAYS
    USE DNS_ARRAYS
    IMPLICIT NONE

    ! ######################################################################
    IF      ( rkm_mode == RKM_IMP3_DIFFUSION ) THEN
      CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(&
          dte, kex(rkm_substep), kim(rkm_substep), kco(rkm_substep),  &
          q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk1d,wrk2d,wrk3d)
      ! pressure-correction algorithm; to be checked
      ! CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(&
      !      dte, kex,kim,kco,  &
      !      q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
      !      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), txc(1,8), &
      !      wrk1d,wrk2d,wrk3d)
    ELSE
      CALL TLAB_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT. Undeveloped formulation.')
      CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

    END IF

    RETURN
  END SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT

  !########################################################################
  !########################################################################
  SUBROUTINE TIME_SUBSTEP_COMPRESSIBLE()
    USE TLAB_VARS, ONLY : iadvection, idiffusion, iviscous, mach
    USE TLAB_ARRAYS
    USE DNS_ARRAYS
    USE THERMO_VARS, ONLY : gama0
    USE BOUNDARY_BUFFER
    USE BOUNDARY_BCS, ONLY : BcsDrift
    USE BOUNDARY_BCS_COMPRESSIBLE

    IMPLICIT NONE

#include "integers.h"

    ! -------------------------------------------------------------------
    TREAL rho_ratio, dt_rho_ratio, prefactor
    TREAL M2_max, dummy
    TINTEGER inb_scal_loc

    ! Pointers to existing allocated space
    TREAL, DIMENSION(:), POINTER :: u, v, w, e, rho, p, T, vis

    ! ###################################################################
    ! Define pointers
    u   => q(:,1)
    v   => q(:,2)
    w   => q(:,3)

    e   => q(:,4)
    rho => q(:,5)
    p   => q(:,6)
    T   => q(:,7)

    vis => q(:,8)

    ! ###################################################################
    ! Evaluate standard RHS of equations
    ! global formulation
    ! ###################################################################
    IF ( imode_eqns == DNS_EQNS_INTERNAL  .AND. &
        iadvection == EQNS_SKEWSYMMETRIC .AND. &
        iviscous   == EQNS_EXPLICIT      .AND. &
        idiffusion == EQNS_EXPLICIT            ) THEN
      CALL RHS_FLOW_GLOBAL_2(rho,u,v,w,p,e,T,s, hq(1,5),hq(1,1),hq(1,2),hq(1,3),hq(1,4),hs,&
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_2(is, rho,u,v,w,s,T, hs, hq(1,4),&
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      END DO

    ELSE
      ! ###################################################################
      ! Evaluate standard RHS of equations
      ! split formulation
      ! ###################################################################
      ! -------------------------------------------------------------------
      ! convective terms
      ! -------------------------------------------------------------------
      IF      ( iadvection == EQNS_DIVERGENCE    ) THEN
        CALL RHS_FLOW_EULER_DIVERGENCE(rho,u,v,w,p,e, hq(1,5),hq(1,1),hq(1,2),hq(1,3),hq(1,4),&
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        DO is = 1,inb_scal
          CALL RHS_SCAL_EULER_DIVERGENCE(rho,u,v,w,s(1,is), hs(1,is),&
              txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        END DO
      ELSE IF ( iadvection == EQNS_SKEWSYMMETRIC ) THEN
        CALL RHS_FLOW_EULER_SKEWSYMMETRIC(rho,u,v,w,p,e,s, hq(1,5),hq(1,1),hq(1,2),hq(1,3),hq(1,4),hs,&
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        DO is = 1,inb_scal
          CALL RHS_SCAL_EULER_SKEWSYMMETRIC(rho,u,v,w,s(1,is), hs(1,is),&
              txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        END DO
      END IF

      ! -------------------------------------------------------------------
      ! viscous terms
      ! -------------------------------------------------------------------
      IF ( itransport .NE. 1 ) THEN
        CALL TLAB_WRITE_ASCII(efile,'TIME_SUBSTEP_COMPRESSIBLE. Section requires to allocate array vis.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
      END IF

      IF      ( iviscous == EQNS_DIVERGENCE ) THEN
        CALL RHS_FLOW_VISCOUS_DIVERGENCE(vis, u,v,w,p, hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8),txc(1,9), wrk2d,wrk3d)
      ELSE IF ( iviscous == EQNS_EXPLICIT   ) THEN
        CALL RHS_FLOW_VISCOUS_EXPLICIT(vis, u,v,w,p, hq(1,1),hq(1,2),hq(1,3),hq(1,4), &
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
      END IF

      ! -------------------------------------------------------------------
      ! diffusion/conduction terms
      ! -------------------------------------------------------------------
      IF      ( idiffusion == EQNS_DIVERGENCE ) THEN
        ! diffusion transport of enthalpy is accumulated in txc and used then in RHS_FLOW_CONDUCTION
        txc(:,1) = C_0_R; txc(:,2) = C_0_R; txc(:,3) = C_0_R
        DO is = 1,inb_scal
          CALL RHS_SCAL_DIFFUSION_DIVERGENCE(is, vis, s, T, hs, &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        END DO
        CALL RHS_FLOW_CONDUCTION_DIVERGENCE(vis, s, T, hq(1,4), &
            txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
      ELSE IF ( idiffusion == EQNS_EXPLICIT   ) THEN
        DO is = 1,inb_scal
          CALL RHS_SCAL_DIFFUSION_EXPLICIT(is, vis, s, T, hs, hq(1,4), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        END DO
        CALL RHS_FLOW_CONDUCTION_EXPLICIT(vis, s, T, hq(1,4), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
      END IF

    END IF

    ! ###################################################################
    ! Evaluate chemical RHS
    ! If LES is on, each of these subroutines should have its
    ! counterpart in LES library
    ! NOT YET DEVELOPED FOR THE ENERGY FORMULATION !!!
    ! ###################################################################
#ifdef CHEMISTRY
    IF ( icalc_scal == 1 ) THEN
      IF      ( imech_type == CHEM_TYPE_PETERS1991  ) THEN
        CALL CHEM_PETERS1991(i0, rho, T, gama, s, hs, hq(1,4), txc)
      ELSE IF ( imech_type == CHEM_TYPE_PETERS1988  ) THEN
        CALL CHEM_PETERS1988(rho, T, gama, s, hs, hq(1,4))
      ELSE IF ( imech_type == CHEM_TYPE_UNIDECOMP   ) THEN
        CALL CHEM_UNIDECOMP(rho, T, gama, s, hs, hq(1,4))
      ELSE IF ( imech_type == CHEM_TYPE_BSZELDOVICH ) THEN
        CALL CHEM_BSZELDOVICH(rho, T, gama, s, hs, hq(1,4))
      ELSE IF ( imech_type == CHEM_TYPE_ONESTEP     ) THEN
        CALL CHEM_ONESTEP(rho, T, gama, s, hs, hq(1,4))
      ELSE IF ( imech_type == CHEM_TYPE_QUASIBS     ) THEN
        CALL CHEM_QUASIBS(rho, T, gama, s, hs, hq(1,4))
      END IF
    END IF
#endif

    ! ###################################################################
    ! Impose boundary conditions
    ! Temperature array T is used as auxiliary array because it is no
    ! longer used until the fields are updated
    ! ###################################################################
#define GAMMA_LOC(i) txc(i,6)
#define AUX_LOC(i)   T(i)

    CALL THERMO_GAMMA(imax, jmax, kmax, s, T, GAMMA_LOC(1))

    ! Maximum Mach for Poinsot & Lele reference pressure BC
    IF ( BcsDrift ) THEN
      M2_max = C_0_R
      DO i = 1,isize_field
        dummy = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))*rho(i)/(GAMMA_LOC(i)*p(i))
        M2_max = MAX(M2_max,dummy)
      END DO
#ifdef USE_MPI
      CALL MPI_ALLREDUCE(M2_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
      M2_max = dummy
#endif
    END IF

    IF ( .NOT. g(2)%periodic ) THEN
      CALL BOUNDARY_BCS_Y(isize_field, M2_max,        rho,u,v,w,p,GAMMA_LOC(1),s, &
          hq(1,5),hq(1,1),hq(1,2),hq(1,3),hq(1,4),hs,&
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), AUX_LOC(:), wrk2d,wrk3d)
    END IF

    IF ( .NOT. g(1)%periodic ) THEN
      CALL BOUNDARY_BCS_X(isize_field, M2_max, etime, rho,u,v,w,p,GAMMA_LOC(1),s, &
          hq(1,5),hq(1,1),hq(1,2),hq(1,3),hq(1,4),hs, txc, AUX_LOC(:), wrk1d,wrk2d,wrk3d)
    END IF

#undef GAMMA_LOC
#undef AUX_LOC

    ! ###################################################################
    ! Impose buffer zone as relaxation terms
    ! ###################################################################
    IF ( BuffType == DNS_BUFFER_RELAX .OR. BuffType == DNS_BUFFER_BOTH ) THEN
      CALL BOUNDARY_BUFFER_RELAX_FLOW(q,hq)
      CALL BOUNDARY_BUFFER_RELAX_SCAL(s,hs, q)
    END IF

    ! ###################################################################
    ! Perform the time stepping
    ! ###################################################################
    rho_ratio = C_1_R
    prefactor = (gama0-C_1_R)*mach*mach

    IF ( icalc_flow == 1 ) THEN
      IF ( icalc_scal == 1 ) THEN; inb_scal_loc = inb_scal
      ELSE;                          inb_scal_loc = 0
      END IF

      ! -------------------------------------------------------------------
      ! Total energy formulation
      ! -------------------------------------------------------------------
      IF ( imode_eqns == DNS_EQNS_TOTAL ) THEN
        !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
        !$omp do
        DO i = 1,isize_field
          rho_ratio    = rho(i)
          rho(i)       = rho(i) + dte*hq(i,5)
          rho_ratio    = rho_ratio/rho(i)
          dt_rho_ratio = dte/rho(i)

          e(i) = rho_ratio*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) ) &
              + dt_rho_ratio*hq(i,4)

          u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i,1)
          v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i,2)
          w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i,3)

          e(i) = e(i) - prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i))

          DO is = 1,inb_scal_loc
            s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
          END DO
        END DO
        !$omp end do
        !$omp end parallel

        ! -------------------------------------------------------------------
        ! Internal energy formulation
        ! -------------------------------------------------------------------
      ELSE IF ( imode_eqns == DNS_EQNS_INTERNAL ) THEN
        !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
        !$omp do
        DO i = 1,isize_field
          rho_ratio    = rho(i)
          rho(i)       = rho(i) + dte*hq(i,5)
          rho_ratio    = rho_ratio/rho(i)
          dt_rho_ratio = dte/rho(i)

          e(i) = rho_ratio*e(i) + dt_rho_ratio*hq(i,4)
          u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i,1)
          v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i,2)
          w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i,3)

          DO is = 1,inb_scal_loc
            s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
          END DO
        END DO
        !$omp end do
        !$omp end parallel

      END IF

    ELSE
      IF ( icalc_scal == 1 ) THEN
        DO is = 1,inb_scal
          !$omp parallel default( shared ) private( i, dt_rho_ratio )
          !$omp do
          DO i = 1,isize_field
            dt_rho_ratio = dte/rho(i)
            s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
          END DO
          !$omp end do
          !$omp end parallel
        END DO
      END IF
    END IF

    ! ###################################################################
    ! Impose buffer zone as filter
    ! ###################################################################
    IF ( BuffType == DNS_BUFFER_FILTER .OR. BuffType == DNS_BUFFER_BOTH ) THEN
      CALL BOUNDARY_BUFFER_FILTER&
          (rho,u,v,w,e,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    END IF

    RETURN
  END SUBROUTINE TIME_SUBSTEP_COMPRESSIBLE

END MODULE TIME

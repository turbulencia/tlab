#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/12 - J.P. Mellado
!#              Created
!# 2011/11/28 - C. Ansorge
!#              Combining BLAS and OMP
!# 2012/10/01 - C. Ansorge
!#              Semi-implicit formulation
!#
!########################################################################
!# DESCRIPTION
!#
!# Runge-Kutta explicit 3th order from Williamson 1980
!# Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
!# Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
!#
!########################################################################
SUBROUTINE TIME_RUNGEKUTTA(q,hq,s,hs, &
     x_inf,y_inf,z_inf,q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
     l_q, l_hq, l_txc, l_tags, l_comm)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : isize_field, inb_flow,inb_scal, icalc_particle
  USE DNS_GLOBAL, ONLY : icalc_flow,icalc_scal, imode_eqns
  USE DNS_GLOBAL, ONLY : isize_particle, inb_particle
  USE DNS_GLOBAL, ONLY : rtime
  USE LAGRANGE_GLOBAL, ONLY : particle_number, particle_bumper, inb_particle_evolution
  USE DNS_LOCAL

#ifdef USE_MPI
  USE DNS_MPI, ONLY : particle_vector, ims_pro,ims_npro,ims_npro_i,ims_npro_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif 

  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(*)             :: x_inf,y_inf,z_inf, q_inf,s_inf
  TREAL, DIMENSION(*)             :: txc, wrk1d,wrk2d,wrk3d, vaux

  TREAL, DIMENSION(isize_particle,*) :: l_q, l_hq
  TREAL, DIMENSION(*)                :: l_comm, l_txc
  INTEGER(8), DIMENSION(*)           :: l_tags

! -------------------------------------------------------------------
  TINTEGER i, is
  TREAL dte, etime, alpha
  TREAL kdt(5), kco(4), ktime(5)
  TREAL kex(3), kim(3) 

  TINTEGER srt,end,siz ! Variables for OpenMP Paritioning  
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

  IF      ( rkm_mode .EQ. RKM_EXP3           ) THEN; rkm_endstep = 3;
  ELSE IF ( rkm_mode .EQ. RKM_EXP4           ) THEN; rkm_endstep = 5;
  ELSE IF ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) THEN; rkm_endstep = 3;
  ENDIF

! -------------------------------------------------------------------
! Scheme coeefficients
! -------------------------------------------------------------------
  IF      ( rkm_mode .EQ. RKM_EXP3 ) THEN 
     kdt(1:3)   = (/ C_1_R/C_3_R, C_15_R/C_16_R,  C_8_R/C_15_R /) 
     ktime(1:3) = (/ C_0_R,       C_1_R/C_3_R,    C_3_R/C_4_R  /) 
     kco(1:2)   = (/-C_5_R/C_9_R,-C_153_R/C_128_R /) 

  ELSE IF ( rkm_mode .EQ. RKM_EXP4 ) THEN 
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

  ELSE IF ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) THEN 
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

  ENDIF

! -------------------------------------------------------------------
! Initialize arrays to zero for the explcit low-storage algorithm
! -------------------------------------------------------------------
  IF ( rkm_mode .EQ. RKM_EXP3 .OR. rkm_mode .EQ. RKM_EXP4 ) THEN 
     IF ( icalc_flow .EQ. 1 ) THEN
        DO is = 1,inb_flow; hq(:,is) = C_0_R; ENDDO
     ENDIF
     IF ( icalc_particle .EQ. 1 ) THEN
        DO is = 1,inb_particle_evolution; l_hq(:,is) = C_0_R; ENDDO
     ENDIF
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal; hs(:,is) = C_0_R; ENDDO
     ENDIF
  ENDIF

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
     IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
        IF    ( rkm_mode .EQ. RKM_EXP3 .OR. rkm_mode .EQ. RKM_EXP4 ) THEN 
           CALL TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT(&
                dte,etime, q,hq,s,hs,txc, vaux, wrk1d,wrk2d,wrk3d, &
                l_q, l_hq, l_txc, l_tags, l_comm)

        ELSE 
           CALL TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT(&
                dte,etime, kex(rkm_substep), kim(rkm_substep), kco(rkm_substep), &
                q,hq,s,hs,txc, vaux, wrk1d,wrk2d,wrk3d)
        ENDIF
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN
        CALL TIME_SUBSTEP_ANELASTIC(&
             dte,etime, q,hq,s,hs,txc, vaux, wrk1d,wrk2d,wrk3d)
     ELSE
        CALL TIME_SUBSTEP_COMPRESSIBLE(&
             dte,etime, q,hq,s,hs, x_inf,y_inf,z_inf, q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d)
     ENDIF

! -------------------------------------------------------------------
! Update RHS hq and hs in the explicit low-storage algorithm
! -------------------------------------------------------------------
     IF ( ( rkm_mode .EQ. RKM_EXP3 .OR. rkm_mode .EQ. RKM_EXP4 ) .AND. &
          rkm_substep .LT. rkm_endstep ) THEN

!$omp parallel default(shared) &
#ifdef USE_BLAS
!$omp private (ilen,srt,end,siz,alpha,is)
#else
!$omp private (i,   srt,end,siz,alpha,is)
#endif 

        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
#ifdef USE_BLAS 
        ilen = siz
#endif 

        alpha = kco(rkm_substep)

        IF ( icalc_flow .EQ. 1 ) THEN
           DO is = 1,inb_flow
#ifdef USE_BLAS
              CALL DSCAL(ilen, alpha, hq(srt,is), 1)
#else
              DO i = srt,end
                 hq(i,is) = alpha*hq(i,is)
              ENDDO
#endif
           ENDDO
        ENDIF

#ifdef USE_MPI
        IF ( icalc_particle .EQ. 1 ) THEN
           DO is = 1,inb_particle_evolution
              DO i = 1,particle_vector(ims_pro+1)
                 l_hq(i,is) = alpha*l_hq(i,is)
              ENDDO
           ENDDO
        ENDIF
#else
        IF ( icalc_particle .EQ. 1 ) THEN
           DO is = 1,inb_particle_evolution
              DO i = 1,particle_number
                 l_hq(i,is) = alpha*l_hq(i,is)
              ENDDO
           ENDDO
        ENDIF
#endif

        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal
#ifdef USE_BLAS
              CALL DSCAL(ilen, alpha, hs(srt,is), 1)
#else
              DO i = srt,end 
                 hs(i,is) = alpha*hs(i,is)
              ENDDO
#endif
           ENDDO
        ENDIF
!$omp end parallel 
     ENDIF

! -------------------------------------------------------------------
! Profiling data
! -------------------------------------------------------------------
#ifdef USE_PROFILE
     CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)
     idummy = t_end-t_srt

#ifdef USE_MPI 
     CALL MPI_REDUCE(idummy,t_dif,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD) 
     IF ( ims_pro .EQ. 0 ) THEN 
        WRITE(time_string,999) ims_npro, ims_npro_i, ims_npro_k, rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro 
999     FORMAT(I5.5,' (ims_npro_i X ims_npro_k:',I4.4,'x',I4.4,1x,') RK-Substep',I1,':', E13.5,'s') 
        CALL IO_WRITE_ASCII(lfile, time_string) 
     ENDIF
#else 
     t_dif = idummy
     WRITE(time_string,999) rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro 
999  FORMAT('RK-Substep',I1,':', E13.5,'s') 
     CALL IO_WRITE_ASCII(lfile,time_string) 
#endif 

#endif 

  ENDDO

  RETURN
END SUBROUTINE TIME_RUNGEKUTTA

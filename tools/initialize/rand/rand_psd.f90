!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Power spectral densities
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE RAND_PSD(nx,ny,nz, ispectrum, spc_param, random,seed, u)

  USE DNS_GLOBAL, ONLY : imax_total, jmax_total, kmax_total, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : scalex, scaley, scalez
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TINTEGER ispectrum, random, seed
  TREAL spc_param(*)
  TREAL, DIMENSION(isize_txc_dimz,nz), INTENT(INOUT) :: u

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,jglobal,kglobal, ip
  TREAL pow_dst, pow_org, phase
  TREAL f, f0,f1, fi,fj,fk
  TREAL RAN0

! #######################################################################
  f0 = spc_param(1); f1 = spc_param(4)

  DO k = 1,nz
#ifdef USE_MPI
     kglobal = k + ims_offset_k
#else
     kglobal = k
#endif
     IF ( kglobal .LE. kmax_total/2+1 ) THEN; fk = M_REAL(kglobal-1)/scalez
     ELSE;                                    fk =-M_REAL(kmax_total+1-kglobal)/scalez; ENDIF

     DO j = 1,ny
        jglobal = j ! No MPI decomposition along Oy
        IF ( jglobal .LE. jmax_total/2+1 ) THEN; fj = M_REAL(jglobal-1)/scaley
        ELSE;                                    fj =-M_REAL(jmax_total+1-jglobal)/scaley; ENDIF
              

        DO i = 1,nx/2+1
#ifdef USE_MPI
           iglobal = i + ims_offset_i/2
#else
           iglobal = i
#endif
           IF ( iglobal .LE. imax_total/2+1 ) THEN; fi = M_REAL(iglobal-1)/scalex
           ELSE;                                    fi =-M_REAL(imax_total+1-iglobal)/scalex; ENDIF
           
           f = SQRT(fi**2 + fj**2 + fk**2)

! -----------------------------------------------------------------------
! target psd
! -----------------------------------------------------------------------
           IF     (ispectrum .EQ. 1) THEN; pow_dst = C_1_R
           ELSEIF (ispectrum .EQ. 2) THEN; pow_dst =(f/f0)**4/(1.+12./5.*(f/f0)**2)**(17./6.)
           ELSEIF (ispectrum .EQ. 3) THEN; pow_dst = f**4     * EXP(-C_2_R*(f/f0)**2)
           ELSEIF (ispectrum .EQ. 4) THEN; pow_dst = f**2     * EXP(-C_2_R* f/f0    )
           ELSEIF (ispectrum .EQ. 5) THEN; pow_dst = f**4/16. * EXP(-C_2_R*(f/f0)**2)
! Gaussian bell (sigma=f1)
           ELSEIF (ispectrum .EQ. 6) THEN; pow_dst = EXP(-C_05_R*((f-f0)/f1)**2)/(f1*SQRT(C_2_R*C_PI_R))
           ENDIF

           IF ( (f-spc_param(2))*(spc_param(3)-f) .GT. C_0_R ) pow_dst = C_0_R ! Clip

           IF ( f .EQ. C_0_R ) THEN
              pow_dst = C_0_R
           ELSE
              IF ( jmax_total .EQ. 1 .OR. kmax_total .EQ. 1 ) THEN ! 2D spectrum
                 pow_dst = pow_dst / (C_PI_R*f)
              ELSE
                 pow_dst = pow_dst / (2*C_PI_R*f**2)
              ENDIF
           ENDIF

! -----------------------------------------------------------------------
! phase and scaling of complex data
! -----------------------------------------------------------------------
           ip = (nx+2)*(j-1) + 2*i 

           IF ( random .EQ. 2 ) THEN
              IF ( iglobal .EQ. 1 .OR. iglobal .EQ. imax_total/2+1 ) THEN; phase = C_0_R
              ELSE;                                                        phase = (RAN0(seed)-C_05_R)*C_2_R*C_PI_R; ENDIF
              u(ip-1,k) = SQRT(pow_dst)*COS(phase)
              u(ip  ,k) = SQRT(pow_dst)*SIN(phase)

           ELSE
              pow_org = u(ip-1,k)**2 + u(ip,k)**2
              
              IF ( pow_org .GT. C_0_R ) pow_dst = SQRT(pow_dst/pow_org)
           
              u(ip-1,k) = u(ip-1,k) * pow_dst
              u(ip  ,k) = u(ip  ,k) * pow_dst

           ENDIF

        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RAND_PSD


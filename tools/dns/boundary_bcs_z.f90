#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2003/06/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!# 
!# BCs at zmin and zmax.
!# Not yet fully developed like X or Y case because we did not need it
!# so far.
!#
!########################################################################
!# ARGUMENTS 
!# 
!########################################################################
SUBROUTINE BOUNDARY_BCS_Z(M2_max, rho, u,v,w,p, gama, z1, h0,h1,h2,h3,h4, zh1,&
     tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture, gama0, THERMO_AI
  USE DNS_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL M2_max

  TREAL, DIMENSION(imax,jmax,kmax)   :: rho, u, v, w, p, gama, h0, h1, h2, h3, h4
  TREAL, DIMENSION(imax,jmax,kmax)   :: tmp1, tmp2, tmp3, tmp4, tmp5
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1, zh1

  TREAL wrk1d(*)
  TREAL wrk2d(imax*jmax,*)
  TREAL wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i, nt, is, iwrk, inb_scal_loc, bcs(2,1)
  TREAL prefactor, pl_const

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_Z' )
#endif

#define hr_loc(i)  wrk2d(i,1)
#define hu_loc(i)  wrk2d(i,2)
#define hv_loc(i)  wrk2d(i,3)
#define hw_loc(i)  wrk2d(i,4)
#define he_loc(i)  wrk2d(i,5)
#define hz1_loc(i) wrk2d(i,6)

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  iwrk = 7
  nt   = imax*jmax
  prefactor = (gama0-C_1_R)*mach*mach

! ###################################################################
! Poinsot & Lele reference pressure BC
! The local value of c is added later at the boundary
! Note that pl_const has dimensions of 1/length
! ###################################################################
  IF ( bcs_euler_drift .EQ. 1 ) THEN
     pl_const = bcs_sigma_out*(C_1_R-M2_max) /g(3)%scale
  ELSE
     pl_const = C_0_R
  ENDIF

! ###################################################################
! Nonreflective BCs at zmin and zmax
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), rho, tmp1, wrk3d, wrk2d(1,iwrk), wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u,   tmp2, wrk3d, wrk2d(1,iwrk), wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v,   tmp3, wrk3d, wrk2d(1,iwrk), wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w,   tmp4, wrk3d, wrk2d(1,iwrk), wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), p,   tmp5, wrk3d, wrk2d(1,iwrk), wrk3d)

! -------------------------------------------------------------------
! BCs at zmin
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        CALL BOUNDARY_BCS_FLOW_NR_2(i0, nt, pl_const, bcs_p_kmin,&
             rho, w, u, v, p, gama, tmp1, tmp4, tmp2, tmp3, tmp5, buoyancy%vector(3),&
             hr_loc(1), hu_loc(1), hv_loc(1), hw_loc(1), he_loc(1))
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        CALL BOUNDARY_BCS_FLOW_NR_1(i0, nt, pl_const, bcs_p_kmin,&
             rho, w, u, v, p, gama, tmp1, tmp4, tmp2, tmp3, tmp5, buoyancy%vector(3),&
             hr_loc(1), hu_loc(1), hv_loc(1), hw_loc(1), he_loc(1))
     ENDIF
     DO i = 1,imax*jmax
        h0(i,1,1) = h0(i,1,1) + hr_loc(i)
        h1(i,1,1) = h1(i,1,1) + hv_loc(i)
        h2(i,1,1) = h2(i,1,1) + hw_loc(i)
        h3(i,1,1) = h3(i,1,1) + hu_loc(i)
        h4(i,1,1) = h4(i,1,1) + he_loc(i)*prefactor
     ENDDO
     IF ( imixture .GT. 0 ) THEN
        DO i = 1,imax*jmax
!           h4(i,1,1) = h4(i,1,1) + hr_loc(i)*THERMO_AI(6,1,NSP)
           h4(i,1,1) = h4(i,1,1) + hr_loc(i)*THERMO_AI(6,1,inb_scal+1)
        ENDDO
     ENDIF
#ifdef USE_MPI
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! BCs at zmax
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. ims_npro-1 ) THEN
#endif
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        CALL BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_const, bcs_p_kmax,&
             rho(1,1,kmax), w(1,1,kmax), u(1,1,kmax), v(1,1,kmax), p(1,1,kmax), gama(1,1,kmax),&
             tmp1(1,1,kmax), tmp4(1,1,kmax), tmp2(1,1,kmax), tmp3(1,1,kmax), tmp5(1,1,kmax), buoyancy%vector(3),&
             hr_loc(1), hu_loc(1), hv_loc(1), hw_loc(1), he_loc(1))
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        CALL BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_const, bcs_p_kmax,&
             rho(1,1,kmax), w(1,1,kmax), u(1,1,kmax), v(1,1,kmax), p(1,1,kmax), gama(1,1,kmax),&
             tmp1(1,1,kmax), tmp4(1,1,kmax), tmp2(1,1,kmax), tmp3(1,1,kmax), tmp5(1,1,kmax), buoyancy%vector(3),&
             hr_loc(1), hu_loc(1), hv_loc(1), hw_loc(1), he_loc(1))
     ENDIF
     DO i = 1,imax*jmax
        h0(i,1,kmax) = h0(i,1,kmax) + hr_loc(i)
        h1(i,1,kmax) = h1(i,1,kmax) + hv_loc(i)
        h2(i,1,kmax) = h2(i,1,kmax) + hw_loc(i)
        h3(i,1,kmax) = h3(i,1,kmax) + hu_loc(i)
        h4(i,1,kmax) = h4(i,1,kmax) + he_loc(i)*prefactor
     ENDDO
     IF ( imixture .GT. 0 ) THEN
        DO i = 1,imax*jmax
!           h4(i,1,kmax) = h4(i,1,kmax) + hr_loc(i)*THERMO_AI(6,1,NSP)
           h4(i,1,kmax) = h4(i,1,kmax) + hr_loc(i)*THERMO_AI(6,1,inb_scal+1)
        ENDDO
     ENDIF
#ifdef USE_MPI
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif

! ###################################################################
! Scalar nonreflective BCs at zmin and zmax
! ###################################################################
  IF ( icalc_scal .EQ. 1 ) THEN
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        inb_scal_loc = inb_scal + 1
     ELSE
        inb_scal_loc = inb_scal 
     ENDIF
     DO is = 1,inb_scal_loc
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1(1,1,1,is), tmp2, wrk3d, wrk2d(1,iwrk),wrk3d)

! -------------------------------------------------------------------
! BCs at zmin
! -------------------------------------------------------------------
#ifdef USE_MPI
        IF ( ims_pro .EQ. 0 ) THEN
#endif
           CALL BOUNDARY_BCS_SCAL_NR(i0, nt, pl_const, bcs_p_kmin, &
                rho, w, z1(1,1,1,is), p, gama, tmp1, tmp4, tmp2, tmp5, buoyancy%vector(3), hz1_loc(1))
! special case affects only energy equation
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
           ELSE
              DO i = 1,imax*jmax
                 zh1(i,1,1,is) = zh1(i,1,1,is) + hz1_loc(i)
              ENDDO
           ENDIF
           IF ( imixture .GT. 0 ) THEN
! special case
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
                 DO i = 1,imax*jmax
                    h4(i,1,1) = h4(i,1,1) + &
                         hz1_loc(i)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
                 ENDDO
! general case
              ELSE
                 DO i = 1,imax*jmax
!                    h4(i,1,1) = h4(i,1,1) + hz1_loc(i)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                    h4(i,1,1) = h4(i,1,1) + hz1_loc(i)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
                 ENDDO
              ENDIF
           ENDIF

#ifdef USE_MPI
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! BCs at zmax
! -------------------------------------------------------------------
#ifdef USE_MPI
        IF ( ims_pro .EQ. ims_npro-1 ) THEN
#endif
           CALL BOUNDARY_BCS_SCAL_NR(i1, nt, pl_const, bcs_p_kmax,&
                rho(1,1,kmax), w(1,1,kmax), z1(1,1,kmax,is), p(1,1,kmax), gama(1,1,kmax),&
                tmp1(1,1,kmax), tmp4(1,1,kmax), tmp2(1,1,kmax), tmp5(1,1,kmax),&
                buoyancy%vector(3), hz1_loc(1) )
! special case affects only energy equation
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
           ELSE
              DO i = 1,imax*jmax
                 zh1(i,1,kmax,is) = zh1(i,1,kmax,is) + hz1_loc(i)
              ENDDO
           ENDIF
           IF ( imixture .GT. 0 ) THEN
! special case
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
                 DO i = 1,imax*jmax
                    h4(i,1,kmax) = h4(i,1,kmax) + &
                         hz1_loc(i)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
                 ENDDO
! general case
              ELSE
                 DO i = 1,imax*jmax
!                    h4(i,1,kmax) = h4(i,1,kmax) + hz1_loc(i)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                    h4(i,1,kmax) = h4(i,1,kmax) + hz1_loc(i)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
                 ENDDO
              ENDIF
           ENDIF

#ifdef USE_MPI
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif
        
     ENDDO
  ENDIF
  
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_Z' )
#endif

  RETURN
END SUBROUTINE BOUNDARY_BCS_Z

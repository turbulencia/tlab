!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Compute the scalar quantity da/dx_i*du_j/dx_i*da/dx_j (strain1) and
!# that normalized by grad(a)*grad(a) (strain2)
!#
!########################################################################
!# ARGUMENTS 
!#
!# tmp1     Out  Returns grad(a)*grad(a)
!#
!########################################################################
      SUBROUTINE FI_STRAIN_A(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
           dx, dy, dz, a, u, v, w, strain1, strain2, &
           normal1, normal2, normal3, tmp1, wrk1d, wrk2d, wrk3d)

      IMPLICIT NONE

#include "types.h"
#include "integers.h"

      TINTEGER imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc

      TREAL dx(imax)
      TREAL dy(jmax)
      TREAL dz(kmax)

      TREAL a(imax, jmax, kmax)
      TREAL u(imax, jmax, kmax)
      TREAL v(imax, jmax, kmax)
      TREAL w(imax, jmax, kmax)
      TREAL strain1(imax, jmax, kmax)
      TREAL strain2(imax, jmax, kmax)

      TREAL normal1(imax,jmax,kmax)
      TREAL normal2(imax,jmax,kmax)
      TREAL normal3(imax,jmax,kmax)
      TREAL tmp1(imax,jmax,kmax)
      TREAL wrk1d(imax,*)
      TREAL wrk2d(imax,*)
      TREAL wrk3d(imax,jmax,kmax)

! -------------------------------------------------------------------
      TINTEGER i

! ###################################################################
! Compute normals
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, a, normal1, i0, i0, wrk1d, wrk2d, wrk3d)
      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, a, normal2, i0, i0, wrk1d, wrk2d, wrk3d)
      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, a, normal3, i0, i0, wrk1d, wrk2d, wrk3d)

! Compute gradient terms with u
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = normal1(i,1,1)*tmp1(i,1,1)*normal1(i,1,1)
      ENDDO
      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal2(i,1,1)*tmp1(i,1,1)*normal1(i,1,1)
      ENDDO
      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal3(i,1,1)*tmp1(i,1,1)*normal1(i,1,1)
      ENDDO
      
! Compute gradient terms with v
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal1(i,1,1)*tmp1(i,1,1)*normal2(i,1,1)
      ENDDO
      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal2(i,1,1)*tmp1(i,1,1)*normal2(i,1,1)
      ENDDO
      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal3(i,1,1)*tmp1(i,1,1)*normal2(i,1,1)
      ENDDO
      
! Compute gradient terms with 3
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal1(i,1,1)*tmp1(i,1,1)*normal3(i,1,1)
      ENDDO
      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal2(i,1,1)*tmp1(i,1,1)*normal3(i,1,1)
      ENDDO
      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      DO i = 1,imax*jmax*kmax
         strain1(i,1,1) = strain1(i,1,1) + normal3(i,1,1)*tmp1(i,1,1)*normal3(i,1,1)
      ENDDO
      
! norm of the gradient of a
      DO i = 1,imax*jmax*kmax
         tmp1(i,1,1) = normal1(i,1,1)**2 + normal2(i,1,1)**2 + normal3(i,1,1)**2 
      ENDDO

! normalize
      DO i = 1,imax*jmax*kmax
         IF ( tmp1(i,1,1) .GT. C_0_R ) THEN
            strain2(i,1,1) = strain1(i,1,1)/tmp1(i,1,1)
         ELSE
            strain2(i,1,1) = strain1(i,1,1)
         ENDIF
      ENDDO
         

      RETURN
      END

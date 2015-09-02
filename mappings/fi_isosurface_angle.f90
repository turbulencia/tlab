!########################################################################
!# Tool/Library FIELDS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/28 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the angle between isosurfaces of a and b at each point;
!# the gradient vectors are calculated first and then the angle between
!# them.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
      SUBROUTINE FI_ISOSURFACE_ANGLE(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
           dx, dy, dz, a, b, result, tmp1, tmp2, tmp3, tmp4, wrk1d, wrk2d, wrk3d)


      IMPLICIT NONE

#include "types.h"
#include "integers.h"

      TINTEGER imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc

      TREAL dx(imax)
      TREAL dy(jmax)
      TREAL dz(kmax)

      TREAL a(*), b(*), result(*)
      TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*)

      TREAL wrk1d(imax,*)
      TREAL wrk2d(imax,*)
      TREAL wrk3d(imax,jmax,kmax)

! -------------------------------------------------------------------
      TINTEGER ij

! ###################################################################
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, a, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
           dx, b, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
      DO ij = 1,imax*jmax*kmax
         result(ij) = tmp1(ij)*tmp2(ij)
         tmp3(ij) = tmp1(ij)*tmp1(ij)         
         tmp4(ij) = tmp2(ij)*tmp2(ij)         
      ENDDO

      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, a, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
           dy, b, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
      DO ij = 1,imax*jmax*kmax
         result(ij) = result(ij) + tmp1(ij)*tmp2(ij)
         tmp3(ij) = tmp3(ij) + tmp1(ij)*tmp1(ij)         
         tmp4(ij) = tmp4(ij) + tmp2(ij)*tmp2(ij)         
      ENDDO

      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, a, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
      CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
           dz, b, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
      DO ij = 1,imax*jmax*kmax
         result(ij) = result(ij) + tmp1(ij)*tmp2(ij)
         tmp3(ij) = tmp3(ij) + tmp1(ij)*tmp1(ij)         
         tmp4(ij) = tmp4(ij) + tmp2(ij)*tmp2(ij)         
         IF ( tmp3(ij) .GT. C_0_R .AND. tmp4(ij) .GT. C_0_R ) THEN
            result(ij) = result(ij)/SQRT(tmp3(ij)*tmp4(ij))
         ENDIF
      ENDDO

      RETURN
      END

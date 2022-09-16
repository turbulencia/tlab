#include "dns_const.h"

subroutine VELOCITY_MEAN(u, v, w, wrk1d)
    use TLAB_TYPES, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imode_sim, imax, jmax, kmax
    use TLAB_VARS, only: qbg
    use TLAB_VARS, only: coriolis
    implicit none

    real(wp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w
    real(wp), dimension(jmax, *), intent(INOUT) :: wrk1d

    ! -------------------------------------------------------------------
    integer(wi) iq, j
    real(wp) PROFILES, calpha, salpha
    external PROFILES

    !########################################################################
    if (imode_sim == DNS_MODE_TEMPORAL) then

        ! Construct reference profiles into array wrk1d
        do iq = 1, 3
            do j = 1, jmax
                wrk1d(j, iq) = PROFILES(qbg(iq), g(2)%nodes(j))
            end do
        end do

        ! Construct velocity field
        if (coriolis%type == EQNS_COR_NORMALIZED) then
            calpha = COS(coriolis%parameters(1)); salpha = SIN(coriolis%parameters(1))
            wrk1d(:, 3) = wrk1d(:, 3)*SIGN(1.0_wp, coriolis%vector(2)) ! right angular velocity vector (Garratt, 1992, p.42)

            do j = 1, jmax
                u(:, j, :) = u(:, j, :) + wrk1d(j, 1)*calpha + wrk1d(j, 3)*salpha
                v(:, j, :) = v(:, j, :) + wrk1d(j, 2)
                w(:, j, :) = w(:, j, :) - wrk1d(j, 1)*salpha + wrk1d(j, 3)*calpha
            end do

        else
            do j = 1, jmax
                u(:, j, :) = u(:, j, :) + wrk1d(j, 1)
                v(:, j, :) = v(:, j, :) + wrk1d(j, 2)
                w(:, j, :) = w(:, j, :) + wrk1d(j, 3)
            end do

        end if

        ! -------------------------------------------------------------------
    else if (imode_sim == DNS_MODE_SPATIAL) then
! I need to pass rho; need to check
! #define rho_vi(j) wrk1d(j,1)
! #define u_vi(j)   wrk1d(j,2)
! #define aux(j)    wrk1d(j,3)
!     ycenter = g(2)%nodes(1) + g(2)%scale *qbg(iq)%ymean_rel
!     DO j = 1,jmax
!       u_vi(j) = PROFILES( qbg(1), ycenter, g(2)%nodes(j) )
!     ENDDO
!     rho_vi(:) = rho(1,:,1)
!
!     ycenter = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean_rel
!     CALL FLOW_SPATIAL_VELOCITY(imax,jmax, &
!     qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%diam, ycenter, &
!     qbg(1)%parameters(2), qbg(1)%parameters(3), qbg(1)%parameters(4), &
!     g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), rho, u, v, aux(1), wrk3d)
!     IF ( g(3)%size .GT. 1 ) THEN
!       DO k = 2,kmax
!         u(:,:,k) = u(:,:,1)
!         v(:,:,k) = v(:,:,1)
!       ENDDO
!       w = w + qbg(3)%mean
!     ENDIF
! #undef rho_vi
! #undef u_vi
! #undef aux

    end if

    ! -------------------------------------------------------------------
    if (g(3)%size == 1) w = 0.0_wp

    return
end subroutine VELOCITY_MEAN

#include "types.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2007/10/05 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the equilibrium T and q_l for given enthalpy and pressure.
!# Iterative method based on (rho,e,q_i)
!#
!########################################################################
subroutine THERMO_AIRWATER_PH_RE(nx, ny, nz, z1, p, h, T)

    use THERMO_VARS, only: CRATIO_INV, THERMO_AI
    use THERMO_THERMAL

    implicit none

    TINTEGER nx, ny, nz
    TREAL z1(nx*ny*nz, *), T(*), h(*), p(*)

! -------------------------------------------------------------------
    TINTEGER ij, iter, niter
    TREAL r_loc(1), e_loc(1), t_loc(1), z1_loc(2), dummy

    integer, parameter :: i1 = 1

! ###################################################################
    niter = 5

    do ij = 1, nx*ny*nz
! -------------------------------------------------------------------
! initialize, q_l=0
! -------------------------------------------------------------------
        z1_loc(1) = z1(ij, 1)
        z1_loc(2) = C_0_R
        t_loc(1) = (h(ij) - THERMO_AI(6, 1, 2) - z1(ij, 1)*(THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 2)))/ &
                   (THERMO_AI(1, 1, 2) + z1(ij, 1)*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2)))

! -------------------------------------------------------------------
! iteration
! -------------------------------------------------------------------
        do iter = 1, niter
! calculate density from temperature/composition
            call THERMO_THERMAL_DENSITY(1, z1_loc, p(ij), t_loc, r_loc)

! calculate energy
            e_loc = h(ij) - CRATIO_INV*p(ij)/r_loc(1)

! solve equilibrium (rho,e,q_i)
            call THERMO_AIRWATER_RE(i1, z1_loc, e_loc, r_loc, t_loc, dummy)

        end do
        z1(ij, 2) = z1_loc(2)
        T(ij) = t_loc(1)

    end do

    return
end subroutine THERMO_AIRWATER_PH_RE

! ###################################################################
! ###################################################################
subroutine THERMO_ANELASTIC_AIRWATER_PH_RE(nx, ny, nz, s, e, p, wrk3d)

    use THERMO_VARS, only: CRATIO_INV, MRATIO
    use THERMO_ANELASTIC

    implicit none

    TINTEGER, intent(IN) :: nx, ny, nz
    TREAL, dimension(nx*ny*nz, *), intent(INOUT) :: s
    TREAL, dimension(*), intent(IN) :: e, p
    TREAL, dimension(nx*ny*nz), intent(OUT) :: wrk3d

! -------------------------------------------------------------------
    TINTEGER ij, jk, is, i, iter, niter
    TREAL r_loc(1), en_loc(1), t_loc(1), z1_loc(2), dummy
    TREAL p_loc(1), e_loc(1)

    integer, parameter :: i1 = 1

! ###################################################################
    niter = 5

    s(:, 3) = C_0_R ! initialize, q_l=0

    do iter = 1, niter ! iteration
! calculate density in wrk3d
        call THERMO_ANELASTIC_DENSITY(nx, ny, nz, s, e, p, wrk3d)

! calculate energy
        ij = 0
        do jk = 0, ny*nz - 1
            is = mod(jk, ny) + 1
            P_LOC = MRATIO*p(is)
            E_LOC = e(is)

            do i = 1, nx
                ij = ij + 1

                r_loc = wrk3d(ij)
                z1_loc(1) = s(ij, 2)
                z1_loc(2) = s(ij, 3)
                en_loc = s(ij, 1) - E_LOC - CRATIO_INV*P_LOC/r_loc

! solve equilibrium (rho,e,q_i)
                call THERMO_AIRWATER_RE(i1, z1_loc, en_loc, r_loc, t_loc, dummy)

                s(ij, 3) = z1_loc(2)

            end do
        end do
    end do

    return
end subroutine THERMO_ANELASTIC_AIRWATER_PH_RE

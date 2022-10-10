#include "types.h"
#include "dns_const.h"

subroutine FI_CHEM(chemistry, nx, ny, nz, is, s, source)

    use TLAB_TYPES, only: term_dt, profiles_dt
    use TLAB_VARS, only: sbg, damkohler
    use TLAB_VARS, only: g
    use PROFILES
    implicit none

    type(term_dt), intent(IN) :: chemistry
    TINTEGER, intent(IN) :: nx, ny, nz, is
    TREAL, dimension(nx*ny*nz, *), intent(IN) :: s
    TREAL, dimension(nx*ny*nz), intent(OUT) :: source

! -----------------------------------------------------------------------
    TREAL dummy, dummy2
    type(profiles_dt) prof_loc
    TINTEGER i, j, k

!########################################################################
    select case (chemistry%type)

    case (EQNS_CHEM_LAYEREDRELAXATION)
        prof_loc%type = PROFILE_TANH
        prof_loc%ymean = sbg(is)%ymean
        prof_loc%thick = -chemistry%parameters(3)*C_05_R
        prof_loc%mean = C_05_R
        prof_loc%delta = C_1_R
        prof_loc%lslope = C_0_R
        prof_loc%uslope = C_0_R
        do i = 1, nx
            do k = 1, nz
                do j = 1, ny
                    source(i + (j - 1)*nx + (k - 1)*nx*ny) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j) - chemistry%parameters(2)) ! strength constant
                end do
            end do
        end do

        dummy = -damkohler(is)/chemistry%parameters(1)
        source = dummy*source*s(:, is)

    case (EQNS_CHEM_QUADRATIC)
        dummy = damkohler(is)*chemistry%parameters(is)
        source = dummy*s(:, 2)*s(:, 3)

    case (EQNS_CHEM_QUADRATIC3)
        dummy = damkohler(is)*chemistry%parameters(is)

        if (is >= 1 .and. is <= 3) then
            source = dummy*s(:, 2)*s(:, 3)
        else if (is >= 4 .and. is <= 6) then
            source = dummy*s(:, 4)*s(:, 5)
        else if (is >= 7 .and. is <= 9) then
            source = dummy*s(:, 7)*s(:, 8)
        end if

    case (EQNS_CHEM_OZONE)
        dummy = damkohler(is)
        if (is == 4) dummy = -dummy

        source = -chemistry%parameters(1)/(C_1_R + chemistry%parameters(2)*s(:, 1))
        source = exp(source)

        if (is == 4) then
            dummy2 = C_1_R + chemistry%parameters(3)
            source = dummy*(dummy2*s(:, 4) - source*s(:, 2)*s(:, 3))
        else
            source = dummy*(s(:, 4) - source*s(:, 2)*s(:, 3))
        end if

    end select

    return
end subroutine FI_CHEM

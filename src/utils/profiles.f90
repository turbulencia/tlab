#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

subroutine PROFILES_READBLOCK(bakfile, inifile, block, tag, var)
    use TLAB_TYPES
    use TLAB_CONSTANTS, only: efile
    use TLAB_PROCS
    implicit none

    character(len=*), intent(in) :: bakfile, inifile, block, tag
    type(profiles_dp), intent(out) :: var

    character(len=512) sRes
    real(cp) derivative

    call SCANINICHAR(bakfile, inifile, block, 'Profile'//TRIM(ADJUSTL(tag)), 'none', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; var%type = PROFILE_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'tanh') then; var%type = PROFILE_TANH
    else if (TRIM(ADJUSTL(sRes)) == 'tanhsymmetric') then; var%type = PROFILE_TANH_SYM
    else if (TRIM(ADJUSTL(sRes)) == 'tanhantisymmetric') then; var%type = PROFILE_TANH_ANTISYM
    else if (TRIM(ADJUSTL(sRes)) == 'linear') then; var%type = PROFILE_LINEAR
    else if (TRIM(ADJUSTL(sRes)) == 'linearcrop') then; var%type = PROFILE_LINEAR_CROP
    else if (TRIM(ADJUSTL(sRes)) == 'erf') then; var%type = PROFILE_ERF
    else if (TRIM(ADJUSTL(sRes)) == 'erfsurface') then; var%type = PROFILE_ERF_SURFACE
    else if (TRIM(ADJUSTL(sRes)) == 'erfantisym') then; var%type = PROFILE_ERF_ANTISYM
    else if (TRIM(ADJUSTL(sRes)) == 'bickley') then; var%type = PROFILE_BICKLEY
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; var%type = PROFILE_GAUSSIAN
    else if (TRIM(ADJUSTL(sRes)) == 'ekman') then; var%type = PROFILE_EKMAN_U
    else if (TRIM(ADJUSTL(sRes)) == 'ekmanp') then; var%type = PROFILE_EKMAN_U_P
    else if (TRIM(ADJUSTL(sRes)) == 'parabolic') then; var%type = PROFILE_PARABOLIC
    else if (TRIM(ADJUSTL(sRes)) == 'mixedlayer') then; var%type = PROFILE_MIXEDLAYER
! the following 2 are used in initialize/flow/pressure_mean; should be cleaned
    else if (TRIM(ADJUSTL(sRes)) == 'enthalpyerf') then; var%type = -PROFILE_ERF
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong '//TRIM(ADJUSTL(tag))//' profile.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, block, 'Mean'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then ! Backwards compatibility
        call SCANINIREAL(bakfile, inifile, block, TRIM(ADJUSTL(tag)), '0.0', var%mean)
    else
        call SCANINIREAL(bakfile, inifile, block, 'Mean'//TRIM(ADJUSTL(tag)), '0.0', var%mean)
    end if
    call SCANINIREAL(bakfile, inifile, block, 'YCoor'//TRIM(ADJUSTL(tag)), '0.5', var%ymean)
    call SCANINIREAL(bakfile, inifile, block, 'Delta'//TRIM(ADJUSTL(tag)), '0.0', var%delta)
    call SCANINIREAL(bakfile, inifile, block, 'Thick'//TRIM(ADJUSTL(tag)), '0.0', var%thick)
    call SCANINIREAL(bakfile, inifile, block, 'BottomSlope'//TRIM(ADJUSTL(tag)), '0.0', var%lslope)
    call SCANINIREAL(bakfile, inifile, block, 'UpperSlope'//TRIM(ADJUSTL(tag)), '0.0', var%uslope)
    ! alternative to provide the variable thick in terms of the maximum derivative
    call SCANINICHAR(bakfile, inifile, block, 'Derivative'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) /= 'void') then
        call SCANINIREAL(bakfile, inifile, block, 'Derivative'//TRIM(ADJUSTL(tag)), '0.0', derivative)
        call PROFILES_DERTOTHICK(derivative, var)
    end if

    call SCANINIREAL(bakfile, inifile, block, 'SurfaceThick'//TRIM(ADJUSTL(tag)), '1.0', var%parameters(3))
    call SCANINIREAL(bakfile, inifile, block, 'SurfaceDelta'//TRIM(ADJUSTL(tag)), '0.0', var%parameters(4))
    ! alternative to provide the variable thick in terms of the maximum derivative
    call SCANINICHAR(bakfile, inifile, block, 'SurfaceDerivative'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) /= 'void') then
        call SCANINIREAL(bakfile, inifile, block, 'SurfaceDerivative'//TRIM(ADJUSTL(tag)), '0.0', derivative)
        call PROFILES_DERTOTHICK(derivative, var)
    end if

    call SCANINIREAL(bakfile, inifile, block, 'ScaleHeight', '0.0', var%parameters(5))

    call SCANINIREAL(bakfile, inifile, block, 'Diam'//TRIM(ADJUSTL(tag)), '0.0', var%diam)
    var%parameters(6) = var%diam

    return
end subroutine PROFILES_READBLOCK

function PROFILES(var, ycenter, y) result(f)
    use TLAB_TYPES, only: profiles_dp
    implicit none

    type(profiles_dp), intent(in) :: var
    TREAL, intent(in) :: ycenter, y
    TREAL f

    ! -------------------------------------------------------------------
    TREAL yrel, xi, amplify, zamp, cnought

    ! ###################################################################
    yrel = y - ycenter ! position relative to ycenter
    amplify = C_0_R    ! default

    ! -------------------------------------------------------------------
    ! base state varying between two constant levels
    ! -------------------------------------------------------------------
    if (var%thick == C_0_R) then
        if (var%type > 0) amplify = C_05_R*sign(C_1_R, yrel)

    else
        xi = yrel/var%thick

        select case (var%type)

        case (PROFILE_LINEAR)
            amplify = -xi

        case (PROFILE_TANH)
            amplify = C_05_R*TANH(-C_05_R*xi)

        case (PROFILE_TANH_SYM)
            amplify = C_05_R*(TANH(-C_05_R*(xi - C_05_R*var%parameters(6)/var%thick)) &
                              + TANH(C_05_R*(xi + C_05_R*var%parameters(6)/var%thick)) - C_1_R)

        case (PROFILE_TANH_ANTISYM)
            amplify = C_025_R*(TANH(-C_05_R*(xi - C_05_R*var%parameters(6)/var%thick)) &
                               - TANH(C_05_R*(xi + C_05_R*var%parameters(6)/var%thick)))

        case (PROFILE_ERF, PROFILE_ERF_ANTISYM, PROFILE_ERF_SURFACE)
            amplify = C_05_R*ERF(-C_05_R*xi)

        case (PROFILE_PARABOLIC, PROFILE_PARABOLIC_SURFACE)
            amplify = (C_1_R + C_05_R*xi)*(C_1_R - C_05_R*xi)

        case (PROFILE_BICKLEY)
            amplify = C_1_R/(COSH(C_05_R*xi))**C_2_R

        case (PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SURFACE)
            amplify = EXP(-C_05_R*xi**C_2_R)

        case (PROFILE_GAUSSIAN_SYM)
            amplify = EXP(-C_05_R*(xi - C_05_R*var%parameters(6)/var%thick)**C_2_R) &
                      + EXP(-C_05_R*(xi + C_05_R*var%parameters(6)/var%thick)**C_2_R)

        case (PROFILE_GAUSSIAN_ANTISYM)
            amplify = EXP(-C_05_R*(xi - C_05_R*var%parameters(6)/var%thick)**C_2_R) &
                      - EXP(-C_05_R*(xi + C_05_R*var%parameters(6)/var%thick)**C_2_R)

        case (PROFILE_EKMAN_U)
            amplify = C_1_R - EXP(-xi)*COS(xi)

        case (PROFILE_EKMAN_U_P)
            amplify = C_1_R - EXP(-xi)*COS(xi) ! + perturbation:

            cnought = C_PI_R*C_PI_R/C_4_R/C_4_R       ! Maximum initial Perturbation is at y=pi/2*var%thick
            zamp = SQRT(C_2_R)*xi*EXP(-xi*xi/C_8_R/cnought)/(var%thick*var%thick*C_4_R*cnought)**C_1_5_R
            amplify = amplify + zamp                  ! Add Perturbations

        case (PROFILE_EKMAN_V)
            amplify = -EXP(-xi)*SIN(xi)

        end select

    end if

    ! var%mean profile plus two linear-varying layers
    f = var%mean + var%delta*amplify &
        + var%lslope*yrel*C_05_R*(C_1_R - sign(C_1_R, yrel)) &
        + var%uslope*yrel*C_05_R*(C_1_R + sign(C_1_R, yrel))

    ! -------------------------------------------------------------------
    ! special profiles
    ! -------------------------------------------------------------------
    select case (var%type)

    case (PROFILE_LINEAR_CROP)
        if (yrel < C_0_R) then
            f = MIN(var%lslope*yrel, var%lslope*var%thick)
        else
            f = MAX(var%uslope*yrel, var%uslope*var%thick)
        end if

    case (PROFILE_MIXEDLAYER)
        if (yrel < C_0_R) then
            f = MIN(var%lslope*yrel, var%lslope*var%thick)
        else
            f = MAX(var%uslope*yrel, var%uslope*var%thick)
        end if
        f = f - C_025_R*var%uslope*var%thick*(C_1_R - SIGN(C_1_R, y - var%thick))

    case (PROFILE_ERF_SURFACE)
        xi = y/var%parameters(3)
        f = f + var%parameters(4)*C_05_R*(C_1_R + ERF(-C_05_R*xi))

    end select

    return
end function PROFILES

subroutine PROFILES_DERTOTHICK(derivative, var)  ! Obtain thick from the value of the maximum derivative
    use TLAB_TYPES, only: profiles_dp, cp
    use TLAB_CONSTANTS, only: efile
    use TLAB_PROCS
    implicit none

    real(cp), intent(in) :: derivative
    type(profiles_dp), intent(inout) :: var

    real(cp) thick_ratio    ! for readibility

    select case (var%type)

    case (PROFILE_TANH, PROFILE_TANH_SYM, PROFILE_TANH_ANTISYM)
        thick_ratio = 4.0_cp
        var%thick = -var%delta/derivative/thick_ratio

    case (PROFILE_ERF, PROFILE_ERF_ANTISYM)
        thick_ratio = 2.0_cp*sqrt(C_PI_R)
        var%thick = -var%delta/(derivative - var%uslope)/thick_ratio

    case (PROFILE_ERF_SURFACE)
        thick_ratio = 2.0_cp*sqrt(C_PI_R)
        var%parameters(3) = -var%parameters(4)/derivative/thick_ratio

    case default
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Undevelop derivative input for this case.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)

    end select

    return
end subroutine PROFILES_DERTOTHICK

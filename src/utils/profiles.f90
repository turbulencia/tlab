#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

function PROFILES(iflag, thick, delta, mean, ycenter, param, y)

    implicit none

    TINTEGER, intent(IN) :: iflag                                 ! type of profile
    TREAL, intent(IN) :: thick, delta, mean, ycenter, param(*) ! parameters defining the profile
    TREAL y, PROFILES

    ! -------------------------------------------------------------------
    TREAL yrel, xi, amplify, zamp, cnought

    ! ###################################################################
    yrel = y - ycenter ! position relative to ycenter
    amplify = C_0_R    ! default

    ! -------------------------------------------------------------------
    ! base state varying between two constant levels
    ! -------------------------------------------------------------------
    if (thick == C_0_R) then
        if (iflag > 0) then
            if (yrel <= C_0_R) then; amplify = C_05_R
            else if (yrel > C_0_R) then; amplify = -C_05_R; end if
        end if

    else
        xi = yrel/thick

        select case (iflag)

        case (PROFILE_LINEAR)
            amplify = -xi

        case (PROFILE_TANH)
            amplify = C_05_R*TANH(-C_05_R*xi)

        case (PROFILE_TANH_SYM)
            amplify = C_05_R*(TANH(-C_05_R*(xi - C_05_R*param(6)/thick)) + TANH(C_05_R*(xi + C_05_R*param(6)/thick)) - C_1_R)

        case (PROFILE_TANH_ANTISYM)
            amplify = C_025_R*(TANH(-C_05_R*(xi - C_05_R*param(6)/thick)) - TANH(C_05_R*(xi + C_05_R*param(6)/thick)))

        case (PROFILE_ERF, PROFILE_LINEAR_ERF, PROFILE_ERF_ANTISYM, PROFILE_ERF_SURFACE, PROFILE_LINEAR_ERF_SURFACE)
            amplify = C_05_R*ERF(-C_05_R*xi)

        case (PROFILE_PARABOLIC, PROFILE_PARABOLIC_SURFACE)
            amplify = (C_1_R + C_05_R*xi)*(C_1_R - C_05_R*xi)

        case (PROFILE_BICKLEY)
            amplify = C_1_R/(COSH(C_05_R*xi))**C_2_R

        case (PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SURFACE)
            amplify = EXP(-C_05_R*xi**C_2_R)

        case (PROFILE_GAUSSIAN_SYM)
            amplify = EXP(-C_05_R*(xi - C_05_R*param(6)/thick)**C_2_R) + EXP(-C_05_R*(xi + C_05_R*param(6)/thick)**C_2_R)

        case (PROFILE_GAUSSIAN_ANTISYM)
            amplify = EXP(-C_05_R*(xi - C_05_R*param(6)/thick)**C_2_R) - EXP(-C_05_R*(xi + C_05_R*param(6)/thick)**C_2_R)

        case (PROFILE_EKMAN_U)
            amplify = C_1_R - EXP(-xi)*COS(xi)

        case (PROFILE_EKMAN_U_P)
            amplify = C_1_R - EXP(-xi)*COS(xi) ! + perturbation:

            cnought = C_PI_R*C_PI_R/C_4_R/C_4_R       ! Maximum initial Perturbation is at y=pi/2*thick
            zamp = SQRT(C_2_R)*xi*EXP(-xi*xi/C_8_R/cnought)/(thick*thick*C_4_R*cnought)**C_1_5_R
            amplify = amplify + zamp                  ! Add Perturbations

        case (PROFILE_EKMAN_V)
            amplify = -EXP(-xi)*SIN(xi)

        end select

    end if

    ! mean profile
    if (ABS(delta) > C_0_R) then
        PROFILES = mean + delta*amplify
    else
        PROFILES = mean
    end if

    ! -------------------------------------------------------------------
    ! special profiles
    ! -------------------------------------------------------------------
    ! two linear-varying layers
    if (iflag == PROFILE_LINEAR_ERF .or. iflag == PROFILE_LINEAR_ERF_SURFACE) then
        if (yrel < C_0_R) then; PROFILES = PROFILES + param(1)*yrel
        else; PROFILES = PROFILES + param(2)*yrel; end if
    end if

    ! cropped linear
    if (iflag == PROFILE_LINEAR_CROP) then
        if (yrel < C_0_R) then; PROFILES = MIN(param(1)*yrel, param(1)*thick)
        else; PROFILES = MAX(param(2)*yrel, param(2)*thick); end if
    end if

    ! mixed layer
    if (iflag == PROFILE_MIXEDLAYER) then
        if (yrel < C_0_R) then; PROFILES = MIN(param(1)*yrel, param(1)*thick)
        else; PROFILES = MAX(param(2)*yrel, param(2)*thick); end if
        PROFILES = PROFILES - C_025_R*param(2)*thick*(C_1_R - SIGN(C_1_R, y - thick))
    end if

    ! adding surface flux
    if (iflag == PROFILE_ERF_SURFACE .or. iflag == PROFILE_LINEAR_ERF_SURFACE) then
        xi = y/param(3)
        PROFILES = PROFILES + param(4)*C_05_R*(C_1_R + ERF(-C_05_R*xi))
    end if

    return
end function PROFILES

subroutine PROFILE_READBLOCK(bakfile, inifile, block, tag, variable)
    use TLAB_TYPES
    use TLAB_CONSTANTS, only: efile
    use TLAB_PROCS
    character(len=*), intent(in) :: bakfile, inifile, block, tag
    type(background_dt), intent(out) :: variable

    character(len=512) sRes

    call SCANINICHAR(bakfile, inifile, block, 'Profile'//TRIM(ADJUSTL(tag)), 'none', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; variable%type = PROFILE_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'tanh') then; variable%type = PROFILE_TANH
    else if (TRIM(ADJUSTL(sRes)) == 'tanhsymmetric') then; variable%type = PROFILE_TANH_SYM
    else if (TRIM(ADJUSTL(sRes)) == 'tanhantisymmetric') then; variable%type = PROFILE_TANH_ANTISYM
    else if (TRIM(ADJUSTL(sRes)) == 'linear') then; variable%type = PROFILE_LINEAR
    else if (TRIM(ADJUSTL(sRes)) == 'linearcrop') then; variable%type = PROFILE_LINEAR_CROP
    else if (TRIM(ADJUSTL(sRes)) == 'linearerf') then; variable%type = PROFILE_LINEAR_ERF
    else if (TRIM(ADJUSTL(sRes)) == 'linearerfsurface') then; variable%type = PROFILE_LINEAR_ERF_SURFACE
    else if (TRIM(ADJUSTL(sRes)) == 'erf') then; variable%type = PROFILE_ERF
    else if (TRIM(ADJUSTL(sRes)) == 'erfsurface') then; variable%type = PROFILE_ERF_SURFACE
    else if (TRIM(ADJUSTL(sRes)) == 'erfantisym') then; variable%type = PROFILE_ERF_ANTISYM
    else if (TRIM(ADJUSTL(sRes)) == 'bickley') then; variable%type = PROFILE_BICKLEY
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; variable%type = PROFILE_GAUSSIAN
    else if (TRIM(ADJUSTL(sRes)) == 'ekman') then; variable%type = PROFILE_EKMAN_U
    else if (TRIM(ADJUSTL(sRes)) == 'ekmanp') then; variable%type = PROFILE_EKMAN_U_P
    else if (TRIM(ADJUSTL(sRes)) == 'parabolic') then; variable%type = PROFILE_PARABOLIC
    else if (TRIM(ADJUSTL(sRes)) == 'mixedlayer') then; variable%type = PROFILE_MIXEDLAYER
! the following 2 are used in initialize/flow/pressure_mean; should be cleaned
    else if (TRIM(ADJUSTL(sRes)) == 'enthalpyerf') then; variable%type = -PROFILE_ERF
    else if (TRIM(ADJUSTL(sRes)) == 'enthalpylinearerf') then; variable%type = -PROFILE_LINEAR_ERF
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong '//TRIM(ADJUSTL(tag))//' profile.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, block, 'Mean'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then ! Backwards compatibility
        call SCANINIREAL(bakfile, inifile, block, TRIM(ADJUSTL(tag)), '0.0', variable%mean)
    else
        call SCANINIREAL(bakfile, inifile, block, 'Mean'//TRIM(ADJUSTL(tag)), '0.0', variable%mean)
    end if
    call SCANINIREAL(bakfile, inifile, block, 'YCoor'//TRIM(ADJUSTL(tag)), '0.5', variable%ymean)
    call SCANINIREAL(bakfile, inifile, block, 'Thick'//TRIM(ADJUSTL(tag)), '0.0', variable%thick)
    call SCANINIREAL(bakfile, inifile, block, 'Delta'//TRIM(ADJUSTL(tag)), '0.0', variable%delta)

    call SCANINIREAL(bakfile, inifile, block, 'Diam'//TRIM(ADJUSTL(tag)), '0.0', variable%diam)
    variable%parameters(6) = variable%diam

    call SCANINICHAR(bakfile, inifile, block, 'BottomSlope'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then ! Backwards compatibility
        call SCANINIREAL(bakfile, inifile, block, 'BottomSlope', '0.0', variable%parameters(1))
    else
        call SCANINIREAL(bakfile, inifile, block, 'BottomSlope'//TRIM(ADJUSTL(tag)), '0.0', variable%parameters(1))
    end if
    call SCANINICHAR(bakfile, inifile, block, 'UpperSlope'//TRIM(ADJUSTL(tag)), 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then ! Backwards compatibility
        call SCANINIREAL(bakfile, inifile, block, 'UpperSlope', '0.0', variable%parameters(2))
    else
        call SCANINIREAL(bakfile, inifile, block, 'UpperSlope'//TRIM(ADJUSTL(tag)), '0.0', variable%parameters(2))
    end if
    call SCANINIREAL(bakfile, inifile, block, 'SurfaceThick'//TRIM(ADJUSTL(tag)), '1.0', variable%parameters(3))
    call SCANINIREAL(bakfile, inifile, block, 'SurfaceDelta'//TRIM(ADJUSTL(tag)), '0.0', variable%parameters(4))

    call SCANINIREAL(bakfile, inifile, block, 'ScaleHeight', '0.0', variable%parameters(5))

end subroutine PROFILE_READBLOCK

! Distributions and density functions

module Distributions
    use TLab_Constants, only: wi, wp, pi_wp, efile, wfile, MAX_PARS
    implicit none
    private

    type, public :: distributions_dt
        sequence
        integer(wi) :: type
        real(wp) :: mean
        real(wp) :: sigma
        real(wp) :: parameters(MAX_PARS) = 0.0_wp
    end type distributions_dt

    integer, parameter, public :: TYPE_DF_NONE = 0
    integer, parameter, public :: TYPE_DF_UNIFORM = 1
    integer, parameter, public :: TYPE_DF_QUARTIC = 3
    integer, parameter, public :: TYPE_DF_QUADRATIC = 4
    integer, parameter, public :: TYPE_DF_GAUSSIAN = 6
    integer, parameter, public :: TYPE_DF_ERF = 7

    public :: Distributions_Compute

contains
    function Distributions_Compute(locVar, x) result(f)
        type(distributions_dt), intent(in) :: locVar
        real(wp), intent(in) :: x
        real(wp) :: f

        ! -----------------------------------------------------------------------
        real(wp) x0, x1

        ! #######################################################################
        x0 = locVar%mean
        x1 = locVar%sigma

        select case (locVar%type)
        case (TYPE_DF_UNIFORM)
            f = 1.0_wp

        case (TYPE_DF_QUARTIC)
            f = x**4*exp(-2.0_wp*(x/x0)**2)

        case (TYPE_DF_QUADRATIC)
            f = x**2*exp(-2.0_wp*x/x0)

        case (TYPE_DF_GAUSSIAN)
            f = exp(-0.5_wp*((x - x0)/x1)**2)/(x1*sqrt(2.0_wp*pi_wp))

            ! case ()
            !     f = (x/x0)**4/(1.0_wp + 12.0_wp/5.0_wp*(x/x0)**2)**(17.0_wp/6.0_wp)

        case (TYPE_DF_ERF)
            f = 0.5_wp*(erf((log(f) - x0)/x1) + 1.0_wp)

        case default

        end select

        ! -----------------------------------------------------------------------
        if ((x - locVar%parameters(1))*(locVar%parameters(2) - x) < 0.0_wp) f = 0.0_wp ! Clip

        return
    end function Distributions_Compute

end module Distributions

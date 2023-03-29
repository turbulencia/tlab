module THERMO_VARS
    use TLAB_CONSTANTS, only: MAX_NSP, MAX_PROF, wp, wi

    implicit none
    save

    ! Compressibility, different combinations of parameters \gamma0 and mach to save calculations
    real(wp) :: gama0                           ! Specific heat ratio, Cp0/Cv0 = Cp0/(Cp0-R0)
    real(wp) :: GRATIO                          ! (gama0-1)/gama0 = R0/Cp0
    real(wp) :: MRATIO                          ! gama0 mach^2 = (U0^2/T0)/R0
    real(wp) :: RRATIO                          ! 1/MRATIO = R0/(U0^2/T0)
    real(wp) :: CRATIO_INV                      ! GRATIO*MRATIO = (gamma0-1)*mach^2 = (U0^2/T0)/Cp0

    ! Anelastic formulation
    real(wp) :: scaleheight                     ! Same as Fr/MRATIO in compressible formulation

    ! Mixture
    integer(wi) :: imixture
    character*128 :: chemkin_file               ! File with thermodynamic data, if used

    ! NSP_MAX is defined in global TLAB_CONSTANTS because it is used as maximum number of scalars
    integer(wi) :: NSP = 0                      ! Number of components (species) in a mixture
    character(len=32) :: THERMO_SPNAME(MAX_NSP) = ''
    real(wp) :: WGHT_INV(MAX_NSP)               ! Inverse of molar masses, i.e., gas constants
    real(wp) :: THERMO_R(MAX_NSP)               ! Normalized gas constants

    ! Caloric data
    integer(wi), parameter :: MAX_NCP = 7       ! Caloric data; cp(T), formation enthalpy, formation entropy
    integer(wi) :: NCP                          ! Number of terms in polynomial fit to cp
    real(wp) :: THERMO_AI(MAX_NCP, 2, MAX_NSP)  ! Polynomial coefficients. Last to terms are formation enthalpy and formation entropy
    !                                           The second index indicates each of the 2 temperature intervals considered in the fit
    real(wp) :: THERMO_TLIM(3, MAX_NSP)         ! Temperature limits of the two temperature intervals in the polynomial fits.

    ! Saturation vapor pressures
    integer(wi), parameter :: MAX_NPSAT = 10    ! Polynomial fit to saturation pressure
    integer(wi) :: NPSAT
    real(wp) :: THERMO_PSAT(MAX_NPSAT), NEWTONRAPHSON_ERROR

    logical nondimensional                      ! Nondimensional formulation
    real(wp) :: WREF, TREF                      ! Reference values; together with gama0, they contain all information
    real(wp), parameter :: RGAS = 8314_wp       ! Universal gas constant, J /kg /K

    real(wp) :: thermo_param(MAX_PROF)          ! Additional data

    real(wp) :: dsmooth                         ! Smoothing factor for derivaative discontinuity in inifinitely fast chemistry and saturation adjustment

end module THERMO_VARS

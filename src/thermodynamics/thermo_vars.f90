module THERMO_VARS
    use TLAB_CONSTANTS, only: MAX_NSP, MAX_PROF, wp, wi

    implicit none
    save

    real(wp) :: gama0                           ! Compressibility. Specific heat ratio
    real(wp) :: GRATIO                          ! (gama0-1)/gama0 to save calculations
    real(wp) :: MRATIO                          ! gama0 mach^2 to save calculations
    real(wp) :: MRATIO_INV                      ! 1/MRATIO
    real(wp) :: scaleheight                     ! for anelastic formulation; Fr/MRATIO in compressible formulation

    integer(wi) :: imixture
    character*128 :: chemkin_file               ! File with thermodynamic data, if used

    ! NSP_MAX is defined in global TLAB_CONSTANTS because it is used as maximum number of scalars
    integer(wi) :: NSP = 0                      ! Number of components (species) in a mixture
    character(len=32) :: THERMO_SPNAME(MAX_NSP) = ''
    real(wp) :: WGHT_INV(MAX_NSP)               ! Inverse of molar masses, i.e., gas constants
    real(wp) :: THERMO_R(MAX_NSP)               ! Normalized gas constants

    integer(wi), parameter :: MAX_NCP = 7       ! Caloric data; cp(T), formation enthalpy, formation entropy
    integer(wi) :: NCP                          ! Number of terms in polynomial fit to cp
    real(wp) :: THERMO_AI(MAX_NCP, 2, MAX_NSP)  ! Polynomial coefficients. Last to terms are formation enthalpy and formation entropy
    !                                         The second index indicates each of the 2 temperature intervals considered in the fit
    real(wp) :: THERMO_TLIM(3, MAX_NSP)         ! Temperature limits of the two temperature intervals in the polynomial fits.

    real(wp) :: dsmooth                         ! Smoothing factor for derivaative discontinuity in inifinitely fast chemistry and saturation adjustment

    integer(wi), parameter :: MAX_NPSAT = 10    ! Polynomial fit to saturation pressure
    integer(wi) :: NPSAT
    real(wp) :: THERMO_PSAT(MAX_NPSAT), NEWTONRAPHSON_ERROR

    logical nondimensional                      ! Nondimensional formulation
    real(wp) :: WREF, TREF                      ! Reference values; together with gama0, they contain all information
    real(wp), parameter :: RGAS = 8314_wp       ! Universal gas constant, J /kg /K

    real(wp) :: thermo_param(MAX_PROF)          ! Additional data

end module THERMO_VARS

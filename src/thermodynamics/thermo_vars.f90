#include "types.h"

module THERMO_VARS
    use TLAB_CONSTANTS, only: MAX_NSP, MAX_PROF

    implicit none
    save

    TREAL :: gama0                          ! Compressibility. Specific heat ratio
    TREAL :: GRATIO                         ! (gama0-1)/gama0 to save calculations
    TREAL :: MRATIO                         ! gama0 mach^2 to save calculations

    TINTEGER :: imixture
    character*128 :: chemkin_file           ! File with thermodynamic data, if used

    ! NSP_MAX is defined in global TLAB_CONSTANTS because it is used as maximum number of scalars
    TINTEGER :: NSP                         ! Number of components (species) in a mixture
    character*16, dimension(MAX_NSP) :: THERMO_SPNAME
    TREAL, dimension(MAX_NSP) :: WGHT_INV   ! Inverse of molar masses, i.e., gas constants

    TINTEGER, parameter :: MAX_NCP = 7      ! Caloric data; cp(T), formation enthalpy, formation entropy
    TINTEGER :: NCP                         ! Number of terms in polynomial fit to cp
    TREAL :: THERMO_AI(MAX_NCP, 2, MAX_NSP) ! Polynomial coefficients. Last to terms are formation enthalpy and formation entropy
    !                                         The second index indicates each of the 2 temperature intervals considered in the fit
    TREAL :: THERMO_TLIM(3, MAX_NSP)        ! Temperature limits of the two temperature intervals in the polynomial fits.

    TREAL :: dsmooth                        ! Smoothing factor for derivaative discontinuity in inifinitely fast chemistry and saturation adjustment

    TINTEGER, parameter :: MAX_NPSAT = 10   ! Polynomial fit to saturation pressure
    TINTEGER :: NPSAT
    TREAL :: THERMO_PSAT(MAX_NPSAT), NEWTONRAPHSON_ERROR

    logical nondimensional                  ! Nondimensional formulation
    TREAL :: WREF, TREF                     ! Reference values; together with gama0, they contain all information
    TREAL, parameter :: RGAS = 8314d0       ! Universal gas constant, J /kg /K

    TREAL :: thermo_param(MAX_PROF)         ! Additional data

end module THERMO_VARS

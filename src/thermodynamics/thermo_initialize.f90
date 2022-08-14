#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# inb_scal          # of scalars transported during the simulation and saved
!# inb_scal_array    # of scalars in array z1 (normally = inb_scal)
!# NSP               # of species in the mixture (NSP>=inb_scal)
!#
!# The code handles reactive/non-reactive, multiple/single species:
!# 1. Reactive => Multispecies.
!# 2. Multispecies admits the case Y_i=f_i(Z), Z conserved scalar, which implies inb_scal=1.
!# 3. Reactive + general multispecies retains NSP-1 species in scalar
!#    array z1 (the last one is obtained by sum Y_i=1) and an
!#    additional conserved scalar, i.e. inb_scal=NSP.
!# 4. Non-reactive + general multispecies retains NSP-1 species, w/o
!#    additional conserved scalar.
!#
!# Multispecies implies that reference T_0 in non-dimensionalization is 298 K.
!#
!# Saturation pressure implies that reference R_0 in non-dimensionalization
!# is such that reference pressure is 1 bar.
!#
!# The new reference value of gamma0 is calculated here based on the
!# reference species
!#
!########################################################################
subroutine THERMO_INITIALIZE

    use TLAB_TYPES
    use THERMO_VARS

    use TLAB_CONSTANTS, only: efile, lfile
    use TLAB_VARS, only: inb_scal, inb_scal_array
    use TLAB_VARS, only: damkohler
    use TLAB_PROCS
    implicit none

! -------------------------------------------------------------------
    TREAL WGHT(MAX_NSP)
    TINTEGER icp, is, im
    TREAL TREF_LOC, HREF(MAX_NSP), SREF(MAX_NSP)
    TINTEGER ISPREF                     ! reference species for CPREF and WREF
    TREAL CPREF
    TREAL WRK1D_LOC(MAX_NPSAT)
    TREAL PREF_LOC
    TINTEGER ipsat, i, j
    TREAL tmp1, tmp2
    character*46 str

! ###################################################################
! Species tags and molar masses in kg/kmol
! ###################################################################
    WGHT = C_1_R        ! Initialize molar masses to 1 kg /kmol

    select case (imixture)

! -------------------------------------------------------------------
! Burke-Schuman case
! Transport just mixture fraction, and then equilibrium
! 4 species + Nitrogen + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_BS, MIXT_TYPE_QUASIBS)
        NSP = 5
        inb_scal = 1
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'CH4'
        THERMO_SPNAME(2) = 'O2'
        THERMO_SPNAME(3) = 'H2O'
        THERMO_SPNAME(4) = 'CO2'
        THERMO_SPNAME(5) = 'N2'

        WGHT(1) = 16.0
        WGHT(2) = 32.0
        WGHT(3) = 18.0
        WGHT(4) = 44.0
        WGHT(5) = 28.0

! -------------------------------------------------------------------
! Peters Mechanism for Methane
! 7 species + Nitrogen + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_PETERS1991, MIXT_TYPE_PETERS1988)
        NSP = 8
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'CH4'
        THERMO_SPNAME(2) = 'O2'
        THERMO_SPNAME(3) = 'H2O'
        THERMO_SPNAME(4) = 'CO2'
        THERMO_SPNAME(5) = 'CO'
        THERMO_SPNAME(6) = 'H2'
        THERMO_SPNAME(7) = 'H'
        THERMO_SPNAME(8) = 'N2'

        WGHT(1) = 16.0
        WGHT(2) = 32.0
        WGHT(3) = 18.0
        WGHT(4) = 44.0
        WGHT(5) = 28.0
        WGHT(6) = 2.0
        WGHT(7) = 1.0
        WGHT(8) = 28.0

! -------------------------------------------------------------------
! Unimolecular decomposition flame
! 1 reactant + 1 product + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_UNIDECOMP)
        NSP = 2
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'R'
        THERMO_SPNAME(2) = 'P'

        WGHT(1) = 32.0
        WGHT(2) = 32.0

! -------------------------------------------------------------------
! Unimolecular decomposition flame
! 2 reactant + 1 product + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_ONESTEP)
        NSP = 4
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'R'   ! Reactant
        THERMO_SPNAME(2) = 'O'   ! Oxidizer
        THERMO_SPNAME(3) = 'P'   ! Product
        THERMO_SPNAME(4) = 'I'   ! Inert

        WGHT(1) = 32.0
        WGHT(2) = 32.0
        WGHT(3) = 32.0
        WGHT(4) = 32.0

! -------------------------------------------------------------------
! Swaminathan & Bilger Mechanism for Methane
! 5 species + Nitrogen + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_BILGER1997)
        NSP = 8
        inb_scal = NSP - 1
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'CH4'
        THERMO_SPNAME(2) = 'O2'
        THERMO_SPNAME(3) = 'H2O'
        THERMO_SPNAME(4) = 'CO2'
        THERMO_SPNAME(5) = 'CO'
        THERMO_SPNAME(6) = 'AR'
        THERMO_SPNAME(7) = 'N2'
        THERMO_SPNAME(8) = 'H2'

        WGHT(1) = 16.0
        WGHT(2) = 32.0
        WGHT(3) = 18.0
        WGHT(4) = 44.0
        WGHT(5) = 28.0
        WGHT(6) = 40.0
        WGHT(7) = 28.0
        WGHT(8) = 2.0

! -------------------------------------------------------------------
! Water vapor and air
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIR)
        NSP = 2
        inb_scal = MAX(inb_scal, 1) ! at least one scalar
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'H2O'
        THERMO_SPNAME(2) = 'AIR'

        WGHT(1) = 18.015_cp    ! from Iribarne and Godson, 1981
        WGHT(2) = 28.9644_cp

! -------------------------------------------------------------------
! Water vapor and air
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIRVAPOR)
        NSP = 2
        inb_scal = MAX(inb_scal, NSP) ! using inb_scal read in the inifile
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'H2O'
        THERMO_SPNAME(2) = 'AIR'

        WGHT(1) = 18.015_cp    ! from Iribarne and Godson, 1981
        WGHT(2) = 28.9644_cp

! -------------------------------------------------------------------
! Water vapor, air and liquid water
! Compressible:   Transport    q_t, and q_l from equilibrium; add space for q_l
! Incompressible: Transport h, q_t, and q_l from equilibrium; add space for q_l
! If non-equilibrium calculation, then inb_scal includes q_l and inb_scal_array = inb_scal (default)
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIRWATER)
        NSP = 3
        if (damkohler(3) <= C_0_R) inb_scal_array = inb_scal + 1 ! using inb_scal read in the inifile

        THERMO_SPNAME(1) = 'H2Ov'
        THERMO_SPNAME(2) = 'AIR'
        THERMO_SPNAME(3) = 'H2Ol'

        WGHT(1) = 18.015_cp    ! from Iribarne and Godson, 1981
        WGHT(2) = 28.9644_cp
        WGHT(3) = 18.015_cp

! -------------------------------------------------------------------
! Linearized thermodynamics for stratocumulus case
! The 1. scalar is the scaled total water (mixing fraction)
! The 2. scalar is the enthalpy deviations from the pure mixing case (described by mixing fraction only)
! The 3. scalar is the normalized concentration of liquid.
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIRWATER_LINEAR)
        inb_scal_array = inb_scal + 1 ! using inb_scal read in the inifile
        NSP = inb_scal_array

        THERMO_SPNAME(1) = 'Chi'      ! Mixture fraction
        THERMO_SPNAME(2) = 'Psi'      ! Deviation in the enthalpy from the mixture fraction
        do is = 3, inb_scal
            write (THERMO_SPNAME(is), *) is; THERMO_SPNAME(is) = 'Scalar'//TRIM(ADJUSTL(THERMO_SPNAME(is)))
        end do
        THERMO_SPNAME(NSP) = 'Liquid' ! Normalized Liquid

        WGHT(1) = 18.015_cp           ! unused, but defined for re-normalization below
        WGHT(2) = 28.9644_cp
        WGHT(3) = 18.015_cp
        WGHT(4:) = C_1_R

    end select

! ###################################################################
! Caloric equations
!
! Specific Heat Cpi, R^0 and SREF in Jules/(Kelvin Kmol), and
! enthalpy of formation HREF in Jules/Kmol.
!
! General formulation is CHEMKIN format: 7-coefficient NASA
! polynomials (see Burcat&Ruscic):
!
! THERMO_AI(i,im,k) = a_i of species k at
!         im=1 high temperature
!         im=2 low   temperature
!
! C_{p,i} = \sum_1^5 a_i T^{i-1}
! h_i     = \sum_1^5 a_i T^i/i + a_6
! s_{T,i} = a_1 ln(T) + \sum_2^5 a_i T^{i-1}/(i-1) + a_7
!
! i.e., dh_i = C_{p,i}dT and ds_{T,i}=C_{p,i}dT/T, where a_6 is
! related to the formation enthalpy, and a_7 to the formation entropy.
! HREF and SREF (at TREF_LOC) are used to fix last Cpi 6-7 coefficients.
!
! Note that Burcat&Ruscic give values devided by R^0
!
! The variables NCP gives the number of coefficients a_i.
! The simplified situation NCP=1 assumes also only one range, which then gives C_p constant.
! This is used to expedite the calculation of T in the energy formulation.
!
! The pressure contribution to the entropy still needs to be added
!
! ###################################################################
    THERMO_AI = C_0_R       ! Initialize to zero
    HREF = C_0_R
    SREF = C_0_R

    THERMO_TLIM(1, :) = 200.0_cp     ! Default T limits for the 2 intervals of polynomial fits
    THERMO_TLIM(2, :) = 5000.0_cp    ! These intervals are currently not used
    THERMO_TLIM(3, :) = 5000.0_cp

    select case (imixture)
! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR METHANE
! -------------------------------------------------------------------
    case (MIXT_TYPE_BS, MIXT_TYPE_QUASIBS)
        TREF_LOC = 298.0d0          ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF(1) = -74.0*C_1E6_R
        HREF(2) = C_0_R
        HREF(3) = -241.82*C_1E6_R
        HREF(4) = -393.51*C_1E6_R
        HREF(5) = C_0_R

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF(1) = 186.37*C_1E3_R
        SREF(2) = 205.15*C_1E3_R
        SREF(3) = 188.83*C_1E3_R
        SREF(4) = 213.78*C_1E3_R
        SREF(5) = 191.61*C_1E3_R

! Heat capacity polynomial. Values taken from fit with T-TREF_LOC instead T
        NCP = 2
        do im = 1, 2
            THERMO_AI(1, im, 1) = 35.70*C_1E3_R - 42.4833*TREF_LOC
            THERMO_AI(1, im, 2) = 28.96*C_1E3_R - 6.21666*TREF_LOC
            THERMO_AI(1, im, 3) = 32.76*C_1E3_R - 11.9570*TREF_LOC
            THERMO_AI(1, im, 4) = 37.22*C_1E3_R - 17.6500*TREF_LOC
            THERMO_AI(1, im, 5) = 28.88*C_1E3_R - 4.70833*TREF_LOC

            THERMO_AI(2, im, 1) = 42.4833
            THERMO_AI(2, im, 2) = 6.21666
            THERMO_AI(2, im, 3) = 11.9570
            THERMO_AI(2, im, 4) = 17.6500
            THERMO_AI(2, im, 5) = 4.70833

! 6th and 7th coefficient are calculated from reference enthalpy
            do is = 1, NSP
                THERMO_AI(6, im, is) = HREF(is) &
                                       - THERMO_AI(1, im, is)*TREF_LOC &
                                       - THERMO_AI(2, im, is)*TREF_LOC*TREF_LOC*C_05_R
                THERMO_AI(7, im, is) = SREF(is) &
                                       - THERMO_AI(2, im, is)*TREF_LOC
            end do
        end do

        do is = 1, NSP  ! Change heat capacities from molar to mass specific
            THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
        end do

! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR SIMPLE REACTION
! -------------------------------------------------------------------
    case (MIXT_TYPE_UNIDECOMP)
        TREF_LOC = 298.0d0          ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF(1) = C_0_R
        HREF(2) = -86.71502*C_1E6_R

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF(1) = 205.15*C_1E3_R
        SREF(2) = 205.15*C_1E3_R

! Heat capacity polynomial
        NCP = 1
        do im = 1, 2
            THERMO_AI(1, im, 1) = 29.099*C_1E3_R
            THERMO_AI(1, im, 2) = 29.099*C_1E3_R
! THERMO_AI(1,im,2) = 59.099*C_1E3_R

! 6th and 7th coefficient are calculated from reference enthalpy
            do is = 1, NSP
                THERMO_AI(6, im, is) = HREF(is) &
                                       - THERMO_AI(1, im, is)*TREF_LOC
                THERMO_AI(7, im, is) = SREF(is)
            end do
        end do

        do is = 1, NSP  ! Change heat capacities from molar to mass specific
            THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
        end do
    
! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR SIMPLE REACTION
! -------------------------------------------------------------------
    case (MIXT_TYPE_ONESTEP)
        TREF_LOC = 298.0d0          ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF(1) = C_0_R
        HREF(2) = C_0_R
        HREF(3) = -86.71502*C_1E6_R
        HREF(4) = C_0_R

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF(1) = 205.15*C_1E3_R
        SREF(2) = 205.15*C_1E3_R
        SREF(3) = 205.15*C_1E3_R
        SREF(4) = 205.15*C_1E3_R

! Heat capacity polynomial
        NCP = 1
        do im = 1, 2
            THERMO_AI(1, im, 1) = 29.099*C_1E3_R
            THERMO_AI(1, im, 2) = 29.099*C_1E3_R
            THERMO_AI(1, im, 3) = 29.099*C_1E3_R
            THERMO_AI(1, im, 4) = 29.099*C_1E3_R

! 6th and 7th coefficient are calculated from reference enthalpy
            do is = 1, NSP
                THERMO_AI(6, im, is) = HREF(is) &
                                       - THERMO_AI(1, im, is)*TREF_LOC
                THERMO_AI(7, im, is) = SREF(is)
            end do
        end do

        do is = 1, NSP  ! Change heat capacities from molar to mass specific
            THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
        end do
    
! -------------------------------------------------------------------
! Water vapor, air and water liquid
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)
        TREF_LOC = 273.15_cp

! Enthalpy of Formation in J /kg
        HREF(1) = 1870.0_cp*TREF_LOC                ! values s.t. THERMO_AI(6,im,1:2) = 0, i.e., liquid-water enthalpy
        HREF(2) = 1007.0_cp*TREF_LOC
        HREF(3) = 1870.0_cp*TREF_LOC -2501600_cp      ! latent heat of vaporization at 273.15 K is 2501.6 kJ/kg

! Entropy of Formation in J /kg /K
        SREF(1) = C_0_R
        SREF(2) = C_0_R
        SREF(3) = -2501600_cp/TREF_LOC

! Heat capacity polynomial in J /kg /K; using only the low temperature range 
        NCP = 1
        do im = 1, 2
            THERMO_AI(1, im, 1) = 1870.0_cp     ! water vapor
            THERMO_AI(1, im, 2) = 1007.0_cp     ! dry air
            THERMO_AI(1, im, 3) = 4217.6_cp     ! liquid water

! 6th and 7th coefficient are calculated from reference enthalpy
            do is = 1, NSP
                THERMO_AI(6, im, is) = HREF(is) - THERMO_AI(1, im, is)*TREF_LOC
                THERMO_AI(7, im, is) = SREF(is)
            end do
        end do

    case (MIXT_TYPE_CHEMKIN) ! Load thermodynamic data from chemkin file
        call THERMO_READ_CHEMKIN(chemkin_file)
        THERMO_AI = THERMO_AI*RGAS
        NCP = 5
    
        do is = 1, NSP  ! Change heat capacities from molar to mass specific
            THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
        end do

    end select

! -------------------------------------------------------------------
! Combination of pure species
! Intermediate species is CO + H2 (species 5 and 8)
! -------------------------------------------------------------------
    if (imixture == MIXT_TYPE_BILGER1997) then
        WGHT(5) = C_05_R*(WGHT(5) + WGHT(8))
        THERMO_AI(:, :, 5) = C_05_R*(THERMO_AI(:, :, 5) + THERMO_AI(:, :, 8))
    end if

! ###################################################################
! Phase change. Polynomial fit to saturation vapor pressure
!
! p_sat(T) = \sum_1^9 a_i T^{i-1}
!
! ###################################################################
    THERMO_PSAT = C_0_R ! Initialize to zero
    NPSAT = 0

    select case (imixture)

    case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)

! Flatau et al., J. Applied Meteorol., 1507-1513, 1992
        NPSAT = 9
        WRK1D_LOC(1) = 0.611213476*C_1E3_R
        WRK1D_LOC(2) = 0.444007856*C_1E2_R
        WRK1D_LOC(3) = 0.143064234*C_1E1_R
        WRK1D_LOC(4) = 0.264461437*C_1EM1_R
        WRK1D_LOC(5) = 0.305930558*C_1EM3_R
        WRK1D_LOC(6) = 0.196237241*C_1EM5_R
        WRK1D_LOC(7) = 0.892344772*C_1EM8_R
        WRK1D_LOC(8) = -0.373208410*C_1EM10_R
        WRK1D_LOC(9) = 0.209339997*C_1EM13_R

! going from powers of (T-T_ref) to T, with T_ref = 273.15 K
        TREF_LOC = 273.15D0
        do ipsat = 1, NPSAT
            THERMO_PSAT(ipsat) = C_0_R
            do i = ipsat, NPSAT
                tmp1 = C_1_R
                do j = i - 1, i - ipsat + 1, -1
                    tmp1 = tmp1*M_REAL(j)
                end do
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat) + &
                                     WRK1D_LOC(i)*TREF_LOC**(i - 1)*tmp1*(-1)**(i - ipsat)
            end do
            tmp2 = C_1_R
            do j = ipsat - 1, 1, -1
                tmp2 = tmp2*M_REAL(j)
            end do
            THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/tmp2/TREF_LOC**(ipsat - 1)
        end do

    case default

    end select

! ###################################################################
! Nondimensionalization
! ###################################################################
    TREF = 298.0d0                  ! K
    ISPREF = 2                      ! Species 2 is taken as reference
    WREF = WGHT(ISPREF)             ! kg /kmol

    THERMO_TLIM = THERMO_TLIM/TREF  ! Temperature limis for polynomial fits to cp

    WGHT_INV = WREF/WGHT            ! Inverse of molar masses, i.e., normalized gas constants

    CPREF = C_0_R                   ! Reference cp
    do icp = NCP, 1, -1
        CPREF = CPREF*TREF + THERMO_AI(icp, 2, ISPREF)
    end do

    do is = 1, NSP
        do im = 1, 2
            THERMO_AI(6, im, is) = THERMO_AI(6, im, is)/(CPREF*TREF)    ! Formation enthalpy
            THERMO_AI(7, im, is) = THERMO_AI(7, im, is)/CPREF           ! Formation entropy
            do icp = 1, NCP
                THERMO_AI(icp, im, is) = THERMO_AI(icp, im, is)*(TREF**(icp - 1))/CPREF
            end do
        end do
    end do

    gama0 = CPREF*WREF/(CPREF*WREF - RGAS)

    PREF_LOC = C_1E5_R      ! Pa
    do ipsat = 1, NPSAT
        THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/PREF_LOC
        THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)*(TREF**(ipsat - 1))
    end do

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
    call TLAB_WRITE_ASCII(lfile, 'Thermodynamic properties have been initialized.')
    do is = 1, NSP
        write (str, *) is; str = 'Setting Species'//TRIM(ADJUSTL(str))//'='//TRIM(ADJUSTL(THERMO_SPNAME(is)))
        call TLAB_WRITE_ASCII(lfile, str)
    end do
    write (str, 1010) 'Setting WREF = ', WREF
    call TLAB_WRITE_ASCII(lfile, str)
    write (str, 1010) 'Setting TREF = ', TREF
    call TLAB_WRITE_ASCII(lfile, str)
    if (NPSAT > 0) then
        write (str, 1010) 'Setting RREF = ', PREF_LOC/(RGAS/WREF*TREF)
        call TLAB_WRITE_ASCII(lfile, str)
    end if
    write (str, 1020) 'Setting CPREF = ', CPREF
    call TLAB_WRITE_ASCII(lfile, str)
    write (str, 1020) 'Setting Gama0 = ', gama0
    call TLAB_WRITE_ASCII(lfile, str)

    return

1010 format(A14, 1X, G_FORMAT_R)
1020 format(A15, 1X, G_FORMAT_R)

end subroutine THERMO_INITIALIZE

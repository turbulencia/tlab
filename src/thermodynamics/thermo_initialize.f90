#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# inb_scal          # of scalars transported during the simulation and saved
!# inb_scal_array    # of scalars in array s (normally = inb_scal)
!# NSP               # of species in the mixture (NSP>=inb_scal)
!#
!# The code handles reactive/non-reactive, multiple/single species:
!# 1. Reactive => Multispecies.
!# 2. Multispecies admits the case Y_i=f_i(Z), Z conserved scalar, which implies inb_scal=1.
!# 3. Reactive + general multispecies retains NSP-1 species in scalar
!#    array s (the last one is obtained by sum Y_i=1) and an
!#    additional conserved scalar, i.e. inb_scal=NSP.
!# 4. Non-reactive + general multispecies retains NSP-1 species, w/o
!#    additional conserved scalar.
!#
!# Multispecies implies that reference T_0 in non-dimensionalization is 298 K.
!#
!# Saturation pressure implies that reference R_0 in non-dimensionalization
!# is such that reference pressure is 1 bar.
!#
!# gama0 has alread been read in dns.ini.
!# If needed, the new reference value of gamma0 is calculated here based on the reference species
!#
!########################################################################
subroutine THERMO_INITIALIZE()
    use TLAB_CONSTANTS, only: efile, lfile, wi, wp
    use TLAB_VARS, only: inb_scal, inb_scal_array, imode_eqns, mach
    use TLAB_VARS, only: damkohler, transport, radiation
    use TLAB_PROCS
    use THERMO_VARS
    implicit none

! -------------------------------------------------------------------
    real(wp) WGHT(MAX_NSP)
    integer(wi) icp, is, im, inb_scal_loc
    real(wp) TREF_LOC, HREF_LOC(MAX_NSP), SREF_LOC(MAX_NSP)
    integer(wi) ISPREF                     ! reference species for CPREF and WREF
    real(wp) CPREF
    real(wp) WRK1D_LOC(MAX_NPSAT)
    real(wp) PREF_LOC
    integer(wi) ipsat, i, j
    real(wp) tmp1, tmp2
    character*46 str
    logical :: molar_data = .true.

! ###################################################################
! Species tags
! Thermal equation, molar masses in kg/kmol
! ###################################################################
    inb_scal_loc = inb_scal     ! Control that inb_scal read in dns.ini is correct

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

        THERMO_SPNAME(1) = 'CH4'; WGHT(1) = 16.0_wp
        THERMO_SPNAME(2) = 'O2 '; WGHT(2) = 32.0_wp
        THERMO_SPNAME(3) = 'H2O'; WGHT(3) = 18.0_wp
        THERMO_SPNAME(4) = 'CO2'; WGHT(4) = 44.0_wp
        THERMO_SPNAME(5) = 'N2 '; WGHT(5) = 28.0_wp

! -------------------------------------------------------------------
! Peters Mechanism for Methane
! 7 species + Nitrogen + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_PETERS1991, MIXT_TYPE_PETERS1988)
        NSP = 8
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'CH4'; WGHT(1) = 16.0_wp
        THERMO_SPNAME(2) = 'O2 '; WGHT(2) = 32.0_wp
        THERMO_SPNAME(3) = 'H2O'; WGHT(3) = 18.0_wp
        THERMO_SPNAME(4) = 'CO2'; WGHT(4) = 44.0_wp
        THERMO_SPNAME(5) = 'CO '; WGHT(5) = 28.0_wp
        THERMO_SPNAME(6) = 'H2 '; WGHT(6) = 2.0_wp
        THERMO_SPNAME(7) = 'H  '; WGHT(7) = 1.0_wp
        THERMO_SPNAME(8) = 'N2 '; WGHT(8) = 28.0_wp

! -------------------------------------------------------------------
! Unimolecular decomposition flame
! 1 reactant + 1 product + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_UNIDECOMP)
        NSP = 2
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'R'; WGHT(1) = 32.0_wp
        THERMO_SPNAME(2) = 'P'; WGHT(2) = 32.0_wp

! -------------------------------------------------------------------
! Unimolecular decomposition flame
! 2 reactant + 1 product + Conserved Scalar
! -------------------------------------------------------------------
    case (MIXT_TYPE_ONESTEP)
        NSP = 4
        inb_scal = NSP
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'R'; WGHT(1) = 32.0_wp      ! Reactant
        THERMO_SPNAME(2) = 'O'; WGHT(2) = 32.0_wp      ! Oxidizer
        THERMO_SPNAME(3) = 'P'; WGHT(3) = 32.0_wp      ! Product
        THERMO_SPNAME(4) = 'I'; WGHT(4) = 32.0_wp      ! Inert

! -------------------------------------------------------------------
! Water vapor and air
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIR)
        NSP = 2
        inb_scal = max(inb_scal, 1) ! at least one scalar
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'H2O'; WGHT(1) = 18.015_wp    ! from Iribarne and Godson, 1981
        THERMO_SPNAME(2) = 'AIR'; WGHT(2) = 28.9644_wp

! -------------------------------------------------------------------
! Water vapor and air
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIRVAPOR)
        NSP = 2
        inb_scal = max(inb_scal, NSP) ! using inb_scal read in the inifile
        inb_scal_array = inb_scal

        THERMO_SPNAME(1) = 'H2O'; WGHT(1) = 18.015_wp    ! from Iribarne and Godson, 1981
        THERMO_SPNAME(2) = 'AIR'; WGHT(2) = 28.9644_wp

! -------------------------------------------------------------------
! Water vapor, air and liquid water
! Compressible:   Transport    q_t, and q_l from equilibrium; add space for q_l
! Incompressible: Transport h, q_t, and q_l from equilibrium; add space for q_l
! If non-equilibrium calculation, then inb_scal includes q_l and inb_scal_array = inb_scal (default)
! -------------------------------------------------------------------
    case (MIXT_TYPE_AIRWATER)
        NSP = 3
        if (damkohler(3) <= 0.0_wp) inb_scal_array = inb_scal + 1 ! using inb_scal read in the inifile

        THERMO_SPNAME(1) = 'H2Ov'; WGHT(1) = 18.015_wp    ! from Iribarne and Godson, 1981
        THERMO_SPNAME(2) = 'AIR '; WGHT(2) = 28.9644_wp
        THERMO_SPNAME(3) = 'H2Ol'; WGHT(1) = 18.015_wp

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
            write (THERMO_SPNAME(is), *) is; THERMO_SPNAME(is) = 'Scalar'//trim(adjustl(THERMO_SPNAME(is)))
        end do
        THERMO_SPNAME(NSP) = 'Liquid' ! Normalized Liquid

        WGHT(1) = 18.015_wp           ! unused, but defined for re-normalization below
        WGHT(2) = 28.9644_wp
        WGHT(3) = 18.015_wp

    end select

    if (inb_scal_loc /= inb_scal) then
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Incorrect number of Schmidt numbers.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! -------------------------------------------------------------------
! Other processes that need correct inb_scal_array
! should go in initialization of each process
! -------------------------------------------------------------------
! By default, transport and radiation are caused by last scalar
! The variable inb_scal_array is only available at the end of this routine
    transport%scalar = inb_scal_array
    radiation%scalar = inb_scal_array

    if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        if (radiation%type /= EQNS_NONE) then
            radiation%active(inb_scal_array) = .true. ! liquid
            radiation%active(inb_scal_array + 1) = .true. ! buoyancy
        end if

    end if

    ! ###################################################################
! Caloric equations
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
! HREF_LOC and SREF_LOC (at TREF_LOC) are used to fix last Cpi 6-7 coefficients.
!
! Note that Burcat&Ruscic give values devided by R^0
!
! NCP gives the number of coefficients a_i.
! The simplified situation NCP=1 assumes also only one range, which then gives C_p constant.
! This is used to expedite the calculation of T in the energy formulation.
!
! The pressure contribution to the entropy still needs to be added
!
! ###################################################################
    THERMO_AI(:,:,:) = 0.0_wp               ! Initialize to zero
    HREF_LOC(:) = 0.0_wp
    SREF_LOC(:) = 0.0_wp

    THERMO_TLIM(1, :) = 200.0_wp    ! Default T limits for the 2 intervals of polynomial fits
    THERMO_TLIM(2, :) = 5000.0_wp   ! These intervals are currently not used
    THERMO_TLIM(3, :) = 5000.0_wp

    NCP = 1                         ! Default order of heat capacity polynomial

    select case (imixture)
    case (MIXT_TYPE_BS, MIXT_TYPE_QUASIBS)      ! Methane
        TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF_LOC(1) = -74.0e+6_wp
        HREF_LOC(2) = 0.0_wp
        HREF_LOC(3) = -241.82e+6_wp
        HREF_LOC(4) = -393.51e+6_wp
        HREF_LOC(5) = 0.0_wp

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF_LOC(1) = 186.37e+3_wp
        SREF_LOC(2) = 205.15e+3_wp
        SREF_LOC(3) = 188.83e+3_wp
        SREF_LOC(4) = 213.78e+3_wp
        SREF_LOC(5) = 191.61e+3_wp

! Heat capacity polynomial. Values taken from fit with T-TREF_LOC instead T
        NCP = 2
        THERMO_AI(1, :, 1) = 35.70e+3_wp - 42.4833_wp*TREF_LOC; THERMO_AI(2, :, 1) = 42.4833_wp
        THERMO_AI(1, :, 2) = 28.96e+3_wp - 6.21666_wp*TREF_LOC; THERMO_AI(2, :, 2) = 6.21666_wp
        THERMO_AI(1, :, 3) = 32.76e+3_wp - 11.9570_wp*TREF_LOC; THERMO_AI(2, :, 3) = 11.9570_wp
        THERMO_AI(1, :, 4) = 37.22e+3_wp - 17.6500_wp*TREF_LOC; THERMO_AI(2, :, 4) = 17.6500_wp
        THERMO_AI(1, :, 5) = 28.88e+3_wp - 4.70833_wp*TREF_LOC; THERMO_AI(2, :, 5) = 4.70833_wp

! -------------------------------------------------------------------
    case (MIXT_TYPE_UNIDECOMP)      ! Simple reaction, R -> P 
        TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF_LOC(1) = 0.0_wp
        HREF_LOC(2) = -86.71502e+6_wp

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF_LOC(1) = 205.15e+3_wp
        SREF_LOC(2) = 205.15e+3_wp

! Heat capacity polynomial
        THERMO_AI(1, :, 1) = 29.099e+3_wp
        THERMO_AI(1, :, 2) = 29.099e+3_wp   ! = 59.099e+3_wp

! -------------------------------------------------------------------
    case (MIXT_TYPE_ONESTEP)        ! Simple reaction, R + O -> P
        TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

! Enthalpy of Formation in Jules/Kmol
        HREF_LOC(1) = 0.0_wp
        HREF_LOC(2) = 0.0_wp
        HREF_LOC(3) = -86.71502e+6_wp
        HREF_LOC(4) = 0.0_wp

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF_LOC(1) = 205.15e+3_wp
        SREF_LOC(2) = 205.15e+3_wp
        SREF_LOC(3) = 205.15e+3_wp
        SREF_LOC(4) = 205.15e+3_wp

! Heat capacity polynomial
        THERMO_AI(1, :, 1) = 29.099e+3_wp
        THERMO_AI(1, :, 2) = 29.099e+3_wp
        THERMO_AI(1, :, 3) = 29.099e+3_wp
        THERMO_AI(1, :, 4) = 29.099e+3_wp

! -------------------------------------------------------------------
    case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)
        TREF_LOC = 273.15_wp

        molar_data = .false. ! this block gives data already in spefic

! Enthalpy of Formation in J /kg
        HREF_LOC(1) = 1870.0_wp*TREF_LOC                 ! values s.t. THERMO_AI(6,im,1:2) = 0, i.e., liquid-water enthalpy
        HREF_LOC(2) = 1007.0_wp*TREF_LOC
        HREF_LOC(3) = 1870.0_wp*TREF_LOC - 2501600_wp    ! latent heat of vaporization at 273.15 K is 2501.6 kJ/kg

! Entropy of Formation in J /kg /K
        SREF_LOC(1) = 0.0_wp
        SREF_LOC(2) = 0.0_wp
        SREF_LOC(3) = -2501600_wp/TREF_LOC

! Heat capacity polynomial in J /kg /K; using only the low temperature range
        THERMO_AI(1, :, 1) = 1870.0_wp     ! water vapor
        THERMO_AI(1, :, 2) = 1007.0_wp     ! dry air
        THERMO_AI(1, :, 3) = 4217.6_wp     ! liquid water

! -------------------------------------------------------------------
    case (MIXT_TYPE_CHEMKIN) ! Load thermodynamic data from chemkin file
        call THERMO_READ_CHEMKIN(chemkin_file)
        THERMO_AI = THERMO_AI*RGAS
        NCP = 5             ! Fifth-order polynomial

    end select

! -------------------------------------------------------------------
! Final calculations

    ! 6th and 7th coefficients of heat capacity polynomial are calculated from reference enthalpy
    do is = 1, NSP
        THERMO_AI(6, :, is) = HREF_LOC(is) - THERMO_AI(1, :, is)*TREF_LOC - THERMO_AI(2, :, is)*TREF_LOC*TREF_LOC*0.5_wp
        THERMO_AI(7, :, is) = SREF_LOC(is) - THERMO_AI(2, :, is)*TREF_LOC
    end do

    if (molar_data) then ! Change heat capacities from molar to mass specific
        do is = 1, NSP
            THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
        end do
    end if

! ###################################################################
! Phase change. Polynomial fit to saturation vapor pressure
!
! p_sat(T) = \sum_1^9 a_i T^{i-1}
!
! ###################################################################
    THERMO_PSAT(:) = 0.0_wp ! Initialize to zero
    NPSAT = 0

    select case (imixture)

    case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)

! Flatau et al., J. Applied Meteorol., 1507-1513, 1992
        NPSAT = 9
        WRK1D_LOC(1) = 0.611213476e+3_wp
        WRK1D_LOC(2) = 0.444007856e+2_wp 
        WRK1D_LOC(3) = 0.143064234e+1_wp 
        WRK1D_LOC(4) = 0.264461437e-1_wp
        WRK1D_LOC(5) = 0.305930558e-3_wp
        WRK1D_LOC(6) = 0.196237241e-5_wp
        WRK1D_LOC(7) = 0.892344772e-8_wp
        WRK1D_LOC(8) = -0.373208410e-10_wp
        WRK1D_LOC(9) = 0.209339997e-13_wp

! going from powers of (T-T_ref) to T, with T_ref = 273.15 K
        TREF_LOC = 273.15_wp
        do ipsat = 1, NPSAT
            THERMO_PSAT(ipsat) = 0.0_wp
            do i = ipsat, NPSAT
                tmp1 = 1.0_wp
                do j = i - 1, i - ipsat + 1, -1
                    tmp1 = tmp1*real(j, wp)
                end do
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat) + &
                                     WRK1D_LOC(i)*TREF_LOC**(i - 1)*tmp1*(-1)**(i - ipsat)
            end do
            tmp2 = 1.0_wp
            do j = ipsat - 1, 1, -1
                tmp2 = tmp2*real(j, wp)
            end do
            THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/tmp2/TREF_LOC**(ipsat - 1)
        end do

    case default

    end select

! ###################################################################
! Final calculations
! ###################################################################
    WGHT_INV(:) = RGAS/WGHT(:)              ! Specific gas constants, J /kg /K

    ! Nondimensionalization
    ISPREF = 2                              ! Species 2 is taken as reference
    WREF = WGHT(ISPREF)                     ! kg /kmol
    TREF = 298.0_wp                          ! K
    CPREF = 0.0_wp                           ! J /kg /K
    do icp = NCP, 1, -1
        CPREF = CPREF*TREF + THERMO_AI(icp, 2, ISPREF)
    end do

    if (imixture /= MIXT_TYPE_NONE) then ! othewise, gama0 is read in dns.ini
        gama0 = CPREF*WREF/(CPREF*WREF - RGAS)  ! Specific heat ratio
    end if
    if (gama0 > 0.0_wp) GRATIO = (gama0 - 1.0_wp)/gama0 ! Value of R_0/(C_{p,0}W_0) is called GRATIO
    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
        MRATIO = 1.0_wp
    else
        MRATIO = gama0*mach*mach
    end if

    if (nondimensional) then
        ! Thermal equation of state
        WGHT_INV(:) = WGHT_INV(:)/WGHT_INV(ISPREF)    ! normalized gas constants (Inverse of molar masses)

        ! Caloric equations of state
        do is = 1, NSP
            do im = 1, 2
                THERMO_AI(6, im, is) = THERMO_AI(6, im, is)/(CPREF*TREF)    ! Formation enthalpy
                THERMO_AI(7, im, is) = THERMO_AI(7, im, is)/CPREF           ! Formation entropy
                do icp = 1, NCP
                    THERMO_AI(icp, im, is) = THERMO_AI(icp, im, is)*(TREF**(icp - 1))/CPREF
                end do
            end do
        end do
        THERMO_TLIM = THERMO_TLIM/TREF          ! Temperature limis for polynomial fits to cp

        ! Saturation vapor pressure
        PREF_LOC = 1e5_wp      ! Pa
        do ipsat = 1, NPSAT
            THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/PREF_LOC
            THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)*(TREF**(ipsat - 1))
        end do

    else
        MRATIO = 1.0_wp
        GRATIO = 1.0_wp

    end if

    if (imixture == MIXT_TYPE_NONE) then
        do is = 1, NSP
            THERMO_SPNAME(is) = ''
        end do
    end if

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
    call TLAB_WRITE_ASCII(lfile, 'Thermodynamic properties have been initialized.')
    do is = 1, NSP
        write (str, *) is; str = 'Setting Species'//trim(adjustl(str))//'='//trim(adjustl(THERMO_SPNAME(is)))
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

1010 format(A14, 1x, G_FORMAT_R)
1020 format(A15, 1x, G_FORMAT_R)

end subroutine THERMO_INITIALIZE

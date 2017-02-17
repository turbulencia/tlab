#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created (reactini.f)
!# 2006/05/01 - J.P. Mellado
!#              Extracted from reactini.f
!# 2007/05/09 - J.P. Mellado
!#              Air-Vapor mixture added
!# 2013/06/12 - A. de Lozar
!#              Linear AirWater added for stratocumulus
!# 
!########################################################################
!# DESCRIPTION
!#
!# Extracted from reactini.f to retain only the thermodynamics
!# initialization of the mixture and make it usuable for 
!# multispecies (e.g. RTI) as well.
!#
!# inb_scal          # of scalars transported during the simulation and saved
!# inb_scal_array    # of scalars in array z1 (normally = inb_scal)
!# NSP               # of species in the mixture (NSP>=inb_scal)
!#
!# The code handles reactive/non-reactive, multiple/single species:
!# 1. Reactive => Multispecies.
!# 2. Multispecies admits the case Y_i=f_i(Z), Z conserved scalar,
!#    which implies inb_scal=1.
!# 3. Reactive + general multispecies retains NSP-1 species in scalar
!#    array z1 (the last one is obtained by sum Y_i=1) and an 
!#    additional conserved scalar, i.e. inb_scal=NSP.
!# 4. Non-reactive + general multispecies retains NSP-1 species, w/o
!#    additional conserved scalar.
!#
!# Multispecies implies that reference T_0 in non-dimensionalization
!# is 298 K.
!# 
!# Saturation pressure implies that reference R_0 in non-dimensionalization
!# is such that reference pressure is 1 bar.
!#
!# The new reference value of gamma0 is calculated here based on the
!# reference species
!#
!########################################################################
SUBROUTINE THERMO_INITIALIZE
  
  USE THERMO_GLOBAL

  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL,    ONLY : inb_scal, inb_scal_array
  USE DNS_GLOBAL,    ONLY : damkohler
  IMPLICIT NONE

! -------------------------------------------------------------------
  TINTEGER icp,is,im,ipsat,i,j
  TINTEGER ISPREF
  TREAL CPREF, RREF
  TREAL tmp1, tmp2
  TREAL HREF(MAX_NSP), SREF(MAX_NSP), WRK1D_LOC(10), tloc
  CHARACTER*46 str

! ###################################################################
! Species 2 is taken as reference
  ISPREF = 2

! ###################################################################
! Species tags and molecular weights
!
! Molecular Weight in kg/kmol
! ###################################################################
  SELECT CASE ( imixture )
     
! -------------------------------------------------------------------
! Burke-Schuman case
! Transport just mixture fraction, and then equilibrium
! 4 species + Nitrogen + Conserved Scalar         
! -------------------------------------------------------------------
  CASE( MIXT_TYPE_BS, MIXT_TYPE_QUASIBS )
     NSP            = 5
     inb_scal       = 1   
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
  CASE( MIXT_TYPE_PETERS1991, MIXT_TYPE_PETERS1988 )
     NSP            = 8
     inb_scal       = NSP
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
  CASE( MIXT_TYPE_UNIDECOMP )
     NSP            = 2
     inb_scal       = NSP
     inb_scal_array = inb_scal

     THERMO_SPNAME(1) = 'R'
     THERMO_SPNAME(2) = 'P'

     WGHT(1) = 32.0
     WGHT(2) = 32.0

! -------------------------------------------------------------------
! Unimolecular decomposition flame
! 2 reactant + 1 product + Conserved Scalar
! -------------------------------------------------------------------
  CASE( MIXT_TYPE_ONESTEP )
     NSP            = 4
     inb_scal       = NSP
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
  CASE( MIXT_TYPE_BILGER1997 )
     NSP            = 8
     inb_scal       = NSP-1
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
     WGHT(8) =  2.0

! -------------------------------------------------------------------
! Water vapor and air
! -------------------------------------------------------------------
  CASE( MIXT_TYPE_AIRVAPOR )
     NSP            = 2
     inb_scal       = NSP
     inb_scal_array = inb_scal

     THERMO_SPNAME(1) = 'H2O'
     THERMO_SPNAME(2) = 'AIR'

!     WGHT(1) = 18.01538 ! from Burcat&Ruscic
!     WGHT(2) = 28.96518
     WGHT(1) = 18.015    ! from B. Stevens
     WGHT(2) = 28.9644

! -------------------------------------------------------------------
! Water vapor, air and liquid water
! Compressible:   Transport    q_t, and q_l from equilibrium; add space for q_l
! Incompressible: Transport h, q_t, and q_l from equilibrium; add space for q_l
! If non-equilibrium calculation, then inb_scal includes q_l and inb_scal_array = inb_scal (default)
! -------------------------------------------------------------------
  CASE( MIXT_TYPE_AIRWATER )
     NSP            = 3
     IF ( damkohler(3) .LE. C_0_R ) inb_scal_array = inb_scal+1 ! using inb_scal read in the inifile

     THERMO_SPNAME(1) = 'H2Ov'
     THERMO_SPNAME(2) = 'AIR'
     THERMO_SPNAME(3) = 'H2Ol'

     WGHT(1) = 18.015
     WGHT(2) = 28.9644
     WGHT(3) = 18.015

! -------------------------------------------------------------------
! Linearized thermodynamics for stratocumulus case
! The 1. scalar is the scaled total water (mixing fraction)
! The 2. scalar is the enthalpy deviations from the pure mixing case (described by mixing fraction only)
! The 3. scalar is the normalized concentration of liquid.
! -------------------------------------------------------------------
  CASE( MIXT_TYPE_AIRWATER_LINEAR )
     inb_scal_array = inb_scal + 1 ! using inb_scal read in the inifile
     NSP            = inb_scal_array 

     THERMO_SPNAME(1) = 'Chi'      ! Mixture fraction 
     THERMO_SPNAME(2) = 'Psi'      ! Deviation in the enthalpy from the mixture fraction
     DO is = 3,inb_scal
        WRITE(THERMO_SPNAME(is),*) is; THERMO_SPNAME(is) = 'Scalar'//TRIM(ADJUSTL(THERMO_SPNAME(is))) 
     ENDDO
     THERMO_SPNAME(NSP) = 'Liquid' ! Normalized Liquid
     
     WGHT(1) = 18.015              ! unused, but defined for re-normalization below
     WGHT(2) = 28.9644
     WGHT(3) = 18.015
     WGHT(4:)= C_1_R
     
  END SELECT

! ###################################################################
! Thermodynamic data
! 
! Specific Heat Cpi, R^0 and SREF in Jules/(Kelvin Kmol), and
! enthalpy of formation HREF in Jules/Kmol.
!
! General formulation is CHEMKIN format: 7-coefficient NASA 
! polynomials (see Burcat&Ruscic):
!
! THERMO_AI(i,im,k) = a_i of species k at 
!         im=1 hight temperature
!         im=2 low   temperature
!
! C_{p,i} = \sum_1^5 a_i T^{i-1}
! h_i     = \sum_1^5 a_i T^i/i + a_6 
! s_{T,i} = a_1 ln(T) + \sum_2^5 a_i T^{i-1}/(i-1) + a_7
!
! i.e., dh_i = C_{p,i}dT and ds_{T,i}=C_{p,i}dT/T, where a_6 is 
! related to the formation enthalpy, and a_7 to the formation entropy.
! HREF and SREF (at TREF) are used to fix last Cpi 6-7 coefficients. 
!
! Note that Burcat&Ruscic give values devided by R^0
!
! The variables NCP_CHEMKIN gives the number of coefficients a_i.
! The simplified situation NCP_CHEMKIN=1 assumes also only one range,
! which then gives C_p constant. This is used to expedite the
! calculation of T in the energy formulation.
! (MAX_NCP is just the maximum number, 7, used to allocate the space
! in all the thermo arrays)
!
! The pressure contribution to the entropy still needs to be added
!
! Saturation pressure expansion is needed so far only on the 
! vapor case. It is
!
! p_sat(T) = \sum_1^9 a_i T^{i-1}
!
! ###################################################################
  TREF = 298.0e0
  RGAS = 8.314*C_1E3_R

  IF ( iuse_chemkin .EQ. 0 ) THEN
     SELECT CASE ( imixture )
! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR METHANE
! -------------------------------------------------------------------
     CASE( MIXT_TYPE_BS, MIXT_TYPE_QUASIBS )

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

! Heat capacity polynomial
! (Values taken from fit with T-TREF instead T)
        DO im = 1, 2
           THERMO_AI(1,im,1) = 35.70*C_1E3_R - 42.4833*TREF
           THERMO_AI(1,im,2) = 28.96*C_1E3_R - 6.21666*TREF
           THERMO_AI(1,im,3) = 32.76*C_1E3_R - 11.9570*TREF
           THERMO_AI(1,im,4) = 37.22*C_1E3_R - 17.6500*TREF
           THERMO_AI(1,im,5) = 28.88*C_1E3_R - 4.70833*TREF

           THERMO_AI(2,im,1) = 42.4833
           THERMO_AI(2,im,2) = 6.21666
           THERMO_AI(2,im,3) = 11.9570
           THERMO_AI(2,im,4) = 17.6500
           THERMO_AI(2,im,5) = 4.70833

! All other coefficients are zero
           DO is=1, NSP
              DO icp=3,MAX_NCP
                 THERMO_AI(icp,im,is) = C_0_R
              ENDDO
           ENDDO

! 6th and 7th coefficient are calculated from reference enthalpy
           DO is=1, NSP
              THERMO_AI(6,im,is) = HREF(is) &
                   - THERMO_AI(1,im,is)*TREF&
                   - THERMO_AI(2,im,is)*TREF*TREF*C_05_R
              THERMO_AI(7,im,is) = SREF(is) &
                   - THERMO_AI(2,im,is)*TREF
           ENDDO
        ENDDO

        DO is=1, NSP
           THERMO_TLIM(1,is) = 200.0e0
           THERMO_TLIM(2,is) = 5000.0e0
           THERMO_TLIM(3,is) = 5000.0e0
        ENDDO

        NCP_CHEMKIN = 2

! saturation pressure no needed
        DO ipsat = 1,MAX_SAT
           THERMO_PSAT(ipsat) = C_0_R 
        ENDDO
        NPSAT = 0

! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR SIMPLE REACTION
! -------------------------------------------------------------------
     CASE( MIXT_TYPE_UNIDECOMP )

! Enthalpy of Formation in Jules/Kmol
        HREF(1) = C_0_R
        HREF(2) = -86.71502*C_1E6_R

! Entropy of Formation in Jules/(Kelvin Kmol)
        SREF(1) = 205.15*C_1E3_R
        SREF(2) = 205.15*C_1E3_R

! Heat capacity polynomial
        DO im = 1, 2
           THERMO_AI(1,im,1) = 29.099*C_1E3_R
           THERMO_AI(1,im,2) = 29.099*C_1E3_R
! THERMO_AI(1,im,2) = 59.099*C_1E3_R

! All other coefficients are zero
           DO is=1, NSP
              DO icp=2,MAX_NCP
                 THERMO_AI(icp,im,is) = C_0_R
              ENDDO
           ENDDO

! 6th and 7th coefficient are calculated from reference enthalpy
           DO is=1, NSP
              THERMO_AI(6,im,is) = HREF(is) &
                   - THERMO_AI(1,im,is)*TREF
              THERMO_AI(7,im,is) = SREF(is) 
           ENDDO
        ENDDO

        DO is=1, NSP
           THERMO_TLIM(1,is) = 200.0e0
           THERMO_TLIM(2,is) = 5000.0e0
           THERMO_TLIM(3,is) = 5000.0e0
        ENDDO

        NCP_CHEMKIN = 1

! saturation pressure no needed
        DO ipsat = 1,MAX_SAT
           THERMO_PSAT(ipsat) = C_0_R 
        ENDDO
        NPSAT = 0

! -------------------------------------------------------------------
! FIXED THERMODYNAMIC DATA FOR SIMPLE REACTION
! -------------------------------------------------------------------
     CASE( MIXT_TYPE_ONESTEP )

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
        DO im = 1, 2
           THERMO_AI(1,im,1) = 29.099*C_1E3_R
           THERMO_AI(1,im,2) = 29.099*C_1E3_R
           THERMO_AI(1,im,3) = 29.099*C_1E3_R
           THERMO_AI(1,im,4) = 29.099*C_1E3_R

! All other coefficients are zero
           DO is=1, NSP
              DO icp=2,MAX_NCP
                 THERMO_AI(icp,im,is) = C_0_R
              ENDDO
           ENDDO

! 6th and 7th coefficient are calculated from reference enthalpy
           DO is=1, NSP
              THERMO_AI(6,im,is) = HREF(is) &
                   - THERMO_AI(1,im,is)*TREF
              THERMO_AI(7,im,is) = SREF(is) 
           ENDDO
        ENDDO

        DO is=1, NSP
           THERMO_TLIM(1,is) = 200.0e0
           THERMO_TLIM(2,is) = 5000.0e0
           THERMO_TLIM(3,is) = 5000.0e0
        ENDDO

        NCP_CHEMKIN = 1

! saturation pressure no needed
        DO ipsat = 1,MAX_SAT
           THERMO_PSAT(ipsat) = C_0_R 
        ENDDO
        NPSAT = 0

! -------------------------------------------------------------------
! Water vapor, air and water liquid
! -------------------------------------------------------------------
     CASE( MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR )

! Enthalpy of Formation in Jules/Kmol
!        HREF(1) =-241.826*C_1E6_R    ! from Burcat&Ruscic
!        HREF(2) =-0.126  *C_1E6_R  
        HREF(1) = 33.688*C_1E3_R*TREF ! from B. Stevens; values s.t. THERMO_AI(6,im,1:2) = 0
        HREF(2) = 29.167*C_1E3_R*TREF ! latent heat of vaporization at 298 K is (2501.6-58.690) kJ/kg
        HREF(3) = 33.688*C_1E3_R*TREF-44.009*C_1E6_R
!        HREF(1) = C_0_R              ! for testing
!        HREF(2) = C_0_R 
!        HREF(3) =-44.009*C_1E6_R

! Entropy of Formation in Jules/(Kelvin Kmol)
!        SREF(1) = 188.829*C_1E3_R ! from Burcat&Ruscic
!        SREF(2) = 198.824*C_1E3_R
        SREF(1) = C_0_R            ! from B. Stevens
        SREF(2) = C_0_R
        SREF(3) =-44.009*C_1E6_R/TREF 

! Heat capacity polynomial; using only the low temperature range
        DO im = 1, 2
!           THERMO_AI(1,im,1) = 34.907*C_1E3_R ! from Burcat&Ruscic
!           THERMO_AI(1,im,2) = 29.668*C_1E3_R
           THERMO_AI(1,im,1) = 33.688*C_1E3_R  ! from B. Stevens; values s.t. THERMO_AI(6,im,1:2) = 0
           THERMO_AI(1,im,2) = 29.167*C_1E3_R
           THERMO_AI(1,im,3) = 75.980*C_1E3_R

! All other coefficients are zero
           DO is=1, NSP
              DO icp=2,MAX_NCP
                 THERMO_AI(icp,im,is) = C_0_R
              ENDDO
           ENDDO

! 6th and 7th coefficient are calculated from reference enthalpy
           DO is=1, NSP
              THERMO_AI(6,im,is) = HREF(is) - THERMO_AI(1,im,is)*TREF
              THERMO_AI(7,im,is) = SREF(is) 
           ENDDO
        ENDDO

        DO is=1, NSP
           THERMO_TLIM(1,is) = 200.0e0
           THERMO_TLIM(2,is) = 5000.0e0
           THERMO_TLIM(3,is) = 5000.0e0
        ENDDO

        NCP_CHEMKIN = 1

! Saturation pressure; Flatau et al., J. Applied Meteorol., 1507-1513, 1992
        WRK1D_LOC(1) = 0.611213476*C_1E3_R
        WRK1D_LOC(2) = 0.444007856*C_1E2_R
        WRK1D_LOC(3) = 0.143064234*C_1E1_R
        WRK1D_LOC(4) = 0.264461437*C_1EM1_R
        WRK1D_LOC(5) = 0.305930558*C_1EM3_R
        WRK1D_LOC(6) = 0.196237241*C_1EM5_R
        WRK1D_LOC(7) = 0.892344772*C_1EM8_R
        WRK1D_LOC(8) =-0.373208410*C_1EM10_R
        WRK1D_LOC(9) = 0.209339997*C_1EM13_R

        NPSAT = 9

        DO ipsat = NPSAT+1,MAX_SAT
           THERMO_PSAT(ipsat) = C_0_R 
        ENDDO

! going from powers of (T-T_ref) to T, with T_ref = 273.15 K
        tloc = 273.15D0
        DO ipsat = 1, NPSAT
           THERMO_PSAT(ipsat) = C_0_R
           DO i = ipsat,NPSAT
              tmp1 = C_1_R
              DO j = i-1, i-ipsat+1,-1
                 tmp1 = tmp1*M_REAL(j)
              ENDDO
              THERMO_PSAT(ipsat) =  THERMO_PSAT(ipsat) +&
                   WRK1D_LOC(i)*tloc**(i-1)*tmp1*(-1)**(i-ipsat)
           ENDDO
           tmp2 = C_1_R
           DO j = ipsat-1,1,-1
              tmp2 = tmp2*M_REAL(j)
           ENDDO
           THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/tmp2/tloc**(ipsat-1)
        ENDDO

! -------------------------------------------------------------------
     CASE DEFAULT
        
        CALL IO_WRITE_ASCII(efile, 'THERMO_INITIALIZE: Must use chemkin data.')
        CALL DNS_STOP(DNS_ERROR_THERMOCONT)

     END SELECT

! -------------------------------------------------------------------
! ALTERNATIVE: LOAD THERMODYNAMIC DATA FROM CHEMKIN FILE
! -------------------------------------------------------------------
  ELSE
     CALL THERMO_READ_CHEMKIN(chemkin_file)

     DO is=1, NSP
        DO im=1,2
           DO icp = 1, MAX_NCP
              THERMO_AI(icp,im,is) = THERMO_AI(icp,im,is)*RGAS
           ENDDO
        ENDDO
     ENDDO

     NCP_CHEMKIN = 5

  ENDIF

! -------------------------------------------------------------------
! Combination of pure species
! Intermediate species is CO + H2 (species 5 and 8)
! -------------------------------------------------------------------
  IF ( imixture .EQ. MIXT_TYPE_BILGER1997 ) THEN
     WGHT(5) = C_05_R*( WGHT(5) + WGHT(8) )
     HREF(5) = C_05_R*( HREF(5) + HREF(8) )
     SREF(5) = C_05_R*( SREF(5) + SREF(8) )
     DO im=1,2
        DO icp=1,MAX_NCP
           THERMO_AI(icp,im,5) = C_05_R*( THERMO_AI(icp,im,5) + THERMO_AI(icp,im,8) )
        ENDDO
     ENDDO

  ENDIF

! ###################################################################
! Final calculations
!
! FROM THIS POINT ON, NO USE OF HREF() and SREF() SHOULD BE PERMITED, AND
! NORMALIZED REFERENCE TEMPERATURE IS 1.0
!
! ###################################################################
! -------------------------------------------------------------------
! Change heat capacities from molar to mass specific
! -------------------------------------------------------------------
  DO is = 1,NSP
     DO im = 1,2
        DO icp = 1,MAX_NCP
           THERMO_AI(icp,im,is) = THERMO_AI(icp,im,is)/WGHT(is)
        ENDDO
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! Compute reference enthalpy (in case chemkin has been used)
! -------------------------------------------------------------------
  DO is = 1,NSP
     HREF(is) = C_0_R
     DO icp = NCP_CHEMKIN, 1, -1
        HREF(is) = HREF(is)*TREF + THERMO_AI(icp,2,is)/M_REAL(icp)
     ENDDO
     HREF(is) = HREF(is)*TREF + THERMO_AI(6,2,is)
  ENDDO

! -------------------------------------------------------------------
! Nondimensional Forms
! -------------------------------------------------------------------
  DO is = 1,NSP
     THERMO_TLIM(1,is) = THERMO_TLIM(1,is)/TREF
     THERMO_TLIM(2,is) = THERMO_TLIM(2,is)/TREF
     THERMO_TLIM(3,is) = THERMO_TLIM(3,is)/TREF
  ENDDO

  WREF = WGHT(ISPREF)
  DO is = 1,NSP
     WGHT_INV(is) = WREF/WGHT(is)
  ENDDO

  CPREF = C_0_R
  DO icp=NCP_CHEMKIN, 1, -1
     CPREF = CPREF*TREF + THERMO_AI(icp,2,ISPREF)
  ENDDO

  gama0 = CPREF*WREF/(CPREF*WREF-RGAS)
! Value of R_0/(C_{p,0}W_0) is called GRATIO
  IF ( gama0 .GT. C_0_R ) GRATIO = (gama0-C_1_R)/gama0

  DO is = 1,NSP
     HREF(is) = HREF(is)/(CPREF*TREF)

     DO im = 1,2
        THERMO_AI(6,im,is) = THERMO_AI(6,im,is)/(CPREF*TREF)
        THERMO_AI(7,im,is) = THERMO_AI(7,im,is)/CPREF
        DO icp = 1,NCP_CHEMKIN
           THERMO_AI(icp,im,is) = THERMO_AI(icp,im,is)*(TREF**(icp-1))/CPREF
        ENDDO
     ENDDO

  ENDDO

! Saturation pressure; RREF is taken to have pressure of 1bar
  IF ( NPSAT .GT. 0 ) THEN
     RREF = C_1E5_R/(RGAS/WREF*TREF)
     DO ipsat = 1,NPSAT
        THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/(RREF*RGAS/WREF*TREF)
        THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)*(TREF**(ipsat-1))
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
  CALL IO_WRITE_ASCII(lfile, 'Thermodynamic properties have been initialized.')
  DO is = 1,NSP
     WRITE(str,*) is; str = 'Species'//TRIM(ADJUSTL(str))//'='//TRIM(ADJUSTL(THERMO_SPNAME(is)))
     CALL IO_WRITE_ASCII(lfile, str)
  ENDDO
  WRITE(str,'(A8,1X,E12.5E3)') 'WREF  = ', WREF
  CALL IO_WRITE_ASCII(lfile, str)
  WRITE(str,'(A8,1X,E12.5E3)') 'TREF  = ', TREF
  CALL IO_WRITE_ASCII(lfile, str)
  IF ( NPSAT .GT. 0 ) THEN
     WRITE(str,'(A8,1X,E12.5E3)') 'RREF  = ', RREF
     CALL IO_WRITE_ASCII(lfile, str)
  ENDIF
  WRITE(str,'(A8,1X,E12.5E3)') 'CPREF = ', CPREF
  CALL IO_WRITE_ASCII(lfile, str)
  WRITE(str,'(A8,1X,E12.5E3)') 'Gama0 = ', gama0
  CALL IO_WRITE_ASCII(lfile, str)
  CALL IO_WRITE_ASCII(lfile, '#')

  RETURN
END SUBROUTINE THERMO_INITIALIZE

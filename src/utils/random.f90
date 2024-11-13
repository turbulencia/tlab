!########################################################################
!#
!# Random number generators.
!# From Press, Flannery, Teukolsky & Vettering, "Numerical Recipes"
!#
!########################################################################

!########################################################################
! Gaussian pdf according to Box-Muller method
!########################################################################
function RANG(mean, sigma, seed)
    use TLab_Constants, only: wp, wi
    implicit none

    real(wp) mean, sigma
    integer(wi) seed
    real(wp) RANG

! -----------------------------------------------------------------------
    real(wp) v1, v2, r
    real(wp) RAN0

! #######################################################################
1   v1 = 2.0_wp*RAN0(seed) - 1.0_wp
    v2 = 2.0_wp*RAN0(seed) - 1.0_wp
    r = v1*v1 + v2*v2
    if (r >= 1.0_wp) GO TO 1
    v2 = v1*sqrt(-2.0_wp*log(r)/r)

    RANG = mean + v2*sigma

    return
end function RANG

! #######################################################################
! Uniform pdf
! #######################################################################
function RAN0(IDUM)
    use TLab_Constants, only: wp, wi

    implicit none

    real(wp) RAN0
    integer(wi) IDUM

    integer(wi) IA, IM, IQ, IR, NTAB, NDIV
    real(wp) AM, EPS, RNMX
    parameter(IA=16807, IM=2147483647, AM=1.0_wp/real(IM, wp), &
              IQ=127773, IR=2836)
    parameter(NTAB=32, NDIV=1 + (IM - 1)/NTAB, EPS=1.2e-7_wp, &
              RNMX=1.0_wp - EPS)

    integer(wi) j, k, iv(NTAB), iy
    save iv, iy
    data iv/NTAB*0/, iy/0/

    if (IDUM <= 0 .or. iy == 0) then
        IDUM = max(-IDUM, 1)
        do j = NTAB + 8, 1, -1
            k = IDUM/IQ
            IDUM = IA*(IDUM - k*IQ) - IR*k
            if (IDUM < 0) IDUM = IDUM + IM
            if (j <= NTAB) iv(j) = IDUM
        end do
        iy = iv(1)
    end if

    k = IDUM/IQ
    IDUM = IA*(IDUM - k*IQ) - IR*k
    if (IDUM < 0) IDUM = IDUM + IM
    j = 1 + iy/NDIV
    iy = iv(j)
    iv(j) = IDUM
    RAN0 = min(AM*iy, RNMX)

    return
end function RAN0

! #######################################################################
! Uniform pdf
! #######################################################################
function RAN1(IDUM)
    use TLab_Constants, only: wp, wi

    implicit none

    real(wp) RAN1
    integer(wi) IDUM

    integer(wi) M1, M2, M3
    integer(wi) IA1, IA2, IA3
    integer(wi) IC1, IC2, IC3
    real(wp) RM1, RM2
    real(wp) R(97)

    integer(wi) IX1, IX2, IX3, IFF, J

    parameter(M1=259200, IA1=7141, IC1=54773, RM1=1.0_wp/real(M1, wp))
    parameter(M2=134456, IA2=8121, IC2=28411, RM2=1.0_wp/real(M2, wp))
    parameter(M3=243000, IA3=4561, IC3=51349)
    data IFF/0/

    if (IDUM < 0 .or. IFF == 0) then
        IFF = 1
        IX1 = mod(IC1 - IDUM, M1)
        IX1 = mod(IA1*IX1 + IC1, M2)
        IX2 = mod(IX1, M2)
        IX1 = mod(IA1*IX1 + IC1, M1)
        IX3 = mod(IX1, M3)
        do J = 1, 97
            IX1 = mod(IA1*IX1 + IC1, M1)
            IX2 = mod(IA2*IX2 + IC2, M2)
            R(J) = (real(IX1, wp) + real(IX2, wp)*RM2)*RM1
        end do
        IDUM = 1
    end if
    IX1 = mod(IA1*IX1 + IC1, M1)
    IX2 = mod(IA2*IX2 + IC2, M2)
    IX3 = mod(IA3*IX3 + IC2, M3)
    J = 1 + (97*IX3)/M3
    if (J > 97 .or. J < 1) then
        write (*, *) 'Error J out of limits in RAN1'
        stop
    end if
    RAN1 = R(J)
    R(J) = (real(IX1, wp) + real(IX2, wp)*RM2)*RM1

    return
end function RAN1

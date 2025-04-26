subroutine SL_CORRELATION_1(ilog, u, v, w, z1, corr, &
                            strain, vorticity, gradient, tmp1, tmp2, wrk1d, wrk2d, wrk3d)
                            use TLab_Constants, only: wp, wi

    use TLAB_VARS
    use OPR_Partial
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN

    implicit none

    integer(wi) ilog
    real(wp) u(imax, jmax, kmax)
    real(wp) v(imax, jmax, kmax)
    real(wp) w(imax, jmax, kmax)
    real(wp) z1(imax, jmax, kmax)
    real(wp) strain(imax, jmax, kmax)
    real(wp) vorticity(imax, jmax, kmax)
    real(wp) gradient(imax, jmax, kmax)

    real(wp) tmp1(imax, jmax, kmax)
    real(wp) tmp2(imax, jmax, kmax)

    real(wp) corr(jmax, *)

    real(wp) wrk1d(jmax, *)
    real(wp) wrk2d(imax, kmax, *)
    real(wp) wrk3d(imax, jmax, kmax)

! -------------------------------------------------------------------
    integer(wi) j, bcs(2, 2)
    real(wp) mean_1, mean_2, var_1, var_2, delta_w
    real(wp) AVG1V2D, COV2V2D

    character*32 fname
    character*250 line1
    character*550 line2

! ###################################################################
    bcs = 0

    if (delta_u == C_0_R) then
        delta_w = C_1_R
    else
        do j = 1, jmax
            wrk1d(j, 1) = AVG1V2D(imax, jmax, kmax, j, i1, u)
        end do
        call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), wrk1d(1, 1), wrk1d(1, 2))
        delta_w = delta_u/maxval(abs(wrk1d(1:jmax, 2)))
    end if

! ###################################################################
! Define fields
! ###################################################################
    call FI_STRAIN(imax, jmax, kmax, u, v, w, strain, tmp1, tmp2)
    call FI_VORTICITY(imax, jmax, kmax, u, v, w, vorticity, tmp1, tmp2)
    call FI_GRADIENT(imax, jmax, kmax, z1, gradient, tmp1)

    if (ilog == 1) then
        strain = log(strain)
        vorticity = log(vorticity)
        gradient = log(gradient)
    end if

! ###################################################################
! Compute plane correlations
! ###################################################################
    do j = 1, jmax

        corr(j, 1) = COV2V2D(imax, jmax, kmax, j, vorticity, strain)
        mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, vorticity)
        var_1 = AVG1V2D(imax, jmax, kmax, j, i2, vorticity) - mean_1*mean_1
        mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, strain)
        var_2 = AVG1V2D(imax, jmax, kmax, j, i2, strain) - mean_2*mean_2
        if (var_1 > C_0_R .and. var_2 > C_0_R) then
            corr(j, 1) = (corr(j, 1) - mean_1*mean_2)/sqrt(var_1*var_2)
        else
            corr(j, 1) = C_2_R
        end if

        corr(j, 2) = COV2V2D(imax, jmax, kmax, j, vorticity, gradient)
        mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, vorticity)
        var_1 = AVG1V2D(imax, jmax, kmax, j, i2, vorticity) - mean_1*mean_1
        mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, gradient)
        var_2 = AVG1V2D(imax, jmax, kmax, j, i2, gradient) - mean_2*mean_2
        if (var_1 > C_0_R .and. var_2 > C_0_R) then
            corr(j, 2) = (corr(j, 2) - mean_1*mean_2)/sqrt(var_1*var_2)
        else
            corr(j, 2) = C_2_R
        end if

        corr(j, 3) = COV2V2D(imax, jmax, kmax, j, gradient, strain)
        mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, gradient)
        var_1 = AVG1V2D(imax, jmax, kmax, j, i2, gradient) - mean_1*mean_1
        mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, strain)
        var_2 = AVG1V2D(imax, jmax, kmax, j, i2, strain) - mean_2*mean_2
        if (var_1 > C_0_R .and. var_2 > C_0_R) then
            corr(j, 3) = (corr(j, 3) - mean_1*mean_2)/sqrt(var_1*var_2)
        else
            corr(j, 3) = C_2_R
        end if

    end do

! ###################################################################
! TkStat Output
! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        write (fname, *) itime; fname = 'slc'//trim(adjustl(fname))

        open (UNIT=23, FILE=fname, STATUS='unknown')
        write (23, '(A8,E14.7E3)') 'RTIME = ', rtime
        write (23, '(A8)') 'IMAX = 1'
        write (23, '(A7,I8)') 'JMAX = ', jmax

! -------------------------------------------------------------------
! Header
! -------------------------------------------------------------------
! Independent variables
        line2 = 'I J Y SW'

! Main fields
        line1 = 'W-S W-G S-G'
        write (i23, 1010) 'GROUP = MainFields '//trim(adjustl(line1))
        line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

        write (i23, 1010) trim(adjustl(line2))

1010    format(A)

! -------------------------------------------------------------------
! Body
! -------------------------------------------------------------------
        do j = 1, jmax
            write (23, 1020) 1, j, g(2)%nodes(j), (g(2)%nodes(j) - qbg(1)%ymean)/delta_w, &
                corr(j, 1), corr(j, 2), corr(j, 3)
        end do

1020    format(I3, 1x, I3, 2(1x, E12.5E3), 3(1x, E12.5E3))

#ifdef USE_MPI
    end if
#endif

    return
end subroutine SL_CORRELATION_1

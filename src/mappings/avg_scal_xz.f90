#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array tmp3 until pressure block
!#
!# Reynolds and Favre averages
!#
!########################################################################

subroutine AVG_SCAL_XZ(is, q, s, s_local, dsdx, dsdy, dsdz, tmp1, tmp2, tmp3, mean2d)
    use TLAB_CONSTANTS, only: MAX_AVG_TEMPORAL
    use TLAB_CONSTANTS, only: efile, lfile, wp, wi
    use TLAB_VARS
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_POINTERS_3D, only: p_wrk3d
    use THERMO_VARS, only: imixture, thermo_param
    use THERMO_ANELASTIC
    use THERMO_AIRWATER
    use IBM_VARS, only: gamma_0, gamma_1, scal_bcs
    use AVGS, only: AVG_IK_V
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif
    use TLAB_PROCS
    use FI_SOURCES, only: bbackground, FI_BUOYANCY, FI_BUOYANCY_SOURCE, FI_TRANSPORT, FI_TRANSPORT_FLUX
    use FI_GRADIENT_EQN
    use OPR_PARTIAL
    
    implicit none

    integer, intent(IN) :: is
    real(wp), intent(IN) :: q(imax, jmax, kmax, inb_flow_array)
    real(wp), intent(IN) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), intent(IN) :: s_local(imax, jmax, kmax)
    real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: dsdx, dsdy, dsdz, tmp1, tmp2, tmp3
    real(wp), intent(INOUT) :: mean2d(jmax,MAX_AVG_TEMPORAL)

    target q, tmp3

    ! -----------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 10
    integer j, k, bcs(2, 2), is_loc
    real(wp) diff, dummy, coefT, coefR, coefQ, c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv, im

    character*32 name, groupname(MAX_VARS_GROUPS)
    character*250 line1, varname(MAX_VARS_GROUPS)

    ! Pointers to existing allocated space
    real(wp), dimension(:, :, :), pointer :: u, v, w, p, rho, vis

    ! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    ! Define pointers
    u => q(:, :, :, 1)
    v => q(:, :, :, 2)
    w => q(:, :, :, 3)
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
        rho => q(:, :, :, 5)
        p => q(:, :, :, 6)
        if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) vis => q(:, :, :, 8)
    else
        p => tmp3
    end if

    if (idiffusion == EQNS_NONE) then; diff = 0.0_wp
    else; diff = visc/schmidt(is)
    end if

    c23 = 2.0_wp/3.0_wp

    ! -----------------------------------------------------------------------
    ! Dependent variables
    ng = 1; ig(ng) = 1
#define rS(j)     mean2d(j,ig(1)  )
#define fS(j)     mean2d(j,ig(1)+1)
#define rS_y(j)   mean2d(j,ig(1)+2)
#define fS_y(j)   mean2d(j,ig(1)+3)
#define rQ(j)     mean2d(j,ig(1)+4)
#define fQ(j)     mean2d(j,ig(1)+5)
    sg(ng) = 6

    groupname(ng) = 'Mean'
    varname(ng) = 'rS fS rS_y fS_y rQ fQ'
    if (imode_ibm == 1) then
        varname(ng) = trim(adjustl(varname(ng)))//' eps_0 eps_1 Sbcs'
#define ep_0(j)   mean2d(j,ig(1)+6)
#define ep_1(j)   mean2d(j,ig(1)+7)
#define Sbcs(j)   mean2d(j,ig(1)+8)
        sg(ng) = sg(ng) + 3
    end if
    if (radiation%active(is)) then
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR ) then
            varname(ng) = trim(adjustl(varname(ng)))//' rQrad rQradC'
        else
            varname(ng) = trim(adjustl(varname(ng)))//' rQrad rFrad'
        end if
        sg(ng) = sg(ng) + 2
    end if
    if (imixture == MIXT_TYPE_AIRWATER_LINEAR .or. imixture == MIXT_TYPE_AIRWATER) then
        varname(ng) = trim(adjustl(varname(ng)))//' rQeva'
        sg(ng) = sg(ng) + 1
    end if
    if (transport%active(is)) then
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR ) then
            varname(ng) = trim(adjustl(varname(ng)))//' rQtra rQtraC'
        else
            varname(ng) = trim(adjustl(varname(ng)))//' rQtra rFtra'
        end if
        sg(ng) = sg(ng) + 2
    end if

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu(j)    mean2d(j,ig(2)  )
#define Rsv(j)    mean2d(j,ig(2)+1)
#define Rsw(j)    mean2d(j,ig(2)+2)
#define fS2(j)    mean2d(j,ig(2)+3)
#define fS3(j)    mean2d(j,ig(2)+4)
#define fS4(j)    mean2d(j,ig(2)+5)
#define rS2(j)    mean2d(j,ig(2)+6)
#define rS3(j)    mean2d(j,ig(2)+7)
#define rS4(j)    mean2d(j,ig(2)+8)
    sg(ng) = 9

    groupname(ng) = 'Fluctuations'
    varname(ng) = 'Rsu Rsv Rsw fS2 fS3 fS4 rS2 rS3 rS4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rss_t(j)  mean2d(j,ig(3)  )
#define Css(j)    mean2d(j,ig(3)+1)
#define Pss(j)    mean2d(j,ig(3)+2)
#define Ess(j)    mean2d(j,ig(3)+3)
#define Tssy1(j)  mean2d(j,ig(3)+4)
#define Tssy2(j)  mean2d(j,ig(3)+5)
#define Tssy_y(j) mean2d(j,ig(3)+6)
#define Dss(j)    mean2d(j,ig(3)+7)
#define Qss(j)    mean2d(j,ig(3)+8)
    sg(ng) = 9

    groupname(ng) = 'RssBudget'
    varname(ng) = 'Rss_t Css Pss Ess Tssy1 Tssy2 Tssy_y Dss Qss'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu_t(j)  mean2d(j,ig(4)  )
#define Csu(j)    mean2d(j,ig(4)+1)
#define Psu(j)    mean2d(j,ig(4)+2)
#define Esu(j)    mean2d(j,ig(4)+3)
#define PIsu(j)   mean2d(j,ig(4)+4)
#define Tsuy1(j)  mean2d(j,ig(4)+5)
#define Tsuy2(j)  mean2d(j,ig(4)+6)
#define Tsuy_y(j) mean2d(j,ig(4)+7)
#define Dsu(j)    mean2d(j,ig(4)+8)
#define Gsu(j)    mean2d(j,ig(4)+9)
#define Bsu(j)    mean2d(j,ig(4)+10)
#define Fsu(j)    mean2d(j,ig(4)+11)
#define Qsu(j)    mean2d(j,ig(4)+12)
    sg(ng) = 13

    groupname(ng) = 'RsuBudget'
    varname(ng) = 'Rsu_t Csu Psu Esu PIsu Tsuy1 Tsuy2 Tsuy_y Dsu Gsu Bsu Fsu Qsu'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsv_t(j)  mean2d(j,ig(5)  )
#define Csv(j)    mean2d(j,ig(5)+1)
#define Psv(j)    mean2d(j,ig(5)+2)
#define Esv(j)    mean2d(j,ig(5)+3)
#define PIsv(j)   mean2d(j,ig(5)+4)
#define Tsvy1(j)  mean2d(j,ig(5)+5)
#define Tsvy2(j)  mean2d(j,ig(5)+6)
#define Tsvy3(j)  mean2d(j,ig(5)+7)
#define Tsvy_y(j) mean2d(j,ig(5)+8)
#define Dsv(j)    mean2d(j,ig(5)+9)
#define Gsv(j)    mean2d(j,ig(5)+10)
#define Bsv(j)    mean2d(j,ig(5)+11)
#define Fsv(j)    mean2d(j,ig(5)+12)
#define Qsv(j)    mean2d(j,ig(5)+13)
    sg(ng) = 14

    groupname(ng) = 'RsvBudget'
    varname(ng) = 'Rsv_t Csv Psv Esv PIsv Tsvy1 Tsvy2 Tsvy3 Tsvy_y Dsv Gsv Bsv Fsv Qsv'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsw_t(j)  mean2d(j,ig(6)  )
#define Csw(j)    mean2d(j,ig(6)+1)
#define Psw(j)    mean2d(j,ig(6)+2)
#define Esw(j)    mean2d(j,ig(6)+3)
#define PIsw(j)   mean2d(j,ig(6)+4)
#define Tswy1(j)  mean2d(j,ig(6)+5)
#define Tswy2(j)  mean2d(j,ig(6)+6)
#define Tswy_y(j) mean2d(j,ig(6)+7)
#define Dsw(j)    mean2d(j,ig(6)+8)
#define Gsw(j)    mean2d(j,ig(6)+9)
#define Bsw(j)    mean2d(j,ig(6)+10)
#define Fsw(j)    mean2d(j,ig(6)+11)
#define Qsw(j)    mean2d(j,ig(6)+12)
    sg(ng) = 13

    groupname(ng) = 'RswBudget'
    varname(ng) = 'Rsw_t Csw Psw Esw PIsw Tswy1 Tswy2 Tswy_y Dsw Gsw Bsw Fsw Qsw'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define S_x2(j)  mean2d(j,ig(7)  )
#define S_y2(j)  mean2d(j,ig(7)+1)
#define S_z2(j)  mean2d(j,ig(7)+2)
#define S_x3(j)  mean2d(j,ig(7)+3)
#define S_y3(j)  mean2d(j,ig(7)+4)
#define S_z3(j)  mean2d(j,ig(7)+5)
#define S_x4(j)  mean2d(j,ig(7)+6)
#define S_y4(j)  mean2d(j,ig(7)+7)
#define S_z4(j)  mean2d(j,ig(7)+8)
    sg(ng) = 9

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'S_x2 S_y2 S_z2 S_x3 S_y3 S_z3 S_x4 S_y4 S_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
    sg(ng) = 2*inb_scal_array

    groupname(ng) = 'CrossScalars'
    varname(ng) = ''
    do is_loc = 1, inb_scal_array
        write (name, *) is_loc
        varname(ng) = trim(adjustl(varname(ng)))//' Cs'//trim(adjustl(name))
        varname(ng) = trim(adjustl(varname(ng)))//' Css'//trim(adjustl(name))
    end do

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on y and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rR(j)       mean2d(j,ig(9))
#define rU(j)       mean2d(j,ig(9)+1)
#define rV(j)       mean2d(j,ig(9)+2)
#define rW(j)       mean2d(j,ig(9)+3)
#define fU(j)       mean2d(j,ig(9)+4)
#define fV(j)       mean2d(j,ig(9)+5)
#define fW(j)       mean2d(j,ig(9)+6)
#define fU_y(j)     mean2d(j,ig(9)+7)
#define fV_y(j)     mean2d(j,ig(9)+8)
#define fW_y(j)     mean2d(j,ig(9)+9)
#define rU_y(j)     mean2d(j,ig(9)+10)
#define rV_y(j)     mean2d(j,ig(9)+11)
#define rW_y(j)     mean2d(j,ig(9)+12)
#define Rvu(j)      mean2d(j,ig(9)+13)
#define Rvv(j)      mean2d(j,ig(9)+14)
#define Rvw(j)      mean2d(j,ig(9)+15)
#define Rss_y(j)    mean2d(j,ig(9)+16)
#define Rsu_y(j)    mean2d(j,ig(9)+17)
#define Rsv_y(j)    mean2d(j,ig(9)+18)
#define Rsw_y(j)    mean2d(j,ig(9)+19)
#define Fy(j)       mean2d(j,ig(9)+20)
#define Fy_y(j)     mean2d(j,ig(9)+21)
#define Tau_yy(j)   mean2d(j,ig(9)+22)
#define Tau_yy_y(j) mean2d(j,ig(9)+23)
#define Tau_yx(j)   mean2d(j,ig(9)+24)
#define Tau_yx_y(j) mean2d(j,ig(9)+25)
#define Tau_yz(j)   mean2d(j,ig(9)+26)
#define Tau_yz_y(j) mean2d(j,ig(9)+27)
#define rP(j)       mean2d(j,ig(9)+28)
#define aux(j)      mean2d(j,ig(9)+29)
    sg(ng) = 30

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLAB_WRITE_ASCII(efile, 'AVG_SCAL_XZ. Not enough space in local arrays.')
        call TLAB_STOP(DNS_ERROR_AVGTMP)
    end if
    mean2d(:, 1:nv) = 0.0_wp

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1 ! the last group is not written out

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating scal statistics at It'//trim(adjustl(line1))//'...'
    call TLAB_WRITE_ASCII(lfile, line1)

    ! #######################################################################
    ! Preliminary for IBM usage
    ! #######################################################################
    ! Asign gammas for conditional averages (c.f. Pope, p.170 [5.305])
    ! write out scalar boundary values applied in solids
    im = 5
    if (imode_ibm == 1) then
        ep_0(:) = gamma_0; ep_1(:) = gamma_1
        Sbcs(:) = scal_bcs(:, is)
        im = im + 3
    end if

    ! #######################################################################
    ! Preliminary data of velocity and density
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, jmax, u, g(1)%jac, g(3)%jac, rU(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, v, g(1)%jac, g(3)%jac, rV(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, w, g(1)%jac, g(3)%jac, rW(1), wrk1d, area)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        rR(:) = 1.0_wp

        fU(:) = rU(:)
        fV(:) = rV(:)
        fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, jmax, rho, g(1)%jac, g(3)%jac, rR(1), wrk1d, area)

        dsdx = rho*u
        dsdy = rho*v
        dsdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, fU(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, fV(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, fW(1), wrk1d, area)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rU(1), rU_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rV(1), rV_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rW(1), rW_y(1))

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fU(1), fU_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fV(1), fV_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fW(1), fW_y(1))

    dsdx = v*u
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) dsdx = dsdx*rho
    dsdy = v*v
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) dsdy = dsdy*rho
    dsdz = v*w
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) dsdz = dsdz*rho
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, Rvu(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, Rvv(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, Rvw(1), wrk1d, area)
    Rvu(:) = Rvu(:)/rR(:) - fV(:)*fU(:)
    Rvv(:) = Rvv(:)/rR(:) - fV(:)*fV(:)
    Rvw(:) = Rvw(:)/rR(:) - fV(:)*fW(:)

    ! #######################################################################
    ! Scalar
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, jmax, s_local, g(1)%jac, g(3)%jac, rS(1), wrk1d, area)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        fS(:) = rS(:)
    else
        p_wrk3d = rho*s_local
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, fS(1), wrk1d, area)
        fS(:) = fS(:)/rR(:)
    end if

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rS(1), rS_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fS(1), fS_y(1))

    ! -----------------------------------------------------------------------
    ! Moments
    do j = 1, jmax
        p_wrk3d(:, j, :) = s_local(:, j, :) - rS(j)
    end do
    tmp1 = p_wrk3d*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, rS2(1), wrk1d, area)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, rS3(1), wrk1d, area)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, rS4(1), wrk1d, area)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        fS2(:) = rS2(:)
        fS3(:) = rS3(:)
        fS4(:) = rS4(:)

    else
        do j = 1, jmax
            p_wrk3d(:, j, :) = s_local(:, j, :) - fS(j)
        end do
        tmp1 = p_wrk3d*p_wrk3d*rho
        call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, fS2(1), wrk1d, area)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, fS3(1), wrk1d, area)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, fS4(1), wrk1d, area)
        fS2(:) = fS2(:)/rR(:)
        fS3(:) = fS3(:)/rR(:)
        fS4(:) = fS4(:)/rR(:)

    end if

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fS2(1), Rss_y(1))

    ! -----------------------------------------------------------------------
    ! Cross terms
    do j = 1, jmax
        p_wrk3d(:, j, :) = s_local(:, j, :) - fS(j)
    end do
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) p_wrk3d = p_wrk3d*rho

    do j = 1, jmax
        dsdx(:, j, :) = p_wrk3d(:, j, :)*(u(:, j, :) - fU(j))
        dsdy(:, j, :) = p_wrk3d(:, j, :)*(v(:, j, :) - fV(j))
        dsdz(:, j, :) = p_wrk3d(:, j, :)*(w(:, j, :) - fW(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, Rsu(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, Rsv(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, Rsw(1), wrk1d, area)
    Rsu(:) = Rsu(:)/rR(:)
    Rsv(:) = Rsv(:)/rR(:)
    Rsw(:) = Rsw(:)/rR(:)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsu(1), Rsu_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsv(1), Rsv_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsw(1), Rsw_y(1))

    ! -----------------------------------------------------------------------
    ! turbulent transport terms
    do j = 1, jmax
        tmp1(:, j, :) = dsdy(:, j, :)*(s_local(:, j, :) - fS(j))
        dsdx(:, j, :) = dsdx(:, j, :)*(v(:, j, :) - fV(j))
        dsdy(:, j, :) = dsdy(:, j, :)*(v(:, j, :) - fV(j))
        dsdz(:, j, :) = dsdz(:, j, :)*(v(:, j, :) - fV(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, Tssy1(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, Tsuy1(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, Tsvy1(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, Tswy1(1), wrk1d, area)

    ! -----------------------------------------------------------------------
    ! Pressure terms in transport equations
    call AVG_IK_V(imax, jmax, kmax, jmax, p, g(1)%jac, g(3)%jac, rP(1), wrk1d, area)

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s_local, dsdx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s_local, dsdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s_local, dsdz)
    do j = 1, jmax
        tmp1(:, j, :) = (p(:, j, :) - rP(j))*(s_local(:, j, :) - fS(j))
        dsdx(:, j, :) = (p(:, j, :) - rP(j))*dsdx(:, j, :)
        dsdy(:, j, :) = (p(:, j, :) - rP(j))*(dsdy(:, j, :) - fS_y(j))
        dsdz(:, j, :) = (p(:, j, :) - rP(j))*dsdz(:, j, :)
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, Tsvy3(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, PIsu(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, PIsv(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, PIsw(1), wrk1d, area)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rP(1), aux(1))
    Gsv(:) = (rS(:) - fS(:))*aux(:)

    ! #######################################################################
    ! Cross-scalar terms
    ! #######################################################################
    k = ig(8) - 1

    do is_loc = 1, inb_scal_array
        call AVG_IK_V(imax, jmax, kmax, jmax, s(:,:,:, is_loc), g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
        do j = 1, jmax
            tmp1(:, j, :) = (s(:, j, :, is_loc) - aux(j))*(s_local(:, j, :) - fS(j))
            tmp2(:, j, :) = tmp1(:, j, :)*(s_local(:, j, :) - fS(j))
        end do
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area)
    end do

    ! #######################################################################
    ! Source terms
    ! #######################################################################
    dsdx = 0.0_wp; dsdy = 0.0_wp; dsdz = 0.0_wp; tmp1 = 0.0_wp; tmp2 = 0.0_wp; tmp3 = 0.0_wp

    if (radiation%active(is)) then ! Radiation in tmp1 and dsdx
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, s(1, 1, 1, radiation%scalar(is)), tmp2)
            call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), tmp2, tmp1)
            call OPR_RADIATION_FLUX(radiation, imax, jmax, kmax, g(2), tmp2, dsdx)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, ribackground, tmp1)
            tmp2 = 0.0_wp

        else
            call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(:,:,:, radiation%scalar(is)), tmp1)
            call OPR_RADIATION_FLUX(radiation, imax, jmax, kmax, g(2), s(:,:,:, radiation%scalar(is)), dsdx)
        end if
    end if

    if (transport%active(is)) then ! Transport in tmp3 and dsdz
        call FI_TRANSPORT(transport, 1, imax, jmax, kmax, is, s, tmp3, dsdy)
        call FI_TRANSPORT_FLUX(transport, imax, jmax, kmax, is, s, dsdz)
    end if

    if (is > inb_scal) then     ! Diagnostic variables; I overwrite tmp1 and dsdx and recalculate them.
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            coefQ = 1.0_wp                        ! Coefficient in the evaporation term
            coefR = 0.0_wp                        ! Coefficient in the radiation term
            coefT = 0.0_wp                        ! Coefficient in the transport term
            if (is == inb_scal_array + 1) then ! Default values are for liquid; defining them for buoyancy
                coefQ = buoyancy%parameters(inb_scal_array)/froude
                coefR = buoyancy%parameters(inb_scal)/froude
                if ( transport%active(is) ) then
                    do is_loc = 1, inb_scal
                        coefT = coefT + transport%parameters(is_loc)/settling*buoyancy%parameters(is_loc)/froude
                    end do
                end if
            end if

            call THERMO_AIRWATER_LINEAR_SOURCE(imax*jmax*kmax, s, dsdx, dsdy, dsdz) ! calculate xi in dsdx
            call FI_GRADIENT(imax, jmax, kmax, dsdx, tmp2, tmp1)

            dummy = -diff*coefQ
            tmp2 = dsdz*tmp2*dummy         ! evaporation source

            if (transport%active(is) .or. radiation%active(is)) then ! preparing correction terms into dsdz
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), dsdx, tmp1)
                dsdz = dsdz*tmp1
            end if

            if (radiation%active(is)) then ! radiation source; needs dsdy
                call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(:,:,:, radiation%scalar(is)), tmp1)
                dummy = thermo_param(2)*coefQ
                tmp1 = tmp1*(coefR + dsdy*dummy)
                ! Correction term needs dsdz
                call OPR_RADIATION_FLUX(radiation, imax, jmax, kmax, g(2), s(:,:,:,radiation%scalar(is)), dsdx)
                dsdx = dsdx*dsdz*dummy
            else
                tmp1 = 0.0_wp; dsdx = 0.0_wp
            end if

            if (transport%active(is)) then ! transport source; needs dsdy
                dummy = coefQ
                tmp3 = tmp3*(coefT + dsdy*dummy)
                ! Correction term needs dsdz
                call FI_TRANSPORT_FLUX(transport, imax, jmax, kmax, is, s, dsdy)
                dsdz = dsdy*dsdz*dummy
            else
                tmp3 = 0.0_wp; dsdz = 0.0_wp
            end if

        else
            if (buoyancy%type /= EQNS_EXPLICIT) then
                call FI_GRADIENT(imax, jmax, kmax, s, dsdx, dsdy)
                call FI_BUOYANCY_SOURCE(buoyancy, imax, jmax, kmax, s, dsdx, tmp1) ! dsdx contains gradient
                tmp1 = tmp1*diff/froude
            end if

        end if

    end if

    ! -----------------------------------------------------------------------
    ! Calculating averages
    k = ig(1) + im
    if (radiation%active(is)) then
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area) ! correction term or flux
    end if
    if (imixture == MIXT_TYPE_AIRWATER_LINEAR .or. imixture == MIXT_TYPE_AIRWATER) then           ! evaporation
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area)
    end if
    if (transport%active(is)) then
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, tmp3, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, mean2d(1, k), wrk1d, area) ! correction term or flux
    end if

    p_wrk3d = tmp1 + tmp2 + tmp3 ! total
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rQ(1), wrk1d, area)
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) p_wrk3d = p_wrk3d*rho
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, fQ(1), wrk1d, area)
    fQ(:) = fQ(:)/rR(:)

    do j = 1, jmax
        tmp1(:, j, :) = (s_local(:, j, :) - fS(j))*p_wrk3d(:, j, :)
        dsdx(:, j, :) = (u(:, j, :) - fU(j))*p_wrk3d(:, j, :)
        dsdy(:, j, :) = (v(:, j, :) - fV(j))*p_wrk3d(:, j, :)
        dsdz(:, j, :) = (w(:, j, :) - fW(j))*p_wrk3d(:, j, :)
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, Qss(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdx, g(1)%jac, g(3)%jac, Qsu(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, Qsv(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdz, g(1)%jac, g(3)%jac, Qsw(1), wrk1d, area)
    Qss(:) = Qss(:)*2.0_wp/rR(:)
    Qsu(:) = Qsu(:)/rR(:)
    Qsv(:) = Qsv(:)/rR(:)
    Qsw(:) = Qsw(:)/rR(:)

    ! #######################################################################
    ! Derivatives
    ! #######################################################################
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s_local, dsdx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s_local, dsdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s_local, dsdz)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*dsdx + dsdy*dsdy + dsdz*dsdz
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Ess(1), wrk1d, area)
    Ess(:) = Ess(:)*diff*2.0_wp

    ! -----------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)

    ! Transport term
    p_wrk3d = (tmp2*2.0_wp - tmp1 - tmp3)*c23*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_yy(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yy(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tsvy2(1), wrk1d, area)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yy(1), Tau_yy_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*((tmp1*2.0_wp - tmp2 - tmp3)*c23*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Esu(1), wrk1d, area)

    p_wrk3d = dsdy*((tmp2*2.0_wp - tmp1 - tmp3)*c23*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Esv(1), wrk1d, area)

    p_wrk3d = dsdz*((tmp3*2.0_wp - tmp1 - tmp2)*c23*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Esw(1), wrk1d, area)

    ! -----------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp2)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp1)

    ! Transport term
    p_wrk3d = (tmp1 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_yx(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yx(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tsuy2(1), wrk1d, area)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yx(1), Tau_yx_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdy*((tmp1 + tmp2)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp1 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esv(:) = Esv(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp2)

    ! Transport term
    p_wrk3d = (tmp3 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_yz(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yz(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tswy2(1), wrk1d, area)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yz(1), Tau_yz_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp3 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esv(:) = Esv(:) + aux(:)

    p_wrk3d = dsdy*((tmp3 + tmp2)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp3)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp1)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp3 + tmp1)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp3 + tmp1)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    ! Moments
    do j = 1, jmax
        p_wrk3d(:, j, :) = dsdy(:, j, :) - rS_y(j)
    end do

    tmp1 = dsdx*dsdx
    tmp2 = p_wrk3d*p_wrk3d
    tmp3 = dsdz*dsdz
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, S_x2(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, S_y2(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp3, g(1)%jac, g(3)%jac, S_z2(1), wrk1d, area)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*p_wrk3d
    tmp3 = tmp3*dsdz
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, S_x3(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, S_y3(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp3, g(1)%jac, g(3)%jac, S_z3(1), wrk1d, area)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*p_wrk3d
    tmp3 = tmp3*dsdz
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, S_x4(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, S_y4(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp3, g(1)%jac, g(3)%jac, S_z4(1), wrk1d, area)

    ! -----------------------------------------------------------------------
    ! Molecular fluxes
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) dsdy = dsdy*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, dsdy, g(1)%jac, g(3)%jac, Fy(1), wrk1d, area)

    ! Contribution to turbulent transport
    do j = 1, jmax
        p_wrk3d(:, j, :) = (dsdy(:, j, :) - Fy(j))*(s_local(:, j, :) - fS(j))
        tmp1(:, j, :) = (dsdy(:, j, :) - Fy(j))*(u(:, j, :) - fU(j))
        tmp2(:, j, :) = (dsdy(:, j, :) - Fy(j))*(v(:, j, :) - fV(j))
        tmp3(:, j, :) = (dsdy(:, j, :) - Fy(j))*(w(:, j, :) - fW(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tssy2(1), wrk1d, area)
    Tssy2(:) = -Tssy2(:)*diff*2.0_wp
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp1, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Tsuy2(:) = Tsuy2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp2, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Tsvy2(:) = Tsvy2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, jmax, tmp3, g(1)%jac, g(3)%jac, aux(1), wrk1d, area)
    Tswy2(:) = Tswy2(:) - aux(:)*diff

    Fy(:) = Fy(:)*diff
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Fy(1), Fy_y(1))

    ! Contribution to dissipation
    Ess(:) = (Ess(:) - Fy(:)*rS_y(:) - Fy(:)*rS_y(:))/rR(:)
    Esu(:) = (Esu(:) - Tau_yx(:)*rS_y(:) - Fy(:)*rU_y(:))/rR(:)
    Esv(:) = (Esv(:) - Tau_yy(:)*rS_y(:) - Fy(:)*rV_y(:))/rR(:)
    Esw(:) = (Esw(:) - Tau_yz(:)*rS_y(:) - Fy(:)*rW_y(:))/rR(:)

    ! #######################################################################
    ! Source terms in transport equations
    ! #######################################################################
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        if (buoyancy%type == EQNS_EXPLICIT) then
            call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, p_wrk3d)
        else
            call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, p_wrk3d, bbackground)
        end if
        dummy = 1.0_wp/froude
        p_wrk3d = p_wrk3d*dummy
    else
        p_wrk3d = rho*buoyancy%vector(2)
    end if
    do j = 1, jmax
        p_wrk3d(:, j, :) = (s_local(:, j, :) - fS(j))*p_wrk3d(:, j, :)
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Bsv(1), wrk1d, area)
    Bsv(:) = Bsv(:)/rR(:)

    ! #######################################################################
    ! Complete budget equations
    ! #######################################################################
    ! Transport terms
    aux(:) = Tssy1(:) + Tssy2(:)
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tssy_y(1))
    aux(:) = Tsuy1(:) + Tsuy2(:)
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tsuy_y(1))
    aux(:) = Tsvy1(:) + Tsvy2(:) + Tsvy3(:)
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tsvy_y(1))
    aux(:) = Tswy1(:) + Tswy2(:)
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tswy_y(1))

    ! Convective terms
    Css(:) = -fV(:)*Rss_y(:)
    Csu(:) = -fV(:)*Rsu_y(:)
    Csv(:) = -fV(:)*Rsv_y(:)
    Csw(:) = -fV(:)*Rsw_y(:)

    ! Production terms
    Pss(:) = -Rsv(:)*fS_y(:)*2.0_wp
    Psu(:) = -Rsv(:)*fU_y(:) - Rvu(:)*fS_y(:)
    Psv(:) = -Rsv(:)*fV_y(:) - Rvv(:)*fS_y(:)
    Psw(:) = -Rsv(:)*fW_y(:) - Rvw(:)*fS_y(:)

    ! Diffusion variable-density terms
    Dss(:) = (rS(:) - fS(:))*Fy_y(:)*2.0_wp
    Dsu(:) = (rS(:) - fS(:))*Tau_yx_y(:) + (rU(:) - fU(:))*Fy_y(:)
    Dsv(:) = (rS(:) - fS(:))*Tau_yy_y(:) + (rV(:) - fV(:))*Fy_y(:)
    Dsw(:) = (rS(:) - fS(:))*Tau_yz_y(:) + (rW(:) - fW(:))*Fy_y(:)

    ! Coriolis terms
    dummy = coriolis%vector(2)
    Fsu(:) = dummy*Rsw(:)
    Fsw(:) = -dummy*Rsu(:)

    ! Transient terms
    Rss_t(:) = Css(:) + Pss(:) - Ess(:) + Qss(:) + (Dss(:) - Tssy_y(:))/rR(:)
    Rsu_t(:) = Csu(:) + Psu(:) - Esu(:) + Bsu(:) - Fsu(:) + Qsu(:) + (PIsu(:) + Dsu(:) - Gsu(:) - Tsuy_y(:))/rR(:)
    Rsv_t(:) = Csv(:) + Psv(:) - Esv(:) + Bsv(:) - Fsv(:) + Qsv(:) + (PIsv(:) + Dsv(:) - Gsv(:) - Tsvy_y(:))/rR(:)
    Rsw_t(:) = Csw(:) + Psw(:) - Esw(:) + Bsw(:) - Fsw(:) + Qsw(:) + (PIsw(:) + Dsw(:) - Gsw(:) - Tswy_y(:))/rR(:)

    ! ###################################################################
    ! Output
    ! #######################################################################
    ! 11 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (line1, *) is; line1 = 'avg'//trim(adjustl(line1))//'s'
    write (name, *) itime; name = trim(adjustl(line1))//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, jmax, nv, ng, g(2)%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_SCAL_XZ

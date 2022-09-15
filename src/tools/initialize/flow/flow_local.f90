#include "types.h"
#include "dns_const.h"

module FLOW_LOCAL
    use TLAB_TYPES, only: profiles_dt, discrete_dt, cp, ci
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g, qbg
    use IO_FIELDS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none
    save
    ! -------------------------------------------------------------------
    integer(ci) :: flag_u, flag_t, flag_dilatation, flag_mixture

    type(profiles_dt) :: Kini                 ! Geometry of perturbation of initial boundary condition
    real(cp) :: norm_ini_u, norm_ini_p          ! Scaling of perturbation
    type(discrete_dt) :: fp                     ! Discrete perturbation

    integer(ci) :: flag_wall ! Boundary conditions: 0  Free-Slip/Free-Slip
    ! 1  No-Slip/Free-Slip
    ! 2  Free-Slip/No-Slip
    ! 3  No-Slip/No-Slip
    ! -------------------------------------------------------------------
    integer(ci) i, j, k

    integer(ci) im, idsp, kdsp
    real(cp) wx, wz, wx_1, wz_1
    type(profiles_dt) prof_loc

    real(cp), dimension(:), pointer :: xn, zn

contains

    ! ###################################################################
    subroutine FLOW_SHAPE(wrk1d)
        real(cp), dimension(jmax, 5), intent(INOUT) :: wrk1d

        ! -------------------------------------------------------------------
        integer(ci) bcs(2, 2)
        real(cp) PROFILES, yr
        external PROFILES

        real(cp), dimension(:), pointer :: yn

        ! ###################################################################
        bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

        yn => g(2)%nodes
        prof_loc= Kini
        prof_loc%delta=C_1_R
        prof_loc%mean=C_0_R
        do j = 1, jmax                               ! Wall-normal velocity
            wrk1d(j, 1) = PROFILES(prof_loc, yn(j))
        end do
        call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
        wrk1d(:, 2) = -wrk1d(:, 2)                     ! Negative of the derivative of f, wall-parallel velocity

        select case (Kini%type)
        case (PROFILE_PARABOLIC_SURFACE)
            ! Zero wall-parallel velocity for no-slip condition, multiply by parabolic again, f=f*f
            wrk1d(:, 2) = C_2_R*wrk1d(:, 2)*wrk1d(:, 1)          ! Wall-parallel velocity
            wrk1d(:, 1) = wrk1d(:, 1)**C_2_R                     ! Wall-normal velocity

        case (PROFILE_GAUSSIAN_SURFACE)
            ! Zero wall-normal derivative of wall-parallel velocity for free-slip and potentialvelocity mode, f=f*tanh
            if (flag_wall == 1 .or. flag_wall == 3) then  ! jmin
                do j = 1, jmax
                    yr = C_05_R*(yn(j) - yn(1))/Kini%thick
                    wrk1d(j, 2) = wrk1d(j, 2)*TANH(yr)**2 - &       ! Wall-parallel velocity
                                  wrk1d(j, 1)*TANH(yr)/COSH(yr)**2/Kini%thick
                    wrk1d(j, 1) = wrk1d(j, 1)*TANH(yr)**2           ! Wall-normal velocity
                end do
            end if

            if (flag_wall == 2 .or. flag_wall == 3) then  ! jmax
                do j = 1, jmax
                    yr = C_05_R*(yn(jmax) - yn(j))/Kini%thick
                    wrk1d(j, 2) = wrk1d(j, 2)*TANH(yr)**2 + &       ! Wall-parallel velocity
                                  wrk1d(j, 1)*TANH(yr)/COSH(yr)**2/Kini%thick
                    wrk1d(j, 1) = wrk1d(j, 1)*TANH(yr)**2           ! Wall-normal velocity
                end do
            end if

        end select

        return
    end subroutine FLOW_SHAPE

    ! ###################################################################
    subroutine VELOCITY_DISCRETE(u, v, w, wrk1d, wrk2d)
        real(cp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w
        real(cp), dimension(jmax, 2), intent(INOUT) :: wrk1d
        real(cp), dimension(imax, kmax, 3), intent(INOUT) :: wrk2d

        ! -------------------------------------------------------------------
        real(cp) factorx, factorz

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        xn => g(1)%nodes
        zn => g(3)%nodes

        call FLOW_SHAPE(wrk1d)

        wx_1 = C_2_R*C_PI_R/g(1)%scale ! Fundamental wavelengths
        wz_1 = C_2_R*C_PI_R/g(3)%scale

        wrk2d = C_0_R
        do im = 1, fp%size
            wx = M_REAL(fp%modex(im))*wx_1
            wz = M_REAL(fp%modez(im))*wz_1

            ! Factor to impose solenoidal constraint
            if (fp%modex(im) == 0 .and. fp%modez(im) == 0) then; exit
            elseif (fp%modez(im) == 0) then; factorx = C_1_R/wx; factorz = C_0_R
            elseif (fp%modex(im) == 0) then; factorx = C_0_R; factorz = C_1_R/wz
            else; factorx = C_05_R/wx; factorz = C_05_R/wz
            end if

            do k = 1, kmax
wrk2d(:,k,2) = wrk2d(:,k,2) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) )
          wrk2d(:,k,1) = wrk2d(:,k,1) + fp%amplitude(im) *SIN( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) ) *factorx
          wrk2d(:,k,3) = wrk2d(:,k,3) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *SIN( wz *zn(kdsp+k) +fp%phasez(im) ) *factorz
            end do

        end do

        do k = 1, kmax
            do j = 1, jmax
                u(:, j, k) = wrk2d(:, k, 1)*wrk1d(j, 2)
                v(:, j, k) = wrk2d(:, k, 2)*wrk1d(j, 1)
                w(:, j, k) = wrk2d(:, k, 3)*wrk1d(j, 2)
            end do
        end do

        if (norm_ini_u >= C_0_R) call FLOW_NORMALIZE(u, v, w)

        return
    end subroutine VELOCITY_DISCRETE

    ! ###################################################################
    subroutine VELOCITY_BROADBAND(u, v, w, ax, ay, az, tmp4, tmp5, wrk1d, wrk2d, wrk3d)
        use TLAB_VARS, only: visc

        real(cp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w
        real(cp), dimension(imax, jmax, kmax), intent(INOUT) :: ax, ay, az, tmp4, tmp5, wrk3d
        real(cp), dimension(imax, kmax, *), intent(INOUT) :: wrk2d
        real(cp), dimension(jmax, *), intent(INOUT) :: wrk1d

        ! -------------------------------------------------------------------
        integer(ci) ibc, bcs(2, 2), bcs2(2, 2)
        real(cp) AVG1V2D, dummy
        external AVG1V2D

        ! ###################################################################
        bcs = 0

        call FLOW_SHAPE(wrk1d)

        dummy = visc
        call IO_READ_FIELDS('flow.rand', IO_FLOW, imax, jmax, kmax, 3, 1, u, wrk3d)
        call IO_READ_FIELDS('flow.rand', IO_FLOW, imax, jmax, kmax, 3, 2, v, wrk3d)
        call IO_READ_FIELDS('flow.rand', IO_FLOW, imax, jmax, kmax, 3, 3, w, wrk3d)
        visc = dummy

        do j = 1, jmax   ! Remove mean
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, u)
            u(:, j, :) = u(:, j, :) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, v)
            v(:, j, :) = v(:, j, :) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, w)
            w(:, j, :) = w(:, j, :) - dummy
        end do

        ! ###################################################################
        select case (flag_u)
        case (2) ! Velocity given
            do j = 1, jmax
                u(:, j, :) = u(:, j, :)*wrk1d(j, 2)
                v(:, j, :) = v(:, j, :)*wrk1d(j, 1)
                w(:, j, :) = w(:, j, :)*wrk1d(j, 2)
            end do

        case (3)  ! Vorticity given, solve lap(u) = - rot(vort), vort = rot(u)
            call FI_CURL(imax, jmax, kmax, u, v, w, ax, ay, az, tmp4, wrk2d, wrk3d)
            do j = 1, jmax
                ax(:, j, :) = -ax(:, j, :)*wrk1d(j, 2)
                ay(:, j, :) = -ay(:, j, :)*wrk1d(j, 1)
                az(:, j, :) = -az(:, j, :)*wrk1d(j, 2)
            end do
            call FI_CURL(imax, jmax, kmax, ax, ay, az, u, v, w, tmp4, wrk2d, wrk3d)

            ! Solve lap(u) = - (rot(vort))_x
            if (flag_wall == 0) then; ibc = 3         ! FreeSlip
            else; ibc = 0
            end if  ! NoSlip
            if (g(1)%periodic .and. g(3)%periodic) then
                wrk2d(:, :, 1:2) = C_0_R                      ! bcs
                call OPR_POISSON_FXZ(.false., imax, jmax, kmax, g, ibc, &
                                     u, wrk3d, tmp4, tmp5, wrk2d(1, 1, 1), wrk2d(1, 1, 2), wrk1d, wrk1d(1, 5), wrk3d)
            else                                          ! General treatment, need global variable ipos,jpos,kpos,ci,cj,ck
#ifdef USE_CGLOC
                call CGPOISSON(1, imax, jmax, kmax, g(3)%size, u, ax, ay, az, ipos, jpos, kpos, ci, cj, ck, wrk2d)
#endif
            end if

            ! Solve lap(v) = - (rot(vort))_y
            ibc = 0                                       ! No penetration
            if (g(1)%periodic .and. g(3)%periodic) then
                wrk2d(:, :, 1:2) = C_0_R                      ! bcs
                call OPR_POISSON_FXZ(.false., imax, jmax, kmax, g, ibc, &
                                     v, wrk3d, tmp4, tmp5, wrk2d(1, 1, 1), wrk2d(1, 1, 2), wrk1d, wrk1d(1, 5), wrk3d)
            else                                          ! General treatment
#ifdef USE_CGLOC
                call CGPOISSON(1, imax, jmax, kmax, g(3)%size, v, ax, ay, az, ipos, jpos, kpos, ci, cj, ck, wrk2d)
#endif
            end if

            ! Solve lap(w) = - (rot(vort))_z
            if (g(3)%size > 1) then
                if (flag_wall == 0) then; ibc = 3         ! FreeSlip
                else; ibc = 0
                end if  ! NoSlip
                if (g(1)%periodic .and. g(3)%periodic) then
                    wrk2d(:, :, 1:2) = C_0_R                      ! bcs
                    call OPR_POISSON_FXZ(.false., imax, jmax, kmax, g, ibc, &
                                         w, wrk3d, tmp4, tmp5, wrk2d(1, 1, 1), wrk2d(1, 1, 2), wrk1d, wrk1d(1, 5), wrk3d)
                else                                          ! General treatment
#ifdef USE_CGLOC
                    call CGPOISSON(1, imax, jmax, kmax, g(3)%size, w, ax, ay, az, ipos, jpos, kpos, ci, cj, ck, wrk2d)
#endif
                end if
            end if

            ! ###################################################################
        case (4) ! Vector potential given
            do j = 1, jmax
                ax(:, j, :) = u(:, j, :)*wrk1d(j, 1) ! Horizontal components of vector potential give vertical velocity
                ay(:, j, :) = v(:, j, :)*wrk1d(j, 2)
                az(:, j, :) = w(:, j, :)*wrk1d(j, 1)
            end do

            bcs2 = 0
            if (flag_wall == 1 .or. flag_wall == 3) bcs2(1, 1) = 1 ! bcs at ymin = 1
            if (flag_wall == 2 .or. flag_wall == 3) bcs2(2, 1) = 1 ! bcs at ymax = 1
            ! Cannot use fi_curl. I need to impose BCs to zero to get zero velocity there
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs2, g(2), az, u, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), ay, tmp4, wrk3d, wrk2d, wrk3d)
            u = u - tmp4
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), ax, v, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), az, tmp4, wrk3d, wrk2d, wrk3d)
            v = v - tmp4
            if (g(3)%size > 1) then
                call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), ay, w, wrk3d, wrk2d, wrk3d)
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs2, g(2), ax, tmp4, wrk3d, wrk2d, wrk3d)
                w = w - tmp4
            end if

        end select

        ! ###################################################################
        ! Remove dilatation (vort was not really a vorticity field because it was not solenoidal)
        if (flag_dilatation == 1) then
            call FI_SOLENOIDAL(flag_wall, imax, jmax, kmax, u, v, w, ax, ay, az, tmp4, tmp5, wrk1d, wrk2d, wrk3d)
        end if

        if (g(3)%size == 1) w = C_0_R       ! Impose zero spanwise velocity in 2D case

        if (norm_ini_u >= C_0_R) call FLOW_NORMALIZE(u, v, w)

        return
    end subroutine VELOCITY_BROADBAND

    ! ###################################################################
    subroutine FLOW_NORMALIZE(u, v, w)
        real(cp), dimension(imax, jmax, kmax) :: u, v, w

        ! -------------------------------------------------------------------
        real(cp) AVG1V2D, dummy, amplify
        external AVG1V2D

        ! ###################################################################
        amplify = C_0_R                                      ! Maximum across the layer
        do j = 1, jmax
            dummy = AVG1V2D(imax, jmax, kmax, j, 2, u) + AVG1V2D(imax, jmax, kmax, j, 2, v) + AVG1V2D(imax, jmax, kmax, j, 2, w)
            amplify = MAX(dummy, amplify)
        end do
        amplify = C_05_R*amplify

        amplify = SQRT(norm_ini_u/amplify)           ! Scaling factor to normalize to maximum TKE

        u = u*amplify
        v = v*amplify
        w = w*amplify

        return
    end subroutine FLOW_NORMALIZE

end module FLOW_LOCAL

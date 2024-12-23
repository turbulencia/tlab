#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Interpolate a 3d field from an original grid into a new one
!# Routines for every direction call the same kernel routine INTERPOLATE_1D,
!# which calls in turn the routines from the library spline
!########################################################################
module OPR_INTERPOLATORS
    use FDM, only: grid_dt
    use TLab_Constants, only: efile, wp, wi
    use TLAB_VARS, only: isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLab_Constants, only: lfile
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLabMPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLabMPI_PROCS
#endif
    implicit none
    private

    public :: OPR_INTERPOLATE

#ifdef USE_MPI
    integer(wi) id
#endif

contains
!########################################################################
!########################################################################
    subroutine OPR_INTERPOLATE(nx, ny, nz, nx_dst, ny_dst, nz_dst, &
                               g, x_org, y_org, z_org, x_dst, y_dst, z_dst, u_org, u_dst, txc)

        integer(wi) nx, ny, nz, nx_dst, ny_dst, nz_dst
        type(grid_dt), intent(IN) :: g(3)
        real(wp), dimension(nx + 1), intent(INOUT) :: x_org
        real(wp), dimension(ny + 1), intent(INOUT) :: y_org
        real(wp), dimension(nz + 1), intent(INOUT) :: z_org
        real(wp), dimension(*), intent(IN) :: x_dst, y_dst, z_dst
        real(wp), dimension(nx*ny*nz), intent(IN) :: u_org
        real(wp), dimension(nx_dst*ny_dst*nz_dst), intent(OUT) :: u_dst
        real(wp), dimension(isize_txc_field, *), intent(INOUT) :: txc

! -------------------------------------------------------------------

#ifdef USE_MPI
        integer(wi) npage
#endif

! ###################################################################
! This should be OPR_INTERPOLATE_INITIALIZE
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLab_Write_ASCII(lfile, 'Initialize MPI type 1 for Ox interpolation.')
            id = TLabMPI_I_AUX1
            npage = nz*ny
            if (MOD(npage, ims_npro_i) /= 0) then ! add space for MPI transposition
                npage = npage/ims_npro_i
                npage = (npage + 1)*ims_npro_i
            end if
            call TLabMPI_TypeI_Create(ims_npro_i, nx, npage, 1, 1, 1, 1, &
                                      ims_size_i(id), ims_ds_i(:, id), ims_dr_i(:, id), ims_ts_i(id), ims_tr_i(id))

            call TLab_Write_ASCII(lfile, 'Initialize MPI type 2 for Ox interpolation.')
            id = TLabMPI_I_AUX2
            npage = nz*ny
            if (MOD(npage, ims_npro_i) /= 0) then ! add space for MPI transposition
                npage = npage/ims_npro_i
                npage = (npage + 1)*ims_npro_i
            end if
            call TLabMPI_TypeI_Create(ims_npro_i, nx_dst, npage, 1, 1, 1, 1, &
                                      ims_size_i(id), ims_ds_i(:, id), ims_dr_i(:, id), ims_ts_i(id), ims_tr_i(id))
        end if

        if (ims_npro_k > 1) then
            call TLab_Write_ASCII(lfile, 'Initialize MPI type 1 for Oz interpolation.')
            id = TLabMPI_K_AUX1
            npage = nx_dst*ny_dst
            call TLabMPI_TypeK_Create(ims_npro_k, nz, npage, 1, 1, 1, 1, &
                                      ims_size_k(id), ims_ds_k(:, id), ims_dr_k(:, id), ims_ts_k(id), ims_tr_k(id))

            call TLab_Write_ASCII(lfile, 'Initialize MPI type 2 for Oz interpolation.')
            id = TLabMPI_K_AUX2
            npage = nx_dst*ny_dst
            call TLabMPI_TypeK_Create(ims_npro_k, nz_dst, npage, 1, 1, 1, 1, &
                                      ims_size_k(id), ims_ds_k(:, id), ims_dr_k(:, id), ims_ts_k(id), ims_tr_k(id))

        end if
#endif

! #######################################################################
! Always interpolating along Ox
        if (g(1)%size > 1) then
            call OPR_INTERPOLATE_X(nx, ny, nz, nx_dst, g(1)%periodic, g(1)%scale, x_org, x_dst, &
                                   u_org, txc(1, 1), txc(1, 2), txc(1, 3))
        else
            txc(1:nx*ny*nz, 1) = u_org(1:nx*ny*nz)
        end if

        if (g(2)%size > 1) then
            call OPR_INTERPOLATE_Y(nx_dst, ny, nz, ny_dst, g(2)%periodic, g(2)%scale, y_org, y_dst, &
                                   txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
        else
            txc(1:nx_dst*ny*nz, 2) = txc(1:nx_dst*ny*nz, 1)
        end if

        if (g(3)%size > 1) then
            call OPR_INTERPOLATE_Z(nx_dst, ny_dst, nz, nz_dst, g(3)%periodic, g(3)%scale, z_org, z_dst, &
                                   txc(1, 2), u_dst, txc(1, 1), txc(1, 3))
        else
            u_dst(1:nx_dst*ny_dst*nz_dst) = txc(1:nx_dst*ny_dst*nz_dst, 2)
        end if

        return
    end subroutine OPR_INTERPOLATE

! #######################################################################
! Interpolation in X
! #######################################################################
    subroutine OPR_INTERPOLATE_X(nx, ny, nz, nx_dst, periodic, scalex, &
                                 x_org, x_dst, u_org, u_dst, u_tmp1, u_tmp2)

        logical periodic
        integer(wi) nx, ny, nz, nx_dst
        real(wp) scalex
        real(wp), dimension(*) :: x_org, x_dst
        real(wp), dimension(nx*ny*nz), target :: u_org
        real(wp), dimension(nx_dst*ny*nz), target :: u_dst
        real(wp), dimension(isize_txc_field), target :: u_tmp1, u_tmp2

        ! -----------------------------------------------------------------------
        integer(wi) nyz, nx_total, nx_total_dst

        real(wp), dimension(:), pointer :: p_a, p_b

        ! #######################################################################
        ! -------------------------------------------------------------------
        ! Transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            id = TLabMPI_I_AUX1
            u_tmp2(1:nx*ny*nz) = u_org(1:nx*ny*nz) ! Need additional space for transposition
            call TLabMPI_TRPF_I(u_tmp2, u_tmp1, id)

            p_a => u_tmp1
            p_b => u_tmp2

            nyz = ims_size_i(id)
            nx_total = nx*ims_npro_i
            nx_total_dst = nx_dst*ims_npro_i

        else
#endif
            p_a => u_org
            p_b => u_dst

            nyz = ny*nz
            nx_total = nx
            nx_total_dst = nx_dst

#ifdef USE_MPI
        end if
#endif

        ! -----------------------------------------------------------------------
        call INTERPOLATE_1D(periodic, nx_total, nyz, nx_total_dst, scalex, x_org, x_dst, p_a, p_b)

        ! -------------------------------------------------------------------
        ! Transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            id = TLabMPI_I_AUX2
            call TLabMPI_TRPB_I(u_tmp2, u_tmp1, id)
            u_dst(1:nx_dst*ny*nz) = u_tmp1(1:nx_dst*ny*nz)
        end if
#endif
        nullify (p_a, p_b)

        return
    end subroutine OPR_INTERPOLATE_X

    ! ###################################################################
    ! Interpolation in Oz direction
    ! ###################################################################
    subroutine OPR_INTERPOLATE_Z(nx, ny, nz, nz_dst, periodic, scalez, &
                                 z_org, z_dst, u_org, u_dst, u_tmp1, u_tmp2)

        logical periodic
        integer(wi) nx, ny, nz, nz_dst
        real(wp) scalez
        real(wp), dimension(*) :: z_org, z_dst
        real(wp), dimension(nx*ny*nz), target :: u_org
        real(wp), dimension(nx*ny*nz_dst), target :: u_dst
        real(wp), dimension(isize_txc_field), target :: u_tmp1, u_tmp2

        ! -----------------------------------------------------------------------
        integer(wi) nxy, nz_total, nz_total_dst

        real(wp), dimension(:), pointer :: p_a, p_b

        ! #######################################################################
        ! -------------------------------------------------------------------
        ! Transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            id = TLabMPI_K_AUX1
            call TLabMPI_TRPF_K(u_org, u_tmp2, id)

            p_a => u_tmp2
            p_b => u_tmp1

            nxy = ims_size_k(id)
            nz_total = nz*ims_npro_k
            nz_total_dst = nz_dst*ims_npro_k

        else
#endif
            p_a => u_org
            p_b => u_dst

            nxy = nx*ny
            nz_total = nz
            nz_total_dst = nz_dst
#ifdef USE_MPI
        end if
#endif

        ! -------------------------------------------------------------------
        ! Make z direction the first one
        ! -------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(p_a, nxy, nxy, nz_total, u_tmp1, nz_total)
#else
        call TLab_Transpose(p_a, nxy, nz_total, nxy, u_tmp1, nz_total)
#endif

        ! -----------------------------------------------------------------------
        call INTERPOLATE_1D(periodic, nz_total, nxy, nz_total_dst, scalez, z_org, z_dst, u_tmp1, u_tmp2)

        ! -------------------------------------------------------------------
        ! Put arrays back in the right order
        ! -------------------------------------------------------------------
#ifdef USE_ESSL
        call DGETMO(u_tmp2, nz_total_dst, nz_total_dst, nxy, p_b, nxy)
#else
        call TLab_Transpose(u_tmp2, nz_total_dst, nxy, nz_total_dst, p_b, nxy)
#endif

        ! -------------------------------------------------------------------
        ! Transposition
        ! -------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            id = TLabMPI_K_AUX2
            call TLabMPI_TRPB_K(u_tmp1, u_dst, id)
        end if
#endif
        nullify (p_a, p_b)

        return
    end subroutine OPR_INTERPOLATE_Z

    ! #######################################################################
    ! Interpolation in Y
    ! #######################################################################
    subroutine OPR_INTERPOLATE_Y(nx, ny, nz, ny_dst, periodic, scaley, &
                                 y_org, y_dst, u_org, u_dst, u_tmp1, u_tmp2)

        logical periodic
        integer(wi) nx, ny, nz, ny_dst
        real(wp) scaley
        real(wp) y_org(ny + 1)
        real(wp) y_dst(ny_dst)
        real(wp), dimension(nx, ny, nz) :: u_org, u_tmp1
        real(wp), dimension(nx, ny_dst, nz) :: u_dst, u_tmp2

        ! -----------------------------------------------------------------------
        integer(wi) ikmax, nyz

        ! #######################################################################
        ! -------------------------------------------------------------------
        ! Make y direction the first one (x direction the last one)
        ! -------------------------------------------------------------------
        nyz = ny*nz
#ifdef USE_ESSL
        call DGETMO(u_org, nx, nx, nyz, u_tmp1, nyz)
#else
        call TLab_Transpose(u_org, nx, nyz, nx, u_tmp1, nyz)
#endif

        ! -----------------------------------------------------------------------
        ikmax = nx*nz
        call INTERPOLATE_1D(periodic, ny, ikmax, ny_dst, scaley, y_org, y_dst, u_tmp1, u_tmp2)

        ! -------------------------------------------------------------------
        ! Put arrays back in the right order
        ! -------------------------------------------------------------------
        nyz = ny_dst*nz
#ifdef USE_ESSL
        call DGETMO(u_tmp2, nyz, nyz, nx, u_dst, nx)
#else
        call TLab_Transpose(u_tmp2, nyz, nx, nyz, u_dst, nx)
#endif

        return
    end subroutine OPR_INTERPOLATE_Y

    ! #######################################################################
    ! #######################################################################

    ! #######################################################################
    ! Interpolation in 1D
    ! #######################################################################
    subroutine INTERPOLATE_1D(periodic, imax, kmax, imax_dst, scalex, x_org, x_dst, u_org, u_dst)
        use TLab_Arrays, only: wrk1d
        logical periodic
        integer(wi) imax, kmax, imax_dst
        real(wp) scalex
        real(wp) x_org(imax + 1)
        real(wp) x_dst(imax_dst)
        real(wp) u_org(imax, *)
        real(wp) u_dst(imax_dst, *)

        real(wp) rdum
        integer(wi) k
        integer(wi), dimension(2) :: CSpline_BCType
        real(wp), dimension(2) :: CSpline_BCVal

        ! #######################################################################

        if (size(wrk1d) < 12*imax + 1) then
            call TLab_Write_ASCII(efile, 'INTERPOLATE_1D. Temporary Array not large enough')
            call TLab_Stop(DNS_ERROR_CURFIT)
        end if
        !------------------------------------------! the periodic case
        if (periodic) then                    !
            Cspline_BCType(1) = CS_BCS_PERIODIC  !
            Cspline_BCType(2) = CSpline_BCType(1)!
            CSpline_BCVal(:) = 0.0_wp             !
            x_org(imax + 1) = x_org(1) + scalex    ! extend x_org for periodic point
            do k = 1, kmax - 1                      !
                rdum = u_org(imax + 1, k)              ! use u_org(imax+1,k) to extend array and
                u_org(imax + 1, k) = u_org(1, k)      ! avoid the copy of the whole line
                call CUBIC_SPLINE(CSpline_BCType, CSpline_BCVal, &
                                  imax + 1, imax_dst, x_org, u_org(1, k), x_dst, u_dst(1, k), &
                                  wrk1d(imax + 2, 1))                 !
                u_org(imax + 1, k) = rdum            ! set u_org back to stored value rdum
            end do                                !
            wrk1d(1:imax, 1) = u_org(1:imax, kmax)     ! cannot avoid the copy for the last line
            wrk1d(imax + 1, 1) = u_org(1, kmax)          ! as u_org(imax+1,kmax) is out of bounds
            call CUBIC_SPLINE(CSpline_BCType, CSpline_BCVal, &
                              imax + 1, imax_dst, x_org, wrk1d, x_dst, u_dst(1, kmax), &
                              wrk1d(imax + 2, 1))
            !---------------------------------------! the aperiodic case
        else
            CSpline_BCType(1) = CS_BCS_NATURAL; CSpline_BCVal(1) = 0.0_wp
            CSpline_BCType(2) = CS_BCS_NATURAL; CSpline_BCVal(2) = 0.0_wp
            do k = 1, kmax
                call CUBIC_SPLINE(CSpline_BCType, CSpline_BCVal, &
                                  imax, imax_dst, x_org, u_org(1, k), x_dst, u_dst(1, k), &
                                  wrk1d)
            end do
        end if

        return
    end subroutine INTERPOLATE_1D

end module OPR_INTERPOLATORS

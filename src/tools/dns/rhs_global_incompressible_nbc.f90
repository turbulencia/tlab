#include "dns_error.h"
#include "dns_const.h"
#include "info_vars.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#ifdef USE_PSFFT
#include "nb3dfft_defines.inc"
#endif

subroutine RHS_GLOBAL_INCOMPRESSIBLE_NBC(u, v, w, s, &
                                         tmpu, tmpw, tmp11, tmp12, tmp21, tmp22, tmp31, tmp32, tmp41, tmp42, &
                                         bt1, bt2, bt3, bt4, &
                                         h1, h2, h3, hs)
    use, intrinsic :: iso_c_binding, only: c_int, c_loc, c_ptr, c_f_pointer

    use OMP_LIB, only: omp_get_thread_num

    use TLAB_CONSTANTS, only: lfile, wfile, efile, tfile
    !
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: inb_flow, inb_scal, inb_scal_array
    use TLAB_VARS, only: isize_field, isize_wrk1d, imax, jmax, kmax
    use TLAB_VARS, only: rbackground, ribackground
    !
    use THERMO_ANELASTIC
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use TIME, only: rkm_substep, rkm_endstep, dte
    use DNS_LOCAL, only: use_tower
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use FI_SOURCES
    use DNS_TOWER
    use PHASEAVG

#ifdef USE_PSFFT
    use DNS_LOCAL, only: nbcsetup
#endif

    use TLAB_MPI_VARS, only: ims_npro, ims_pro, ims_err, ims_size_i, ims_size_k

    use NB3DFFT, only: nb3dfft_nbc_prepare, nb3dfft_nbc_finish, nb3dfft_infoType
    use NB3DFFT, only: nb3dfft_nbc_schedl_start, nb3dfft_nbc_worker_start
    use NB3DFFT, only: nb3dfft_nbc_schedl_stop, nb3dfft_schedlType
    use NB3DFFT, only: nb3dfft_r2r_yxcomm, nb3dfft_r2r_yzcomm
    use NB3DFFT, only: nb3dfft_r2r_xycomm, nb3dfft_r2r_zycomm
    use NB3DFFT, only: nb3dfft_r2r_xunpack, nb3dfft_r2r_zunpack
    use NB3DFFT, only: nb3dfft_r2r_y1unpack, nb3dfft_r2r_y2unpack
    use NB3DFFT, only: nb3dfft_r2r_ready, mytype

    use MPI

    implicit none

    !
    ! PARAMETERS
    !
    integer(wi), parameter :: nmeasure = 3

    real(wp), dimension(isize_field), intent(IN) :: u, v, w
    real(wp), dimension(isize_field, inb_scal_array), intent(IN) :: s

    real(wp), dimension(isize_field), intent(INOUT) :: h1, h2, h3
    real(wp), dimension(isize_field, inb_scal), target, intent(OUT) :: hs
    real(wp), dimension(isize_field), intent(INOUT) :: tmpu, tmpw, tmp11, tmp12, tmp21, tmp22, tmp31, tmp32, tmp41, tmp42
    real(wp), dimension(isize_field) :: bt1, bt2, bt3, bt4
    real(wp), dimension(isize_wrk1d, *) :: wrk1d
    real(wp), dimension(*) :: wrk2d, wrk3d
    real(wp), dimension(:), pointer :: p_h
    !
    ! LOCAL VARIABLES
    !
    integer(wi) :: nxy_trans, nyz_trans, nxy, id, imeasure, ij, k, is, commID, iq
    integer(wi) :: finished, ip_b, ip_t, ibc, bcs(2, 2)
    real(wp) tdummy
    real(wp), dimension(:), pointer :: p_bcs
    !
    type(nb3dfft_infoType), dimension(24) :: info
    integer(KIND=4), dimension(2) :: cur_time
    real(KIND=mytype) :: t_comp, t_test, t_ser, t_init, t_wait, t_tmp
    real(KIND=8), dimension(6) :: t_snd, t_rcv
    type(nb3dfft_schedlType) :: nbcsetup_
    integer :: pkg_cnt
    real(wp) rtime, rtime_loc, t_run, ptime, ctime_loc

    target h1, h2, h3

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING SUBROUTINE, RHS_GLOBAL_INCOMPRESSIBLE_NBC')
#endif

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero
    BcsScalJmin%ref(:, :, :) = 0.0_wp
    BcsScalJmax%ref(:, :, :) = 0.0_wp

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        do is = 1, inb_scal
            if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) then
                p_bcs => hs(ip_b:, is); BcsScalJmin%ref(1:imax, k, is) = p_bcs(1:imax); end if
            if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) then
                p_bcs => hs(ip_t:, is); BcsScalJmax%ref(1:imax, k, is) = p_bcs(1:imax); end if
        end do
        ip_b = ip_b + nxy ! bottom BC address
        ip_t = ip_t + nxy ! top BC address
    end do

    nbcsetup_ = nbcsetup
    pkg_cnt = 24*ims_npro
    t_comp = 0
    t_test = 0
    t_ser = 0
    ptime = 0

    nxy = imax*jmax

    rtime = -MPI_WTime()
    t_run = rtime
    t_init = rtime

! #######################################################################
! Advection-diffusion terms
! #######################################################################
    call NB3DFFT_NBC_PREPARE(24, .false.)

!$omp parallel num_threads(2)  &
!$omp default(shared)
    if (omp_get_thread_num() == 0) then
        call NB3DFFT_NBC_SCHEDL_START(nbcsetup_)
    else if (omp_get_thread_num() == 1) then
        call NB3DFFT_NBC_WORKER_START()
        ; rtime_loc = MPI_WTime(); commID = 1; 
        info(FUYX)%timer_start = rtime_loc; info(FUYX)%id = commID; commID = commID + 1     ! STEP1
        info(BUXY)%timer_start = rtime_loc; info(BUXY)%id = commID; commID = commID + 1     ! STEP1
        info(FWYZ)%timer_start = rtime_loc; info(FWYZ)%id = commID; commID = commID + 1     ! STEP1
        info(BWZY)%timer_Start = rtime_loc; info(BWZY)%id = commID; commID = commID + 1     ! STEP1
        !
        if (inb_scal > 0) then
            info(F1YX)%timer_start = rtime_loc; info(F1YX)%id = commID; commID = commID + 1  ! STEP2
            info(B1XY)%timer_start = rtime_loc; info(B1XY)%id = commID; commID = commID + 1  ! STEP2
            info(F1YZ)%timer_start = rtime_loc; info(F1YZ)%id = commID; commID = commID + 1  ! STEP2
            info(B1ZY)%timer_start = rtime_loc; info(B1ZY)%id = commID; commID = commID + 1  ! STEP2
        end if

        info(FUYZ)%timer_start = rtime_loc; info(FUYZ)%id = commID; commID = commID + 1    ! STEP3
        info(BUZY)%timer_start = rtime_loc; info(BUZY)%id = commID; commID = commID + 1    ! STEP3
        info(FVYZ)%timer_start = rtime_loc; info(FVYZ)%id = commID; commID = commID + 1    ! STEP3
        info(BVZY)%timer_start = rtime_loc; info(BVZY)%id = commID; commID = commID + 1    ! STEP3

        if (inb_scal == 2) then
            info(F2YX)%timer_start = rtime_loc; info(F2YX)%id = commID; commID = commID + 1 ! STEP4
            info(B2XY)%timer_start = rtime_loc; info(B2XY)%id = commID; commID = commID + 1 ! STEP4
            info(F2YZ)%timer_start = rtime_loc; info(F2YZ)%id = commID; commID = commID + 1 ! STEP4
            info(B2ZY)%timer_start = rtime_loc; info(B2ZY)%id = commID; commID = commID + 1 ! STEP4
        end if

        info(FVYX)%timer_start = rtime_loc; info(FVYX)%id = commID; commID = commID + 1    ! STEP5
        info(BVXY)%timer_start = rtime_loc; info(BVXY)%id = commID; commID = commID + 1    ! STEP5
        info(FWYX)%timer_start = rtime_loc; info(FWYX)%id = commID; commID = commID + 1    ! STEP5
        info(BWXY)%timer_start = rtime_loc; info(BWXY)%id = commID; commID = commID + 1    ! STEP5

        info(FPYX)%timer_start = rtime_loc; info(FPYX)%id = commID; commID = commID + 1    ! STEP6
        info(BPXY)%timer_start = rtime_loc; info(BPXY)%id = commID; commID = commID + 1    ! STEP6
        info(FPYZ)%timer_start = rtime_loc; info(FPYZ)%id = commID; commID = commID + 1    ! STEP6
        info(BPZY)%timer_start = rtime_loc; info(BPZY)%id = commID; ! STEP6

        t_init = t_init + MPI_WTime()

        id = TLAB_MPI_I_PARTIAL; nyz_trans = ims_size_i(id)
        id = TLAB_MPI_K_PARTIAL; nxy_trans = ims_size_k(id)
        !
        ! kick off transpose U y->x and W y->z
        call NB3DFFT_R2R_YXCOMM(u, bt1, bt1, tmp11, info(FUYX), t_tmp); t_comp = t_comp + t_tmp
        call NB3DFFT_R2R_YZCOMM(w, tmpw, tmpw, bt2, info(FWYZ), t_tmp); t_comp = t_comp + t_tmp
        call NB3DFFT_R2R_YXCOMM(w, bt3, bt3, tmp31, info(FWYX), t_tmp); t_comp = t_comp + t_tmp
        !
        ! Vertical derivatives, and Vertical advection
        !
        t_tmp = -MPI_WTime()
        call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), v, v, tmp21, tmp22) ! store v transposed in tmp22
        h2 = h2 + tmp21
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), u, v, tmp21, tmpu, tmp22) ! using tmp22
        h1 = h1 + tmp21
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), w, v, tmp21, tmpu, tmp22) ! using tmp22
        h3 = h3 + tmp21
        t_ser = t_ser + (t_tmp + MPI_WTime())

        call NB3DFFT_R2R_YZCOMM(u, tmp41, tmp41, bt4, info(FUYZ), t_tmp); t_comp = t_comp + t_tmp

        t_tmp = -MPI_WTime()
        do is = 1, inb_scal
            call OPR_BURGERS_Y(OPR_B_U_IN, is, imax, jmax, kmax, bcs, g(2), s(1, is), v, tmp21, tmpu, tmp22) ! using tmp22
            hs(:, is) = hs(:, is) + tmp21
        end do
        t_ser = t_ser + (t_tmp + MPI_WTime())

        finished = 0
        do while (finished /= 2)
            if (NB3DFFT_R2R_READY(info(FUYX), t_tmp)) then
                t_test = t_test + t_tmp
                ! u du/dx + 1/Re d2u/dx2
                call NB3DFFT_R2R_XUNPACK(bt1, tmp11, info(FUYX), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                call DNS_TRANSPOSE(bt1, g(1)%size, nyz_trans, g(1)%size, tmpu, nyz_trans)
                call OPR_BURGERS_1D(0, nyz_trans, bcs, g(1), tmpu, tmpu, tmp11)
                call DNS_TRANSPOSE(tmp11, nyz_trans, g(1)%size, nyz_trans, bt1, g(1)%size)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_XYCOMM(bt1, bt1, tmp12, tmp11, info(BUXY), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if

            if (NB3DFFT_R2R_READY(info(FWYZ), t_tmp)) then
                t_test = t_test + t_tmp
                ! w dw/dz + 1/Re d2w/dz2
                call NB3DFFT_R2R_ZUNPACK(tmpw, bt2, info(FWYZ), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                call OPR_BURGERS_1D(0, nxy_trans, bcs, g(3), tmpw, tmpw, bt2)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_ZYCOMM(bt2, bt2, tmp22, tmp21, info(BWZY), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
        end do
        !
        do while (finished /= 4)
            if (NB3DFFT_R2R_READY(info(FWYX), t_tmp)) then
                t_test = t_test + t_tmp
                ! u dw/dx + 1/Re d2w/dx2
                call NB3DFFT_R2R_XUNPACK(bt3, tmp31, info(FWYX), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call DNS_TRANSPOSE(bt3, g(1)%size, nyz_trans, g(1)%size, tmp31, nyz_trans)
                call OPR_BURGERS_1D(0, nyz_trans, bcs, g(1), tmp31, tmpu, tmp32)
                call DNS_TRANSPOSE(tmp32, nyz_trans, g(1)%size, nyz_trans, bt3, g(1)%size)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_XYCOMM(bt3, bt3, tmp32, tmp31, info(BWXY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(FUYZ), t_tmp)) then
                t_test = t_test + t_tmp
                ! w du/dz + 1/Re d2u/dz2
                call NB3DFFT_R2R_ZUNPACK(tmp41, bt4, info(FUYZ), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call OPR_BURGERS_1D(0, nxy_trans, bcs, g(3), tmp41, tmpw, bt4)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_ZYCOMM(bt4, bt4, tmp42, tmp41, info(BUZY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
        end do
        !
        !
        do while (finished /= 8)
            if (NB3DFFT_R2R_READY(info(BUXY), t_tmp)) then
                t_test = t_test + t_tmp; 
                call NB3DFFT_R2R_Y1UNPACK(bt1, tmp11, info(BUXY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h1 = h1 + bt1
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_YXCOMM(v, bt1, bt1, tmp11, info(FVYX), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(BWZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt2, tmp21, info(BWZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h3 = h3 + bt2
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_YZCOMM(v, tmp21, tmp21, bt2, info(FVYZ), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(BWXY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y1UNPACK(bt3, tmp31, info(BWXY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h3 = h3 + bt3
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_YXCOMM(s(:, 1), bt3, bt3, tmp31, info(F1YX), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(BUZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt4, tmp41, info(BUZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h1 = h1 + bt4
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_YZCOMM(s(:, 1), tmp41, tmp41, bt4, info(F1YZ), t_tmp); t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
        end do
        ! u and w are finished and v and s1 are initiated
        !
        do while (finished /= 12)
            if (NB3DFFT_R2R_READY(info(FVYX), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_XUNPACK(bt1, tmp11, info(FVYX), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call DNS_TRANSPOSE(bt1, g(1)%size, nyz_trans, g(1)%size, tmp11, nyz_trans)
                call OPR_BURGERS_1D(0, nyz_trans, bcs, g(1), tmp11, tmpu, tmp12)
                call DNS_TRANSPOSE(tmp12, nyz_trans, g(1)%size, nyz_trans, bt1, g(1)%size)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_XYCOMM(bt1, bt1, tmp12, tmp11, info(BVXY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(FVYZ), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_ZUNPACK(tmp21, bt2, info(FVYZ), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call OPR_BURGERS_1D(0, nxy_trans, bcs, g(3), tmp21, tmpw, bt2)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_ZYCOMM(bt2, bt2, tmp22, tmp21, info(BVZY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(F1YX), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_XUNPACK(bt3, tmp31, info(F1YX), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call DNS_TRANSPOSE(bt3, g(1)%size, nyz_trans, g(1)%size, tmp31, nyz_trans)
                call OPR_BURGERS_1D(0, nyz_trans, bcs, g(1), tmp31, tmpu, tmp32)
                call DNS_TRANSPOSE(tmp32, nyz_trans, g(1)%size, nyz_trans, bt3, g(1)%size)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_XYCOMM(bt3, bt3, tmp32, tmp31, info(B1XY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(F1YZ), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_ZUNPACK(tmp41, bt4, info(F1YZ), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call OPR_BURGERS_1D(0, nxy_trans, bcs, g(3), tmp41, tmpw, bt4)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_ZYCOMM(bt4, bt4, tmp42, tmp41, info(B1ZY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
        end do
        !
        do while (finished /= 16)
            if (NB3DFFT_R2R_READY(info(BVXY), t_tmp)) then
                t_test = t_test + t_tmp; 
                call NB3DFFT_R2R_Y1UNPACK(bt1, tmp11, info(BVXY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h2 = h2 + bt1
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                if (inb_scal > 1) &
                    call NB3DFFT_R2R_YXCOMM(s(:, 2), bt1, bt1, tmp11, info(F2YX), t_tmp); 
                t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(BVZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt2, tmp21, info(BVZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                h2 = h2 + bt2
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                if (inb_scal > 1) &
                    call NB3DFFT_R2R_YZCOMM(s(:, 2), tmp21, tmp21, bt2, info(F2YZ), t_tmp); 
                t_comp = t_comp + t_tmp
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(B1XY), t_tmp)) then
                t_test = t_test + t_tmp; 
                call NB3DFFT_R2R_Y1UNPACK(bt3, tmp31, info(B1XY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                hs(:, 1) = hs(:, 1) + bt3
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(B1ZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt4, tmp41, info(B1ZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                hs(:, 1) = hs(:, 1) + bt4
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                finished = finished + 1
            end if
        end do
        ! u,v,w,s1 are finished and s2 initiated
        ! we can prepare the pressure solver, and update tendencies

        t_tmp = -MPI_WTime()
        !
        ! Source terms
        !
        call FI_SOURCES_FLOW(u, s, h1, tmp31)
        call FI_SOURCES_SCAL(s, hs, tmp31, tmp32, tmp42)
        !
        ! Impose buffer zone as relaxation terms
        !
        if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
            call BOUNDARY_BUFFER_RELAX_FLOW()
        end if
        !
        ! Calculate divergence for pressure solver
        !
        tdummy = 1.0_wp/dte
        tmp32 = h1 + u*tdummy
        tmp42 = h3 + w*tdummy
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp32)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp42)
        end if
        t_ser = t_ser + (t_tmp + MPI_WTime())

        call NB3DFFT_R2R_YXCOMM(tmp32, bt3, bt3, tmp31, info(FPYX), t_tmp); t_comp = t_comp + t_tmp
        call NB3DFFT_R2R_YZCOMM(tmp42, tmp41, tmp41, bt4, info(FPYZ), t_tmp); t_comp = t_comp + t_tmp

        ! Oy source term for pressure solver below
        !
        if (inb_scal < 2) &
            finished = finished + 2

        do while (finished /= 20)
            if (inb_scal > 1) then
                if (NB3DFFT_R2R_READY(info(F2YX), t_tmp)) then
                    t_test = t_test + t_tmp
                    call NB3DFFT_R2R_XUNPACK(bt1, tmp11, info(F2YX), t_tmp); t_comp = t_comp + t_tmp; 
                    !
                    t_tmp = -MPI_WTime()
                    call DNS_TRANSPOSE(bt1, g(1)%size, nyz_trans, g(1)%size, tmp11, nyz_trans)
                    call OPR_BURGERS_1D(0, nyz_trans, bcs, g(1), tmp11, tmpu, tmp12)
                    call DNS_TRANSPOSE(tmp12, nyz_trans, g(1)%size, nyz_trans, bt1, g(1)%size)
                    t_ser = t_ser + (t_tmp + MPI_WTime())
                    !
                    call NB3DFFT_R2R_XYCOMM(bt1, bt1, tmp12, tmp11, info(B2XY), t_tmp); t_comp = t_comp + t_tmp; 
                    finished = finished + 1
                end if
                if (NB3DFFT_R2R_READY(info(F2YZ), t_tmp)) then
                    t_test = t_test + t_tmp
                    call NB3DFFT_R2R_ZUNPACK(tmp21, bt2, info(F2YZ), t_tmp); t_comp = t_comp + t_tmp; 
                    !
                    t_tmp = -MPI_WTime()
                    call OPR_BURGERS_1D(0, nxy_trans, bcs, g(3), tmp21, tmpw, bt2)
                    t_ser = t_ser + (t_tmp + MPI_WTime())
                    !
                    call NB3DFFT_R2R_ZYCOMM(bt2, bt2, tmp22, tmp21, info(B2ZY), t_tmp); t_comp = t_comp + t_tmp; 
                    finished = finished + 1
                end if
            end if
            !
            if (NB3DFFT_R2R_READY(info(FPYX), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_XUNPACK(bt3, tmp31, info(FPYX), t_tmp); t_comp = t_comp + t_tmp; 
                !
                t_tmp = -MPI_WTime()
                call DNS_TRANSPOSE(bt3, g(1)%size, nyz_trans, g(1)%size, tmp31, nyz_trans)
                call OPR_PARTIAL1(nyz_trans, bcs, g(1), tmp31, tmp32, wrk2d)
                call DNS_TRANSPOSE(tmp32, nyz_trans, g(1)%size, nyz_trans, bt3, g(1)%size)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_XYCOMM(bt3, bt3, tmp32, tmp31, info(BPXY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(FPYZ), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_ZUNPACK(tmp41, bt4, info(FPYZ), t_tmp); t_comp = t_comp + t_tmp; 
                t_tmp = -MPI_WTime()
                call OPR_PARTIAL1(nxy_trans, bcs, g(3), tmp41, bt4, wrk2d)
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                call NB3DFFT_R2R_ZYCOMM(bt4, bt4, tmp42, tmp41, info(BPZY), t_tmp); t_comp = t_comp + t_tmp; 
                finished = finished + 1
            end if
        end do
        !
        do while (finished /= 22)
            if (NB3DFFT_R2R_READY(info(B2XY), t_tmp)) then
                t_test = t_test + t_tmp; 
                call NB3DFFT_R2R_Y1UNPACK(bt1, tmp11, info(B2XY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                hs(:, 2) = hs(:, 2) + bt1
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(B2ZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt2, tmp21, info(B2ZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                hs(:, 2) = hs(:, 2) + bt2
                t_ser = t_ser + t_tmp + MPI_WTime()
                !
                finished = finished + 1
            end if
        end do
        !
        t_tmp = -MPI_WTime()
        tdummy = 1.0_wp/dte
        tmp11 = h2 + v*tdummy
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp11)
        end if
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp11, tmp12)
        t_ser = t_ser + (t_tmp + MPI_WTime())
        !
        do while (finished /= 24)
            if (NB3DFFT_R2R_READY(info(BPXY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y1UNPACK(bt3, tmp31, info(BPXY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                tmp12 = tmp12 + bt3
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                finished = finished + 1
            end if
            if (NB3DFFT_R2R_READY(info(BPZY), t_tmp)) then
                t_test = t_test + t_tmp
                call NB3DFFT_R2R_Y2UNPACK(bt4, tmp41, info(BPZY), t_tmp); t_comp = t_comp + t_tmp
                !
                t_tmp = -MPI_WTime()
                tmp12 = tmp12 + bt4
                t_ser = t_ser + (t_tmp + MPI_WTime())
                !
                finished = finished + 1
            end if
        end do
        call NB3DFFT_NBC_SCHEDL_STOP(nbcsetup_)
    else
        PRTERR1("there MUST NOT be more threads than 2 on first level")
    end if
!$omp end parallel
    t_run = t_run + MPI_WTime()
    !               rhs-stuff    testing    packing     init
    t_wait = t_run - t_ser - t_test - t_comp - t_init
    call NB3DFFT_NBC_FINISH(nbcsetup_, t_run, pkg_cnt)
    ! CALL NB3DFFT_NBC_FINISH()
    nbcsetup = nbcsetup_

    t_snd = (/t_run, t_ser, t_comp, t_test, t_init, t_wait/)
    call MPI_Reduce(t_snd, t_rcv, 6, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)

#ifdef USE_PROFILE
    if (ims_pro == 0) then
        t_rcv = t_rcv/ims_npro
        t_run = t_rcv(1); t_ser = t_rcv(2); t_comp = t_rcv(3); 
        t_test = t_rcv(4); t_init = t_rcv(5); t_wait = t_rcv(6); 
        write (*, 903) ims_npro, t_run, &
            t_ser/t_run, t_comp/t_run, t_test/t_run, t_init/t_run, t_wait/t_run
903     format( &
            'NBC-time nproc', i6, &
            '(total [s], serial[1], comp[1], test[1], init[1], wait [1]);', &
            F10.3, 5(F6.3, ';'))
    end if
#endif

    ptime = -MPI_WTime()

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        p_bcs => h2(ip_b:); BcsFlowJmin%ref(1:imax, k, 2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        p_bcs => h2(ip_t:); BcsFlowJmax%ref(1:imax, k, 2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
    end do

! Adding density in BCs
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        BcsFlowJmin%ref(:, :, 2) = BcsFlowJmin%ref(:, :, 2)*rbackground(1)
        BcsFlowJmax%ref(:, :, 2) = BcsFlowJmax%ref(:, :, 2)*rbackground(g(2)%size)
    end if

! pressure in tmp12, Oy derivative in tmp11
    call OPR_POISSON_FXZ(imax, jmax, kmax, g, i3, tmp12, tmp41, tmp42, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp11)

    if (use_tower .and. rkm_substep == rkm_endstep) then
        call DNS_TOWER_ACCUMULATE(tmp12, i4, wrk1d)
    end if

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp12, tmp41)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp12, tmp42)

    if (mod(itime+1, phaseAvg%stride) == 0) then
        call SPACE_AVG(tmp12, avg_p, 1, wrk2d, (itime+1)/phaseAvg%stride, nitera_first, nitera_save/phaseAvg%stride, 4)
    end if

    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp41, h1)
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp11, h2)
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp42, h3)

    else
        h1 = h1 - tmp41
        h2 = h2 - tmp11
        h3 = h3 - tmp42

    end if

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
    BcsFlowJmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowJmax%ref = 0.0_wp

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (iq == 1) p_h => h1(1:)
        if (iq == 2) p_h => h2(1:)
        if (iq == 3) p_h => h3(1:)
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), p_h, &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp11)
        end if
    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hs(1, is), &
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), tmp11)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalJmax%type(is) /= DNS_SFC_STATIC) then
            call BOUNDARY_BCS_SURFACE_Y(is, bcs, s, hs, tmp11, tmp12)
        end if
    end do

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
    ip_b = 1
    do k = 1, kmax
        h1(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 1)
        h2(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 2)
        h3(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 3)
        do is = 1, inb_scal
            hs(ip_b:ip_b + imax - 1, is) = BcsScalJmin%ref(1:imax, k, is)
        end do
        ip_b = ip_b + nxy
    end do

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        h1(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 1)
        h2(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 2)
        h3(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 3)
        do is = 1, inb_scal
            hs(ip_t:ip_t + imax - 1, is) = BcsScalJmax%ref(1:imax, k, is)
        end do
        ip_t = ip_t + nxy
    end do

    ptime = ptime + MPI_WTime()
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC')
#endif
    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_NBC

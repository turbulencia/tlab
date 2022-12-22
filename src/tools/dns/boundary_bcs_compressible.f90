#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!  Nonperiodic characteristic BCs at xmin and xmax
module BOUNDARY_BCS_COMPRESSIBLE
    use TLAB_CONSTANTS, only: efile
    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS, only: imixture, gama0, THERMO_AI
    use BOUNDARY_INFLOW
    use BOUNDARY_BCS
    use OPR_PARTIAL
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif

    implicit none
    private

    public :: BOUNDARY_BCS_X, BOUNDARY_BCS_Y

contains
    subroutine BOUNDARY_BCS_X(iaux, M2_max, etime, rho, u, v, w, p, gama, z1, &
                              h0, h1, h2, h3, h4, zh1, txc, aux2d, wrk1d, wrk2d, wrk3d)

#include "integers.h"

        TINTEGER iaux

        TREAL M2_max, etime

        TREAL, dimension(imax, jmax, kmax) :: rho, u, v, w, p, gama, h0, h1, h2, h3, h4
        TREAL, dimension(imax, jmax, kmax, *) :: z1, zh1, txc
        TREAL, dimension(jmax, kmax, *) :: aux2d
        TREAL, dimension(*) :: wrk1d, wrk2d, wrk3d

        target aux2d

! -------------------------------------------------------------------
        TINTEGER j, k, is, nt, inb_scal_loc, isize, iflag_min, iflag_max, idir, ip0, bcs(2, 1)
        TREAL prefactor, pl_out_min, pl_out_max, pl_inf_min, pl_inf_max, pl_aux

        TREAL, dimension(:, :, :), pointer :: tmin, mmin, tmax, mmax, inf_rhs

! ###################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_X')
#endif

#define hr_loc(j,k)  aux2d(j,k,1)
#define hu_loc(j,k)  aux2d(j,k,2)
#define hv_loc(j,k)  aux2d(j,k,3)
#define hw_loc(j,k)  aux2d(j,k,4)
#define he_loc(j,k)  aux2d(j,k,5)
#define hz1_loc(j,k) aux2d(j,k,6)

#define r_loc(j,k)   aux2d(j,k,7)
#define u_loc(j,k)   aux2d(j,k,8)
#define v_loc(j,k)   aux2d(j,k,9)
#define w_loc(j,k)   aux2d(j,k,10)
#define p_loc(j,k)   aux2d(j,k,11)
#define g_loc(j,k)   aux2d(j,k,12)
#define z1_loc(j,k)  aux2d(j,k,13)

#define drdn_loc(j,k)  aux2d(j,k,14)
#define dudn_loc(j,k)  aux2d(j,k,15)
#define dvdn_loc(j,k)  aux2d(j,k,16)
#define dwdn_loc(j,k)  aux2d(j,k,17)
#define dpdn_loc(j,k)  aux2d(j,k,18)
#define dz1dn_loc(j,k) aux2d(j,k,19)

        bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

        ip0 = 19

        nt = jmax*kmax
        prefactor = (gama0 - C_1_R)*mach*mach

        if (iaux < nt*(19 + 5*(inb_flow + inb_scal_array))) then
            call TLAB_WRITE_ASCII(efile, 'BOUNDARY_BCS_X. Not enough space in txc.')
            call TLAB_STOP(DNS_ERROR_IBC)
        end if

! Define pointers
        inf_rhs => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        tmin => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        mmin => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        tmax => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        mmax => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)

! -------------------------------------------------------------------
! Type of characteristic BCs
! 1. only nonreflective
! 2. add fluctuation
! 3. add mean
! 4. add fluctuation+mean
!
! Relaxation towards a mean profile (Poinsot & Lele term)
! The local value of c is added later at the boundary
! Note that pl_??? has dimensions of 1/length
! -------------------------------------------------------------------
        idir = 1

        if (imode_sim == DNS_MODE_TEMPORAL) then ! not used
        else if (imode_sim == DNS_MODE_SPATIAL) then; iflag_min = -4; iflag_max = 3; end if

        pl_out_min = C_0_R ! default is only nonreflective
        if (BcsFlowImin%cout > 0) then
            pl_out_min = BcsFlowImin%cout*(C_1_R - M2_max)/g(1)%scale
        end if

        pl_inf_min = C_0_R ! jet inflow region (dimensions 1/time)
        if (BcsFlowImin%cinf > 0) then
            pl_inf_min = BcsFlowImin%cinf*qbg(1)%mean/qbg(1)%diam
        end if

        pl_aux = C_0_R     ! far from jet inflow region
        if (BcsFlowJmin%cinf > 0) then
            pl_aux = BcsFlowJmin%cinf/g(2)%scale
        end if

        pl_out_max = C_0_R ! default is only nonreflective
        if (BcsFlowImax%cout > 0) then
            pl_out_max = BcsFlowImax%cout*(C_1_R - M2_max)/g(1)%scale
        end if

        pl_inf_max = C_0_R
        if (BcsFlowImax%cinf > 0) then
            pl_inf_max = BcsFlowImax%cinf/g(1)%scale
        end if

! ###################################################################
! forcing terms in array inf_rhs
        if (inflow_mode /= 0) then
            isize = inb_flow + inb_scal_array
            inf_rhs(:, :, isize) = C_0_R

            if (inflow_mode == 1 .or. inflow_mode == 4) then
                call BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs, wrk2d, wrk3d)
            elseif (inflow_mode == 2 .or. inflow_mode == 3) then
                call BOUNDARY_INFLOW_BROADBAND(etime, inf_rhs, txc, wrk1d, wrk2d, wrk3d)
            end if
        end if

! ###################################################################
! Transverse terms
! ###################################################################
        call BOUNDARY_BCS_TRANSVERSE_X(u, v, w, p, rho, gama, z1, &
                                       tmin, mmin, tmax, mmax, txc(1, 1, 1, 1), txc(1, 1, 1, 2), txc(1, 1, 1, 3), wrk2d, wrk3d)

! ###################################################################
! Flow
! ###################################################################
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, txc(1, 1, 1, 2), wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), p, txc(1, 1, 1, 5), wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), rho, txc(1, 1, 1, 1), wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, txc(1, 1, 1, 3), wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, txc(1, 1, 1, 4), wrk3d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! Nonreflective BCs at xmin
! -------------------------------------------------------------------
        do k = 1, kmax
            do j = 1, jmax
                r_loc(j, k) = rho(1, j, k)
                u_loc(j, k) = u(1, j, k)
                v_loc(j, k) = v(1, j, k)
                w_loc(j, k) = w(1, j, k)
                p_loc(j, k) = p(1, j, k)
                g_loc(j, k) = gama(1, j, k)
                drdn_loc(j, k) = txc(1, j, k, 1)
                dudn_loc(j, k) = txc(1, j, k, 2)
                dvdn_loc(j, k) = txc(1, j, k, 3)
                dwdn_loc(j, k) = txc(1, j, k, 4)
                dpdn_loc(j, k) = txc(1, j, k, 5)
            end do
        end do
        if (imode_eqns == DNS_EQNS_TOTAL) then
            call BOUNDARY_BCS_FLOW_NR_2(i0, nt, pl_out_min, BcsFlowImin%ref(1, 1, 5), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(1), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        else if (imode_eqns == DNS_EQNS_INTERNAL) then
            call BOUNDARY_BCS_FLOW_NR_3(iflag_min, idir, nt, pl_aux, pl_inf_min, inf_rhs, BcsFlowImin%ref, &
                                        BcsFlowImin%ref(1, 1, inb_flow + 1), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(1), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! add transverse terms
            call BOUNDARY_BCS_FLOW_NR_4(iflag_min, idir, nt, BcsFlowImin%ctan, &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        tmin(:, :, 1), tmin(:, :, 2), tmin(:, :, 3), tmin(:, :, 4), tmin(:, :, 5), &
                                        mmin(:, :, 1), mmin(:, :, 5), &
                                        hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! edge corrections
            call BOUNDARY_BCS_FLOW_NR_EDGE(iflag_min, jmax, kmax, BcsFlowImin%ctan, &
                                           r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                           mmin(:, :, 1), mmin(:, :, 2), mmin(:, :, 3), mmin(:, :, 4), mmin(:, :, 5), &
                                           hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        end if
        do k = 1, kmax
            do j = 1, jmax
                h0(1, j, k) = h0(1, j, k) + hr_loc(j, k)
                h1(1, j, k) = h1(1, j, k) + hu_loc(j, k)
                h2(1, j, k) = h2(1, j, k) + hv_loc(j, k)
                h3(1, j, k) = h3(1, j, k) + hw_loc(j, k)
                h4(1, j, k) = h4(1, j, k) + he_loc(j, k)*prefactor
            end do
        end do
        if (imixture > 0) then
            do k = 1, kmax
                do j = 1, jmax
!           h4(1,j,k) = h4(1,j,k) + hr_loc(j,k)*THERMO_AI(6,1,NSP)
                    h4(1, j, k) = h4(1, j, k) + hr_loc(j, k)*THERMO_AI(6, 1, inb_scal + 1)
                end do
            end do
        end if

! -------------------------------------------------------------------
! Nonreflective BCs at xmax
! -------------------------------------------------------------------
        do k = 1, kmax
            do j = 1, jmax
                r_loc(j, k) = rho(imax, j, k)
                u_loc(j, k) = u(imax, j, k)
                v_loc(j, k) = v(imax, j, k)
                w_loc(j, k) = w(imax, j, k)
                p_loc(j, k) = p(imax, j, k)
                g_loc(j, k) = gama(imax, j, k)
                drdn_loc(j, k) = txc(imax, j, k, 1)
                dudn_loc(j, k) = txc(imax, j, k, 2)
                dvdn_loc(j, k) = txc(imax, j, k, 3)
                dwdn_loc(j, k) = txc(imax, j, k, 4)
                dpdn_loc(j, k) = txc(imax, j, k, 5)
            end do
        end do
        if (imode_eqns == DNS_EQNS_TOTAL) then
            call BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_out_max, BcsFlowImax%ref(1, 1, 5), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(1), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        else if (imode_eqns == DNS_EQNS_INTERNAL) then
            call BOUNDARY_BCS_FLOW_NR_3(iflag_max, idir, nt, pl_out_max, pl_inf_max, inf_rhs, BcsFlowImax%ref, &
                                        BcsFlowImax%ref(1, 1, inb_flow + 1), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(1), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! add transverse terms
            call BOUNDARY_BCS_FLOW_NR_4(iflag_max, idir, nt, BcsFlowImax%ctan, &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        tmax(:, :, 1), tmax(:, :, 2), tmax(:, :, 3), tmax(:, :, 4), tmax(:, :, 5), &
                                        mmax(:, :, 1), mmax(:, :, 5), &
                                        hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! edge corrections
            call BOUNDARY_BCS_FLOW_NR_EDGE(iflag_max, jmax, kmax, BcsFlowImax%ctan, &
                                           r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                           mmax(:, :, 1), mmax(:, :, 2), mmax(:, :, 3), mmax(:, :, 4), mmax(:, :, 5), &
                                           hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        end if
        do k = 1, kmax
            do j = 1, jmax
                h0(imax, j, k) = h0(imax, j, k) + hr_loc(j, k)
                h1(imax, j, k) = h1(imax, j, k) + hu_loc(j, k)
                h2(imax, j, k) = h2(imax, j, k) + hv_loc(j, k)
                h3(imax, j, k) = h3(imax, j, k) + hw_loc(j, k)
                h4(imax, j, k) = h4(imax, j, k) + he_loc(j, k)*prefactor
            end do
        end do
        if (imixture > 0) then
            do k = 1, kmax
                do j = 1, jmax
!           h4(imax,j,k) = h4(imax,j,k) + hr_loc(j,k)*THERMO_AI(6,1,NSP)
                    h4(imax, j, k) = h4(imax, j, k) + hr_loc(j, k)*THERMO_AI(6, 1, inb_scal + 1)
                end do
            end do
        end if

! ###################################################################
! Scalar
! ###################################################################
        if (icalc_scal == 1) then
            if (imixture == MIXT_TYPE_AIRWATER) then; inb_scal_loc = inb_scal + 1
            else; inb_scal_loc = inb_scal; end if
            do is = 1, inb_scal_loc
                call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), z1(1, 1, 1, is), txc(1, 1, 1, 3), wrk3d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! Nonreflective BCs at xmin
! -------------------------------------------------------------------
                do k = 1, kmax
                    do j = 1, jmax
                        r_loc(j, k) = rho(1, j, k)
                        u_loc(j, k) = u(1, j, k)
                        z1_loc(j, k) = z1(1, j, k, is)
                        p_loc(j, k) = p(1, j, k)
                        g_loc(j, k) = gama(1, j, k)
                        drdn_loc(j, k) = txc(1, j, k, 1)
                        dudn_loc(j, k) = txc(1, j, k, 2)
                        dz1dn_loc(j, k) = txc(1, j, k, 3)
                        dpdn_loc(j, k) = txc(1, j, k, 5)
                    end do
                end do
                call BOUNDARY_BCS_SCAL_NR_3(iflag_min, idir, nt, pl_aux, pl_inf_min, &
                            inf_rhs, inf_rhs(:, :, 5 + is), BcsFlowImin%ref, BcsScalImin%ref, BcsScalImin%ref(1, 1, inb_scal + 1), &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            drdn_loc(1, 1), dudn_loc(1, 1), dz1dn_loc(1, 1), dpdn_loc(1, 1), &
                                            buoyancy%vector(1), hz1_loc(1, 1))
! add transverse terms
                call BOUNDARY_BCS_SCAL_NR_4(iflag_min, nt, BcsScalImin%ctan, &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            tmin(:, :, 1), tmin(:, :, 2), tmin(:, :, 5), tmin(:, :, 5 + is), &
                                            hz1_loc(1, 1))
! edge corrections
                call BOUNDARY_BCS_SCAL_NR_EDGE(iflag_min, jmax, kmax, BcsScalImin%ctan, &
                                               r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                      mmin(:, :, 1), mmin(:, :, 2), mmin(:, :, 3), mmin(:, :, 5), mmin(:, :, 5 + is), hz1_loc(1, 1))
! special case affects only energy equation
                if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                else
                    do k = 1, kmax
                        do j = 1, jmax
                            zh1(1, j, k, is) = zh1(1, j, k, is) + hz1_loc(j, k)
                        end do
                    end do
                end if
                if (imixture > 0) then
! special case
                    if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                        do k = 1, kmax
                            do j = 1, jmax
                                h4(1, j, k) = h4(1, j, k) + hz1_loc(j, k)*(THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 1))
                            end do
                        end do
! general case
                    else
                        do k = 1, kmax
                            do j = 1, jmax
!                    h4(1,j,k) = h4(1,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                                h4(1, j, k) = h4(1, j, k) + hz1_loc(j, k)*(THERMO_AI(6, 1, is) - THERMO_AI(6, 1, inb_scal + 1))
                            end do
                        end do
                    end if
                end if

! -------------------------------------------------------------------
! Nonreflective BCs at xmax
! -------------------------------------------------------------------
                do k = 1, kmax
                    do j = 1, jmax
                        r_loc(j, k) = rho(imax, j, k)
                        u_loc(j, k) = u(imax, j, k)
                        z1_loc(j, k) = z1(imax, j, k, 1)
                        p_loc(j, k) = p(imax, j, k)
                        g_loc(j, k) = gama(imax, j, k)
                        drdn_loc(j, k) = txc(imax, j, k, 1)
                        dudn_loc(j, k) = txc(imax, j, k, 2)
                        dz1dn_loc(j, k) = txc(imax, j, k, 3)
                        dpdn_loc(j, k) = txc(imax, j, k, 5)
                    end do
                end do
                call BOUNDARY_BCS_SCAL_NR_3(iflag_max, idir, nt, pl_out_max, pl_inf_max, &
                            inf_rhs, inf_rhs(:, :, 5 + is), BcsFlowImax%ref, BcsScalImax%ref, BcsScalImax%ref(1, 1, inb_scal + 1), &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            drdn_loc(1, 1), dudn_loc(1, 1), dz1dn_loc(1, 1), dpdn_loc(1, 1), &
                                            buoyancy%vector(1), hz1_loc(1, 1))
! add transverse terms
                call BOUNDARY_BCS_SCAL_NR_4(iflag_max, nt, BcsScalImax%ctan, &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            tmax(:, :, 1), tmax(:, :, 2), tmax(:, :, 5), tmax(:, :, 5 + is), &
                                            hz1_loc(1, 1))
! edge corrections
                call BOUNDARY_BCS_SCAL_NR_EDGE(iflag_max, jmax, kmax, BcsScalImax%ctan, &
                                               r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                      mmax(:, :, 1), mmax(:, :, 2), mmax(:, :, 3), mmax(:, :, 5), mmax(:, :, 5 + is), hz1_loc(1, 1))
! special case affects only energy equation
                if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                else
                    do k = 1, kmax
                        do j = 1, jmax
                            zh1(imax, j, k, is) = zh1(imax, j, k, is) + hz1_loc(j, k)
                        end do
                    end do
                end if
                if (imixture > 0) then
! special case
                    if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                        do k = 1, kmax
                            do j = 1, jmax
                                h4(imax, j, k) = h4(imax, j, k) + hz1_loc(j, k)*(THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 1))
                            end do
                        end do
! general case
                    else
                        do k = 1, kmax
                            do j = 1, jmax
!                    h4(imax,j,k) = h4(imax,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                               h4(imax, j, k) = h4(imax, j, k) + hz1_loc(j, k)*(THERMO_AI(6, 1, is) - THERMO_AI(6, 1, inb_scal + 1))
                            end do
                        end do
                    end if
                end if

            end do

        end if

#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_X')
#endif

#undef hr_loc
#undef hu_loc
#undef hv_loc
#undef hw_loc
#undef he_loc
#undef hz1_loc

#undef r_loc
#undef u_loc
#undef v_loc
#undef w_loc
#undef p_loc
#undef g_loc
#undef z1_loc

#undef drdn_loc
#undef dudn_loc
#undef dvdn_loc
#undef dwdn_loc
#undef dpdn_loc
#undef dz1dn_loc

        return
    end subroutine BOUNDARY_BCS_X

!########################################################################
!#
!# Non-periodic characteristic BCs at ymin and ymax.
!# The flunctuating inflow forcing has not yet been implemented like
!# in BOUNDARY_BCS_X.
!#
!########################################################################
    subroutine BOUNDARY_BCS_Y(iaux, M2_max, rho, u, v, w, p, gama, z1, &
                              h0, h1, h2, h3, h4, zh1, tmp1, tmp2, tmp3, tmp4, tmp5, aux2d, wrk2d, wrk3d)

#include "integers.h"

        TINTEGER iaux
        TREAL M2_max

        TREAL, dimension(imax, jmax, kmax) :: rho, u, v, w, p, gama, h0, h1, h2, h3, h4
        TREAL, dimension(imax, jmax, kmax) :: tmp1, tmp2, tmp3, tmp4, tmp5
        TREAL, dimension(imax, jmax, kmax, *) :: z1, zh1
        TREAL, dimension(imax, kmax, *) :: aux2d
        TREAL, dimension(*) :: wrk2d, wrk3d

        target aux2d

! -------------------------------------------------------------------
        TINTEGER i, k, is, nt, inb_scal_loc, iflag_min, iflag_max, idir, ip0, bcs(2, 1)
        TINTEGER imin_loc, imax_loc
        TREAL prefactor, pl_out_min, pl_inf_min, pl_out_max, pl_inf_max

        TREAL, dimension(:, :, :), pointer :: tmin, lmin, tmax, lmax, inf_rhs

! ###################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_Y')
#endif

#define hr_loc(i,k)  aux2d(i,k,1)
#define hu_loc(i,k)  aux2d(i,k,2)
#define hv_loc(i,k)  aux2d(i,k,3)
#define hw_loc(i,k)  aux2d(i,k,4)
#define he_loc(i,k)  aux2d(i,k,5)
#define hz1_loc(i,k) aux2d(i,k,6)

#define r_loc(i,k)   aux2d(i,k,7)
#define u_loc(i,k)   aux2d(i,k,8)
#define v_loc(i,k)   aux2d(i,k,9)
#define w_loc(i,k)   aux2d(i,k,10)
#define p_loc(i,k)   aux2d(i,k,11)
#define g_loc(i,k)   aux2d(i,k,12)
#define z1_loc(i,k)  aux2d(i,k,13)

#define drdn_loc(i,k)  aux2d(i,k,14)
#define dudn_loc(i,k)  aux2d(i,k,15)
#define dvdn_loc(i,k)  aux2d(i,k,16)
#define dwdn_loc(i,k)  aux2d(i,k,17)
#define dpdn_loc(i,k)  aux2d(i,k,18)
#define dz1dn_loc(i,k) aux2d(i,k,19)

        bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

        ip0 = 19

        nt = imax*kmax
        prefactor = (gama0 - C_1_R)*mach*mach

        if (iaux < nt*(19 + 5*(inb_flow + inb_scal_array))) then
            call TLAB_WRITE_ASCII(efile, 'RHS_BCS_Y. Not enough space.')
            call TLAB_STOP(DNS_ERROR_JBC)
        end if

! Define pointers
        inf_rhs => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        tmin => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        lmin => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        tmax => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)
        ip0 = ip0 + inb_flow + inb_scal_array
        lmax => aux2d(:, :, ip0 + 1:ip0 + inb_flow + inb_scal_array)

! -------------------------------------------------------------------
! Type of characteristic BCs
! 1. only nonreflective
! 2. add fluctuation
! 3. add mean
! 4. add fluctuation+mean
!
! Relaxation towards a mean profile (Poinsot & Lele term)
! The local value of c is added later at the boundary
! Note that pl_??? has dimensions of 1/length
! -------------------------------------------------------------------
        idir = 2

        pl_out_min = C_0_R ! default is only nonreflective
        iflag_min = -1
        if (BcsFlowJmin%cout > 0) then
            pl_out_min = BcsFlowJmin%cout*(C_1_R - M2_max)/g(2)%scale
            iflag_min = -3
        end if

        pl_inf_min = C_0_R
        if (BcsFlowJmin%cinf > 0) then
            pl_inf_min = BcsFlowJmin%cinf/g(2)%scale
            iflag_min = -3
        end if

        pl_out_max = C_0_R ! default is only nonreflective
        iflag_max = -1
        if (BcsFlowJmax%cout > 0) then
            pl_out_max = BcsFlowJmax%cout*(C_1_R - M2_max)/g(2)%scale
            iflag_max = 3
        end if

        pl_inf_max = C_0_R
        if (BcsFlowJmax%cinf > 0) then
            pl_inf_max = BcsFlowJmax%cinf/g(2)%scale
            iflag_max = 3
        end if

        if (imode_sim == DNS_MODE_TEMPORAL) then; imin_loc = 1; imax_loc = imax
        else if (imode_sim == DNS_MODE_SPATIAL) then; imin_loc = 2; imax_loc = imax - 1; end if

! ###################################################################
! Transverse terms
! ###################################################################
        call BOUNDARY_BCS_TRANSVERSE_Y(u, v, w, p, rho, gama, z1, &
                                       tmin, lmin, tmax, lmax, tmp1, tmp2, tmp3, wrk2d, wrk3d)

! ###################################################################
! Flow
! ###################################################################
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), rho, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp3, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp4, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), p, tmp5, wrk3d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! BCs at ymin
! -------------------------------------------------------------------
        do k = 1, kmax; do i = 1, imax
                r_loc(i, k) = rho(i, 1, k)
                u_loc(i, k) = v(i, 1, k)
                v_loc(i, k) = u(i, 1, k)
                w_loc(i, k) = w(i, 1, k)
                p_loc(i, k) = p(i, 1, k)
                g_loc(i, k) = gama(i, 1, k)
                drdn_loc(i, k) = tmp1(i, 1, k)
                dudn_loc(i, k) = tmp3(i, 1, k)
                dvdn_loc(i, k) = tmp2(i, 1, k)
                dwdn_loc(i, k) = tmp4(i, 1, k)
                dpdn_loc(i, k) = tmp5(i, 1, k)
            end do; end do
        if (imode_eqns == DNS_EQNS_TOTAL) then
            call BOUNDARY_BCS_FLOW_NR_2(i0, nt, pl_out_min, BcsFlowJmin%ref(1, 1, 5), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(2), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        else if (imode_eqns == DNS_EQNS_INTERNAL) then
            call BOUNDARY_BCS_FLOW_NR_3(iflag_min, idir, nt, pl_out_min, pl_inf_min, inf_rhs, BcsFlowJmin%ref, &
                                        BcsFlowJmin%ref(1, 1, inb_flow + 1), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(2), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! add transverse terms
            call BOUNDARY_BCS_FLOW_NR_4(iflag_min, idir, nt, BcsFlowJmin%ctan, &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        tmin(:, :, 1), tmin(:, :, 3), tmin(:, :, 2), tmin(:, :, 4), tmin(:, :, 5), &
                                        lmin(:, :, 1), lmin(:, :, 5), &
                                        hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        end if
        do k = 1, kmax; do i = imin_loc, imax_loc
                h0(i, 1, k) = h0(i, 1, k) + hr_loc(i, k)
                h1(i, 1, k) = h1(i, 1, k) + hv_loc(i, k)
                h2(i, 1, k) = h2(i, 1, k) + hu_loc(i, k)
                h3(i, 1, k) = h3(i, 1, k) + hw_loc(i, k)
                h4(i, 1, k) = h4(i, 1, k) + he_loc(i, k)*prefactor
            end do; end do
        if (imixture > 0) then
            do k = 1, kmax; do i = imin_loc, imax_loc
!        h4(i,1,k) = h4(i,1,k) + hr_loc(i,k)*THERMO_AI(6,1,NSP)
                    h4(i, 1, k) = h4(i, 1, k) + hr_loc(i, k)*THERMO_AI(6, 1, inb_scal + 1)
                end do; end do
        end if

! -------------------------------------------------------------------
! BCs at ymax
! -------------------------------------------------------------------
        do k = 1, kmax; do i = 1, imax
                r_loc(i, k) = rho(i, jmax, k)
                u_loc(i, k) = v(i, jmax, k)
                v_loc(i, k) = u(i, jmax, k)
                w_loc(i, k) = w(i, jmax, k)
                p_loc(i, k) = p(i, jmax, k)
                g_loc(i, k) = gama(i, jmax, k)
                drdn_loc(i, k) = tmp1(i, jmax, k)
                dudn_loc(i, k) = tmp3(i, jmax, k)
                dvdn_loc(i, k) = tmp2(i, jmax, k)
                dwdn_loc(i, k) = tmp4(i, jmax, k)
                dpdn_loc(i, k) = tmp5(i, jmax, k)
            end do; end do
        if (imode_eqns == DNS_EQNS_TOTAL) then
            call BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_out_max, BcsFlowJmax%ref(1, 1, 5), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(2), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        else if (imode_eqns == DNS_EQNS_INTERNAL) then
            call BOUNDARY_BCS_FLOW_NR_3(iflag_max, idir, nt, pl_out_max, pl_inf_max, inf_rhs, BcsFlowJmax%ref, &
                                        BcsFlowJmax%ref(1, 1, inb_flow + 1), &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        drdn_loc(1, 1), dudn_loc(1, 1), dvdn_loc(1, 1), dwdn_loc(1, 1), dpdn_loc(1, 1), &
                                        buoyancy%vector(2), hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
! add transverse terms
            call BOUNDARY_BCS_FLOW_NR_4(iflag_max, idir, nt, BcsFlowJmax%ctan, &
                                        r_loc(1, 1), u_loc(1, 1), v_loc(1, 1), w_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                        tmax(:, :, 1), tmax(:, :, 3), tmax(:, :, 2), tmax(:, :, 4), tmax(:, :, 5), &
                                        lmax(:, :, 1), lmax(:, :, 5), &
                                        hr_loc(1, 1), hu_loc(1, 1), hv_loc(1, 1), hw_loc(1, 1), he_loc(1, 1))
        end if
        do k = 1, kmax; do i = imin_loc, imax_loc
                h0(i, jmax, k) = h0(i, jmax, k) + hr_loc(i, k)
                h1(i, jmax, k) = h1(i, jmax, k) + hv_loc(i, k)
                h2(i, jmax, k) = h2(i, jmax, k) + hu_loc(i, k)
                h3(i, jmax, k) = h3(i, jmax, k) + hw_loc(i, k)
                h4(i, jmax, k) = h4(i, jmax, k) + he_loc(i, k)*prefactor
            end do; end do
        if (imixture > 0) then
            do k = 1, kmax; do i = imin_loc, imax_loc
!        h4(i,jmax,k) = h4(i,jmax,k) + hr_loc(i,k)*THERMO_AI(6,1,NSP)
                    h4(i, jmax, k) = h4(i, jmax, k) + hr_loc(i, k)*THERMO_AI(6, 1, inb_scal + 1)
                end do; end do
        end if

! ###################################################################
! Scalar
! ###################################################################
        if (icalc_scal == 1) then
            if (imixture == MIXT_TYPE_AIRWATER) then; inb_scal_loc = inb_scal + 1
            else; inb_scal_loc = inb_scal; end if

            do is = 1, inb_scal_loc
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), z1(1, 1, 1, is), tmp2, wrk3d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! BCs at ymin
! -------------------------------------------------------------------
                do k = 1, kmax; do i = 1, imax
                        r_loc(i, k) = rho(i, 1, k)
                        u_loc(i, k) = v(i, 1, k)
                        z1_loc(i, k) = z1(i, 1, k, is)
                        p_loc(i, k) = p(i, 1, k)
                        g_loc(i, k) = gama(i, 1, k)
                        drdn_loc(i, k) = tmp1(i, 1, k)
                        dudn_loc(i, k) = tmp3(i, 1, k)
                        dz1dn_loc(i, k) = tmp2(i, 1, k)
                        dpdn_loc(i, k) = tmp5(i, 1, k)
                    end do; end do
                call BOUNDARY_BCS_SCAL_NR_3(iflag_min, idir, nt, pl_out_min, pl_inf_min, &
                            inf_rhs, inf_rhs(:, :, 5 + is), BcsFlowJmin%ref, BcsScalJmin%ref, BcsScalJmin%ref(1, 1, inb_scal + 1), &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            drdn_loc(1, 1), dudn_loc(1, 1), dz1dn_loc(1, 1), dpdn_loc(1, 1), &
                                            buoyancy%vector(2), hz1_loc(1, 1))
! add transverse terms
                call BOUNDARY_BCS_SCAL_NR_4(iflag_min, nt, BcsScalJmin%ctan, &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            tmin(:, :, 1), tmin(:, :, 3), tmin(:, :, 5), tmin(:, :, 5 + is), &
                                            hz1_loc(1, 1))
                if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                else
                    do k = 1, kmax; do i = imin_loc, imax_loc
                            zh1(i, 1, k, is) = zh1(i, 1, k, is) + hz1_loc(i, k)
                        end do; end do
                end if
                if (imixture > 0) then
! special case
                    if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                        do k = 1, kmax; do i = imin_loc, imax_loc
                                h4(i, 1, k) = h4(i, 1, k) + hz1_loc(i, k)*(THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 1))
                            end do; end do
! general case
                    else
                        do k = 1, kmax; do i = imin_loc, imax_loc
!                 h4(i,1,k) = h4(i,1,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                                h4(i, 1, k) = h4(i, 1, k) + hz1_loc(i, k)*(THERMO_AI(6, 1, is) - THERMO_AI(6, 1, inb_scal + 1))
                            end do; end do
                    end if
                end if

! -------------------------------------------------------------------
! BCs at ymax
! -------------------------------------------------------------------
                do k = 1, kmax; do i = 1, imax
                        r_loc(i, k) = rho(i, jmax, k)
                        u_loc(i, k) = v(i, jmax, k)
                        z1_loc(i, k) = z1(i, jmax, k, is)
                        p_loc(i, k) = p(i, jmax, k)
                        g_loc(i, k) = gama(i, jmax, k)
                        drdn_loc(i, k) = tmp1(i, jmax, k)
                        dudn_loc(i, k) = tmp3(i, jmax, k)
                        dz1dn_loc(i, k) = tmp2(i, jmax, k)
                        dpdn_loc(i, k) = tmp5(i, jmax, k)
                    end do; end do
                call BOUNDARY_BCS_SCAL_NR_3(iflag_max, idir, nt, pl_out_max, pl_inf_max, &
                            inf_rhs, inf_rhs(:, :, 5 + is), BcsFlowJmax%ref, BcsScalJmax%ref, BcsScalJmax%ref(1, 1, inb_scal + 1), &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            drdn_loc(1, 1), dudn_loc(1, 1), dz1dn_loc(1, 1), dpdn_loc(1, 1), &
                                            buoyancy%vector(2), hz1_loc(1, 1))
! add transverse terms
                call BOUNDARY_BCS_SCAL_NR_4(iflag_max, nt, BcsScalJmax%ctan, &
                                            r_loc(1, 1), u_loc(1, 1), z1_loc(1, 1), p_loc(1, 1), g_loc(1, 1), &
                                            tmax(:, :, 1), tmax(:, :, 3), tmax(:, :, 5), tmax(:, :, 5 + is), &
                                            hz1_loc(1, 1))
! special case affects only energy equation
                if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                else
                    do k = 1, kmax; do i = imin_loc, imax_loc
                            zh1(i, jmax, k, is) = zh1(i, jmax, k, is) + hz1_loc(i, k)
                        end do; end do
                end if
                if (imixture > 0) then
! special case
                    if (imixture == MIXT_TYPE_AIRWATER .and. is == 2) then
                        do k = 1, kmax; do i = imin_loc, imax_loc
                                h4(i, jmax, k) = h4(i, jmax, k) + hz1_loc(i, k)*(THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 1))
                            end do; end do
! general case
                    else
                        do k = 1, kmax; do i = imin_loc, imax_loc
!                 h4(i,jmax,k) = h4(i,jmax,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                               h4(i, jmax, k) = h4(i, jmax, k) + hz1_loc(i, k)*(THERMO_AI(6, 1, is) - THERMO_AI(6, 1, inb_scal + 1))
                            end do; end do
                    end if
                end if

            end do

        end if

#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_Y')
#endif

#undef hr_loc
#undef hu_loc
#undef hv_loc
#undef hw_loc
#undef he_loc
#undef hz1_loc

#undef r_loc
#undef u_loc
#undef v_loc
#undef w_loc
#undef p_loc
#undef g_loc
#undef z1_loc

#undef drdn_loc
#undef dudn_loc
#undef dvdn_loc
#undef dwdn_loc
#undef dpdn_loc
#undef dz1dn_loc

        return
    end subroutine BOUNDARY_BCS_Y

!########################################################################
!#
!# Implementation of the nonreflective boundary conditions for the
!# transport equations of density, momentum and total energy.
!#
!# Poinsot&Lele correction term is included, with pl_const.
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag   In   Flag selecting BC at xmin (0) or at xmax (1)
!#
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!# gn        In   Constant body force normal to the boundary
!#
!########################################################################
    subroutine BOUNDARY_BCS_FLOW_NR_2 &
        (iflag, nt, pl_const, pl_pref, &
         r, un, v1, v2, p, gama, drdn, dundn, dv1dn, dv2dn, dpdn, gn, &
         hr, hun, hv1, hv2, he)

        TINTEGER iflag, nt
        TREAL pl_const, pl_pref(*)
        TREAL r(*), un(*), v1(*), v2(*), p(*), gama(*)
        TREAL drdn(*), dundn(*), dv1dn(*), dv2dn(*), dpdn(*), gn
        TREAL hr(*), hun(*), hv1(*), hv2(*), he(*)

! -------------------------------------------------------------------
        TINTEGER i
        TREAL c, Mn, M2, dummy

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag == 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c
                M2 = C_05_R*(un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i))/c/c

                if (un(i) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) > C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R - Mn)/c*dpdn(i) &
                                        - r(i)*gn/c)

                        hr(i) = un(i)*drdn(i) + dummy
                        hun(i) = un(i)*un(i)*drdn(i) + dummy*c*(1 + Mn) + Mn*dpdn(i)
                        hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
                        hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
                        he(i) = un(i)*M2*c*c*drdn(i) + r(i)*un(i)*(v1(i)*dv1dn(i) + v2(i)*dv2dn(i)) &
                                + dummy*c*c*(C_1_R/(gama(i) - C_1_R) + M2 + Mn) &
                                + un(i)*(C_1_R/(gama(i) - C_1_R) + Mn)*dpdn(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R + Mn)/c*dpdn(i) &
                                        - r(i)*gn/c - pl_const*(p(i) - pl_pref(i))/c)

                        hr(i) = dummy
                        hun(i) = dummy*c*(C_1_R + Mn)
                        hv1(i) = dummy*v1(i)
                        hv2(i) = dummy*v2(i)
                        he(i) = dummy*c*c*(C_1_R/(gama(i) - C_1_R) + M2 + Mn)

                    end if

                end if ! subsonic branch

            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag == 1) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c
                M2 = C_05_R*(un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i))/c/c

                if (un(i) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) < C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R + Mn)/c*dpdn(i) &
                                        + r(i)*gn/c)

                        hr(i) = un(i)*drdn(i) + dummy
                        hun(i) = un(i)*un(i)*drdn(i) - (1 - Mn)*c*dummy - Mn*dpdn(i)
                        hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
                        hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
                        he(i) = un(i)*M2*c*c*drdn(i) + r(i)*un(i)*(v1(i)*dv1dn(i) + v2(i)*dv2dn(i)) &
                                + dummy*c*c*(C_1_R/(gama(i) - C_1_R) + M2 - Mn) &
                                + un(i)*(C_1_R/(gama(i) - C_1_R) - Mn)*dpdn(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R - Mn)/c*dpdn(i) &
                                        + r(i)*gn/c - pl_const*(p(i) - pl_pref(i))/c)

                        hr(i) = dummy
                        hun(i) = -dummy*c*(C_1_R - Mn)
                        hv1(i) = dummy*v1(i)
                        hv2(i) = dummy*v2(i)
                        he(i) = dummy*c*c*(C_1_R/(gama(i) - C_1_R) + M2 - Mn)

                    end if

                end if ! subsonic branch

            end do

        end if

        return
    end subroutine BOUNDARY_BCS_FLOW_NR_2

    !########################################################################
!#
!# Implementation of the nonreflective boundary conditions for the
!# transport equations of density, momentum and internal energy.
!#
!# Poinsot&Lele correction term is included, with pl_out.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR to use in the formulations using the
!# internal energy instead of the total energy.
!#
!# The case of internal energy is like that of p, devided by (\gamma-1),
!# and with the formation part coming from the scalar equations.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR_1 adding forcing terms for inflow
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!#                1. nonreflective
!#                2. fluctuation
!#                3. mean
!#                4. fluctuation+mean
!# idir      In   Flag to predefined models for relaxation terms
!#                1. OX direction
!#                2. OY direction
!#                3. all terms
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!# gn        In   Constant body force normal to the boundary
!#
!########################################################################
    subroutine BOUNDARY_BCS_FLOW_NR_3(iflag, idir, nt, pl_out, pl_inf, inf_rhs, bf, bf_shape, &
                                      r, un, v1, v2, p, gama, drdn, dundn, dv1dn, dv2dn, dpdn, gn, hr, hun, hv1, hv2, he)

        TINTEGER iflag, idir, nt
        TREAL pl_out, pl_inf
        TREAL inf_rhs(nt, *), bf(nt, *), bf_shape(*)
        TREAL r(*), un(*), v1(*), v2(*), p(*), gama(*)
        TREAL drdn(*), dundn(*), dv1dn(*), dv2dn(*), dpdn(*), gn
        TREAL hr(*), hun(*), hv1(*), hv2(*), he(*)

! -------------------------------------------------------------------
        TINTEGER i
        TREAL c, Mn, M2, dummy
        TREAL F1, F2, F3, F4, F5

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c
                M2 = C_05_R*(un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i))/c/c

                if (un(i) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) > C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R - Mn)/c*dpdn(i) - r(i)*gn/c)

                        hr(i) = un(i)*drdn(i) + dummy
                        hun(i) = un(i)*un(i)*drdn(i) + dummy*c*(1 + Mn) + Mn*dpdn(i)
                        hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
                        hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
                        he(i) = (un(i)*dpdn(i) + dummy*c*c)/(gama(i) - C_1_R)

! add forcing terms
                        F2 = C_0_R
                        F3 = C_0_R
                        F4 = C_0_R
                        F5 = C_0_R
                        if (abs(iflag) == 2 .or. abs(iflag) == 4) then ! fluctuation
                            F2 = F2 + (inf_rhs(i, 1) - inf_rhs(i, 5)/c/c)*bf_shape(i)
                            F3 = F3 + (inf_rhs(i, 3))*bf_shape(i)
                            F4 = F4 + (inf_rhs(i, 4))*bf_shape(i)
                            F5 = F5 + (inf_rhs(i, 5) + inf_rhs(i, 2)*r(i)*c)*bf_shape(i)
                        end if
                        if (abs(iflag) == 3 .or. abs(iflag) == 4) then ! mean
                            if (idir == 1) then      ! OX direction; possible v1 forcing w/ F3
!                    F2 = F2 - pl_inf*( (r(i)-p(i)/c/c) - (bf(i,1)-bf(i,5)/c/c) )
                                F2 = F2 - pl_inf*(r(i) - bf(i, 1))*bf_shape(i) &
                                     - pl_out*c*(r(i) - bf(i, 1))*(C_1_R - bf_shape(i))
                                F3 = F3 - pl_inf*(v1(i) - bf(i, 3))*bf_shape(i)
                                F4 = F4 - pl_inf*(v2(i) - bf(i, 4))*bf_shape(i) &
                                     - pl_out*c*(v2(i) - bf(i, 4))*(C_1_R - bf_shape(i))
                                F5 = F5 - pl_inf*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*bf_shape(i) &
                                     - pl_out*c*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*(C_1_R - bf_shape(i))
                            else if (idir == 2) then ! OY direction; no un forcing w/ F5
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                F3 = F3 - pl_inf*c*(v1(i) - bf(i, 3))
                                F4 = F4 - pl_inf*c*(v2(i) - bf(i, 4))
                                F5 = F5 - pl_inf*c*(p(i) - bf(i, 5))
                            end if
                        end if

                        dummy = F2 + C_05_R*F5/c/c

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) + un(i)*F2 + C_05_R*(Mn + C_1_R)*F5/c
                        hv1(i) = hv1(i) + r(i)*F3 + v1(i)*dummy
                        hv2(i) = hv2(i) + r(i)*F4 + v2(i)*dummy
                        he(i) = he(i) + C_05_R*F5/(gama(i) - C_1_R)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        F5 = C_0_R
! Poinsot & Lele term. Multiplication by c done in the definition of dummy below
                        if (idir == 1) then      ! OX treatment at xmin; complete forcing
                            F5 = -pl_inf/c*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*bf_shape(i) &
                                 - pl_out*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*(C_1_R - bf_shape(i))
                        else if (idir == 2) then ! OY treatment; no un forcing w/ F5
                            F5 = -pl_out*(p(i) - bf(i, 5))
                        end if

                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R + Mn)/c*dpdn(i) - r(i)*gn/c + F5/c)

                        hr(i) = dummy
                        hun(i) = dummy*c*(C_1_R + Mn)
                        hv1(i) = dummy*v1(i)
                        hv2(i) = dummy*v2(i)
                        he(i) = dummy*c*c/(gama(i) - C_1_R)

                    end if

                end if ! subsonic branch

            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c
                M2 = C_05_R*(un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i))/c/c

                if (un(i) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) < C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R + Mn)/c*dpdn(i) + r(i)*gn/c)

                        hr(i) = un(i)*drdn(i) + dummy
                        hun(i) = un(i)*un(i)*drdn(i) - (1 - Mn)*c*dummy - Mn*dpdn(i)
                        hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
                        hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
                        he(i) = (un(i)*dpdn(i) + dummy*c*c)/(gama(i) - C_1_R)

! add forcing terms: only mean
                        F1 = C_0_R
                        F2 = C_0_R
                        F3 = C_0_R
                        F4 = C_0_R
                        if (abs(iflag) == 3 .or. abs(iflag) == 4) then ! mean
                            if (idir == 1) then      ! OX direction
                                F1 = F1 - pl_inf*c*((p(i) - r(i)*c*un(i)) - (bf(i, 5) - r(i)*c*bf(i, 2)))
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                F3 = F3 - pl_inf*c*(v1(i) - bf(i, 3))
                                F4 = F4 - pl_inf*c*(v2(i) - bf(i, 4))
                            else if (idir == 2) then ! OY direction; no un forcing w/ F1
                                F1 = F1 - pl_inf*c*(p(i) - bf(i, 5))
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                F3 = F3 - pl_inf*c*(v1(i) - bf(i, 3))
                                F4 = F4 - pl_inf*c*(v2(i) - bf(i, 4))
                            end if
                        end if

                        dummy = F2 + C_05_R*F1/c/c

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) + un(i)*F2 + C_05_R*(Mn - C_1_R)*F1/c
                        hv1(i) = hv1(i) + r(i)*F3 + v1(i)*dummy
                        hv2(i) = hv2(i) + r(i)*F4 + v2(i)*dummy
                        he(i) = he(i) + C_05_R*F1/(gama(i) - C_1_R)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        F1 = C_0_R
! Poinsot & Lele term. Multiplication by c done in the definition of dummy below
                        if (idir == 1) then      ! OX treatment; no un forcing w/ F1
                            F1 = -pl_out*(p(i) - bf(i, 5))
                        else if (idir == 2) then ! OY treatment; no un forcing w/ F1
                            F1 = -pl_out*(p(i) - bf(i, 5))
                        end if
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R - Mn)/c*dpdn(i) + r(i)*gn/c + F1/c)

                        hr(i) = dummy
                        hun(i) = -dummy*c*(C_1_R - Mn)
                        hv1(i) = dummy*v1(i)
                        hv2(i) = dummy*v2(i)
                        he(i) = dummy*c*c/(gama(i) - C_1_R)

                    end if

                end if ! subsonic branch

            end do

        end if

        return
    end subroutine BOUNDARY_BCS_FLOW_NR_3

!########################################################################
!#
!# Implementation of the nonreflective boundary conditions for the
!# transport equations of species mass fractions.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR_1 adding forcing terms for inflow
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag   In   Flag selecting BC at xmin (<0) or at xmax (>0)
!#              1. nonreflective
!#              2. fluctuation
!#              3. mean
!#              4. fluctuation+mean
!# idir    In   Flag to predefined models, like OX or OY directions
!# un      In   Velocity normal to the boundary
!# gn      In   Constant body force normal to the boundary
!#
!########################################################################
    subroutine BOUNDARY_BCS_SCAL_NR_3(iflag, idir, nt, pl_out, pl_inf, inf_rhs, inf_rhs_z, bf, bf_z, &
                                      bf_shape, r, un, z1, p, gama, drdn, dundn, dz1dn, dpdn, gn, hz1)

        TINTEGER iflag, idir, nt
        TREAL pl_out, pl_inf
        TREAL inf_rhs(nt, *), inf_rhs_z(nt), bf(nt, *), bf_z(nt), bf_shape(*)
        TREAL r(*), un(*), z1(*), p(*), gama(*)
        TREAL drdn(*), dundn(*), dz1dn(*), dpdn(*), gn
        TREAL hz1(*)

! -------------------------------------------------------------------
        TINTEGER i
        TREAL c, Mn, dummy
        TREAL F1, F2, F5, FZ

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) > C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R - Mn)/c*dpdn(i) - r(i)*gn/c)

                        hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)

! add forcing terms
                        F2 = C_0_R
                        F5 = C_0_R
                        FZ = C_0_R
                        if (abs(iflag) == 2 .or. abs(iflag) == 4) then ! fluctuation
                            F2 = F2 + (inf_rhs(i, 1) - inf_rhs(i, 5)/c/c)*bf_shape(i)
                            F5 = F5 + (inf_rhs(i, 5) + inf_rhs(i, 2)*r(i)*c)*bf_shape(i)
                            FZ = FZ + (inf_rhs_z(i))*bf_shape(i)
                        end if
                        if (abs(iflag) == 3 .or. abs(iflag) == 4) then ! mean
                            if (idir == 1) then      ! OX direction
!                    F2 = F2 - pl_inf*( (r(i)-p(i)/c/c) - (bf(i,1)-bf(i,5)/c/c) )
                                F2 = F2 - pl_inf*(r(i) - bf(i, 1))*bf_shape(i) &
                                     - pl_out*c*(r(i) - bf(i, 1))*(C_1_R - bf_shape(i))
                                F5 = F5 - pl_inf*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*bf_shape(i) &
                                     - pl_out*c*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*(C_1_R - bf_shape(i))
                                FZ = FZ - pl_inf*(z1(i) - bf_z(i))*bf_shape(i) &
                                     - pl_out*c*(z1(i) - bf_z(i))*(C_1_R - bf_shape(i))
                            else if (idir == 2) then ! OY direction; no un forcing w/ F5
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                F5 = F5 - pl_inf*c*(p(i) - bf(i, 5))
                                FZ = FZ - pl_inf*c*(z1(i) - bf_z(i))
                            end if
                        end if

                        dummy = F2 + C_05_R*F5/c/c

                        hz1(i) = hz1(i) + r(i)*FZ + z1(i)*dummy

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        F5 = C_0_R
! Poinsot & Lele term. Multiplication by c is done in the definition of dummy below
                        if (idir == 1) then      ! OX treatment at xmin; complete forcing
                            F5 = -pl_inf/c*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*bf_shape(i) &
                                 - pl_out*(p(i) + r(i)*c*un(i) - (bf(i, 5) + r(i)*c*bf(i, 2)))*(C_1_R - bf_shape(i))
                        else if (idir == 2) then ! OY treatment; no un forcing w/ F5
                            F5 = -pl_out*(p(i) - bf(i, 5))
                        end if
                        dummy = C_05_R*(r(i)*(C_1_R + Mn)*dundn(i) + (C_1_R + Mn)/c*dpdn(i) - r(i)*gn/c + F5/c)

                        hz1(i) = dummy*z1(i)

                    end if

                end if ! subsonic branch

            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) < C_0_R) then
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R + Mn)/c*dpdn(i) + r(i)*gn/c)

                        hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)

! add forcing terms: only mean
                        F1 = C_0_R
                        F2 = C_0_R
                        FZ = C_0_R
                        if (abs(iflag) == 3 .or. abs(iflag) == 4) then ! mean
                            if (idir == 1) then      ! OX direction; only acoustic, with un
                                F1 = F1 - pl_inf*c*((p(i) - r(i)*c*un(i)) - (bf(i, 5) - r(i)*c*bf(i, 2)))
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                FZ = FZ - pl_inf*c*(z1(i) - bf_z(i))
                            else if (idir == 2) then ! OY direction
                                F1 = F1 - pl_inf*c*(p(i) - bf(i, 5))
                                F2 = F2 - pl_inf*c*(r(i) - bf(i, 1))
                                FZ = FZ - pl_inf*c*(z1(i) - bf_z(i))
                            end if
                        end if

                        dummy = F2 + C_05_R*F1/c/c

                        hz1(i) = hz1(i) + r(i)*FZ + z1(i)*dummy

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        F1 = C_0_R
! Poinsot & Lele term. Multiplication by c is done in the definition
! of dummy below
                        if (idir == 1) then      ! OX treatment; no un forcing w/ F1
                            F1 = -pl_out*(p(i) - bf(i, 5))
                        else if (idir == 2) then ! OY treatment; no un forcing w/ F1
                            F1 = -pl_out*(p(i) - bf(i, 5))
                        end if
                        dummy = C_05_R*(r(i)*(C_1_R - Mn)*dundn(i) - (C_1_R - Mn)/c*dpdn(i) + r(i)*gn/c + F1/c)

                        hz1(i) = dummy*z1(i)

                    end if

                end if ! subsonic branch

            end do

        end if

        return
    end subroutine BOUNDARY_BCS_SCAL_NR_3
!########################################################################
!#
!# Implementation of the transverse terms of nonreflective boundary
!# conditions for the transport equations of density, momentum and
!# internal energy.
!#
!# After Lodato et al, JCP 227 (2008), 5105-5143
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!#
!########################################################################
    subroutine BOUNDARY_BCS_FLOW_NR_4(iflag, idir, nt, beta, &
                                      r, un, v1, v2, p, gama, t1, t2, t3, t4, t5, m1, m5, hr, hun, hv1, hv2, he)

        TINTEGER iflag, idir, nt
        TREAL beta

        TREAL, dimension(*) :: r, un, v1, v2, p, gama
        TREAL, dimension(*) :: t1, t2, t3, t4, t5, m1, m5
        TREAL, dimension(*) :: hr, hun, hv1, hv2, he

! -------------------------------------------------------------------
        TINTEGER i
        TREAL c, Mn, dummy

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) > C_0_R) then
                        dummy = C_05_R*t5(i)/c/c - C_05_R*r(i)*t2(i)/c - t1(i)

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) + C_05_R*(Mn - C_1_R)*t5(i)/c - C_05_R*r(i)*(Mn + C_1_R)*t2(i) - t1(i)*un(i)
                        hv1(i) = hv1(i) + dummy*v1(i) - r(i)*t3(i)
                        hv2(i) = hv2(i) + dummy*v2(i) - r(i)*t4(i)
                        he(i) = he(i) - C_05_R*(t5(i) + r(i)*c*t2(i))/(gama(i) - C_1_R)

! recover lateral term for v1 velocity at inflow in Ox
                        if (idir == 1 .or. idir == 2) then
                            hv1(i) = hv1(i) - C_05_R*(m5(i) - m1(i))/c
                        end if

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = -C_05_R*(C_1_R - beta)*(r(i)*c*t2(i) + t5(i))/c/c

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) + dummy*c*(C_1_R + Mn)
                        hv1(i) = hv1(i) + dummy*v1(i)
                        hv2(i) = hv2(i) + dummy*v2(i)
                        he(i) = he(i) + dummy*c*c/(gama(i) - C_1_R)

                    end if

                end if ! subsonic branch

            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) < C_0_R) then
                        dummy = C_05_R*t5(i)/c/c + C_05_R*r(i)*t2(i)/c - t1(i)

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) + C_05_R*(Mn + C_1_R)*t5(i)/c + C_05_R*r(i)*(Mn - C_1_R)*t2(i) - t1(i)*un(i)
                        hv1(i) = hv1(i) + dummy*v1(i) - r(i)*t3(i)
                        hv2(i) = hv2(i) + dummy*v2(i) - r(i)*t4(i)
                        he(i) = he(i) - C_05_R*(t5(i) - r(i)*c*t2(i))/(gama(i) - C_1_R)

! recover lateral term for v1 velocity at inflow in Ox
                        if (idir == 1 .or. idir == 2) then
                            hv1(i) = hv1(i) - C_05_R*(m5(i) - m1(i))/c
                        end if

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = C_05_R*(C_1_R - beta)*(r(i)*c*t2(i) - t5(i))/c/c

                        hr(i) = hr(i) + dummy
                        hun(i) = hun(i) - dummy*c*(C_1_R - Mn)
                        hv1(i) = hv1(i) + dummy*v1(i)
                        hv2(i) = hv2(i) + dummy*v2(i)
                        he(i) = he(i) + dummy*c*c/(gama(i) - C_1_R)
                    end if

                end if ! subsonic branch

            end do

        end if

        return
    end subroutine BOUNDARY_BCS_FLOW_NR_4

    !########################################################################
!#
!# Implementation of the transverse terms of nonreflective boundary
!# conditions for the transport equations of density, momentum and
!# internal energy.
!#
!# After Lodato et al, JCP 227 (2008), 5105-5143
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!#
!########################################################################
    subroutine BOUNDARY_BCS_SCAL_NR_4(iflag, nt, beta, r, un, z1, p, gama, t1, t2, t5, tz1, hz1)

        TINTEGER iflag, nt
        TREAL beta

        TREAL, dimension(*) :: r, un, z1, p, gama
        TREAL, dimension(*) :: t1, t2, t5, tz1
        TREAL, dimension(*) :: hz1

! -------------------------------------------------------------------
        TINTEGER i
        TREAL c, Mn, dummy

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) > C_0_R) then
                        dummy = C_05_R*t5(i)/c/c - C_05_R*r(i)*t2(i)/c - t1(i)

                        hz1(i) = hz1(i) + dummy*z1(i) - r(i)*tz1(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = -C_05_R*(C_1_R - beta)*(r(i)*c*t2(i) + t5(i))/c/c

                        hz1(i) = hz1(i) + dummy*z1(i)

                    end if

                end if ! subsonic branch

            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do i = 1, nt
                c = sqrt(gama(i)*p(i)/r(i))
                Mn = un(i)/c

                if (un(i) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
                    if (un(i) < C_0_R) then
                        dummy = C_05_R*t5(i)/c/c + C_05_R*r(i)*t2(i)/c - t1(i)

                        hz1(i) = hz1(i) + dummy*z1(i) - r(i)*tz1(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
                    else
                        dummy = C_05_R*(C_1_R - beta)*(r(i)*c*t2(i) - t5(i))/c/c

                        hz1(i) = hz1(i) + dummy*z1(i)

                    end if

                end if ! subsonic branch

            end do

        end if

        return
    end subroutine BOUNDARY_BCS_SCAL_NR_4

!########################################################################
!#
!# Implementation of the transverse terms of nonreflective boundary
!# conditions for the transport equations of density, momentum and
!# internal energy.
!#
!# After Lodato et al, JCP 227 (2008), 5105-5143
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!#
!########################################################################
    subroutine BOUNDARY_BCS_FLOW_NR_EDGE(iflag, jmax, kmax, beta, &
                                         r, un, v1, v2, p, gama, m1, m2, m3, m4, m5, hr, hun, hv1, hv2, he)

        TINTEGER iflag, jmax, kmax
        TREAL beta

        TREAL, dimension(jmax, *) :: r, un, v1, v2, p, gama
        TREAL, dimension(jmax, *) :: m1, m2, m3, m4, m5
        TREAL, dimension(jmax, *) :: hr, hun, hv1, hv2, he

! -------------------------------------------------------------------
        TINTEGER j, k
        TREAL c, Mn, dummy, F1, F2, F3, F4, F5

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do k = 1, kmax
                do j = 1, jmax, jmax - 1
                    c = sqrt(gama(j, k)*p(j, k)/r(j, k))
                    Mn = un(j, k)/c

                    if (un(j, k) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow in Ox
! -------------------------------------------------------------------
                        if (un(j, k) > C_0_R) then
                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m5(j, k)
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_0_R
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m5(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_0_R
                                end if
                            end if
                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m1(j, k)
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_0_R
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m1(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_0_R
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hr(j, k) = hr(j, k) + dummy
                            hun(j, k) = hun(j, k) + dummy*un(j, k) + (F5 - F1)*C_05_R/c
                            hv1(j, k) = hv1(j, k) + dummy*v1(j, k) + r(j, k)*F3
                            hv2(j, k) = hv2(j, k) + dummy*v2(j, k) + r(j, k)*F4
                            he(j, k) = he(j, k) + C_05_R*(F1 + F5)/(gama(j, k) - C_1_R)

! -------------------------------------------------------------------
! Outflow in Ox
! -------------------------------------------------------------------
                        else
                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m5(j, k)
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = beta*C_05_R*m5(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m5(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = m3(j, k)
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = m4(j, k)
                                    F5 = beta*(C_05_R*m5(j, k) + r(j, k)*c*m2(j, k))
                                end if
                            end if

                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m1(j, k)
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = beta*C_05_R*m1(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m1(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = m3(j, k)
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = m4(j, k)
                                    F5 = beta*(C_05_R*m1(j, k) + r(j, k)*c*m2(j, k))
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hr(j, k) = hr(j, k) + dummy
                            hun(j, k) = hun(j, k) + dummy*un(j, k) + (F5 - F1)*C_05_R/c
                            hv1(j, k) = hv1(j, k) + dummy*v1(j, k) + r(j, k)*F3
                            hv2(j, k) = hv2(j, k) + dummy*v2(j, k) + r(j, k)*F4
                            he(j, k) = he(j, k) + C_05_R*(F1 + F5)/(gama(j, k) - C_1_R)

                        end if

                    end if ! subsonic branch

                end do
            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do k = 1, kmax
                do j = 1, jmax, jmax - 1
                    c = sqrt(gama(j, k)*p(j, k)/r(j, k))
                    Mn = un(j, k)/c

                    if (un(j, k) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow in Ox
! -------------------------------------------------------------------
                        if (un(j, k) < C_0_R) then
                            dummy = C_0_R

                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F3 = C_0_R
                                    F4 = C_0_R
                                    F5 = C_05_R*m5(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F3 = C_0_R
                                    F4 = C_0_R
                                    F5 = C_05_R*m5(j, k) + r(j, k)*c*m2(j, k)
                                end if
                            end if
                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F3 = C_0_R
                                    F4 = C_0_R
                                    F5 = C_05_R*m1(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F3 = C_0_R
                                    F4 = C_0_R
                                    F5 = C_05_R*m1(j, k) + r(j, k)*c*m2(j, k)
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hr(j, k) = hr(j, k) + dummy
                            hun(j, k) = hun(j, k) + dummy*un(j, k) + (F5 - F1)*C_05_R/c
                            hv1(j, k) = hv1(j, k) + dummy*v1(j, k) + r(j, k)*F3
                            hv2(j, k) = hv2(j, k) + dummy*v2(j, k) + r(j, k)*F4
                            he(j, k) = he(j, k) + C_05_R*(F1 + F5)/(gama(j, k) - C_1_R)

! -------------------------------------------------------------------
! Outflow in Ox
! -------------------------------------------------------------------
                        else
                            if (j == 1) then
                                if (v1(j, k) > C_0_R) then ! Inflow in Oy
                                    F1 = beta*(C_05_R*m5(j, k) - r(j, k)*c*m2(j, k))
                                    F2 = m3(j, k)
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = m4(j, k)
                                    F5 = C_05_R*m5(j, k) + r(j, k)*c*m2(j, k)
                                else                           ! Outflow in Oy
                                    F1 = beta*C_05_R*m5(j, k)
                                    F2 = C_0_R
                                    F3 = C_05_R*m5(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_05_R*m5(j, k)
                                end if
                            end if

                            if (j == jmax) then
                                if (v1(j, k) < C_0_R) then ! Inflow in Oy
                                    F1 = beta*(C_05_R*m1(j, k) - r(j, k)*c*m2(j, k))
                                    F2 = m3(j, k)
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = m4(j, k)
                                    F5 = C_05_R*m1(j, k) + r(j, k)*c*m2(j, k)
                                else                           ! Outflow in Oy
                                    F1 = beta*C_05_R*m1(j, k)
                                    F2 = C_0_R
                                    F3 = -C_05_R*m1(j, k)/r(j, k)/c
                                    F4 = C_0_R
                                    F5 = C_05_R*m1(j, k)
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hr(j, k) = hr(j, k) + dummy
                            hun(j, k) = hun(j, k) + dummy*un(j, k) + (F5 - F1)*C_05_R/c
                            hv1(j, k) = hv1(j, k) + dummy*v1(j, k) + r(j, k)*F3
                            hv2(j, k) = hv2(j, k) + dummy*v2(j, k) + r(j, k)*F4
                            he(j, k) = he(j, k) + C_05_R*(F1 + F5)/(gama(j, k) - C_1_R)

                        end if

                    end if ! subsonic branch

                end do
            end do

        end if

        return
    end subroutine BOUNDARY_BCS_FLOW_NR_EDGE

    !########################################################################
!#
!# Implementation of the transverse terms of nonreflective boundary
!# conditions for the transport equations of density, momentum and
!# internal energy.
!#
!# After Lodato et al, JCP 227 (2008), 5105-5143
!#
!########################################################################
!# ARGUMENTS
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!#
!########################################################################
    subroutine BOUNDARY_BCS_SCAL_NR_EDGE(iflag, jmax, kmax, beta, &
                                         r, un, v1, z1, p, gama, m1, m2, m3, m5, m6, hz1)

        TINTEGER iflag, jmax, kmax
        TREAL beta

        TREAL, dimension(jmax, *) :: r, un, v1, z1, p, gama
        TREAL, dimension(jmax, *) :: m1, m2, m3, m5, m6
        TREAL, dimension(jmax, *) :: hz1

! -------------------------------------------------------------------
        TINTEGER j, k
        TREAL c, Mn, dummy, F1, F2, F5, F6

! ###################################################################
! BCs at x_min
! ###################################################################
        if (iflag < 0) then

            do k = 1, kmax
                do j = 1, jmax, jmax - 1
                    c = sqrt(gama(j, k)*p(j, k)/r(j, k))
                    Mn = un(j, k)/c

                    if (un(j, k) + c > C_0_R) then

! -------------------------------------------------------------------
! Inflow in Ox
! -------------------------------------------------------------------
                        if (un(j, k) > C_0_R) then
                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m5(j, k)
                                    F2 = C_0_R
                                    F5 = C_0_R
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m5(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = C_0_R
                                    F5 = C_0_R
                                end if
                            end if
                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                    F1 = C_05_R*m1(j, k)
                                    F2 = C_0_R
                                    F5 = C_0_R
                                else                           ! & Inflow in Oy
                                    F1 = C_05_R*m1(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = C_0_R
                                    F5 = C_0_R
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hz1(j, k) = hz1(j, k) + dummy*z1(j, k)

! -------------------------------------------------------------------
! Outflow in Ox
! -------------------------------------------------------------------
                        else
                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                else ! & Inflow in Oy
                                    F1 = C_05_R*m5(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = m3(j, k)
                                    F5 = beta*(C_05_R*m5(j, k) + r(j, k)*c*m2(j, k))
                                    F6 = m6(j, k)
                                end if
                            end if

                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                else ! & Inflow in Oy
                                    F1 = C_05_R*m1(j, k) - r(j, k)*c*m2(j, k)
                                    F2 = m3(j, k)
                                    F5 = beta*(C_05_R*m1(j, k) + r(j, k)*c*m2(j, k))
                                    F6 = m6(j, k)
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hz1(j, k) = hz1(j, k) + dummy*z1(j, k) + r(j, k)*F6

                        end if

                    end if ! subsonic branch

                end do
            end do

! ###################################################################
! BCs at x_max
! ###################################################################
        else if (iflag > 0) then

            do k = 1, kmax
                do j = 1, jmax, jmax - 1
                    c = sqrt(gama(j, k)*p(j, k)/r(j, k))
                    Mn = un(j, k)/c

                    if (un(j, k) - c < C_0_R) then

! -------------------------------------------------------------------
! Inflow in Ox
! -------------------------------------------------------------------
                        if (un(j, k) < C_0_R) then
                            dummy = C_0_R

                            if (j == 1) then
                                if (v1(j, k) < C_0_R) then ! & Outflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F5 = C_05_R*m5(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F5 = C_05_R*m5(j, k) + r(j, k)*c*m2(j, k)
                                end if
                            end if
                            if (j == jmax) then
                                if (v1(j, k) > C_0_R) then ! & Outflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F5 = C_05_R*m1(j, k)
                                else                           ! & Inflow in Oy
                                    F1 = C_0_R
                                    F2 = C_0_R
                                    F5 = C_05_R*m1(j, k) + r(j, k)*c*m2(j, k)
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hz1(j, k) = hz1(j, k) + dummy*z1(j, k)

! -------------------------------------------------------------------
! Outflow in Ox
! -------------------------------------------------------------------
                        else
                            if (j == 1) then
                                if (v1(j, k) > C_0_R) then ! Inflow in Oy
                                    F1 = beta*(C_05_R*m5(j, k) - r(j, k)*c*m2(j, k))
                                    F2 = m3(j, k)
                                    F5 = C_05_R*m5(j, k) + r(j, k)*c*m2(j, k)
                                    F6 = m6(j, k)
                                else
                                end if
                            end if

                            if (j == jmax) then
                                if (v1(j, k) < C_0_R) then ! Inflow in Oy
                                    F1 = beta*(C_05_R*m1(j, k) - r(j, k)*c*m2(j, k))
                                    F2 = m3(j, k)
                                    F5 = C_05_R*m1(j, k) + r(j, k)*c*m2(j, k)
                                    F6 = m6(j, k)
                                else
                                end if
                            end if
                            dummy = (F2 + C_05_R*(F1 + F5))/c/c

                            hz1(j, k) = hz1(j, k) + dummy*z1(j, k) + r(j, k)*F6

                        end if

                    end if ! subsonic branch

                end do
            end do

        end if

        return
    end subroutine BOUNDARY_BCS_SCAL_NR_EDGE

!########################################################################
!#
!# Calculate transverse terms at Ox_min (array tmin) and at Ox_max (array tmax)
!# according to Lodato et al, JCP 227 (2008), 5105-5143
!# The sign is the opposite to that paper
!#
!########################################################################
!# ARGUMENTS
!#
!# tmin    In    Transverse term at OxMin
!# tmax    In    Transverse term at OxMax
!#
!########################################################################
    subroutine BOUNDARY_BCS_TRANSVERSE_X(u, v, w, p, r, gamma, z1, &
                                         tmin, mmin, tmax, mmax, tmp1, ddy, ddz, wrk2d, wrk3d)

#include "integers.h"

        TREAL, dimension(imax, jmax, kmax) :: u, v, w, p, r, gamma
#ifdef USE_MPI
        TREAL, dimension(ims_bcs_imax, jmax, kmax) :: tmp1, ddy, ddz
#else
        TREAL, dimension(2*(inb_flow + inb_scal_array), jmax, kmax) :: tmp1, ddy, ddz
#endif
        TREAL, dimension(imax, jmax, kmax, *) :: z1
        TREAL, dimension(jmax, kmax, *) :: tmin, tmax, mmin, mmax

        TREAL, dimension(*) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
        TINTEGER ip, j, k, is, bcs(2, 2)
        TREAL c

! #######################################################################
        bcs = 0

! -------------------------------------------------------------------
! Arrange data
! -------------------------------------------------------------------
        ip = 0

! BCs at x_min
        do k = 1, kmax
            do j = 1, jmax
                tmp1(ip + 1, j, k) = u(1, j, k)
                tmp1(ip + 2, j, k) = v(1, j, k)
                tmp1(ip + 3, j, k) = w(1, j, k)
                tmp1(ip + 4, j, k) = p(1, j, k)
                tmp1(ip + 5, j, k) = r(1, j, k)
                do is = 1, inb_scal_array
                    tmp1(ip + 5 + is, j, k) = z1(1, j, k, is)
                end do
            end do
        end do
        ip = ip + inb_flow + inb_scal_array

! BCs at x_max
        do k = 1, kmax
            do j = 1, jmax
                tmp1(ip + 1, j, k) = u(imax, j, k)
                tmp1(ip + 2, j, k) = v(imax, j, k)
                tmp1(ip + 3, j, k) = w(imax, j, k)
                tmp1(ip + 4, j, k) = p(imax, j, k)
                tmp1(ip + 5, j, k) = r(imax, j, k)
                do is = 1, inb_scal_array
                    tmp1(ip + 5 + is, j, k) = z1(imax, j, k, is)
                end do
            end do
        end do
        ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Construct t1-t5
! -------------------------------------------------------------------
#ifdef USE_MPI
        call OPR_PARTIAL_Y(OPR_P1, ims_bcs_imax, jmax, kmax, bcs, g(2), tmp1, ddy, wrk3d, wrk2d, wrk3d)
! Needs to be checked
        call TLAB_WRITE_ASCII(efile, 'BOUNDARY_BCS_TRANSVERSE_X. To be checked')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
!  imode_fdm_loc = imode_fdm + (TLAB_MPI_K_NRBCX-1)*100
        call OPR_PARTIAL_Z(OPR_P1_BCS, ims_bcs_imax, jmax, kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d, wrk3d)
#else
        call OPR_PARTIAL_Y(OPR_P1, ip, jmax, kmax, bcs, g(2), tmp1, ddy, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, ip, jmax, kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d, wrk3d)
#endif
        ip = 0

! BCs at x_min
        do k = 1, kmax
            do j = 1, jmax
                tmin(j, k, 1) = r(1, j, k)*ddy(ip + 2, j, k) + v(1, j, k)*ddy(ip + 5, j, k) &
                                + r(1, j, k)*ddz(ip + 3, j, k) + w(1, j, k)*ddz(ip + 5, j, k)
                tmin(j, k, 2) = v(1, j, k)*ddy(ip + 1, j, k) + w(1, j, k)*ddz(ip + 1, j, k)
                tmin(j, k, 3) = v(1, j, k)*ddy(ip + 2, j, k) + w(1, j, k)*ddz(ip + 2, j, k) &
                                + ddy(ip + 4, j, k)/r(1, j, k) - buoyancy%vector(2)
                tmin(j, k, 4) = v(1, j, k)*ddy(ip + 3, j, k) + w(1, j, k)*ddz(ip + 3, j, k) &
                                + ddz(ip + 4, j, k)/r(1, j, k) - buoyancy%vector(3)
                tmin(j, k, 5) = v(1, j, k)*ddy(ip + 4, j, k) + w(1, j, k)*ddz(ip + 4, j, k) &
                                + gamma(1, j, k)*p(1, j, k)*(ddy(ip + 2, j, k) + ddz(ip + 3, j, k))
                do is = 1, inb_scal_array
                    tmin(j, k, 5 + is) = v(1, j, k)*ddy(ip + 5 + is, j, k) + w(1, j, k)*ddz(ip + 5 + is, j, k)
                end do
! constructing M1-M5
                c = sqrt(gamma(1, j, k)*p(1, j, k)/r(1, j, k))
                mmin(j, k, 1) = (v(1, j, k) - c)*(ddy(ip + 4, j, k) - ddy(ip + 2, j, k)*r(1, j, k)*c)
                mmin(j, k, 2) = v(1, j, k)*(ddy(ip + 1, j, k))
                mmin(j, k, 3) = v(1, j, k)*(ddy(ip + 5, j, k)*c*c - ddy(ip + 4, j, k))
                mmin(j, k, 4) = v(1, j, k)*(ddy(ip + 3, j, k))
                mmin(j, k, 5) = (v(1, j, k) + c)*(ddy(ip + 4, j, k) + ddy(ip + 2, j, k)*r(1, j, k)*c)
                do is = 1, inb_scal_array
                    mmin(j, k, 5 + is) = v(1, j, k)*ddy(ip + 5 + is, j, k)
                end do
            end do
        end do
        ip = ip + inb_flow + inb_scal_array

! BCs at x_max
        do k = 1, kmax
            do j = 1, jmax
                tmax(j, k, 1) = r(imax, j, k)*ddy(ip + 2, j, k) + v(imax, j, k)*ddy(ip + 5, j, k) &
                                + r(imax, j, k)*ddz(ip + 3, j, k) + w(imax, j, k)*ddz(ip + 5, j, k)
                tmax(j, k, 2) = v(imax, j, k)*ddy(ip + 1, j, k) + w(imax, j, k)*ddz(ip + 1, j, k)
                tmax(j, k, 3) = v(imax, j, k)*ddy(ip + 2, j, k) + w(imax, j, k)*ddz(ip + 2, j, k) &
                                + ddy(ip + 4, j, k)/r(imax, j, k) - buoyancy%vector(2)
                tmax(j, k, 4) = v(imax, j, k)*ddy(ip + 3, j, k) + w(imax, j, k)*ddz(ip + 3, j, k) &
                                + ddz(ip + 4, j, k)/r(imax, j, k) - buoyancy%vector(3)
                tmax(j, k, 5) = v(imax, j, k)*ddy(ip + 4, j, k) + w(imax, j, k)*ddz(ip + 4, j, k) &
                                + gamma(imax, j, k)*p(imax, j, k)*(ddy(ip + 2, j, k) + ddz(ip + 3, j, k))
                do is = 1, inb_scal_array
                    tmax(j, k, 5 + is) = v(imax, j, k)*ddy(ip + 5 + is, j, k) + w(imax, j, k)*ddz(ip + 5 + is, j, k)
                end do
! constructing M1-M5
                c = sqrt(gamma(imax, j, k)*p(imax, j, k)/r(imax, j, k))
                mmax(j, k, 1) = (v(imax, j, k) - c)*(ddy(ip + 4, j, k) - ddy(ip + 2, j, k)*r(imax, j, k)*c)
                mmax(j, k, 2) = v(imax, j, k)*(ddy(ip + 1, j, k))
                mmax(j, k, 3) = v(imax, j, k)*(ddy(ip + 5, j, k)*c*c - ddy(ip + 4, j, k))
                mmax(j, k, 4) = v(imax, j, k)*(ddy(ip + 3, j, k))
                mmax(j, k, 5) = (v(imax, j, k) + c)*(ddy(ip + 4, j, k) + ddy(ip + 2, j, k)*r(imax, j, k)*c)
                do is = 1, inb_scal_array
                    mmax(j, k, 5 + is) = v(imax, j, k)*ddy(ip + 5 + is, j, k)
                end do
            end do
        end do
        ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Change sign
! -------------------------------------------------------------------
        do is = 1, inb_flow + inb_scal_array
            do k = 1, kmax
                do j = 1, jmax
                    tmin(j, k, is) = -tmin(j, k, is)
                    tmax(j, k, is) = -tmax(j, k, is)
                end do
            end do
        end do

        return
    end subroutine BOUNDARY_BCS_TRANSVERSE_X

!########################################################################
!#
!# Calculate transverse terms at Oy_min (array tmin) and at Oy_max (array tmax)
!# according to Lodato et al, JCP 227 (2008), 5105-5143.
!# The sign is the opposite to that paper
!#
!########################################################################
!# ARGUMENTS
!#
!# tmin    In    Transverse term at OyMin
!# tmin    In    Transverse term at OyMin
!#
!########################################################################
    subroutine BOUNDARY_BCS_TRANSVERSE_Y(u, v, w, p, r, gamma, z1, &
                                         tmin, lmin, tmax, lmax, tmp1, ddx, ddz, wrk2d, wrk3d)

#include "integers.h"

        TREAL, dimension(imax, jmax, kmax) :: u, v, w, p, r, gamma
#ifdef USE_MPI
        TREAL, dimension(imax, ims_bcs_jmax, kmax) :: tmp1, ddx, ddz
#else
        TREAL, dimension(imax, 2*(inb_flow + inb_scal_array), kmax) :: tmp1, ddx, ddz
#endif
        TREAL, dimension(imax, jmax, kmax, *) :: z1
        TREAL, dimension(imax, kmax, *) :: tmin, lmin, tmax, lmax

        TREAL, dimension(*) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
        TINTEGER ip, i, k, is, bcs(2, 2)
        TREAL c

! #######################################################################
        bcs = 0

! -------------------------------------------------------------------
! Arrange data
! -------------------------------------------------------------------
        ip = 0

! BCs at y_min
        do k = 1, kmax; do i = 1, imax
                tmp1(i, ip + 1, k) = u(i, 1, k)
                tmp1(i, ip + 2, k) = v(i, 1, k)
                tmp1(i, ip + 3, k) = w(i, 1, k)
                tmp1(i, ip + 4, k) = p(i, 1, k)
                tmp1(i, ip + 5, k) = r(i, 1, k)
                do is = 1, inb_scal_array
                    tmp1(i, ip + 5 + is, k) = z1(i, 1, k, is)
                end do
            end do; end do
        ip = ip + inb_flow + inb_scal_array

! BCs at y_max
        do k = 1, kmax; do i = 1, imax
                tmp1(i, ip + 1, k) = u(i, jmax, k)
                tmp1(i, ip + 2, k) = v(i, jmax, k)
                tmp1(i, ip + 3, k) = w(i, jmax, k)
                tmp1(i, ip + 4, k) = p(i, jmax, k)
                tmp1(i, ip + 5, k) = r(i, jmax, k)
                do is = 1, inb_scal_array
                    tmp1(i, ip + 5 + is, k) = z1(i, jmax, k, is)
                end do
            end do; end do
        ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Construct t1-t5
! -------------------------------------------------------------------
#ifdef USE_MPI
        call OPR_PARTIAL_X(OPR_P1, imax, ims_bcs_jmax, kmax, bcs, g(1), tmp1, ddx, wrk3d, wrk2d, wrk3d)
! Needs to be checked
        call TLAB_WRITE_ASCII(efile, 'BOUNDARY_BCS_TRANSVERSE_Y. To be checked')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
!  imode_fdm_loc = imode_fdm + (TLAB_MPI_K_NRBCY-1)*100
        call OPR_PARTIAL_Z(OPR_P1_BCS, imax, ims_bcs_jmax, kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d, wrk3d)
#else
        call OPR_PARTIAL_X(OPR_P1, imax, ip, kmax, bcs, g(1), tmp1, ddx, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, imax, ip, kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d, wrk3d)
#endif
        ip = 0

! BCs at y_min
        do k = 1, kmax; do i = 1, imax
                tmin(i, k, 1) = r(i, 1, k)*ddx(i, ip + 1, k) + u(i, 1, k)*ddx(i, ip + 5, k) &
                                + r(i, 1, k)*ddz(i, ip + 3, k) + w(i, 1, k)*ddz(i, ip + 5, k)
                tmin(i, k, 2) = u(i, 1, k)*ddx(i, ip + 1, k) + w(i, 1, k)*ddz(i, ip + 1, k) &
                                + ddx(i, ip + 4, k)/r(i, 1, k) - buoyancy%vector(1)
                tmin(i, k, 3) = u(i, 1, k)*ddx(i, ip + 2, k) + w(i, 1, k)*ddz(i, ip + 2, k)
                tmin(i, k, 4) = u(i, 1, k)*ddx(i, ip + 3, k) + w(i, 1, k)*ddz(i, ip + 3, k) &
                                + ddz(i, ip + 4, k)/r(i, 1, k) - buoyancy%vector(3)
                tmin(i, k, 5) = u(i, 1, k)*ddx(i, ip + 4, k) + w(i, 1, k)*ddz(i, ip + 4, k) &
                                + gamma(i, 1, k)*p(i, 1, k)*(ddx(i, ip + 1, k) + ddz(i, ip + 3, k))
                do is = 1, inb_scal_array
                    tmin(i, k, 5 + is) = u(i, 1, k)*ddx(i, ip + 5 + is, k) + w(i, 1, k)*ddz(i, ip + 5 + is, k)
                end do
! constructing L1-L5
                c = sqrt(gamma(i, 1, k)*p(i, 1, k)/r(i, 1, k))
                lmin(i, k, 1) = (u(i, 1, k) - c)*(ddx(i, ip + 4, k) - ddx(i, ip + 1, k)*r(i, 1, k)*c)
                lmin(i, k, 2) = u(i, 1, k)*(ddx(i, ip + 5, k)*c*c - ddx(i, ip + 4, k))
                lmin(i, k, 3) = u(i, 1, k)*(ddx(i, ip + 2, k))
                lmin(i, k, 4) = u(i, 1, k)*(ddx(i, ip + 3, k))
                lmin(i, k, 5) = (u(i, 1, k) + c)*(ddx(i, ip + 4, k) + ddx(i, ip + 1, k)*r(i, 1, k)*c)
                do is = 1, inb_scal_array
                    lmin(i, k, 5 + is) = u(i, 1, k)*ddx(i, ip + 5 + is, k)
                end do
            end do; end do
        ip = ip + inb_flow + inb_scal_array

! BCs at y_max
        do k = 1, kmax; do i = 1, imax
                tmax(i, k, 1) = r(i, jmax, k)*ddx(i, ip + 1, k) + u(i, jmax, k)*ddx(i, ip + 5, k) &
                                + r(i, jmax, k)*ddz(i, ip + 3, k) + w(i, jmax, k)*ddz(i, ip + 5, k)
                tmax(i, k, 2) = u(i, jmax, k)*ddx(i, ip + 1, k) + w(i, jmax, k)*ddz(i, ip + 1, k) &
                                + ddx(i, ip + 4, k)/r(i, jmax, k) - buoyancy%vector(1)
                tmax(i, k, 3) = u(i, jmax, k)*ddx(i, ip + 2, k) + w(i, jmax, k)*ddz(i, ip + 2, k)
                tmax(i, k, 4) = u(i, jmax, k)*ddx(i, ip + 3, k) + w(i, jmax, k)*ddz(i, ip + 3, k) &
                                + ddz(i, ip + 4, k)/r(i, jmax, k) - buoyancy%vector(3)
                tmax(i, k, 5) = u(i, jmax, k)*ddx(i, ip + 4, k) + w(i, jmax, k)*ddz(i, ip + 4, k) &
                                + gamma(i, jmax, k)*p(i, jmax, k)*(ddx(i, ip + 1, k) + ddz(i, ip + 3, k))
                do is = 1, inb_scal_array
                    tmax(i, k, 5 + is) = u(i, jmax, k)*ddx(i, ip + 5 + is, k) + w(i, jmax, k)*ddz(i, ip + 5 + is, k)
                end do
! constructing L1-L5
                c = sqrt(gamma(i, jmax, k)*p(i, jmax, k)/r(i, jmax, k))
                lmax(i, k, 1) = (u(i, jmax, k) - c)*(ddx(i, ip + 4, k) - ddx(i, ip + 1, k)*r(i, jmax, k)*c)
                lmax(i, k, 2) = u(i, jmax, k)*(ddx(i, ip + 5, k)*c*c - ddx(i, ip + 4, k))
                lmax(i, k, 3) = u(i, jmax, k)*(ddx(i, ip + 2, k))
                lmax(i, k, 4) = u(i, jmax, k)*(ddx(i, ip + 3, k))
                lmax(i, k, 5) = (u(i, jmax, k) + c)*(ddx(i, ip + 4, k) + ddx(i, ip + 1, k)*r(i, jmax, k)*c)
                do is = 1, inb_scal_array
                    lmax(i, k, 5 + is) = u(i, jmax, k)*ddx(i, ip + 5 + is, k)
                end do
            end do; end do
        ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Change sign
! -------------------------------------------------------------------
        do is = 1, inb_flow + inb_scal_array
            do k = 1, kmax; do i = 1, imax
                    tmin(i, k, is) = -tmin(i, k, is)
                    tmax(i, k, is) = -tmax(i, k, is)
                end do; end do
        end do

        return
    end subroutine BOUNDARY_BCS_TRANSVERSE_Y

end module BOUNDARY_BCS_COMPRESSIBLE

#include "dns_const.h"
#include "dns_error.h"
#include "types.h"

!########################################################################
!#
!# Runge-Kutta explicit 3th order from Williamson 1980
!# Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
!# Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
!#
!########################################################################
module TIME

#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLAB_CONSTANTS, only: efile, wp, wi
    use TLAB_VARS
    use THERMO_VARS, only: gama0
    use TLAB_PROCS
    use PARTICLE_VARS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
#endif
    implicit none
    save
    private

    integer(wi), public :: rkm_mode             ! Type of Runge-Kutta scheme
    integer(wi), public :: rkm_endstep          ! number of substeps
    integer(wi), public :: rkm_substep          ! substep counter

    real(wp), public :: cfla, cfld, cflr        ! CFL numbers
    real(wp), public :: dtime                   ! time step
    real(wp), public :: dte                     ! time step of each substep
    real(wp) etime                              ! time at each substep

    real(wp) kdt(5), kco(4), ktime(5)           ! explicit scheme coefficients
    real(wp) kex(3), kim(3)                     ! implicit scheme coefficients

    real(wp) schmidtfactor, dx2i
    integer(wi) i, j, k, kdsp, idsp, is
    real(wp) dummy

    public :: TIME_INITIALIZE
    public :: TIME_RUNGEKUTTA
    public :: TIME_COURANT

contains

! ###################################################################
! ###################################################################
    subroutine TIME_INITIALIZE()

        ! ###################################################################
        ! RK coefficients
        select case (rkm_mode)
        case (RKM_EXP3)             ! Runge-Kutta explicit 3th order from Williamson 1980
            rkm_endstep = 3

            kdt(1:3) = [C_1_R/C_3_R, C_15_R/C_16_R, C_8_R/C_15_R]
            ktime(1:3) = [C_0_R, C_1_R/C_3_R, C_3_R/C_4_R]
            kco(1:2) = [-C_5_R/C_9_R, -C_153_R/C_128_R]

        case (RKM_EXP4)             ! Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
            rkm_endstep = 5

            kdt(1) = C_1432997174477_R/C_9575080441755_R
            kdt(2) = C_5161836677717_R/C_13612068292357_R
            kdt(3) = C_1720146321549_R/C_2090206949498_R
            kdt(4) = C_3134564353537_R/C_4481467310338_R
            kdt(5) = C_2277821191437_R/C_14882151754819_R

            ktime(1) = C_0_R
            ktime(2) = C_1432997174477_R/C_9575080441755_R
            ktime(3) = C_2526269341429_R/C_6820363962896_R
            ktime(4) = C_2006345519317_R/C_3224310063776_R
            ktime(5) = C_2802321613138_R/C_2924317926251_R

            kco(1) = -C_567301805773_R/C_1357537059087_R
            kco(2) = -C_2404267990393_R/C_2016746695238_R
            kco(3) = -C_3550918686646_R/C_2091501179385_R
            kco(4) = -C_1275806237668_R/C_842570457699_R

        case (RKM_IMP3_DIFFUSION)   ! Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
            rkm_endstep = 3

            kdt(1:3) = (/C_8_R/C_15_R, C_5_R/C_12_R, C_3_R/C_4_R/)

            kim(1:3) = (/C_111_R/C_256_R, C_1_R/C_2_R, C_2_R/C_9_R/)
            kex(1:3) = (/C_145_R/C_256_R, -C_9_R/C_50_R, C_2_R/C_9_R/)
            kco(1:3) = (/C_0_R, -C_17_R/C_25_R, -C_5_R/C_9_R/)
            ! TO DO - calculate ktime from coefficients  ktime
            ktime(1:3) = (/C_0_R, C_0_R, C_0_R/)

            ! Coefficients from Spalart, Moser, Rogers (1991)
            ! kim = beta/gamma
            ! kex = alpha/gamma
            ! kco = zeta/gamma
            !
            ! alpha = (/ 29./96.,   -3./40,    1./6. /)
            ! beta  = (/ 37./160.,   5./24.,   1./6. /)
            ! gamma = (/  8./15.,    5./12.,   3./4. /)
            ! zeta  = (/  0.,      -17./60.,  -5./12./)

        end select

        ! ###################################################################
        ! maximum diffusivities for TIME_COURANT
        schmidtfactor = C_1_R
        dummy = C_1_R/prandtl
        schmidtfactor = max(schmidtfactor, dummy)
        dummy = C_1_R/minval(schmidt(1:inb_scal))
        schmidtfactor = max(schmidtfactor, dummy)
        schmidtfactor = schmidtfactor*visc

        ! -------------------------------------------------------------------
        ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TIME_COURANT
#ifdef USE_MPI
        idsp = ims_offset_i
        kdsp = ims_offset_k
#else
        idsp = 0
        kdsp = 0
#endif

        dx2i = C_0_R
        if (g(3)%size > 1) then
            do k = 1, kmax; do j = 1, jmax; do i = 1, imax
                    dummy = g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4) + g(3)%jac(k + kdsp, 4)
                    dx2i = max(dx2i, dummy)
                end do; end do; end do
        else
            do k = 1, kmax; do j = 1, jmax; do i = 1, imax
                    dummy = g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4)
                    dx2i = max(dx2i, dummy)
                end do; end do; end do
        end if

        return
    end subroutine TIME_INITIALIZE

! ###################################################################
! ###################################################################
    subroutine TIME_RUNGEKUTTA()
        use TLAB_ARRAYS
        use PARTICLE_ARRAYS
        use DNS_LOCAL
        use DNS_ARRAYS

        ! -------------------------------------------------------------------
        real(wp) alpha

        integer(wi) ij_srt, ij_end, ij_siz ! Variables for OpenMP Paritioning
#ifdef USE_PROFILE
        integer(wi) t_srt, t_end, t_dif, idummy, PROC_CYCLES, MAX_CYCLES
        character*256 time_string
#endif

#ifdef USE_BLAS
        integer ij_len
#endif

        !########################################################################
#ifdef USE_BLAS
        ij_len = isize_field
#endif

        ! -------------------------------------------------------------------
        ! Initialize arrays to zero for the explcit low-storage algorithm
        ! -------------------------------------------------------------------
        if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
            if (icalc_flow == 1) hq = C_0_R
            if (icalc_scal == 1) hs = C_0_R
            if (part%type /= PART_TYPE_NONE) l_hq = C_0_R
        end if

        !########################################################################
        ! Loop over the sub-stages
        !########################################################################
        do rkm_substep = 1, rkm_endstep

            ! -------------------------------------------------------------------
            ! Update transported (or prognostic) variables q and s
            ! -------------------------------------------------------------------
            dte = dtime*kdt(rkm_substep)
            etime = rtime + dtime*ktime(rkm_substep)

#ifdef USE_PROFILE
            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
#endif
            if (part%type /= PART_TYPE_NONE) then
                call TIME_SUBSTEP_PARTICLE()
            end if

            select case (imode_eqns)
            case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
                if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                    call TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT()
                else
                    call TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT()
                end if

            case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
                call TIME_SUBSTEP_COMPRESSIBLE()

            end select

            call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

            call DNS_BOUNDS_LIMIT()
!            if (int(logs_data(1)) /= 0) return ! Error detected

            if (part%type == PART_TYPE_BIL_CLOUD_4) then
                call PARTICLE_TIME_RESIDENCE(dtime, l_g%np, l_q)
                call PARTICLE_TIME_LIQUID_CLIPPING(s, l_q, l_txc)
            end if

            ! -------------------------------------------------------------------
            ! Update RHS hq and hs in the explicit low-storage algorithm
            ! -------------------------------------------------------------------
            if ((rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) .and. &
                rkm_substep < rkm_endstep) then

#ifdef USE_BLAS
!$omp parallel default(shared) &
!$omp private (ij_len,ij_srt,ij_end,ij_siz,alpha,is)
#else
!$omp parallel default(shared) &
!$omp private (i,   ij_srt,ij_end,ij_siz,alpha,is)
#endif

                call DNS_OMP_PARTITION(isize_field, ij_srt, ij_end, ij_siz)
#ifdef USE_BLAS
                ij_len = ij_siz
#endif

                alpha = kco(rkm_substep)

                if (icalc_flow == 1) then
                    do is = 1, inb_flow
#ifdef USE_BLAS
                        call DSCAL(ij_len, alpha, hq(ij_srt, is), 1)
#else
                        hq(ij_srt:ij_end, is) = alpha*hq(ij_srt:ij_end, is)
#endif
                    end do
                end if

                if (icalc_scal == 1) then
                    do is = 1, inb_scal
#ifdef USE_BLAS
                        call DSCAL(ij_len, alpha, hs(ij_srt, is), 1)
#else
                        hs(ij_srt:ij_end, is) = alpha*hs(ij_srt:ij_end, is)
#endif
                    end do
                end if
!$omp end parallel

                if (part%type /= PART_TYPE_NONE) then
                    do is = 1, inb_part
                        l_hq(1:l_g%np, is) = alpha*l_hq(1:l_g%np, is)
                    end do
                end if

            end if

            ! -------------------------------------------------------------------
            ! Profiling data
            ! -------------------------------------------------------------------
#ifdef USE_PROFILE
            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)
            idummy = t_end - t_srt

#ifdef USE_MPI
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
            if (ims_pro == 0) then
                write (time_string, 999) ims_npro, ims_npro_i, ims_npro_k, rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro
999             format(I5.5, ' (ims_npro_i X ims_npro_k:', I4.4, 'x', I4.4, 1x, ') RK-Substep', I1, ':', E13.5, 's')
                call TLAB_WRITE_ASCII(lfile, time_string)
            end if
#else
            t_dif = idummy
            write (time_string, 999) rkm_substep, t_dif/1.0d0/PROC_CYCLES/ims_npro
999         format('RK-Substep', I1, ':', E13.5, 's')
            call TLAB_WRITE_ASCII(lfile, time_string)
#endif

#endif

        end do

        return
    end subroutine TIME_RUNGEKUTTA

!########################################################################
!#
!# Determine the variable time step is a positive cfl is given
!# If negative cfl is prescribed, then constant dtime and this routine
!# calculates CFL and diffustion numbers for log files
!#
!# The reacting case should be reviewed in the new formulation in terms
!# of energy.
!#
!# The diffusion number is fixed in terms of the CFL, which is the input.
!# This depends on the scheme used. From Lele (1992), page 32, we have that, if
!# the sixth order tridiagonal scheme is used, then the maximum CFL number
!# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994
!# used here we have 3.36/1.989, about 1.69.
!# This holds for periodic case. A safety margin leads to the common value of 1.2.
!#
!# If second order finite different operator is used, then the maximum
!# diffusion number is 2.9/6.857, about 0.42.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
!# If the extension by Lamballais et al is used, then the maximum
!# diffusion number is 2.9/pi^2, about 0.29.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
!#
!# If twice the first order finite difference operator is used, then the
!# maximum diffusion number is 2.9/1.989^2, about 0.73.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
!#
!# In incompressible mode the arrays rho, p and vis are not used
!#
!########################################################################
    subroutine TIME_COURANT()
        use DNS_LOCAL, only: logs_data
        use TLAB_POINTERS_3D, only: u, v, w, p_wrk3d, p, rho, vis

        ! -------------------------------------------------------------------
        integer(wi) ipmax, k_glo
        real(wp) dt_loc
        real(wp) pmax(3), dtc, dtd, dtr
#ifdef USE_MPI
        real(wp) pmax_aux(3)
#endif

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i
        kdsp = ims_offset_k
#else
        idsp = 0
        kdsp = 0
#endif

        dtc = C_BIG_R   ! So that the minimum non-zero determines dt at the end
        dtd = C_BIG_R
        dtr = C_BIG_R

        ipmax = 0       ! Initialize counter of time constraints

        ! ###################################################################
        ! CFL number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
        ! -------------------------------------------------------------------
        select case (imode_eqns)
        case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
            if (g(3)%size > 1) then
                do k = 1, kmax
                    k_glo = k + kdsp
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*g(1)%jac(i + idsp, 3) &
                                               + abs(v(i, j, k))*g(2)%jac(j, 3) &
                                               + abs(w(i, j, k))*g(3)%jac(k_glo, 3)
                        end do
                    end do
                end do
            else    ! do I need this?
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*g(1)%jac(i + idsp, 3) &
                                               + abs(v(i, j, k))*g(2)%jac(j, 3)
                        end do
                    end do
                end do
            end if

            ! -------------------------------------------------------------------
            ! Compressible: Calculate global maximum of (u+c)/dx + (v+c)/dy + (w+c)/dz
            ! -------------------------------------------------------------------
        case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
            p_wrk3d = sqrt(gama0*p(:, :, :)/rho(:, :, :)) ! sound speed; positiveness of p and rho checked in routine DNS_CONTROL
            if (g(3)%size > 1) then
                do k = 1, kmax
                    k_glo = k + kdsp
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = (abs(u(i, j, k)) + p_wrk3d(i, j, k))*g(1)%jac(i + idsp, 3) &
                                               + (abs(v(i, j, k)) + p_wrk3d(i, j, k))*g(2)%jac(j, 3) &
                                               + (abs(w(i, j, k)) + p_wrk3d(i, j, k))*g(3)%jac(k_glo, 3)
                        end do
                    end do
                end do
            else
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = (abs(u(i, j, k)) + p_wrk3d(i, j, k))*g(1)%jac(i + idsp, 3) &
                                               + (abs(v(i, j, k)) + p_wrk3d(i, j, k))*g(2)%jac(j, 3)
                        end do
                    end do
                end do
            end if

        end select

        pmax(1) = maxval(p_wrk3d)

        ! ###################################################################
        ! Diffusion number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of \nu*(1/dx^2 + 1/dy^2 + 1/dz^2)
        ! -------------------------------------------------------------------
        select case (imode_eqns)
        case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
            pmax(2) = schmidtfactor*dx2i

            ! -------------------------------------------------------------------
            ! Compressible: Calculate global maximum of \mu/rho*(1/dx^2 + 1/dy^2 + 1/dz^2)
            ! -------------------------------------------------------------------
        case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
            if (itransport == EQNS_TRANS_POWERLAW) then
                if (g(3)%size > 1) then
                    do k = 1, kmax
                        k_glo = k + kdsp
                        do j = 1, jmax
                            do i = 1, imax
                       p_wrk3d(i, j, k) = (g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4) + g(3)%jac(k_glo, 4))*vis(i, j, k)/rho(i, j, k)
                            end do
                        end do
                    end do
                else
                    do k = 1, kmax
                        do j = 1, jmax
                            do i = 1, imax
                                p_wrk3d(i, j, k) = (g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4))*vis(i, j, k)/rho(i, j, k)
                            end do
                        end do
                    end do
                end if

            else ! constant dynamic viscosity
                if (g(3)%size > 1) then
                    do k = 1, kmax
                        k_glo = k + kdsp
                        do j = 1, jmax
                            do i = 1, imax
                                p_wrk3d(i, j, k) = (g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4) + g(3)%jac(k_glo, 4))/rho(i, j, k)
                            end do
                        end do
                    end do
                else
                    do k = 1, kmax
                        do j = 1, jmax
                            do i = 1, imax
                                p_wrk3d(i, j, k) = (g(1)%jac(i + idsp, 4) + g(2)%jac(j, 4))/rho(i, j, k)
                            end do
                        end do
                    end do
                end if

            end if

            pmax(2) = schmidtfactor*maxval(p_wrk3d)

        end select

        ! ###################################################################
        ! Final operations
        ! ###################################################################
#ifdef USE_MPI
        call MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

        if (pmax(1) > C_0_R) dtc = cfla/pmax(1) ! Set time step for the given CFL number
        if (pmax(2) > C_0_R) dtd = cfld/pmax(2) ! Set time step for the given diffusion number

        ! -------------------------------------------------------------------
        if (cfla > C_0_R) then
            if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                dt_loc = min(dtc, dtd)
            else
                dt_loc = dtc
            end if
            dt_loc = min(dt_loc, dtr)

            dtime = dt_loc

        end if

        ! Real CFL and diffusion numbers being used, for the logfile
        logs_data(2) = dtime*pmax(1)
        logs_data(3) = dtime*pmax(2)

        return

    end subroutine TIME_COURANT

!########################################################################
!#
!# Branching among different formulations of the RHS.
!#
!# Be careful to define here the pointers and to enter the RHS routines
!# with individual fields. Definition of pointer inside of RHS routines
!# decreased performance considerably (at least in JUGENE)
!#
!########################################################################
    subroutine TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT()
        use TLAB_ARRAYS
        use PARTICLE_ARRAYS
        use DNS_ARRAYS
        use DNS_LOCAL, only: imode_rhs
        use BOUNDARY_BUFFER
        use FI_SOURCES
        
        ! -----------------------------------------------------------------------
        integer(wi) ij_srt, ij_end, ij_siz    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
        integer ij_len
#endif

        ! #######################################################################
#ifdef USE_BLAS
        ij_len = isize_field
#endif

        select case (iadvection)
        case (EQNS_DIVERGENCE)
            call FI_SOURCES_FLOW(q, s, hq, txc(1, 1))
            call RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3()

            call FI_SOURCES_SCAL(s, hs, txc(1, 1), txc(1, 2))
            do is = 1, inb_scal
                call RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3(is)
            end do

        case (EQNS_SKEWSYMMETRIC)
            call FI_SOURCES_FLOW(q, s, hq, txc(1, 1))
            call RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2()

            call FI_SOURCES_SCAL(s, hs, txc(1, 1), txc(1, 2))
            do is = 1, inb_scal
                call RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2(is)
            end do

        case (EQNS_CONVECTIVE)
            select case (imode_rhs)
            case (EQNS_RHS_SPLIT)
                call FI_SOURCES_FLOW(q, s, hq, txc(1, 1))
                call RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1()

                call FI_SOURCES_SCAL(s, hs, txc(1, 1), txc(1, 2))
                do is = 1, inb_scal
                    call RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1(is)
                end do

            case (EQNS_RHS_COMBINED)
                call FI_SOURCES_FLOW(q, s, hq, txc(1, 1))
                call FI_SOURCES_SCAL(s, hs, txc(1, 1), txc(1, 2))
                call RHS_GLOBAL_INCOMPRESSIBLE_1()

            case (EQNS_RHS_NONBLOCKING)
#ifdef USE_PSFFT
                call RHS_GLOBAL_INCOMPRESSIBLE_NBC(q(1, 1), q(1, 2), q(1, 3), s(1, 1), &
                                                   txc(1, 1), txc(1, 2), &
                                          txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8), txc(1, 9), txc(1, 10), &
                                                   txc(1, 11), txc(1, 12), txc(1, 13), txc(1, 14), &
                                                   hq(1, 1), hq(1, 2), hq(1, 3), hs(1, 1), &
                                                   wrk1d, wrk2d, wrk3d)
#else
                call TLAB_WRITE_ASCII(efile, 'TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT. Need compiling flag -DUSE_PSFFT.')
                call TLAB_STOP(DNS_ERROR_PSFFT)
#endif
            end select
        end select

        if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
            call BOUNDARY_BUFFER_RELAX_SCAL() ! Flow part needs to be taken into account in the pressure
        end if

        ! #######################################################################
        ! Perform the time stepping for incompressible equations
        ! #######################################################################
#ifdef USE_OPENMP
#ifdef USE_BLAS
!$omp parallel default(shared) &
!$omp private (ij_len,is,ij_srt,ij_end,ij_siz)
#else
!$omp parallel default(shared) &
!$omp private (ij,  is,ij_srt,ij_end,ij_siz)
#endif
#endif

        call DNS_OMP_PARTITION(isize_field, ij_srt, ij_end, ij_siz)
#ifdef USE_BLAS
        ij_len = ij_siz
#endif
        do is = 1, inb_flow

#ifdef USE_BLAS
            call DAXPY(ij_len, dte, hq(ij_srt, is), 1, q(ij_srt, is), 1)
#else
            q(ij_srt:ij_end, is) = q(ij_srt:ij_end, is) + dte*hq(ij_srt:ij_end, is)
#endif
        end do

        do is = 1, inb_scal
#ifdef BLAS
            call DAXPY(ij_len, dte, hs(ij_srt, is), 1, s(ij_srt, is), 1)
#else
            s(ij_srt:ij_end, is) = s(ij_srt:ij_end, is) + dte*hs(ij_srt:ij_end, is)
#endif
        end do
#ifdef USE_OPENMP
!$omp end parallel
#endif

        return
    end subroutine TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT

!########################################################################
!########################################################################
    subroutine TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT()
        use TLAB_ARRAYS
        use DNS_ARRAYS

        ! ######################################################################
        if (rkm_mode == RKM_IMP3_DIFFUSION) then
            call RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(kex(rkm_substep), kim(rkm_substep), kco(rkm_substep), &
                                                      q, hq, q(:, 1), q(:, 2), q(:, 3), hq(1, 1), hq(1, 2), hq(1, 3), s, hs, &
                                   txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), wrk1d, wrk2d, wrk3d)
            ! pressure-correction algorithm; to be checked
            ! CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(&
            !      kex,kim,kco,  &
            !      q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
            !      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), txc(1,8), &
            !      wrk1d,wrk2d,wrk3d)
        else
            call TLAB_WRITE_ASCII(efile, 'TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT. Undeveloped formulation.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)

        end if

        return
    end subroutine TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT

!########################################################################
!########################################################################
    subroutine TIME_SUBSTEP_COMPRESSIBLE()
        use TLAB_ARRAYS
        use TLAB_POINTERS
        use DNS_ARRAYS
        use BOUNDARY_BUFFER
        use BOUNDARY_BCS, only: BcsDrift
        use BOUNDARY_BCS_COMPRESSIBLE

        ! -------------------------------------------------------------------
        real(wp) rho_ratio, dt_rho_ratio, prefactor
        real(wp) M2_max, dummy
        integer(wi) inb_scal_loc

        ! ###################################################################
        ! Evaluate standard RHS of equations
        ! global formulation
        ! ###################################################################
        if (imode_eqns == DNS_EQNS_INTERNAL .and. &
            iadvection == EQNS_SKEWSYMMETRIC .and. &
            iviscous == EQNS_EXPLICIT .and. &
            idiffusion == EQNS_EXPLICIT) then
            call RHS_FLOW_GLOBAL_2()

            do is = 1, inb_scal
                call RHS_SCAL_GLOBAL_2(is)
            end do

        else
            ! ###################################################################
            ! Evaluate standard RHS of equations
            ! split formulation
            ! ###################################################################
            ! -------------------------------------------------------------------
            ! convective terms
            ! -------------------------------------------------------------------
            if (iadvection == EQNS_DIVERGENCE) then
                call RHS_FLOW_EULER_DIVERGENCE(rho, u, v, w, p, e, hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), &
                                               txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
                do is = 1, inb_scal
                    call RHS_SCAL_EULER_DIVERGENCE(rho, u, v, w, s(1, is), hs(1, is), &
                                                   txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
                end do
            else if (iadvection == EQNS_SKEWSYMMETRIC) then
                call RHS_FLOW_EULER_SKEWSYMMETRIC(rho, u, v, w, p, e, s, hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), hs, &
                                                  txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
                do is = 1, inb_scal
                    call RHS_SCAL_EULER_SKEWSYMMETRIC(rho, u, v, w, s(1, is), hs(1, is), &
                                                      txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
                end do
            end if

            ! -------------------------------------------------------------------
            ! viscous terms
            ! -------------------------------------------------------------------
            if (itransport /= 1) then
                call TLAB_WRITE_ASCII(efile, 'TIME_SUBSTEP_COMPRESSIBLE. Section requires to allocate array vis.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            if (iviscous == EQNS_DIVERGENCE) then
                call RHS_FLOW_VISCOUS_DIVERGENCE(vis, u, v, w, p, hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), &
                    txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8), txc(1, 9), wrk2d, wrk3d)
            else if (iviscous == EQNS_EXPLICIT) then
                call RHS_FLOW_VISCOUS_EXPLICIT(vis, u, v, w, p, hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), &
                                               txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
            end if

            ! -------------------------------------------------------------------
            ! diffusion/conduction terms
            ! -------------------------------------------------------------------
            if (idiffusion == EQNS_DIVERGENCE) then
                ! diffusion transport of enthalpy is accumulated in txc and used then in RHS_FLOW_CONDUCTION
                txc(:, 1) = C_0_R; txc(:, 2) = C_0_R; txc(:, 3) = C_0_R
                do is = 1, inb_scal
                    call RHS_SCAL_DIFFUSION_DIVERGENCE(is, vis, s, T, hs, &
                                          txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), wrk2d, wrk3d)
                end do
                call RHS_FLOW_CONDUCTION_DIVERGENCE(vis, s, T, hq(1, 4), &
                                          txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), wrk2d, wrk3d)
            else if (idiffusion == EQNS_EXPLICIT) then
                do is = 1, inb_scal
                    call RHS_SCAL_DIFFUSION_EXPLICIT(is, vis, s, T, hs, hq(1, 4), &
                                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
                end do
         call RHS_FLOW_CONDUCTION_EXPLICIT(vis, s, T, hq(1, 4), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)
            end if

        end if

        ! ###################################################################
        ! Impose boundary conditions
        ! Temperature array T is used as auxiliary array because it is no
        ! longer used until the fields are updated
        ! ###################################################################
#define GAMMA_LOC(i) txc(i,6)
#define AUX_LOC(i)   T(i)

        call THERMO_GAMMA(imax, jmax, kmax, s, T, GAMMA_LOC(1))

        ! Maximum Mach for Poinsot & Lele reference pressure BC
        if (BcsDrift) then
            M2_max = C_0_R
            do i = 1, isize_field
                dummy = (u(i)*u(i) + v(i)*v(i) + w(i)*w(i))*rho(i)/(GAMMA_LOC(i)*p(i))
                M2_max = max(M2_max, dummy)
            end do
#ifdef USE_MPI
            call MPI_ALLREDUCE(M2_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
            M2_max = dummy
#endif
        end if

        if (.not. g(2)%periodic) then
            call BOUNDARY_BCS_Y(isize_field, M2_max, rho, u, v, w, p, GAMMA_LOC(1), s, &
                                hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), hs, &
                                txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), AUX_LOC(:), wrk2d, wrk3d)
        end if

        if (.not. g(1)%periodic) then
            call BOUNDARY_BCS_X(isize_field, M2_max, etime, rho, u, v, w, p, GAMMA_LOC(1), s, &
                                hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), hs, txc, AUX_LOC(:), wrk1d, wrk2d, wrk3d)
        end if

#undef GAMMA_LOC
#undef AUX_LOC

        ! ###################################################################
        ! Impose buffer zone as relaxation terms
        ! ###################################################################
        if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
            call BOUNDARY_BUFFER_RELAX_FLOW()
            call BOUNDARY_BUFFER_RELAX_SCAL()
        end if

        ! ###################################################################
        ! Perform the time stepping
        ! ###################################################################
        rho_ratio = C_1_R
        prefactor = (gama0 - C_1_R)*mach*mach

        if (icalc_flow == 1) then
            if (icalc_scal == 1) then; inb_scal_loc = inb_scal
            else; inb_scal_loc = 0
            end if

            ! -------------------------------------------------------------------
            ! Total energy formulation
            ! -------------------------------------------------------------------
            if (imode_eqns == DNS_EQNS_TOTAL) then
                !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
                !$omp do
                do i = 1, isize_field
                    rho_ratio = rho(i)
                    rho(i) = rho(i) + dte*hq(i, 5)
                    rho_ratio = rho_ratio/rho(i)
                    dt_rho_ratio = dte/rho(i)

                    e(i) = rho_ratio*(e(i) + prefactor*C_05_R*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) &
                           + dt_rho_ratio*hq(i, 4)

                    u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i, 1)
                    v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i, 2)
                    w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i, 3)

                    e(i) = e(i) - prefactor*C_05_R*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))

                    do is = 1, inb_scal_loc
                        s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
                    end do
                end do
                !$omp end do
                !$omp end parallel

                ! -------------------------------------------------------------------
                ! Internal energy formulation
                ! -------------------------------------------------------------------
            else if (imode_eqns == DNS_EQNS_INTERNAL) then
                !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
                !$omp do
                do i = 1, isize_field
                    rho_ratio = rho(i)
                    rho(i) = rho(i) + dte*hq(i, 5)
                    rho_ratio = rho_ratio/rho(i)
                    dt_rho_ratio = dte/rho(i)

                    e(i) = rho_ratio*e(i) + dt_rho_ratio*hq(i, 4)
                    u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i, 1)
                    v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i, 2)
                    w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i, 3)

                    do is = 1, inb_scal_loc
                        s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
                    end do
                end do
                !$omp end do
                !$omp end parallel

            end if

        else
            if (icalc_scal == 1) then
                do is = 1, inb_scal
                    !$omp parallel default( shared ) private( i, dt_rho_ratio )
                    !$omp do
                    do i = 1, isize_field
                        dt_rho_ratio = dte/rho(i)
                        s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
                    end do
                    !$omp end do
                    !$omp end parallel
                end do
            end if
        end if

        ! ###################################################################
        ! Impose buffer zone as filter
        ! ###################################################################
        if (BuffType == DNS_BUFFER_FILTER .or. BuffType == DNS_BUFFER_BOTH) then
            call BOUNDARY_BUFFER_FILTER &
                (rho, u, v, w, e, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        end if

        return
    end subroutine TIME_SUBSTEP_COMPRESSIBLE

    !########################################################################
    !########################################################################
    subroutine TIME_SUBSTEP_PARTICLE()
        use DNS_ARRAYS
        use PARTICLE_VARS
        use PARTICLE_ARRAYS

        ! -------------------------------------------------------------------
        integer(wi) is, i

#ifdef USE_MPI
        integer(wi) nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
#else
        real(wp) x_right, z_right
#endif

        !#####################################################################
        call RHS_PART_1()

        !#######################################################################
        ! Update particle properties
        !#######################################################################
        do is = 1, inb_part
            l_q(1:l_g%np, is) = l_q(1:l_g%np, is) + dte*l_hq(1:l_g%np, is)
        end do

        !#####################################################################
        ! Boundary control to see if particles leave processor
        !#####################################################################
#ifdef USE_MPI

        ! -------------------------------------------------------------------
        ! Particle sorting for Send/Recv X-Direction
        ! -------------------------------------------------------------------
        call PARTICLE_MPI_SORT(1, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

        if (ims_pro_i == 0) then !Take care of periodic boundary conditions west
            if (nzone_west /= 0) then
                l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) = &
                    l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) + g(1)%scale
            end if
        end if

        if (ims_pro_i == (ims_npro_i - 1)) then !Take care of periodic boundary conditions east
            if (nzone_east /= 0) then
                l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) = &
                    l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) - g(1)%scale
            end if
        end if

        call PARTICLE_MPI_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, &
                                      l_q, l_hq, l_g%tags, l_g%np)

        ! -------------------------------------------------------------------
        ! Particle sorting for Send/Recv Z-Direction
        ! -------------------------------------------------------------------
        call PARTICLE_MPI_SORT(3, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

        if (ims_pro_k == 0) then !Take care of periodic boundary conditions south
            if (nzone_south /= 0) then
                l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) = &
                    l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) + g(3)%scale
            end if
        end if

        if (ims_pro_k == (ims_npro_k - 1)) then !Take care of periodic boundary conditions north
            if (nzone_north /= 0) then
                l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) = &
                    l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) - g(3)%scale
            end if
        end if

        call PARTICLE_MPI_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, &
                                      l_q, l_hq, l_g%tags, l_g%np)

#else
        !#######################################################################
        ! Serial; would it be faster to use MOD?
        x_right = g(1)%nodes(1) + g(1)%scale
        z_right = g(3)%nodes(1) + g(3)%scale

        do i = 1, l_g%np
            if (l_q(i, 1) > x_right) then
                l_q(i, 1) = l_q(i, 1) - g(1)%scale

            elseif (l_q(i, 1) < g(1)%nodes(1)) then
                l_q(i, 1) = l_q(i, 1) + g(1)%scale

            end if

            if (l_q(i, 3) > z_right) then
                l_q(i, 3) = l_q(i, 3) - g(3)%scale

            elseif (l_q(i, 3) < g(3)%nodes(1)) then
                l_q(i, 3) = l_q(i, 3) + g(3)%scale

            end if
        end do

#endif

        !#######################################################################
        ! Recalculating closest node below in Y direction
        !#######################################################################
        call PARTICLE_LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, g(2)%size, g(2)%nodes)

        return
    end subroutine TIME_SUBSTEP_PARTICLE

end module TIME

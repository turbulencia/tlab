#include "dns_error.h"
#include "dns_const.h"

module BOUNDARY_BCS
    use TLab_Constants, only: wp, wi, BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH, MAX_VARS
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk3d
    use FDM_MatMul
    use FDM_Com1_Jacobian
    implicit none
    save
    private

    type bcs_dt
        sequence
        integer type(MAX_VARS)                              ! dirichlet, neumann for incompressible
        integer SfcType(MAX_VARS)                           ! Type of Surface Model
        real(wp) cpl(MAX_VARS)                              ! Coupling parameter for surface model
        real(wp) cinf, cout, ctan                           ! characteristic formulation for compressible
        real(wp), allocatable, dimension(:, :, :) :: ref    ! reference fields
    end type bcs_dt

    type(bcs_dt), public :: BcsFlowImin, BcsFlowImax, BcsFlowJmin, BcsFlowJmax, BcsFlowKmin, BcsFlowKmax
    type(bcs_dt), public :: BcsScalImin, BcsScalImax, BcsScalJmin, BcsScalJmax, BcsScalKmin, BcsScalKmax

    logical, public :: BcsDrift

    integer(wi) nmin, nmax, nsize
    public :: BOUNDARY_BCS_NEUMANN_Y

    public :: BOUNDARY_BCS_SURFACE_Y

    ! Compressible viscous
    integer, public :: bcs_inf(2, 2, 3), bcs_out(2, 2, 3) ! 1. index: lower and upper values
    !                                                       2. index: derivative order
    !                                                       3. index: direction
    public :: BOUNDARY_BCS_SCAL_READBLOCK, BOUNDARY_BCS_FLOW_READBLOCK
    public :: BOUNDARY_BCS_INITIALIZE

    ! Boundary conditions
    integer, parameter, public :: DNS_BCS_NONE = 0
    integer, parameter, public :: DNS_BCS_NR = 1
    integer, parameter, public :: DNS_BCS_INFLOW = 2
    integer, parameter, public :: DNS_BCS_DIRICHLET = 3
    integer, parameter, public :: DNS_BCS_NEUMANN = 4

    ! Surface Models
    integer, parameter, public :: DNS_SFC_STATIC = 0
    integer, parameter, public :: DNS_SFC_LINEAR = 1

contains
! ###################################################################
! ###################################################################
    subroutine BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, tag, var)
        use TLab_Memory, only: inb_scal
        character(len=*), intent(in) :: bakfile, inifile, tag
        type(bcs_dt), intent(out) :: var

        character(len=512) sRes
        character(len=20) lstr
        integer is

        do is = 1, inb_scal
            write (lstr, *) is; lstr = 'Scalar'//trim(adjustl(lstr))

            call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'none') then; var%type(is) = DNS_BCS_NONE
            else if (trim(adjustl(sRes)) == 'dirichlet') then; var%type(is) = DNS_BCS_DIRICHLET
            else if (trim(adjustl(sRes)) == 'neumann') then; var%type(is) = DNS_BCS_NEUMANN
            else
                call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr)))
                call TLab_Stop(DNS_ERROR_JBC)
            end if

            call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)), 'static', sRes)
            if (sRes == 'static') then
                var%SfcType(is) = DNS_SFC_STATIC
            elseif (sRes == 'linear') then
                var%SfcType(is) = DNS_SFC_LINEAR
            else
                call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)))
                call TLab_Stop(DNS_ERROR_JBC)
            end if

            call ScanFile_Real(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//'Coupling'//trim(adjustl(tag)), '0.0', var%cpl(is))

        end do

        return
    end subroutine BOUNDARY_BCS_SCAL_READBLOCK

! ###################################################################
! ###################################################################
    subroutine BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, tag, var)
        character(len=*), intent(in) :: bakfile, inifile, tag
        type(bcs_dt), intent(out) :: var

        character(len=512) sRes
        integer inormal, itangential(2)

        select case (trim(adjustl(tag)))
        case ('Imin', 'Imax')
            inormal = 1
            itangential = [2, 3]
        case ('Jmin', 'Jmax')
            inormal = 2
            itangential = [1, 3]
        case ('Kmin', 'Kmax')
            inormal = 3
            itangential = [1, 2]
        end select

        call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', 'Velocity'//trim(adjustl(tag)), 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; var%type(1:3) = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; var%type(1:3) = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; var%type(inormal) = DNS_BCS_DIRICHLET
            var%type(itangential) = DNS_BCS_NEUMANN
        else
            call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.Velocity'//trim(adjustl(tag)))
            call TLab_Stop(DNS_ERROR_IBC)
        end if

    end subroutine BOUNDARY_BCS_FLOW_READBLOCK

! ###################################################################
! ###################################################################
    subroutine BOUNDARY_BCS_INITIALIZE()
        use TLab_Constants, only: tag_flow, tag_scal, lfile, efile
#ifdef TRACE_ON
        use TLab_Constants, only: tfile
#endif
        use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL
        use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal, inb_flow_array, inb_scal_array
        use FDM, only: g
        use Tlab_Background, only: pbg, qbg
        use Thermodynamics, only: CRATIO_INV
        use THERMO_THERMAL
        use THERMO_CALORIC
        use BOUNDARY_BUFFER
        use Profiles, only: profiles_dt, Profiles_Calculate, PROFILE_TANH
#ifdef USE_MPI
        use mpi_f08
        use TLab_Memory, only: inb_scal_array
        use TLabMPI_VARS, only: ims_npro_k
        use TLabMPI_VARS, only: ims_bcs_imax, ims_bcs_jmax
        use TLabMPI_Transpose, only: TLabMPI_Trp_PlanK
#endif

! -------------------------------------------------------------------
        integer j, is !, ny
        real(wp) prefactor !, coef(5)
        type(profiles_dt) prof_loc
        integer, parameter :: i1 = 1

#ifdef USE_MPI
        character*32 str
        integer(wi) isize_loc
#endif

! ###################################################################
#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING BOUNDARY_BCS_INITIALLIZE')
#endif

! ###################################################################
! Allocate memory space
! ###################################################################
! we use inb_flow_array and inb_scal_array to have space for aux fields
! we add space at the end for the shape factor
        allocate (BcsFlowImin%ref(jmax, kmax, inb_flow_array + 1))
        allocate (BcsFlowImax%ref(jmax, kmax, inb_flow_array + 1))
        allocate (BcsFlowJmin%ref(imax, kmax, inb_flow_array + 1))
        allocate (BcsFlowJmax%ref(imax, kmax, inb_flow_array + 1))
        allocate (BcsFlowKmin%ref(imax, jmax, inb_flow_array + 1)) ! not yet used
        allocate (BcsFlowKmax%ref(imax, jmax, inb_flow_array + 1))

        allocate (BcsScalImin%ref(jmax, kmax, inb_scal_array + 1))
        allocate (BcsScalImax%ref(jmax, kmax, inb_scal_array + 1))
        allocate (BcsScalJmin%ref(imax, kmax, inb_scal_array + 1))
        allocate (BcsScalJmax%ref(imax, kmax, inb_scal_array + 1))
        allocate (BcsScalKmin%ref(imax, jmax, inb_scal_array + 1)) ! not yet used
        allocate (BcsScalKmax%ref(imax, jmax, inb_scal_array + 1))

! #######################################################################
! Compressible mode
! #######################################################################
        if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
#ifdef USE_MPI
! -------------------------------------------------------------------
! Characteristic BCs
! -------------------------------------------------------------------
            if (.not. g(1)%periodic) then ! Required for NRBCs in Ox
                isize_loc = mod(jmax, ims_npro_k)
                ims_bcs_imax = 2*(inb_flow + inb_scal_array)
                do while (mod(isize_loc*ims_bcs_imax, ims_npro_k) > 0)
                    ims_bcs_imax = ims_bcs_imax + 1
                end do
                write (str, *) ims_bcs_imax

                ! to be checked
                ! isize_loc = ims_bcs_imax*jmax
                ! ims_trp_plan_k(TLAB_MPI_TRP_K_NRBCX) = TLabMPI_Trp_PlanK(kmax, isize_loc, message ='Ox BCs transverse terms. '//trim(adjustl(str))//' planes.')
            end if

            if (.not. g(2)%periodic) then ! Required for NRBCs in Oy
                isize_loc = mod(imax, ims_npro_k)
                ims_bcs_jmax = 2*(inb_flow + inb_scal_array)
                do while (mod(isize_loc*ims_bcs_jmax, ims_npro_k) > 0)
                    ims_bcs_jmax = ims_bcs_jmax + 1
                end do
                write (str, *) ims_bcs_jmax

                ! to be checked
                ! isize_loc = imax*ims_bcs_jmax
                ! ims_trp_plan_k(TLAB_MPI_TRP_K_NRBCY) = TLabMPI_Trp_PlanK(kmax, isize_loc, message ='Oy BCs transverse terms. '//trim(adjustl(str))//' planes.')
            end if
#endif

! ###################################################################
! Compute reference pressure for Poinsot&Lele term in NR Bcs.
! Note that buffer_?(1,1,1,4) is now the energy, and we need p
! ###################################################################
            BcsFlowImin%ref(:, :, 5) = pbg%mean; BcsFlowImax%ref(:, :, 5) = pbg%mean ! default values
            BcsFlowJmin%ref(:, :, 5) = pbg%mean; BcsFlowJmax%ref(:, :, 5) = pbg%mean
            BcsFlowKmin%ref(:, :, 5) = pbg%mean; BcsFlowKmax%ref(:, :, 5) = pbg%mean

            prefactor = 0.5_wp*CRATIO_INV

! -------------------------------------------------------------------
! Using buffer fields; bottom
! -------------------------------------------------------------------
            if (BuffFlowJmin%size > 0) then
                BcsFlowJmin%ref(:, :, 1) = BuffFlowJmin%Ref(:, 1, :, 5)                          ! density
                BcsFlowJmin%ref(:, :, 2) = BuffFlowJmin%Ref(:, 1, :, 2)/BcsFlowJmin%ref(:, :, 1) ! normal velocity, in this case v
                BcsFlowJmin%ref(:, :, 3) = BuffFlowJmin%Ref(:, 1, :, 1)/BcsFlowJmin%ref(:, :, 1)
                BcsFlowJmin%ref(:, :, 4) = BuffFlowJmin%Ref(:, 1, :, 3)/BcsFlowJmin%ref(:, :, 1)
                BcsFlowJmin%ref(:, :, 6) = BuffFlowJmin%Ref(:, 1, :, 4)/BcsFlowJmin%ref(:, :, 1) ! energy, to get pressure into ref5 below
                if (nse_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowJmin%ref(:, :, 6) = BcsFlowJmin%ref(:, :, 6) &
                                               - prefactor*(BcsFlowJmin%ref(:, :, 2)*BcsFlowJmin%ref(:, :, 2) &
                                                            + BcsFlowJmin%ref(:, :, 3)*BcsFlowJmin%ref(:, :, 3) &
                                                            + BcsFlowJmin%ref(:, :, 4)*BcsFlowJmin%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalJmin%ref(:, :, is) = BuffScalJmin%Ref(:, 1, :, is)/BcsFlowJmin%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax*kmax, BcsScalJmin%Ref, &
                                                BcsFlowJmin%ref(1, 1, 6), BcsFlowJmin%Ref(1, 1, 1), BcsFlowJmin%Ref(1, 1, 7), wrk3d)
                call THERMO_THERMAL_PRESSURE(imax*kmax, BcsScalJmin%Ref, &
                                             BcsFlowJmin%Ref(1, 1, 1), BcsFlowJmin%Ref(1, 1, 7), BcsFlowJmin%Ref(1, 1, 5))

! shape factor
                BcsFlowJmin%ref(:, :, inb_flow + 1) = 1.0_wp
                BcsScalJmin%ref(:, :, inb_scal + 1) = 1.0_wp

            end if

! -------------------------------------------------------------------
! Using buffer fields; top
! -------------------------------------------------------------------
            if (BuffFlowJmax%size > 0) then
                j = BuffFlowJmax%size
                BcsFlowJmax%ref(:, :, 1) = BuffFlowJmax%Ref(:, j, :, 5)                          ! density
                BcsFlowJmax%ref(:, :, 2) = BuffFlowJmax%Ref(:, j, :, 2)/BcsFlowJmax%ref(:, :, 1) ! normal velocity, in this case v
                BcsFlowJmax%ref(:, :, 3) = BuffFlowJmax%Ref(:, j, :, 1)/BcsFlowJmax%ref(:, :, 1)
                BcsFlowJmax%ref(:, :, 4) = BuffFlowJmax%Ref(:, j, :, 3)/BcsFlowJmax%ref(:, :, 1)
                BcsFlowJmax%ref(:, :, 6) = BuffFlowJmax%Ref(:, j, :, 4)/BcsFlowJmax%ref(:, :, 1) ! energy, to get pressure into ref5 below
                if (nse_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowJmax%ref(:, :, 6) = BcsFlowJmax%ref(:, :, 6) &
                                               - prefactor*(BcsFlowJmax%ref(:, :, 2)*BcsFlowJmax%ref(:, :, 2) &
                                                            + BcsFlowJmax%ref(:, :, 3)*BcsFlowJmax%ref(:, :, 3) &
                                                            + BcsFlowJmax%ref(:, :, 4)*BcsFlowJmax%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalJmax%ref(:, :, is) = BuffScalJmax%Ref(:, j, :, is)/BcsFlowJmax%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax*kmax, BcsScalJmax%Ref, &
                                                BcsFlowJmax%ref(1, 1, 6), BcsFlowJmax%Ref(1, 1, 1), BcsFlowJmax%Ref(1, 1, 7), wrk3d)
                call THERMO_THERMAL_PRESSURE(imax*kmax, BcsScalJmax%Ref, &
                                             BcsFlowJmax%Ref(1, 1, 1), BcsFlowJmax%Ref(1, 1, 7), BcsFlowJmax%Ref(1, 1, 5))

! shape factor
                BcsFlowJmax%ref(:, :, inb_flow + 1) = 1.0_wp
                BcsScalJmax%ref(:, :, inb_scal + 1) = 1.0_wp
            end if

! -------------------------------------------------------------------
! Using buffer fields; left
! -------------------------------------------------------------------
            if (BuffFlowImin%size > 0) then
                BcsFlowImin%ref(:, :, 1) = BuffFlowImin%Ref(:, :, 1, 5)                          ! density
                BcsFlowImin%ref(:, :, 2) = BuffFlowImin%Ref(:, :, 1, 1)/BcsFlowImin%ref(:, :, 1) ! normal velocity, in this case u
                BcsFlowImin%ref(:, :, 3) = BuffFlowImin%Ref(:, :, 1, 2)/BcsFlowImin%ref(:, :, 1)
                BcsFlowImin%ref(:, :, 4) = BuffFlowImin%Ref(:, :, 1, 3)/BcsFlowImin%ref(:, :, 1)
                BcsFlowImin%ref(:, :, 6) = BuffFlowImin%Ref(:, :, 1, 4)/BcsFlowImin%ref(:, :, 1) ! energy, to get pressure into ref5 below
                if (nse_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowImin%ref(:, :, 6) = BcsFlowImin%ref(:, :, 6) &
                                               - prefactor*(BcsFlowImin%ref(:, :, 2)*BcsFlowImin%ref(:, :, 2) &
                                                            + BcsFlowImin%ref(:, :, 3)*BcsFlowImin%ref(:, :, 3) &
                                                            + BcsFlowImin%ref(:, :, 4)*BcsFlowImin%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalImin%ref(:, :, is) = BuffScalImin%Ref(:, :, 1, is)/BcsFlowImin%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax*kmax, BcsScalImin%Ref, &
                                                BcsFlowImin%ref(1, 1, 6), BcsFlowImin%Ref(1, 1, 1), BcsFlowImin%Ref(1, 1, 7), wrk3d)
                call THERMO_THERMAL_PRESSURE(imax*kmax, BcsScalImin%Ref, &
                                             BcsFlowImin%Ref(1, 1, 1), BcsFlowImin%Ref(1, 1, 7), BcsFlowImin%Ref(1, 1, 5))

! shape factor
                prof_loc = profiles_dt(type=PROFILE_TANH)
                prof_loc%mean = 0.5_wp
                prof_loc%ymean = qbg(1)%ymean
                prof_loc%thick = qbg(1)%diam/8.0_wp
                prof_loc%parameters(5) = 3.0_wp*qbg(1)%diam
                BcsFlowImin%ref(:, :, inb_flow + 1) = Profiles_Calculate(prof_loc, g(2)%nodes(j))
                BcsScalImin%ref(:, :, inb_scal + 1) = Profiles_Calculate(prof_loc, g(2)%nodes(j))

            end if

! -------------------------------------------------------------------
! Using buffer fields; right
! -------------------------------------------------------------------
            if (BuffFlowImax%size > 0) then
                BcsFlowImax%ref(:, :, 1) = BuffFlowImax%Ref(:, :, 1, 5)                          ! density
                BcsFlowImax%ref(:, :, 2) = BuffFlowImax%Ref(:, :, 1, 1)/BcsFlowImax%ref(:, :, 1) ! normal velocity, in this case u
                BcsFlowImax%ref(:, :, 3) = BuffFlowImax%Ref(:, :, 1, 2)/BcsFlowImax%ref(:, :, 1)
                BcsFlowImax%ref(:, :, 4) = BuffFlowImax%Ref(:, :, 1, 3)/BcsFlowImax%ref(:, :, 1)
                BcsFlowImax%ref(:, :, 6) = BuffFlowImax%Ref(:, :, 1, 4)/BcsFlowImax%ref(:, :, 1) ! energy, to get pressure into ref5 below
                if (nse_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowImax%ref(:, :, 6) = BcsFlowImax%ref(:, :, 6) &
                                               - prefactor*(BcsFlowImax%ref(:, :, 2)*BcsFlowImax%ref(:, :, 2) &
                                                            + BcsFlowImax%ref(:, :, 3)*BcsFlowImax%ref(:, :, 3) &
                                                            + BcsFlowImax%ref(:, :, 4)*BcsFlowImax%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalImax%ref(:, :, is) = BuffScalImax%Ref(:, :, 1, is)/BcsFlowImax%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax*kmax, BcsScalImax%Ref, &
                                                BcsFlowImax%ref(1, 1, 6), BcsFlowImax%Ref(1, 1, 1), BcsFlowImax%Ref(1, 1, 7), wrk3d)
                call THERMO_THERMAL_PRESSURE(imax*kmax, BcsScalImax%Ref, &
                                             BcsFlowImax%Ref(1, 1, 1), BcsFlowImax%Ref(1, 1, 7), BcsFlowImax%Ref(1, 1, 5))

! try to use only the coflow values
                BcsFlowImax%ref(:, :, 3) = 0.0_wp

! shape factor
                BcsFlowImax%ref(:, :, inb_flow + 1) = 1.0_wp
                BcsScalImax%ref(:, :, inb_scal + 1) = 1.0_wp

            end if

        end if

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'LEAVING BOUNDARY_BCS_INITIALLIZE')
#endif

        return

    end subroutine BOUNDARY_BCS_INITIALIZE

!########################################################################
!# Calculate the boundary values of a field s.t. the normal derivative is zero
!# Routine format extracted from OPR_Partial_Y
!########################################################################
    subroutine BOUNDARY_BCS_NEUMANN_Y(ibc, nx, ny, nz, g, u, bcs_hb, bcs_ht, tmp1)
        use FDM, only: fdm_dt

        integer(wi), intent(in) :: ibc     ! BCs at jmin/jmax: 1, for Neumann/-
        !                                                   2, for -      /Neumann
        !                                                   3, for Neumann/Neumann
        integer(wi) nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*nz, ny)         ! they are transposed below
        real(wp), intent(inout) :: tmp1(nx*nz, ny)      ! they are transposed below
        real(wp), intent(out) :: bcs_hb(nx*nz), bcs_ht(nx*nz)

        target u, bcs_hb, bcs_ht, tmp1

        ! -------------------------------------------------------------------
        integer(wi) nxz, nxy, ip, idl, ic

        real(wp), pointer :: p_org(:, :), p_dst(:, :)
        real(wp), pointer :: p_bcs_hb(:), p_bcs_ht(:)

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            bcs_hb = 0.0_wp; bcs_ht = 0.0_wp

        else
            ! ###################################################################
            nxy = nx*ny
            nxz = nx*nz

            ! -------------------------------------------------------------------
            ! Make y  direction the last one
            ! -------------------------------------------------------------------
            if (nz > 1) then
#ifdef USE_ESSL
                call DGETMO(u, nxy, nxy, nz, tmp1, nz)
#else
                call TLab_Transpose(u, nxy, nz, nxy, tmp1, nz)
#endif
                p_org => tmp1
                p_dst(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
                p_bcs_hb => tmp1(:, 1)
                p_bcs_ht => tmp1(:, 2)
            else
                p_org => u
                p_dst => tmp1
                p_bcs_hb => bcs_hb
                p_bcs_ht => bcs_ht
            end if

            ! ###################################################################
            ip = ibc*5

            nmin = 1; nmax = g%size
            if (any([BCS_ND, BCS_NN] == ibc)) then
                p_dst(:, 1) = 0.0_wp      ! homogeneous bcs
                nmin = nmin + 1
            end if
            if (any([BCS_DN, BCS_NN] == ibc)) then
                p_dst(:, ny) = 0.0_wp
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            ! -------------------------------------------------------------------
            ! Calculate RHS in system of equations A u' = B u
            call g%der1%matmul(g%der1%rhs, p_org, p_dst, ibc, g%der1%rhs_b, g%der1%rhs_t, p_bcs_hb, p_bcs_ht)

            ! -------------------------------------------------------------------
            ! Solve for u' in system of equations A u' = B u
            select case (g%der1%nb_diag(1))
            case (3)
                call TRIDSS(nsize, nxz, g%der1%lu(nmin:nmax, ip + 1), g%der1%lu(nmin:nmax, ip + 2), g%der1%lu(nmin:nmax, ip + 3), p_dst(:, nmin:nmax))
            case (5)
                call PENTADSS2(nsize, nxz, g%der1%lu(nmin:nmax, ip + 1), g%der1%lu(nmin:nmax, ip + 2), g%der1%lu(nmin:nmax, ip + 3), g%der1%lu(nmin:nmax, ip + 4), g%der1%lu(nmin:nmax, ip + 5), p_dst(:, nmin:nmax))
            end select

            idl = g%der1%nb_diag(1)/2 + 1
            if (any([BCS_ND, BCS_NN] == ibc)) then
                do ic = 1, idl - 1
                    p_bcs_hb(:) = p_bcs_hb(:) + g%der1%lu(1, ip + idl + ic)*p_dst(:, 1 + ic)
                end do
            end if
            if (any([BCS_DN, BCS_NN] == ibc)) then
                do ic = 1, idl - 1
                    p_bcs_ht(:) = p_bcs_ht(:) + g%der1%lu(ny, ip + idl - ic)*p_dst(:, ny - ic)
                end do
            end if

            ! ###################################################################
            ! -------------------------------------------------------------------
            ! Put bcs arrays in correct order
            ! -------------------------------------------------------------------
            if (nz > 1) then
#ifdef USE_ESSL
                if (any([BCS_ND, BCS_NN] == ibc)) call DGETMO(p_bcs_hb, nz, nz, nx, bcs_hb, nx)
                if (any([BCS_DN, BCS_NN] == ibc)) call DGETMO(p_bcs_ht, nz, nz, nx, bcs_ht, nx)
#else
                if (any([BCS_ND, BCS_NN] == ibc)) call TLab_Transpose(p_bcs_hb, nz, nx, nz, bcs_hb, nx)
                if (any([BCS_DN, BCS_NN] == ibc)) call TLab_Transpose(p_bcs_ht, nz, nx, nz, bcs_ht, nx)
#endif
            end if

        end if

        return
    end subroutine BOUNDARY_BCS_NEUMANN_Y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates and updates interactive surface boundary condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine BOUNDARY_BCS_SURFACE_Y(is, bcs, s, hs, tmp1, aux)
#ifdef TRACE_ON
        use TLab_Constants, only: tfile
        use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
        use TLab_Constants, only: lfile
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: isize_field
        use NavierStokes, only: visc, schmidt
        use FDM, only: g
        use Averages, only: AVG1V2D
        use OPR_Partial

        integer(wi) is
        integer(wi), dimension(2, 2), intent(IN) :: bcs          ! Boundary conditions from derivative operator
        real(wp), dimension(isize_field, *) :: s, hs
        real(wp), dimension(isize_field) :: tmp1
        real(wp), dimension(imax, kmax, 6), target :: aux

        integer(wi) nxy, ip, k
        real(wp), dimension(:, :), pointer :: hfx, hfx_anom
        real(wp) :: diff, hfx_avg

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING SUBROUTINE BOUNDARY_BCS_SURFACE_Y')
#endif
        diff = visc/schmidt(is)
        nxy = imax*jmax

        ! vertical derivative of scalar for flux at the boundaries
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp1)

        ! ------------------------------------------------------------
        ! Bottom Boundary
        ! ------------------------------------------------------------
        if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) then
            hfx => aux(:, :, 1)
            hfx_anom => aux(:, :, 2)
            ip = 1
            do k = 1, kmax    ! Calculate the surface flux
                hfx(:, k) = diff*tmp1(ip:ip + imax - 1); ip = ip + nxy
            end do
            hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
            hfx_anom = hfx - hfx_avg
            BcsScalJmin%ref(:, :, is) = BcsScalJmin%ref(:, :, is) + BcsScalJmin%cpl(is)*hfx_anom
        end if

        ! ------------------------------------------------------------
        ! Top Boundary
        ! ------------------------------------------------------------
        if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) then
            hfx => aux(:, :, 3)
            hfx_anom => aux(:, :, 4)
            ip = imax*(jmax - 1) + 1
            do k = 1, kmax; ! Calculate the surface flux
                hfx(:, k) = -diff*tmp1(ip:ip + imax - 1); ip = ip + nxy; 
            end do
            hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
            hfx_anom = hfx - hfx_avg
            BcsScalJmax%ref(:, :, is) = BcsScalJmax%ref(:, :, is) + BcsScalJmax%cpl(is)*hfx_anom
        end if

#ifdef TRACE_ON
        call TLab_Write_ASCII(TFILE, 'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

        return

    end subroutine BOUNDARY_BCS_SURFACE_Y

end module BOUNDARY_BCS

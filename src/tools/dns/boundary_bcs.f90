#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

module BOUNDARY_BCS
    use TLAB_CONSTANTS
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk3d
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

    real(wp), allocatable :: luND(:, :), luDN(:, :), luNN(:, :)
    public :: BOUNDARY_BCS_NEUMANN_Y

    public :: BOUNDARY_BCS_SURFACE_Y

    ! Compressible viscous
    integer, public :: bcs_inf(2, 2, 3), bcs_out(2, 2, 3) ! 1. index: lower and upper values
    !                                                       2. index: derivative order
    !                                                       3. index: direction
    public :: BOUNDARY_BCS_SCAL_READBLOCK, BOUNDARY_BCS_FLOW_READBLOCK
    public :: BOUNDARY_BCS_INITIALIZE

contains
! ###################################################################
! ###################################################################
    subroutine BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, tag, var)
        use TLAB_VARS, only: inb_scal
        character(len=*), intent(in) :: bakfile, inifile, tag
        type(bcs_dt), intent(out) :: var

        character(len=512) sRes
        character(len=20) lstr
        integer is

        do is = 1, inb_scal
            write (lstr, *) is; lstr = 'Scalar'//trim(adjustl(lstr))

            call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'none') then; var%type(is) = DNS_BCS_NONE
            else if (trim(adjustl(sRes)) == 'dirichlet') then; var%type(is) = DNS_BCS_DIRICHLET
            else if (trim(adjustl(sRes)) == 'neumann') then; var%type(is) = DNS_BCS_NEUMANN
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr)))
                call TLAB_STOP(DNS_ERROR_JBC)
            end if

        call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)), 'static', sRes)
            if (sRes == 'static') then
                var%SfcType(is) = DNS_SFC_STATIC
            elseif (sRes == 'linear') then
                var%SfcType(is) = DNS_SFC_LINEAR
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)))
                call TLAB_STOP(DNS_ERROR_JBC)
            end if

   call SCANINIREAL(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr))//'Coupling'//trim(adjustl(tag)), '0.0', var%cpl(is))

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

        call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'Velocity'//trim(adjustl(tag)), 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; var%type(1:3) = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; var%type(1:3) = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; var%type(inormal) = DNS_BCS_DIRICHLET
            var%type(itangential) = DNS_BCS_NEUMANN
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. BoundaryConditions.Velocity'//trim(adjustl(tag)))
            call TLAB_STOP(DNS_ERROR_IBC)
        end if

    end subroutine BOUNDARY_BCS_FLOW_READBLOCK

! ###################################################################
! ###################################################################
    subroutine BOUNDARY_BCS_INITIALIZE()
        use TLAB_TYPES, only: profiles_dt
        use TLAB_CONSTANTS, only: tag_flow, tag_scal, lfile, efile
#ifdef TRACE_ON
        use TLAB_CONSTANTS, only: tfile
#endif
        use TLAB_VARS, only: imode_eqns
        use TLAB_VARS, only: imax, jmax, kmax, inb_flow, inb_scal, inb_flow_array, inb_scal_array
        use TLAB_VARS, only: g
        use TLAB_VARS, only: mach, pbg, qbg
        use THERMO_VARS, only: gama0
        use THERMO_THERMAL
        use BOUNDARY_BUFFER
        use PROFILES
#ifdef USE_MPI
        use MPI
        use TLAB_VARS, only: inb_scal_array
        use TLAB_MPI_VARS, only: ims_npro_k
        use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
        use TLAB_MPI_VARS, only: ims_bcs_imax, ims_bcs_jmax
        use TLAB_MPI_PROCS
#endif

! -------------------------------------------------------------------
        integer j, is, ny
        real(wp) prefactor
        type(profiles_dt) prof_loc
        integer, parameter :: i1 = 1

#ifdef USE_MPI
        character*32 str
        integer(wi) isize_loc
        integer id
#endif

! ###################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_INITIALLIZE')
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
! Incompressible mode
! #######################################################################
        if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then

            ! LU decomposition for the 3 possible neumann boudary conditions
            ny = g(2)%size

            allocate(luND(ny, 5))
            allocate(luDN(ny, 5))
            allocate(luNN(ny, 5))

            select case (g(2)%mode_fdm)
            case (FDM_COM6_JACPENTA)
                call FDM_C1N6M_BCS_LHS(ny, BCS_ND, g(2)%jac, luND(:, 1), luND(:, 2), luND(:, 3), luND(:, 4), luND(:, 5))
                call PENTADFS2(ny - 1, luND(2:, 1), luND(2:, 2), luND(2:, 3), luND(2:, 4), luND(2:, 5))

                call FDM_C1N6M_BCS_LHS(ny, BCS_DN, g(2)%jac, luDN(:, 1), luDN(:, 2), luDN(:, 3), luDN(:, 4), luDN(:, 5))
                call PENTADFS2(ny - 1, luDN(:ny - 1, 1), luDN(:ny - 1, 2), luDN(:ny - 1, 3), luDN(:ny - 1, 4), luDN(:ny - 1, 5))

                call FDM_C1N6M_BCS_LHS(ny, BCS_NN, g(2)%jac, luNN(:, 1), luNN(:, 2), luNN(:, 3), luNN(:, 4), luNN(:, 5))
               call PENTADFS2(ny - 2, luNN(2:ny - 1, 1), luNN(2:ny - 1, 2), luNN(2:ny - 1, 3), luNN(2:ny - 1, 4), luNN(2:ny - 1, 5))

            case default ! some cases still need to be developed
                call FDM_C1N6_BCS_LHS(ny, BCS_ND, g(2)%jac, luND(:, 1), luND(:, 2), luND(:, 3))
                call TRIDFS(ny - 1, luND(2:, 1), luND(2:, 2), luND(2:, 3))

                call FDM_C1N6_BCS_LHS(ny, BCS_DN, g(2)%jac, luDN(:, 1), luDN(:, 2), luDN(:, 3))
                call TRIDFS(ny - 1, luDN(:ny - 1, 1), luDN(:ny - 1, 2), luDN(:ny - 1, 3))

                call FDM_C1N6_BCS_LHS(ny, BCS_NN, g(2)%jac, luNN(:, 1), luNN(:, 2), luNN(:, 3))
                call TRIDFS(ny - 2, luNN(2:ny - 1, 1), luNN(2:ny - 1, 2), luNN(2:ny - 1, 3))

            end select

        else
! #######################################################################
! Compressible mode
! #######################################################################
#ifdef USE_MPI
! -------------------------------------------------------------------
! Characteristic BCs
! -------------------------------------------------------------------
            if (.not. g(1)%periodic) then ! Required for NRBCs in Ox
                id = TLAB_MPI_K_NRBCX
                isize_loc = mod(jmax, ims_npro_k)
                ims_bcs_imax = 2*(inb_flow + inb_scal_array)
                do while (mod(isize_loc*ims_bcs_imax, ims_npro_k) > 0)
                    ims_bcs_imax = ims_bcs_imax + 1
                end do
                write (str, *) ims_bcs_imax
                str = 'Initialize MPI types for Ox BCs transverse terms. '//trim(adjustl(str))//' planes.'
                call TLAB_WRITE_ASCII(lfile, str)
                isize_loc = ims_bcs_imax*jmax
                call TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, 1, 1, 1, 1, &
                                     ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            end if

            if (.not. g(2)%periodic) then ! Required for NRBCs in Oy
                id = TLAB_MPI_K_NRBCY
                isize_loc = mod(imax, ims_npro_k)
                ims_bcs_jmax = 2*(inb_flow + inb_scal_array)
                do while (mod(isize_loc*ims_bcs_jmax, ims_npro_k) > 0)
                    ims_bcs_jmax = ims_bcs_jmax + 1
                end do
                write (str, *) ims_bcs_jmax
                str = 'Initialize MPI types for Oy BCs transverse terms. '//trim(adjustl(str))//' planes.'
                call TLAB_WRITE_ASCII(lfile, str)
                isize_loc = imax*ims_bcs_jmax
                call TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, 1, 1, 1, 1, &
                                     ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            end if
#endif

! ###################################################################
! Compute reference pressure for Poinsot&Lele term in NR Bcs.
! Note that buffer_?(1,1,1,4) is now the energy, and we need p
! ###################################################################
            BcsFlowImin%ref(:, :, 5) = pbg%mean; BcsFlowImax%ref(:, :, 5) = pbg%mean ! default values
            BcsFlowJmin%ref(:, :, 5) = pbg%mean; BcsFlowJmax%ref(:, :, 5) = pbg%mean
            BcsFlowKmin%ref(:, :, 5) = pbg%mean; BcsFlowKmax%ref(:, :, 5) = pbg%mean

            prefactor = 0.5_wp*(gama0 - 1.0_wp)*mach*mach

! -------------------------------------------------------------------
! Using buffer fields; bottom
! -------------------------------------------------------------------
            if (BuffFlowJmin%size > 0) then
                BcsFlowJmin%ref(:, :, 1) = BuffFlowJmin%Ref(:, 1, :, 5)                          ! density
                BcsFlowJmin%ref(:, :, 2) = BuffFlowJmin%Ref(:, 1, :, 2)/BcsFlowJmin%ref(:, :, 1) ! normal velocity, in this case v
                BcsFlowJmin%ref(:, :, 3) = BuffFlowJmin%Ref(:, 1, :, 1)/BcsFlowJmin%ref(:, :, 1)
                BcsFlowJmin%ref(:, :, 4) = BuffFlowJmin%Ref(:, 1, :, 3)/BcsFlowJmin%ref(:, :, 1)
                BcsFlowJmin%ref(:, :, 6) = BuffFlowJmin%Ref(:, 1, :, 4)/BcsFlowJmin%ref(:, :, 1) ! energy, to get pressure into ref5 below
                if (imode_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowJmin%ref(:, :, 6) = BcsFlowJmin%ref(:, :, 6) &
                                               - prefactor*(BcsFlowJmin%ref(:, :, 2)*BcsFlowJmin%ref(:, :, 2) &
                                                            + BcsFlowJmin%ref(:, :, 3)*BcsFlowJmin%ref(:, :, 3) &
                                                            + BcsFlowJmin%ref(:, :, 4)*BcsFlowJmin%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalJmin%ref(:, :, is) = BuffScalJmin%Ref(:, 1, :, is)/BcsFlowJmin%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax, i1, kmax, BcsScalJmin%Ref, &
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
                if (imode_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowJmax%ref(:, :, 6) = BcsFlowJmax%ref(:, :, 6) &
                                               - prefactor*(BcsFlowJmax%ref(:, :, 2)*BcsFlowJmax%ref(:, :, 2) &
                                                            + BcsFlowJmax%ref(:, :, 3)*BcsFlowJmax%ref(:, :, 3) &
                                                            + BcsFlowJmax%ref(:, :, 4)*BcsFlowJmax%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalJmax%ref(:, :, is) = BuffScalJmax%Ref(:, j, :, is)/BcsFlowJmax%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax, i1, kmax, BcsScalJmax%Ref, &
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
                if (imode_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowImin%ref(:, :, 6) = BcsFlowImin%ref(:, :, 6) &
                                               - prefactor*(BcsFlowImin%ref(:, :, 2)*BcsFlowImin%ref(:, :, 2) &
                                                            + BcsFlowImin%ref(:, :, 3)*BcsFlowImin%ref(:, :, 3) &
                                                            + BcsFlowImin%ref(:, :, 4)*BcsFlowImin%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalImin%ref(:, :, is) = BuffScalImin%Ref(:, :, 1, is)/BcsFlowImin%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax, i1, kmax, BcsScalImin%Ref, &
                                                BcsFlowImin%ref(1, 1, 6), BcsFlowImin%Ref(1, 1, 1), BcsFlowImin%Ref(1, 1, 7), wrk3d)
                call THERMO_THERMAL_PRESSURE(imax*kmax, BcsScalImin%Ref, &
                                             BcsFlowImin%Ref(1, 1, 1), BcsFlowImin%Ref(1, 1, 7), BcsFlowImin%Ref(1, 1, 5))

! shape factor
                prof_loc%type = PROFILE_TANH
                prof_loc%thick = qbg(1)%diam/8.0_wp
                prof_loc%delta = 1.0_wp
                prof_loc%ymean = qbg(1)%ymean
                prof_loc%mean = 0.5_wp
                prof_loc%lslope = 0.0_wp
                prof_loc%uslope = 0.0_wp
                prof_loc%parameters = 0.0_wp
                prof_loc%parameters(5) = 3.0_wp*qbg(1)%diam
                BcsFlowImin%ref(:, :, inb_flow + 1) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j))
                BcsScalImin%ref(:, :, inb_scal + 1) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j))

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
                if (imode_eqns == DNS_EQNS_TOTAL) then
                    BcsFlowImax%ref(:, :, 6) = BcsFlowImax%ref(:, :, 6) &
                                               - prefactor*(BcsFlowImax%ref(:, :, 2)*BcsFlowImax%ref(:, :, 2) &
                                                            + BcsFlowImax%ref(:, :, 3)*BcsFlowImax%ref(:, :, 3) &
                                                            + BcsFlowImax%ref(:, :, 4)*BcsFlowImax%ref(:, :, 4))
                end if
                do is = 1, inb_scal
                    BcsScalImax%ref(:, :, is) = BuffScalImax%Ref(:, :, 1, is)/BcsFlowImax%ref(:, :, 1)
                end do
                call THERMO_CALORIC_TEMPERATURE(imax, i1, kmax, BcsScalImax%Ref, &
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
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_INITIALLIZE')
#endif

        return

    end subroutine BOUNDARY_BCS_INITIALIZE

!########################################################################
!# Calculate the boundary values of a field s.t. the normal derivative is zero
!# Routine format extracted from OPR_PARTIAL_Y
!########################################################################
    subroutine BOUNDARY_BCS_NEUMANN_Y(ibc, nx, ny, nz, g, u, bcs_hb, bcs_ht, tmp1)
        use TLAB_TYPES, only: grid_dt

        integer(wi), intent(in) :: ibc     ! BCs at jmin/jmax: 1, for Neumann/-
        !                                                   2, for -      /Neumann
        !                                                   3, for Neumann/Neumann
        integer(wi) nx, ny, nz
        type(grid_dt), intent(in)    :: g
        real(wp),      intent(in)    :: u(nx*nz, ny)         ! they are transposed below
        real(wp),      intent(inout) :: tmp1(nx*nz, ny)      ! they are transposed below
        real(wp),      intent(out)   :: bcs_hb(nx*nz), bcs_ht(nx*nz)

        target u, bcs_hb, bcs_ht, tmp1

        ! -------------------------------------------------------------------
        integer(wi) nxz, nxy

        real(wp), pointer :: p_org(:,:), p_dst(:,:)
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
                call DNS_TRANSPOSE(u, nxy, nz, nxy, tmp1, nz)
#endif
                p_org => tmp1
                p_dst(1:nx*nz,1:ny) => wrk3d(1:nx*ny*nz)
                p_bcs_hb => tmp1(:, 1)
                p_bcs_ht => tmp1(:, 2)
            else
                p_org => u
                p_dst => tmp1
                p_bcs_hb => bcs_hb
                p_bcs_ht => bcs_ht
            end if

            ! ###################################################################
            select case (g%mode_fdm)
            case (FDM_COM6_JACPENTA)
                call FDM_C1N6M_BCS_RHS(ny, nxz, ibc, p_org, p_dst)

                select case (ibc)
                case (BCS_ND)
                    call PENTADSS2(ny - 1, nxz, luND(2:, 1), luND(2:, 2), luND(2:, 3), luND(2:, 4), luND(2:, 5), p_dst(:, 2))
                    p_bcs_hb(:) = p_dst(:, 1) + luND(1,4)*p_dst(:, 2)

                case (BCS_DN)
                    call PENTADSS2(ny - 1, nxz, luDN(:ny - 1, 1), luDN(:ny - 1, 2), luDN(:ny - 1, 3), luDN(:ny - 1, 4), luDN(:ny - 1, 5), p_dst)
                    p_bcs_ht(:) = p_dst(:, ny) + luDN(ny, 2)*p_dst(:, ny - 1)

                case (BCS_NN)
                    call PENTADSS2(ny - 2, nxz,  luNN(2:ny-1,1), luNN(2:ny-1,2), luNN(2:ny-1,3), luNN(2:ny-1,4), luNN(2:ny-1,5), p_dst(:, 2))
                    p_bcs_hb(:) = p_dst(:, 1) + luNN(1, 4)*p_dst(:, 2)
                    p_bcs_ht(:) = p_dst(:, ny) + luNN(ny, 2)*p_dst(:, ny - 1)

                end select

            case default
                call FDM_C1N6_BCS_RHS(ny, nxz, ibc, p_org, p_dst)

                select case (ibc)
                case (BCS_ND)
                    call TRIDSS(ny - 1, nxz, luND(2:, 1), luND(2:, 2), luND(2:, 3), p_dst(:, 2))
                    p_bcs_hb(:) = p_dst(:, 1) + luND(1, 3)*p_dst(:, 2)

                case (BCS_DN)
                    call TRIDSS(ny - 1, nxz, luDN(:ny - 1, 1), luDN(:ny - 1, 2), luDN(:ny - 1, 3), p_dst)
                    p_bcs_ht(:) = p_dst(:, ny) + luDN(ny, 1)*p_dst(:, ny - 1)

                case (BCS_NN)
                    call TRIDSS(ny - 2, nxz, luNN(2:ny - 1, 1), luNN(2:ny - 1, 2), luNN(2:ny - 1, 3), p_dst(:, 2))
                    p_bcs_hb(:) = p_dst(:, 1) + luNN(1, 3)*p_dst(:, 2)
                    p_bcs_ht(:) = p_dst(:, ny) + luNN(ny, 1)*p_dst(:, ny - 1)

                end select

            end select

            ! ###################################################################
            ! -------------------------------------------------------------------
            ! Put bcs arrays in correct order
            ! -------------------------------------------------------------------
            if (nz > 1) then
#ifdef USE_ESSL
                if (ibc == 1 .or. ibc == 3) call DGETMO(p_bcs_hb, nz, nz, nx, bcs_hb, nx)
                if (ibc == 2 .or. ibc == 3) call DGETMO(p_bcs_ht, nz, nz, nx, bcs_ht, nx)
#else
                if (ibc == 1 .or. ibc == 3) call DNS_TRANSPOSE(p_bcs_hb, nz, nx, nz, bcs_hb, nx)
                if (ibc == 2 .or. ibc == 3) call DNS_TRANSPOSE(p_bcs_ht, nz, nx, nz, bcs_ht, nx)
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
        use TLAB_CONSTANTS, only: tfile
        use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
        use TLAB_CONSTANTS, only: lfile
        use TLAB_VARS, only: imax, jmax, kmax, g
        use TLAB_VARS, only: isize_field
        use TLAB_VARS, only: visc, schmidt
        use AVGS, only: AVG1V2D
        use OPR_PARTIAL

        integer(wi) is
        integer(wi), dimension(2, 2), intent(IN) :: bcs          ! Boundary conditions from derivative operator
        real(wp), dimension(isize_field, *) :: s, hs
        real(wp), dimension(isize_field) :: tmp1
        real(wp), dimension(imax, kmax, 6), target :: aux

        integer(wi) nxy, ip, k
        real(wp), dimension(:, :), pointer :: hfx, hfx_anom
        real(wp) :: diff, hfx_avg

#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING SUBROUTINE BOUNDARY_BCS_SURFACE_Y')
#endif
        diff = visc/schmidt(is)
        nxy = imax*jmax

        ! vertical derivative of scalar for flux at the boundaries
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp1)

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
        call TLAB_WRITE_ASCII(TFILE, 'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

        return

    end subroutine BOUNDARY_BCS_SURFACE_Y

end module BOUNDARY_BCS

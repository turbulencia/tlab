#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

subroutine PARTICLE_RANDOM_POSITION(l_q, l_txc, txc)

    use TLAB_TYPES, only: pointers_dt, pointers3d_dt, wp, wi, longi
    use TLAB_CONSTANTS
    use TLAB_VARS
    use PARTICLE_TYPES, only: particle_dt
    use PARTICLE_VARS
    use PARTICLE_ARRAYS, only: l_g ! but this is also varying, like l_q...
    use PARTICLE_INTERPOLATE
    use THERMO_VARS, only: imixture
    use THERMO_AIRWATER
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
    use PARTICLE_ARRAYS, only: ims_np_all
#endif
    ! use PARTICLE_TINIA
use IO_FIELDS
    implicit none

    real(wp), target :: l_q(isize_part, inb_part_array)
    real(wp), target :: l_txc(isize_part, 2)
    real(wp), target :: txc(imax, jmax, kmax, inb_scal)

! -------------------------------------------------------------------
    integer(wi) i, is
    integer(wi), allocatable :: x_seed(:)
    integer(wi) size_seed

    real(wp) xref, zref, xscale, zscale, dy_loc, dummy, dy_frac
    real(wp) y_limits(2)
    integer(wi) j_limits(2)

    real(wp) rnd_number(3), rnd_number_second
    integer(wi) rnd_scal(3)

    integer(wi) nvar
    type(pointers3d_dt), dimension(2) :: data
    type(pointers_dt), dimension(2) :: data_out

    integer(longi) count

!########################################################################
#ifdef USE_MPI
    l_g%np = int(isize_part_total/int(ims_npro, KIND=8))
    if (ims_pro < int(mod(isize_part_total, int(ims_npro, KIND=8)))) then
        l_g%np = l_g%np + 1
    end if
    call MPI_ALLGATHER(l_g%np, 1, MPI_INTEGER4, ims_np_all, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

#else
    l_g%np = int(isize_part_total)

#endif

! Create tags
    count = 0
#ifdef USE_MPI
    do i = 1, ims_pro
        count = count + int(ims_np_all(i), KIND=8)
    end do
#endif
    do i = 1, l_g%np
        l_g%tags(i) = int(i, KIND=8) + count
    end do

! Generate seed - different seed for each processor
    call random_seed(SIZE=size_seed)
    allocate (x_seed(size_seed))
#ifdef USE_MPI
    x_seed = [(i, i=1 + ims_pro, size_seed + ims_pro)]
#else
    x_seed = [(i, i=1, size_seed)]
#endif
    call random_seed(PUT=x_seed)

#ifdef USE_MPI
    xref = g(1)%nodes(ims_offset_i + 1)
    xscale = g(1)%scale/real(ims_npro_i, wp)
    zref = g(3)%nodes(ims_offset_k + 1)
    zscale = g(3)%scale/real(ims_npro_k, wp)
#else
    xref = g(1)%nodes(1)
    xscale = g(1)%scale
    zref = g(3)%nodes(1)
    zscale = g(3)%scale
#endif
    if (g(3)%size == 1) zscale = 0.0_wp ! 2D case

!########################################################################
! Particle position
!########################################################################
    select case (IniP%type)

    case default

        call random_number(l_q(1:l_g%np, 1))
        l_q(1:l_g%np, 1) = xref + l_q(1:l_g%np, 1)*xscale

        call random_number(l_q(1:l_g%np, 3))
        l_q(1:l_g%np, 3) = zref + l_q(1:l_g%np, 3)*zscale

        call random_number(l_q(1:l_g%np, 2))
        l_q(1:l_g%np, 2) = IniP%ymean + (l_q(1:l_g%np, 2) - 0.5_wp)*IniP%diam

    case (PART_INITYPE_HARDCODED)       ! For testing
        l_q(1:l_g%np, 1) = 0.0_wp
        l_q(1:l_g%np, 2) = IniP%ymean
        l_q(1:l_g%np, 3) = 0.0_wp

    case (PART_INITYPE_SCALAR)          ! Use the scalar field to create the particle distribution
        call IO_READ_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, txc)
        is = 1 ! Reference scalar

        y_limits(1) = IniP%ymean - 0.5_wp*IniP%diam
        y_limits(2) = IniP%ymean + 0.5_wp*IniP%diam
        call LOCATE_Y(2, y_limits, j_limits, g(2)%size, g(2)%nodes)
        dy_loc = g(2)%nodes(j_limits(2)) - g(2)%nodes(j_limits(1))

        i = 1
        do while (i <= l_g%np)
            call random_number(rnd_number(1:3))

            rnd_scal(1) = 1 + floor(rnd_number(1)*imax)
            rnd_scal(3) = 1 + floor(rnd_number(3)*kmax)
            rnd_scal(2) = j_limits(1) + floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))
            dy_frac = rnd_number(2)*(j_limits(2) - j_limits(1) + 1) &
                      - floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))

            dummy = (txc(rnd_scal(1), rnd_scal(2), rnd_scal(3), is) - sbg(is)%mean)/sbg(is)%delta
            dummy = abs(dummy + 0.5_wp)

            call random_number(rnd_number_second)

            if (rnd_number_second <= dummy) then

                l_q(i, 1) = xref + rnd_number(1)*xscale
                l_q(i, 3) = zref + rnd_number(3)*zscale
                l_q(i, 2) = g(2)%nodes(rnd_scal(2)) + dy_frac*dy_loc

                i = i + 1

            end if

        end do

    end select

! Calculating closest node below in Y direction
    call LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, g(2)%size, g(2)%nodes)

!########################################################################
! Remaining particle properties
!########################################################################
    select case (part%type)

    case (PART_TYPE_INERTIA)                            ! velocity
        l_q(:, 4:6) = 0.0_wp

    case (PART_TYPE_BIL_CLOUD_3, PART_TYPE_BIL_CLOUD_4) ! scalar fields
        call IO_READ_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, txc)

        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            nvar = 0
            nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 1); data_out(nvar)%field => l_txc(:, 1)
            nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 2); data_out(nvar)%field => l_txc(:, 2)
            l_txc(:, 1:2) = 0.0_wp
            call FIELD_TO_PARTICLE(data(1:nvar), data_out(1:nvar), l_g, l_q)

            l_q(:, 4) = 0.0_wp
            call THERMO_AIRWATER_LINEAR(l_g%np, l_txc(1, 1), l_q(1, 4))

            l_q(:, 5) = l_q(:, 4) ! l_q(:,6) for bil_cloud_4 is set =0 in dns_main at initialization

        end if

    ! case (PART_TYPE_NEW_CASES)
    case (PART_TYPE_TINIA_1)
        ! call PARTICLE_TINIA_INITIALIZE()

    end select

    return
end subroutine PARTICLE_RANDOM_POSITION

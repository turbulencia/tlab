#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

subroutine PARTICLE_RANDOM_POSITION(l_q, l_txc, l_comm, txc, wrk3d)

    use TLAB_TYPES, only: pointers_dt, pointers3d_dt, wp, wi, longi
    use TLAB_CONSTANTS
    use TLAB_VARS
    use PARTICLE_TYPES, only: particle_dt
    use PARTICLE_VARS
    use THERMO_VARS, only: imixture
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
#endif
    use IO_FIELDS
    implicit none

    real(wp), dimension(isize_part, *), target :: l_q, l_txc
    real(wp), dimension(*), target :: l_comm
    real(wp), dimension(imax, jmax, kmax, *), target :: txc
    real(wp), dimension(*) :: wrk3d

! -------------------------------------------------------------------
    integer(wi) i, is
    integer(wi), allocatable :: x_seed(:)
    integer(wi) size_seed

    real(wp) xref, yref, zref, xscale, yscale, zscale, dy_loc, dummy, dy_frac
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
    l_g%np = INT(isize_part_total/INT(ims_npro, KIND=8))
    if (ims_pro < INT(MOD(isize_part_total, INT(ims_npro, KIND=8)))) then
        l_g%np = l_g%np + 1
    end if
    call MPI_ALLGATHER(l_g%np, 1, MPI_INTEGER4, ims_size_p, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

#else
    l_g%np = INT(isize_part_total)

#endif

! Create tags
    count = 0
#ifdef USE_MPI
    do i = 1, ims_pro
        count = count + INT(ims_size_p(i), KIND=8)
    end do
#endif
    do i = 1, l_g%np
        l_g%tags(i) = INT(i, KIND=8) + count
    end do

! Generate seed - different seed for each processor
    call RANDOM_SEED(SIZE=size_seed)
    allocate (x_seed(size_seed))
#ifdef USE_MPI
    x_seed = [(i, i=1 + ims_pro, size_seed + ims_pro)]
#else
    x_seed = [(i, i=1, size_seed)]
#endif
    call RANDOM_SEED(PUT=x_seed)

#ifdef USE_MPI
    xref = g(1)%nodes(ims_offset_i + 1)
    xscale = g(1)%scale/real(ims_npro_i,wp)
    zref = g(3)%nodes(ims_offset_k + 1)
    zscale = g(3)%scale/real(ims_npro_k,wp)
#else
    xref = g(1)%nodes(1)
    xscale = g(1)%scale
    zref = g(3)%nodes(1)
    zscale = g(3)%scale
#endif
    if (g(3)%size == 1) zscale = 0.0_wp ! 2D case

    ! ########################################################################
    select case (part_ini_mode)

    case (1)
        yref = part_ini_ymean - 0.5_wp*part_ini_thick
        yscale = part_ini_thick

        do i = 1, l_g%np
            call RANDOM_NUMBER(rnd_number(1:3))

            l_q(i, 1) = xref + rnd_number(1)*xscale
            l_q(i, 3) = zref + rnd_number(3)*zscale
            l_q(i, 2) = yref + rnd_number(2)*yscale

        end do

    case (2) ! Use the scalar field to create the particle distribution
        call IO_READ_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, txc, wrk3d)
        is = 1 ! Reference scalar

        y_limits(1) = part_ini_ymean - 0.5_wp*part_ini_thick
        y_limits(2) = part_ini_ymean + 0.5_wp*part_ini_thick
        call PARTICLE_LOCATE_Y(2, y_limits, j_limits, g(2)%size, g(2)%nodes)
        dy_loc = g(2)%nodes(j_limits(2)) - g(2)%nodes(j_limits(1))

        i = 1
        do while (i <= l_g%np)
            call RANDOM_NUMBER(rnd_number(1:3))

            rnd_scal(1) = 1 + floor(rnd_number(1)*imax)
            rnd_scal(3) = 1 + floor(rnd_number(3)*kmax)
            rnd_scal(2) = j_limits(1) + floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))
            dy_frac = rnd_number(2)*(j_limits(2) - j_limits(1) + 1) &
                      - floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))

            dummy = (txc(rnd_scal(1), rnd_scal(2), rnd_scal(3), is) - sbg(is)%mean)/sbg(is)%delta
            dummy = abs(dummy + 0.5_wp)

            call RANDOM_NUMBER(rnd_number_second)

            if (rnd_number_second <= dummy) then

                l_q(i, 1) = xref + rnd_number(1)*xscale
                l_q(i, 3) = zref + rnd_number(3)*zscale
                l_q(i, 2) = g(2)%nodes(rnd_scal(2)) + dy_frac*dy_loc

                i = i + 1

            end if

        end do

    end select

! Calculating closest node below in Y direction
    call PARTICLE_LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, g(2)%size, g(2)%nodes)

!########################################################################
! Remaining scalar properties of the lagrangian field
!########################################################################

    if (imode_part == PART_TYPE_BIL_CLOUD_3 .or. imode_part == PART_TYPE_BIL_CLOUD_4) then

        call IO_READ_FIELDS('scal.ics', IO_SCAL, imax, jmax, kmax, inb_scal, 0, txc, wrk3d)

        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            nvar = 0
            nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 1); data_out(nvar)%field => l_txc(:, 1)
            nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 2); data_out(nvar)%field => l_txc(:, 2)
            l_txc(:, 1:2) = 0.0_wp
            call FIELD_TO_PARTICLE(nvar, data, data_out, l_g, l_q, l_comm, wrk3d)

            l_q(:, 4) = 0.0_wp
            call THERMO_AIRWATER_LINEAR(l_g%np, 1, 1, l_txc(1, 1), l_q(1, 4))

            l_q(:, 5) = l_q(:, 4) ! l_q(:,6) for bil_cloud_4 is set =0 in dns_main at initialization

        end if

    end if

    return
end subroutine PARTICLE_RANDOM_POSITION

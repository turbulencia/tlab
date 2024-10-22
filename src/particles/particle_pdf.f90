#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
subroutine PARTICLE_PDF(fname, s, l_g, l_q, l_txc)

    use TLab_Types, only: pointers_dt, pointers3d_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal_array
    use TLAB_VARS, only: g
    use PARTICLE_TYPES, only: particle_dt
    use PARTICLE_VARS, only: isize_part
    use PARTICLE_VARS, only: particle_pdf_subdomain, particle_pdf_max, particle_pdf_interval
    use PARTICLE_INTERPOLATE
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS
#endif

    implicit none

    character*(*) fname
    TREAL, dimension(isize_field, *), target :: s

    type(particle_dt) :: l_g
    TREAL, dimension(isize_part, *) :: l_q
    TREAL, dimension(isize_part, 1), target :: l_txc

! -------------------------------------------------------------------
    TINTEGER nvar, number_of_bins
    type(pointers3d_dt), dimension(1) :: data
    type(pointers_dt), dimension(1) :: data_out

    TREAL, dimension(:), allocatable :: counter_interval
    TLONGINTEGER, dimension(:, :), allocatable :: particle_bins
#ifdef USE_MPI
    TLONGINTEGER, dimension(:, :), allocatable :: particle_bins_local
#endif

    TINTEGER i, j, is
    TREAL particle_pdf_min

!########################################################################
    number_of_bins = int(particle_pdf_max/particle_pdf_interval)

    allocate (particle_bins(number_of_bins, 3))
    allocate (counter_interval(number_of_bins))
#ifdef USE_MPI
    allocate (particle_bins_local(number_of_bins, 3))
#endif

    particle_bins = int(0, KIND=8)

    nvar = 0
    nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => s(:, inb_scal_array); data_out(nvar)%field => l_txc(:, 1)
    l_txc(:, 1) = C_0_R
    call FIELD_TO_PARTICLE(data(1:nvar), data_out(1:nvar), l_g, l_q)

!########################################################################
! Calculating
!########################################################################
    particle_pdf_min = C_0_R  !if needed for future

    particle_bins = 0

    do i = 1, l_g%np

        if (l_q(i, 1)/g(1)%scale >= particle_pdf_subdomain(1) .and. l_q(i, 1)/g(1)%scale <= particle_pdf_subdomain(2)) then
            if (l_q(i, 2)/g(2)%scale >= particle_pdf_subdomain(3) .and. l_q(i, 2)/g(2)%scale <= particle_pdf_subdomain(4)) then
                if (l_q(i, 3)/g(3)%scale >= particle_pdf_subdomain(5) .and. l_q(i, 3)/g(3)%scale <= particle_pdf_subdomain(6)) then

                    j = 1 + int((l_txc(i, 1) - particle_pdf_min)/particle_pdf_interval)
                    particle_bins(j, 1) = particle_bins(j, 1) + 1

                    do is = 4, 5
                        j = 1 + int((l_q(i, is) - particle_pdf_min)/particle_pdf_interval)
                        particle_bins(j, is - 2) = particle_bins(j, is - 2) + 1
                    end do

                end if
            end if
        end if

    end do

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    call MPI_REDUCE(particle_bins, particle_bins_local, number_of_bins*3, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    particle_bins = particle_bins_local
#endif

!#######################################################################
! Writing
!#######################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        counter_interval(1) = 0
        counter_interval(2) = particle_pdf_interval
        do i = 3, number_of_bins
            counter_interval(i) = counter_interval(i - 1) + particle_pdf_interval
        end do

        open (unit=116, file=fname)
        do i = 1, number_of_bins
       write (116, '(F6.3, I20.1, I20.1, I20.1)') counter_interval(i), particle_bins(i, 1), particle_bins(i, 2), particle_bins(i, 3)
        end do
        close (116)
#ifdef USE_MPI
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    deallocate (particle_bins_local)
#endif

    deallocate (particle_bins)
    deallocate (counter_interval)

    return
end subroutine PARTICLE_PDF

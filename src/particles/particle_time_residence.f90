#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# DESCRIPTION
!#
!# Calculates the residence times for the lagrangian particles
!#
!########################################################################
subroutine PARTICLE_TIME_RESIDENCE(dtime, particle_number, l_q)

    use PARTICLE_VARS, only: isize_part, inb_part_array
    use PARTICLE_VARS, only: l_y_lambda, l_y_base

    implicit none

    TREAL dtime
    TINTEGER particle_number
    TREAL, dimension(isize_part, *) :: l_q

    TINTEGER i

    do i = 1, particle_number
        if (l_q(i, 2) > l_y_lambda) then
            l_q(i, inb_part_array - 1) = l_q(i, inb_part_array - 1) + dtime   !time cloud droplets spend on cloud-top
        end if
        if (l_q(i, 2) > l_y_base) then
            l_q(i, inb_part_array) = l_q(i, inb_part_array) + dtime   !time cloud droplets spend in intermediate 2/3 of cloud
        elseif (l_q(i, 2) <= l_y_base) then
            l_q(i, inb_part_array - 1) = C_0_R   !cloud droplets loose memory when "leaving" cloud
            l_q(i, inb_part_array) = C_0_R   !cloud droplets loose memory when "leaving" cloud
        end if
    end do

    return
end subroutine PARTICLE_TIME_RESIDENCE

!########################################################################
!########################################################################
subroutine PARTICLE_RESIDENCE_PDF(fname, particle_number, l_q)

    use PARTICLE_VARS, only: isize_part, inb_part_array
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS
#endif

    implicit none

    character*(*) fname
    TINTEGER particle_number
    TREAL, dimension(isize_part, *) :: l_q

! -------------------------------------------------------------------
    TLONGINTEGER, dimension(:, :), allocatable :: residence_bins
    TREAL, dimension(:), allocatable :: residence_counter_interval
    TINTEGER i, j
    TINTEGER residence_tmax, residence_nbins
    TREAL residence_pdf_interval
#ifdef USE_MPI
    TLONGINTEGER, dimension(:, :), allocatable :: residence_bins_local
#endif

! #####################################################################
    residence_tmax = 100
    residence_nbins = 1000
    residence_pdf_interval = real(residence_tmax)/real(residence_nbins)

!#######################################################################
    allocate (residence_bins(residence_nbins, 2))
    allocate (residence_counter_interval(residence_nbins))
#ifdef USE_MPI
    allocate (residence_bins_local(residence_nbins, 2))
#endif

    residence_bins = int(0, KIND=8)

!#######################################################################
!Start counting of particles
!#######################################################################
    do i = 1, particle_number
        j = 1 + int(l_q(i, inb_part_array - 1)/residence_pdf_interval)
        residence_bins(j, 1) = residence_bins(j, 1) + 1
        j = 1 + int(l_q(i, inb_part_array)/residence_pdf_interval)
        residence_bins(j, 2) = residence_bins(j, 2) + 1
    end do

#ifdef USE_MPI
!#######################################################################
!Reduce all information to root
!#######################################################################
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    call MPI_REDUCE(residence_bins, residence_bins_local, residence_nbins*2, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    residence_bins = residence_bins_local

!#######################################################################
!Create interval for writing
!#######################################################################
    if (ims_pro == 0) then
#endif
        residence_counter_interval(1) = 0
        residence_counter_interval(2) = residence_pdf_interval
        do i = 3, residence_nbins
            residence_counter_interval(i) = residence_counter_interval(i - 1) + residence_pdf_interval
        end do

!#######################################################################
!Write data to file
!#######################################################################
        open (unit=116, file=fname)
        do i = 1, residence_nbins
            write (116, '(F6.3, I20.1, I20.1)') residence_counter_interval(i), residence_bins(i, 1), residence_bins(i, 2)
        end do
        close (116)

#ifdef USE_MPI
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    deallocate (residence_bins_local)
#endif

    deallocate (residence_bins)
    deallocate (residence_counter_interval)

    return
end subroutine PARTICLE_RESIDENCE_PDF

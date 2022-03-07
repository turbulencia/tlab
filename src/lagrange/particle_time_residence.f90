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
SUBROUTINE PARTICLE_TIME_RESIDENCE(dtime, particle_number, l_q)

  USE TLAB_VARS,      ONLY : isize_particle, inb_part_array
  USE LAGRANGE_VARS, ONLY : l_y_lambda, l_y_base

  IMPLICIT NONE

  TREAL dtime
  TINTEGER particle_number
  TREAL, DIMENSION(isize_particle,*) :: l_q

  TINTEGER i

  DO i=1,particle_number
     IF (l_q(i,2) .GT. l_y_lambda)THEN
        l_q(i,inb_part_array-1) = l_q(i,inb_part_array-1) + dtime   !time cloud droplets spend on cloud-top
     ENDIF
     IF (l_q(i,2) .GT. l_y_base)THEN
        l_q(i,inb_part_array  ) = l_q(i,inb_part_array  ) + dtime   !time cloud droplets spend in intermediate 2/3 of cloud
     ELSEIF (l_q(i,2) .LE. l_y_base)THEN
        l_q(i,inb_part_array-1)=C_0_R   !cloud droplets loose memory when "leaving" cloud
        l_q(i,inb_part_array  )=C_0_R   !cloud droplets loose memory when "leaving" cloud
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE PARTICLE_TIME_RESIDENCE

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_RESIDENCE_PDF(fname, particle_number, l_q)

  USE TLAB_VARS, ONLY : isize_particle, inb_part_array
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

  CHARACTER*(*) fname
  TINTEGER particle_number
  TREAL, DIMENSION(isize_particle,*) :: l_q

! -------------------------------------------------------------------
  TLONGINTEGER, DIMENSION(:,:), ALLOCATABLE :: residence_bins
  TREAL,        DIMENSION(:),   ALLOCATABLE :: residence_counter_interval
  TINTEGER i,j
  TINTEGER residence_tmax, residence_nbins
  TREAL residence_pdf_interval
#ifdef USE_MPI
  TLONGINTEGER, DIMENSION(:,:), ALLOCATABLE :: residence_bins_local
#endif

! #####################################################################
  residence_tmax = 100
  residence_nbins = 1000
  residence_pdf_interval = real(residence_tmax)/real(residence_nbins)

!#######################################################################
  ALLOCATE(residence_bins(residence_nbins,2))
  ALLOCATE(residence_counter_interval(residence_nbins))
#ifdef USE_MPI
  ALLOCATE(residence_bins_local(residence_nbins,2))
#endif

  residence_bins=int(0,KIND=8)

!#######################################################################
!Start counting of particles
!#######################################################################
  DO i=1,particle_number
     j = 1 + int( l_q(i,inb_part_array-1) / residence_pdf_interval )
     residence_bins(j,1)=residence_bins(j,1)+1
     j = 1 + int( l_q(i,inb_part_array  ) / residence_pdf_interval )
     residence_bins(j,2)=residence_bins(j,2)+1
  ENDDO

#ifdef USE_MPI
!#######################################################################
!Reduce all information to root
!#######################################################################
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  CALL MPI_REDUCE(residence_bins, residence_bins_local, residence_nbins*2, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
  residence_bins = residence_bins_local

!#######################################################################
!Create interval for writing
!#######################################################################
  IF(ims_pro .EQ. 0) THEN
#endif
     residence_counter_interval(1)=0
     residence_counter_interval(2)=residence_pdf_interval
     DO i=3,residence_nbins
        residence_counter_interval(i)=residence_counter_interval(i-1)+residence_pdf_interval
     ENDDO

!#######################################################################
!Write data to file
!#######################################################################
     OPEN(unit=116, file=fname)
     DO i=1,residence_nbins
        WRITE (116,'(F6.3, I20.1, I20.1)') residence_counter_interval(i), residence_bins(i,1), residence_bins(i,2)
     END DO
     CLOSE(116)

#ifdef USE_MPI
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  DEALLOCATE(residence_bins_local)
#endif

  DEALLOCATE(residence_bins)
  DEALLOCATE(residence_counter_interval)

  RETURN
END SUBROUTINE PARTICLE_RESIDENCE_PDF

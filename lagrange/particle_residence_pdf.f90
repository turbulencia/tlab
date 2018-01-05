#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the pdf for a certain region
!# 
!########################################################################
SUBROUTINE PARTICLE_RESIDENCE_PDF(fname,l_hq,l_q)

  USE DNS_GLOBAL,     ONLY : isize_particle, inb_particle
  USE LAGRANGE_GLOBAL,ONLY : particle_number_local
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q, l_hq 

  TLONGINTEGER, DIMENSION(:,:),   ALLOCATABLE         :: residence_bins
  TREAL, DIMENSION(:),   ALLOCATABLE         :: residence_counter_interval
  CHARACTER(*)  ::  fname
  TINTEGER i,j
  TINTEGER residence_tmax, residence_nbins
  TREAL residence_pdf_interval
#ifdef USE_MPI
  TLONGINTEGER, DIMENSION(:,:),   ALLOCATABLE         :: residence_bins_local
#endif 

! #####################################################################
  residence_tmax = 100
  residence_nbins = 1000
  residence_pdf_interval = real(residence_tmax)/real(residence_nbins)
#ifdef USE_MPI
  ALLOCATE(residence_bins_local(residence_nbins,2)) 
  residence_bins_local = 0.0
#endif 

  ALLOCATE(residence_bins(residence_nbins,2)) 
  ALLOCATE(residence_counter_interval(residence_nbins)) 
  residence_bins=int(0,KIND=8)

#ifdef USE_MPI
!#######################################################################
!Start counting of particles in bins per processor
!#######################################################################

  DO i=1,particle_number_local
     j = 1 + int( l_q(i,inb_particle) / residence_pdf_interval ) !if residence is calles information is on inb_particle=6
     residence_bins_local(j,1)=residence_bins_local(j,1)+1

     j = 1 + int( l_hq(i,inb_particle) / residence_pdf_interval ) !if residence is calles information is on inb_particle=6
     residence_bins_local(j,2)=residence_bins_local(j,2)+1
  ENDDO
!#######################################################################
!Reduce all information to root
!#######################################################################

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  CALL MPI_REDUCE(residence_bins_local, residence_bins, residence_nbins*2, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)

!#######################################################################
!Create interval for writing
!#######################################################################
  IF(ims_pro .EQ. 0) THEN
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
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  DEALLOCATE(residence_bins_local)
#else

  DO i=1,particle_number_local
     j = 1 + int( l_q(i,inb_particle) / residence_pdf_interval ) !if residence is calles information is on inb_particle=6
     residence_bins(j,1)=residence_bins(j,1)+1
     j = 1 + int( l_hq(i,inb_particle) / residence_pdf_interval ) !if residence is calles information is on inb_particle=6
     residence_bins(j,2)=residence_bins(j,2)+1
  ENDDO

!#######################################################################
!Create interval for writing
!#######################################################################
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

#endif 

  DEALLOCATE(residence_bins)
  DEALLOCATE(residence_counter_interval)

  RETURN
END SUBROUTINE PARTICLE_RESIDENCE_PDF

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_PDF(fname,s, l_g,l_q,l_txc,l_comm, wrk3d)

  USE TLAB_TYPES,  ONLY: pointers_dt, pointers3d_dt
  USE TLAB_VARS, ONLY: imax,jmax,kmax, isize_field,isize_particle, inb_scal_array
  USE TLAB_VARS, ONLY: g
  USE LAGRANGE_VARS, ONLY : particle_dt
  USE LAGRANGE_VARS, ONLY : particle_pdf_subdomain, particle_pdf_max, particle_pdf_interval
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  

  CHARACTER*(*) fname
  TREAL, DIMENSION(isize_field,*), TARGET :: s
  TREAL, DIMENSION(*)                     :: wrk3d

  TYPE(particle_dt)                          :: l_g
  TREAL, DIMENSION(isize_particle,*)         :: l_q
  TREAL, DIMENSION(isize_particle,1), TARGET :: l_txc
  TREAL, DIMENSION(*)                        :: l_comm

! -------------------------------------------------------------------  
  TINTEGER nvar, number_of_bins
  TYPE(pointers3d_dt), DIMENSION(1) :: data
  TYPE(pointers_dt),   DIMENSION(1) :: data_out

  TREAL,        DIMENSION(:),   ALLOCATABLE :: counter_interval
  TLONGINTEGER, DIMENSION(:,:), ALLOCATABLE :: particle_bins
#ifdef USE_MPI
  TLONGINTEGER, DIMENSION(:,:), ALLOCATABLE :: particle_bins_local
#endif 

  TINTEGER i,j, is
  TREAL particle_pdf_min

!########################################################################
  number_of_bins = INT(particle_pdf_max/particle_pdf_interval)
  
  ALLOCATE(particle_bins(number_of_bins,3)) 
  ALLOCATE(counter_interval(number_of_bins)) 
#ifdef USE_MPI
  ALLOCATE(particle_bins_local(number_of_bins,3)) 
#endif 

  particle_bins=int(0,KIND=8)

  nvar = 0
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => s(:,inb_scal_array); data_out(nvar)%field => l_txc(:,1)
  l_txc(:,1) = C_0_R
  CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_g,l_q,l_comm, wrk3d)

!########################################################################
! Calculating
!########################################################################
  particle_pdf_min = C_0_R  !if needed for future 

  particle_bins = 0

  DO i=1,l_g%np
     
     IF ( l_q(i,1)/g(1)%scale .GE. particle_pdf_subdomain(1) .AND. l_q(i,1)/g(1)%scale .LE. particle_pdf_subdomain(2)) THEN
        IF ( l_q(i,2)/g(2)%scale .GE. particle_pdf_subdomain(3) .AND. l_q(i,2)/g(2)%scale .LE. particle_pdf_subdomain(4)) THEN
           IF ( l_q(i,3)/g(3)%scale .GE. particle_pdf_subdomain(5) .AND. l_q(i,3)/g(3)%scale .LE. particle_pdf_subdomain(6)) THEN

              j = 1 + int( (l_txc(i,1) - particle_pdf_min) / particle_pdf_interval )
              particle_bins(j,1)=particle_bins(j,1)+1
              
              DO is=4,5
                 j = 1 + int( (l_q(i,is) - particle_pdf_min) / particle_pdf_interval )
                 particle_bins(j,is-2)=particle_bins(j,is-2)+1
              ENDDO
              
           ENDIF
        ENDIF
     ENDIF
     
  ENDDO

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  CALL MPI_REDUCE(particle_bins, particle_bins_local, number_of_bins*3, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
  particle_bins = particle_bins_local
#endif
  
!#######################################################################
! Writing
!#######################################################################
#ifdef USE_MPI 
  IF( ims_pro .EQ. 0) THEN
#endif
     counter_interval(1)=0
     counter_interval(2)=particle_pdf_interval
     DO i=3,number_of_bins
        counter_interval(i)=counter_interval(i-1)+particle_pdf_interval
     ENDDO

     OPEN(unit=116, file=fname)
     DO i=1,number_of_bins
        WRITE (116,'(F6.3, I20.1, I20.1, I20.1)') counter_interval(i), particle_bins(i,1), particle_bins(i,2), particle_bins(i,3)
     END DO
     CLOSE(116)
#ifdef USE_MPI 
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  DEALLOCATE(particle_bins_local)
#endif

  DEALLOCATE(particle_bins)
  DEALLOCATE(counter_interval)

  RETURN
END SUBROUTINE PARTICLE_PDF

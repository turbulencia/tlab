#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_PDF(fname,s, wrk2d,wrk3d, l_txc,l_tags,l_q,l_comm)

  USE DNS_TYPES,  ONLY: pointers_dt, pointers3d_dt
  USE DNS_GLOBAL, ONLY: imax,jmax,kmax, isize_field,isize_particle, inb_particle, inb_scal_array
  USE DNS_GLOBAL, ONLY: g
  USE LAGRANGE_GLOBAL, ONLY :  particle_number_local
  USE LAGRANGE_GLOBAL, ONLY :  number_of_bins, y_particle_pdf_pos, y_particle_pdf_width
  USE LAGRANGE_GLOBAL, ONLY :  y_particle_pdf_pos, y_particle_pdf_width
  USE LAGRANGE_GLOBAL, ONLY :  x_particle_pdf_pos, x_particle_pdf_width
  USE LAGRANGE_GLOBAL, ONLY :  z_particle_pdf_pos, z_particle_pdf_width
  USE LAGRANGE_GLOBAL, ONLY :  particle_pdf_interval
#ifdef USE_MPI
  USE DNS_MPI
#endif


  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
  TREAL, DIMENSION(isize_field,*), TARGET :: s
  TREAL, DIMENSION(*)                     :: wrk2d, wrk3d

  TREAL, DIMENSION(isize_particle,inb_particle)        :: l_q
  TREAL, DIMENSION(isize_particle,1),           TARGET :: l_txc
  TREAL, DIMENSION(*)                :: l_comm
  INTEGER(8), DIMENSION(*)           :: l_tags

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(1) :: data
  TYPE(pointers_dt),   DIMENSION(1) :: data_out

  TLONGINTEGER, DIMENSION(:,:),   ALLOCATABLE         :: particle_bins
  TREAL, DIMENSION(:),   ALLOCATABLE         :: counter_interval
  CHARACTER(*)  ::  fname
  TREAL y_pdf_max, y_pdf_min
  TREAL x_pdf_max, x_pdf_min
  TREAL z_pdf_max, z_pdf_min
  TINTEGER i,j, is, particle_pdf_min
#ifdef USE_MPI
  TLONGINTEGER, DIMENSION(:,:),   ALLOCATABLE         :: particle_bins_local
  ALLOCATE(particle_bins_local(number_of_bins,3)) 
#endif 

  ALLOCATE(particle_bins(number_of_bins,3)) 
  ALLOCATE(counter_interval(number_of_bins)) 

  y_pdf_max=y_particle_pdf_pos+0.5*y_particle_pdf_width
  y_pdf_min=y_particle_pdf_pos-0.5*y_particle_pdf_width
  particle_bins=int(0,KIND=8)

  IF (x_particle_pdf_width .NE. 0) THEN
     x_pdf_max=x_particle_pdf_pos+0.5*x_particle_pdf_width
     x_pdf_min=x_particle_pdf_pos-0.5*x_particle_pdf_width
  ENDIF
  IF (z_particle_pdf_width .NE. 0) THEN
     z_pdf_max=z_particle_pdf_pos+0.5*z_particle_pdf_width
     z_pdf_min=z_particle_pdf_pos-0.5*z_particle_pdf_width
  ENDIF

  nvar = 0
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => s(:,inb_scal_array); data_out(nvar)%field => l_txc(:,1)
  CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_q,l_tags,l_comm, wrk2d,wrk3d)

#ifdef USE_MPI

  particle_bins_local=0.0

!#######################################################################
!Start counting of particles in bins per processor
!#######################################################################
  particle_pdf_min = 0  !if needed for future 

  IF (x_particle_pdf_width .EQ. 0) THEN !ONLY PART OF Y
     DO i=1,particle_number_local
        IF ( l_q(i,2)/g(2)%scale .GE. y_pdf_min .AND. l_q(i,2)/g(2)%scale .LE. y_pdf_max) THEN
           j = 1 + int( (l_txc(i,1) - particle_pdf_min) / particle_pdf_interval )
           particle_bins_local(j,1)=particle_bins_local(j,1)+1

           DO is=4,5
              j = 1 + int( (l_q(i,is) - particle_pdf_min) / particle_pdf_interval )
              particle_bins_local(j,is-2)=particle_bins_local(j,is-2)+1
           ENDDO
        ENDIF
     ENDDO
  ELSE !3D BOX
     DO i=1,particle_number_local
        IF ( l_q(i,1)/g(1)%scale .GE. x_pdf_min .AND. l_q(i,1)/g(1)%scale .LE. x_pdf_max) THEN
           IF ( l_q(i,2)/g(2)%scale .GE. y_pdf_min .AND. l_q(i,2)/g(2)%scale .LE. y_pdf_max) THEN
              IF ( l_q(i,3)/g(3)%scale .GE. z_pdf_min .AND. l_q(i,3)/g(3)%scale .LE. z_pdf_max) THEN
                 j = 1 + int( (l_txc(i,1) - particle_pdf_min) / particle_pdf_interval )
                 particle_bins_local(j,1)=particle_bins_local(j,1)+1

                 DO is=4,5
                    j = 1 + int( (l_q(i,is) - particle_pdf_min) / particle_pdf_interval )
                    particle_bins_local(j,is-2)=particle_bins_local(j,is-2)+1
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

!#######################################################################
!Reduce all information to root
!#######################################################################

  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  CALL MPI_REDUCE(particle_bins_local, particle_bins, number_of_bins*3, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)

!#######################################################################
!Create interval for writing
!#######################################################################
  IF(ims_pro .EQ. 0) THEN
     counter_interval(1)=0
     counter_interval(2)=particle_pdf_interval
     DO i=3,number_of_bins
        counter_interval(i)=counter_interval(i-1)+particle_pdf_interval
     ENDDO

!#######################################################################
!Write data to file
!#######################################################################
     OPEN(unit=116, file=fname)
     DO i=1,number_of_bins
        WRITE (116,'(F6.3, I20.1, I20.1, I20.1)') counter_interval(i), particle_bins(i,1), particle_bins(i,2), particle_bins(i,3)
     END DO
     CLOSE(116)
  END IF


  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
  DEALLOCATE(particle_bins_local)
#else

  particle_pdf_min = 0

  DO i=1,particle_number_local
     IF ( l_q(i,2)/g(2)%scale .GE. y_pdf_min .AND. l_q(i,2)/g(2)%scale .LE. y_pdf_max) THEN
        j = 1 + int( (l_txc(i,1) - particle_pdf_min) / particle_pdf_interval )
        particle_bins(j,1)=particle_bins(j,1)+1

        DO is=4,5
           j = 1 + int( (l_q(i,is) - particle_pdf_min) / particle_pdf_interval )
           particle_bins(j,is-2)=particle_bins(j,is-2)+1
        ENDDO

     ENDIF
  ENDDO

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


#endif 

  DEALLOCATE(particle_bins)
  DEALLOCATE(counter_interval)

  RETURN
END SUBROUTINE PARTICLE_PDF

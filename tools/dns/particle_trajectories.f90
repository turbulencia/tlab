#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

MODULE PARTICLE_TRAJECTORIES

  USE DNS_CONSTANTS,  ONLY : efile, lfile
  USE DNS_GLOBAL,     ONLY : inb_flow_array, inb_scal_array
  USE DNS_GLOBAL,     ONLY : isize_particle, inb_particle
  USE LAGRANGE_GLOBAL,ONLY : particle_dt
  USE LAGRANGE_GLOBAL,ONLY : isize_trajectory, inb_trajectory, isize_l_comm, itrajectory
#ifdef USE_MPI
  USE DNS_MPI,        ONLY : ims_pro, ims_err
#endif

  IMPLICIT NONE
  SAVE
  
  TREAL,      DIMENSION(:,:,:), ALLOCATABLE :: l_trajectories
  INTEGER(8), DIMENSION(:),     ALLOCATABLE :: l_trajectories_tags
  TINTEGER                                  :: counter, isize_time

CONTAINS
  
!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_TRAJECTORIES_INITIALIZE(nitera_save, nitera_last)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nitera_save, nitera_last

! -------------------------------------------------------------------
  CHARACTER*32 name
  CHARACTER*128 str, line
  TINTEGER ims_npro_loc

  TINTEGER j, ierr

!#######################################################################
  isize_time = nitera_save

! Adding space for saving the time
  WRITE(str,*) (isize_trajectory+1)*isize_time; line = 'Allocating array l_trajectories of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) inb_trajectory; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_trajectories(isize_trajectory+1,isize_time,inb_trajectory),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  
  WRITE(str,*) isize_trajectory; line = 'Allocating array l_trajectory_tags of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_trajectories_tags(isize_trajectory),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories_tags.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  
! Initialize
  IF ( itrajectory .EQ. LAG_TRAJECTORY_LARGEST ) THEN ! Read file with tags of largest particles, the ones to track
     WRITE(name,*) nitera_last; name='largest_particle.'//TRIM(ADJUSTL(name))
#ifdef USE_MPI
     IF (ims_pro .EQ. 0) THEN
#endif
#define LOC_UNIT_ID 117
#define LOC_STATUS 'old'
#include "dns_open_file.h"
        READ(LOC_UNIT_ID) ims_npro_loc
!        READ(LOC_UNIT_ID) ims_size_p(1:ims_npro_loc)
        READ(LOC_UNIT_ID, POS=SIZEOFINT*(ims_npro_loc+1) +1) l_trajectories_tags
        CLOSE(LOC_UNIT_ID)
#undef LOC_UNIT_ID
#undef LOC_STATUS
#ifdef USE_MPI
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
     CALL MPI_BCAST(l_trajectories_tags,isize_trajectory,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)
#endif

  ELSE ! default track only isize_trajectory particles 
     DO j=1,isize_trajectory              
        l_trajectories_tags(j) = INT(j, KIND=8)
     ENDDO
          
  ENDIF

  l_trajectories = C_0_R
  counter        = 0
  
  RETURN
END SUBROUTINE PARTICLE_TRAJECTORIES_INITIALIZE

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_TRAJECTORIES_ACCUMULATE(q,s, txc, l_g,l_q,l_hq,l_txc,l_comm, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL,ONLY : isize_field, imax,jmax,kmax
  USE DNS_GLOBAL,ONLY : rtime
  
  IMPLICIT NONE

  TREAL,      DIMENSION(isize_field,*),    TARGET :: q, s, txc
  TYPE(particle_dt)                               :: l_g
  TREAL,      DIMENSION(isize_particle,*), TARGET :: l_q, l_hq, l_txc ! l_hq as aux array
  TREAL,      DIMENSION(isize_l_comm)             :: l_comm
  TREAL,      DIMENSION(*)                        :: wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER i, j
  TINTEGER iv, nvar
  TYPE(pointers3d_dt), DIMENSION(inb_trajectory) :: data_in
  TYPE(pointers_dt),   DIMENSION(inb_trajectory) :: data

!#######################################################################
  counter = counter +1
  
! -------------------------------------------------------------------
! Setting pointers to position
  nvar = 0
  nvar = nvar+1; data(nvar)%field => l_q(:,1)
  nvar = nvar+1; data(nvar)%field => l_q(:,2)
  nvar = nvar+1; data(nvar)%field => l_q(:,3)

! Primitive variables
  DO iv = 4,inb_flow_array
     nvar = nvar+1; data_in(nvar)%field(1:imax,1:jmax,1:kmax) => q(:,iv);   data(nvar)%field => l_txc(:,iv-3)
     l_txc(:,iv-3) = C_0_R ! Field to particle is additive
  ENDDO
  
  DO iv = 1,inb_scal_array
     nvar = nvar+1; data_in(nvar)%field(1:imax,1:jmax,1:kmax) => s(:,iv);   data(nvar)%field => l_txc(:,iv-3+inb_flow_array)
     l_txc(:,iv-3+inb_flow_array) = C_0_R ! Field to particle is additive
  ENDDO

! -------------------------------------------------------------------
! Additional information
  IF ( itrajectory .EQ. LAG_TRAJECTORY_VORTICITY ) THEN
     CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), txc(1,4), wrk2d,wrk3d)
     nvar = nvar+1; data_in(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,1); data(nvar)%field => l_hq(:,1)
     nvar = nvar+1; data_in(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,2); data(nvar)%field => l_hq(:,2)
     nvar = nvar+1; data_in(nvar)%field(1:imax,1:jmax,1:kmax) => txc(:,3); data(nvar)%field => l_hq(:,3)
     l_hq(:,1:3) = C_0_R ! Field to particle is additive
  ENDIF

! -------------------------------------------------------------------
! Interpolation
  IF ( nvar-3 .GT. 0 ) THEN
     iv = nvar -3
     CALL FIELD_TO_PARTICLE(iv, data_in(4), data(4), l_g,l_q,l_comm, wrk2d,wrk3d)
  ENDIF
  
! -------------------------------------------------------------------
! Accumulate time
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     l_trajectories(1,counter,1:nvar) = rtime
#ifdef USE_MPI
  ENDIF
#endif
  
! Accumulating the data
  DO i = 1,l_g%np
     DO j = 1,isize_trajectory
        IF ( l_g%tags(i) .EQ. l_trajectories_tags(j) ) THEN
           DO iv = 1,nvar
              l_trajectories(1+j,counter,iv) = data(iv)%field(i)
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE PARTICLE_TRAJECTORIES_ACCUMULATE

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_TRAJECTORIES_WRITE(fname, wrk3d)

  IMPLICIT NONE
  
#ifdef USE_MPI
#include "mpif.h"
#endif
  
  CHARACTER*(*) fname
  TREAL, DIMENSION(1+isize_trajectory,isize_time) :: wrk3d

! -------------------------------------------------------------------
  CHARACTER(len=32) name
  TINTEGER iv
  
!#######################################################################
  DO iv = 1,inb_trajectory
#ifdef USE_MPI
     wrk3d = C_0_R
     CALL MPI_REDUCE(l_trajectories(1,1,iv), wrk3d, (1+isize_trajectory)*isize_time, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
     CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
     IF(ims_pro .EQ. 0) THEN
#endif
        
        WRITE(name,*) iv; name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))
#define LOC_UNIT_ID 115
#define LOC_STATUS 'unknown'        
#include "dns_open_file.h"
        REWIND(LOC_UNIT_ID)
#ifdef USE_MPI
        WRITE(LOC_UNIT_ID) SNGL(wrk3d)
#else        
        WRITE(LOC_UNIT_ID) SNGL(l_trajectories(:,:,iv))
#endif
        CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
     END IF
#endif
  ENDDO
  
  l_trajectories = C_0_R
  counter        = 0
  
  RETURN
END SUBROUTINE PARTICLE_TRAJECTORIES_WRITE

END MODULE PARTICLE_TRAJECTORIES

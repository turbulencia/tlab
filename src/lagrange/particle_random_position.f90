#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE PARTICLE_RANDOM_POSITION(l_g,l_q,l_txc,l_comm, txc, wrk3d)
  
  USE TLAB_TYPES,  ONLY : pointers_dt, pointers3d_dt 
  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE LAGRANGE_GLOBAL, ONLY : particle_dt, particle_number_total
  USE LAGRANGE_GLOBAL, ONLY : particle_rnd_mode, y_particle_pos, y_particle_width, ilagrange
  USE THERMO_VARS,   ONLY : imixture
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TYPE(particle_dt)                               :: l_g
  TREAL,      DIMENSION(isize_particle,*), TARGET :: l_q, l_txc
  TREAL,      DIMENSION(*),                TARGET :: l_comm
  TREAL,      DIMENSION(imax,jmax,kmax,*), TARGET :: txc
  TREAL,      DIMENSION(*)                        :: wrk3d

! -------------------------------------------------------------------
  TINTEGER  i, j, is
  TINTEGER, ALLOCATABLE :: x_seed(:)
  TINTEGER size_seed
  
  TREAL xref,yref,zref, xscale,yscale,zscale, dy_loc, dummy, dy_frac
  TREAL y_limits(2)
  TINTEGER j_limits(2)
  
  TREAL rnd_number(4), rnd_number_second
  TINTEGER rnd_scal(3)

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(2) :: data
  TYPE(pointers_dt),   DIMENSION(2) :: data_out

  TLONGINTEGER count
  
!########################################################################
#ifdef USE_MPI
  l_g%np = INT( particle_number_total /INT(ims_npro, KIND=8) )
  IF ( ims_pro .LT. INT( MOD(particle_number_total, INT(ims_npro, KIND=8)) ) ) THEN
     l_g%np = l_g%np +1
  ENDIF
  CALL MPI_ALLGATHER(l_g%np,1,MPI_INTEGER4,ims_size_p,1,MPI_INTEGER4,MPI_COMM_WORLD,ims_err)

#else
  l_g%np = INT(particle_number_total)

#endif

! Create tags
  count = 0
#ifdef USE_MPI
  DO i = 1,ims_pro
     count = count +INT(ims_size_p(i),KIND=8)
  ENDDO
#endif
  DO i = 1,l_g%np
     l_g%tags(i) = INT(i, KIND=8) +count
  END DO
  
! Generate seed - different seed for each processor
  CALL RANDOM_SEED(SIZE=size_seed)
  ALLOCATE(x_seed(size_seed))
#ifdef USE_MPI
  x_seed = (/ (j,j=1+ims_pro, size_seed+ims_pro) /)
#else
  x_seed = (/ (j,j=1, size_seed) /)
#endif
  CALL RANDOM_SEED(PUT=x_seed)

#ifdef USE_MPI
  xref   = g(1)%nodes(ims_offset_i +1)
  xscale = g(1)%scale /M_REAL( ims_npro_i )
  zref   = g(3)%nodes(ims_offset_k +1)
  zscale = g(3)%scale /M_REAL( ims_npro_k )
#else
  xref   = g(1)%nodes(1)
  xscale = g(1)%scale
  zref   = g(3)%nodes(1)
  zscale = g(3)%scale
#endif
  IF ( g(3)%size .EQ. 1 ) zscale = C_0_R ! 2D case
     
!########################################################################
  IF ( particle_rnd_mode .EQ. 1 ) THEN
     
     yref   = y_particle_pos -C_05_R *y_particle_width
     yscale = y_particle_width
     
     DO i = 1,l_g%np
        DO j = 1,3 ! three directions
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO
        
        l_q(i,1) = xref + rnd_number(1) *xscale
        l_q(i,3) = zref + rnd_number(3) *zscale
        l_q(i,2)=  yref + rnd_number(2) *yscale
        
     END DO

!########################################################################
! Use the scalar field to create the particle distribution
  ELSE IF ( particle_rnd_mode .EQ. 2 ) THEN
     
     CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d)
     is = 1 ! Reference scalar

     y_limits(1) = y_particle_pos -C_05_R *y_particle_width
     y_limits(2) = y_particle_pos +C_05_R *y_particle_width
     CALL PARTICLE_LOCATE_Y( 2, y_limits, j_limits, g(2)%size, g(2)%nodes )
     dy_loc= g(2)%nodes(j_limits(2)) -g(2)%nodes(j_limits(1))
     
     i = 1
     DO WHILE ( i .LE. l_g%np ) 
        DO j=1,3
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO

        rnd_scal(1) = 1           +floor(rnd_number(1)*imax)
        rnd_scal(3) = 1           +floor(rnd_number(3)*kmax)
        rnd_scal(2) = j_limits(1) +floor(rnd_number(2)*(j_limits(2)-j_limits(1)+1))        
        dy_frac     =                    rnd_number(2)*(j_limits(2)-j_limits(1)+1) &
                                 -floor(rnd_number(2)*(j_limits(2)-j_limits(1)+1))
        
        dummy = ( txc(rnd_scal(1),rnd_scal(2),rnd_scal(3),is) -sbg(is)%mean )/sbg(is)%delta
        dummy = abs( dummy + C_05_R )

        CALL RANDOM_NUMBER(rnd_number_second)
        
        IF ( rnd_number_second .LE. dummy ) THEN
           
           l_q(i,1) = xref + rnd_number(1) *xscale
           l_q(i,3) = zref + rnd_number(3) *zscale
           l_q(i,2) = g(2)%nodes(rnd_scal(2)) + dy_frac *dy_loc
        
           i=i+1

        END IF
        
     END DO

  ENDIF
  
!########################################################################
! Remaining scalar properties of the lagrangian field
!########################################################################
! Calculating closest node below in Y direction
  CALL PARTICLE_LOCATE_Y( l_g%np, l_q(1,2), l_g%nodes, g(2)%size, g(2)%nodes )
  
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN

     CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d)
     
     IF ( imixture .EQ.  MIXT_TYPE_AIRWATER_LINEAR ) THEN
        nvar = 0
        nvar = nvar+1; data(nvar)%field => txc(:,:,:,1); data_out(nvar)%field => l_txc(:,1)
        nvar = nvar+1; data(nvar)%field => txc(:,:,:,2); data_out(nvar)%field => l_txc(:,2)
        l_txc(:,1:2) = C_0_R
        CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_g,l_q,l_comm, wrk3d)
        
        l_q(:,4) = C_0_R
        CALL THERMO_AIRWATER_LINEAR(l_g%np,1,1,l_txc(1,1),l_q(1,4))
        
        l_q(:,5) = l_q(:,4) ! l_q(:,6) for bil_cloud_4 is set =0 in dns_main at initialization
        
     ENDIF
     
  ENDIF
  
  RETURN
END SUBROUTINE PARTICLE_RANDOM_POSITION

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE  PARTICLE_RANDOM_POSITION(l_q,l_hq,l_txc,l_tags,l_comm, txc, wrk1d,wrk2d,wrk3d)
  
  USE DNS_TYPES,  ONLY : pointers_dt
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL,      DIMENSION(isize_particle,*), TARGET :: l_q, l_hq, l_txc
  INTEGER(8), DIMENSION(isize_particle)           :: l_tags
  TREAL,      DIMENSION(isize_l_comm),     TARGET :: l_comm
  TREAL,      DIMENSION(isize_field,*),    TARGET :: txc
  TREAL,      DIMENSION(*)                        :: wrk1d,wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER  i, j, is
  TINTEGER  particle_number_local
  TINTEGER, ALLOCATABLE :: x_seed(:)
  TINTEGER size_seed
  
  TREAL xref,yref,zref, xscale,yscale,zscale, dy_loc, dummy, real_buffer_frac

  TREAL rnd_number(4), rnd_number_second
  TINTEGER rnd_scal(3)

  TINTEGER nvar, npar
  TYPE(pointers_dt), DIMENSION(inb_lag_total_interp) :: data, data_out
  TREAL, DIMENSION(:,:,:), POINTER :: txc_3d

!########################################################################
#ifdef USE_MPI
  particle_number_local = INT( particle_number /INT(ims_npro, KIND=8) )
  IF ( ims_pro .LT. INT( MOD(particle_number, INT(ims_npro, KIND=8)) ) ) THEN
     particle_number_local = particle_number_local +1
  ENDIF
  CALL MPI_ALLGATHER(particle_number_local,1,MPI_INTEGER4,ims_size_p,1,MPI_INTEGER4,MPI_COMM_WORLD,ims_err)

#else
  particle_number_local = INT(particle_number)

#endif
  
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
     
     DO i = 1,particle_number_local
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
     txc_3d(1:imax,1:jmax,1:kmax) => txc(1:isize_field,is)

     ! IF ( jmin_part /sbg(is)%ymean .GT. g(2)%size ) THEN
     !    CALL IO_WRITE_ASCII(efile,'PARTICLE_RANDOM_POSITION. JMIN_PART exceeds YCorrScalar value')
     !    CALL DNS_STOP(DNS_ERROR_PARTICLE)
     ! END IF

     dy_loc= g(2)%nodes(jmin_part+1) -g(2)%nodes(jmin_part)
     
     i = 1
     DO WHILE ( i .LE. particle_number_local ) 
        DO j=1,3
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO
        
        rnd_scal(1) = 1         +floor(rnd_number(1)*imax)
        rnd_scal(3) = 1         +floor(rnd_number(3)*kmax)
        rnd_scal(2) = jmin_part +floor(rnd_number(2)*(jmax_part-jmin_part+1))
        real_buffer_frac =             rnd_number(2)*(jmax_part-jmin_part+1) &
                         -       floor(rnd_number(2)*(jmax_part-jmin_part+1))
        
        dummy = ( txc_3d(rnd_scal(1),rnd_scal(2),rnd_scal(3)) -sbg(is)%mean )/sbg(is)%delta
        dummy = abs( dummy + C_05_R )

        CALL RANDOM_NUMBER(rnd_number_second)
        
        IF ( rnd_number_second .LE. dummy ) THEN
           
           l_q(i,1) = xref + rnd_number(1) *xscale
           l_q(i,3) = zref + rnd_number(3) *zscale
           l_q(i,2) = g(2)%nodes(rnd_scal(2)) + real_buffer_frac *dy_loc
        
           i=i+1

        END IF
        
     END DO

  ENDIF
  
!########################################################################
! Remaining scalar properties of the lagrangian field
!########################################################################
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN

     CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d)
     
     IF ( imixture .EQ.  MIXT_TYPE_AIRWATER_LINEAR ) THEN
        nvar = 0
        nvar = nvar+1; data(nvar)%field => txc(:,1); data_out(nvar)%field => l_txc(:,1)
        nvar = nvar+1; data(nvar)%field => txc(:,2); data_out(nvar)%field => l_txc(:,2)        
        CALL FIELD_TO_PARTICLE(nvar, data, npar, data_out, l_q,l_hq,l_tags,l_comm, wrk1d,wrk2d,wrk3d)
        
        CALL THERMO_AIRWATER_LINEAR(isize_particle,1,1,l_txc(1,1),l_q(1,4))
        
        l_q(:,5) = l_q(:,4) ! l_hq(:,6) for bil_cloud_4 is set =0 in dns_main at initialization
        
     ENDIF
     
  ENDIF
  
  RETURN
END SUBROUTINE PARTICLE_RANDOM_POSITION

!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# Generate random particle position 
!# Three options are included
!# 1) Fill hole domain
!# 2) Fill only certain part of y-domain -> y_particle_width, y_particle_pos
!# 3) Use scalar profile -> rnd_mode=2
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
!# COMMENTS
!# In case of needing memoy consider using l_txc instead of l_buffer
!#
!########################################################################

SUBROUTINE  PARTICLE_RANDOM_POSITION_OLD(l_q,l_tags,l_hq, txc, wrk1d,wrk2d,wrk3d)
  
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(isize_particle,inb_particle)   :: l_q, l_hq
  TREAL, DIMENSION(isize_particle)                :: l_tags
  TREAL, DIMENSION(*)                             :: wrk1d,wrk2d,wrk3d
  TREAL, DIMENSION(imax,jmax,kmax,*)              :: txc

  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE        :: l_buffer
  TINTEGER  i,j, is
  TINTEGER  particle_number_local
 ! TINTEGER, DIMENSION(1)::x_seed
  TINTEGER, ALLOCATABLE                            ::x_seed(:)
  TINTEGER size_seed

  TREAL real_buffer_calc, real_buffer_mean, real_buffer_delta, real_buffer_frac
  TREAL rnd_number(4) 
  TINTEGER rnd_scal(3)
  TREAL rnd_number_second

#ifdef USE_MPI
  TINTEGER particle_number_each
  TREAL real_buffer_i, real_nbuffer_i, real_buffer_k, real_nbuffer_k
  TREAL, DIMENSION(jmax) :: dummy
  TINTEGER, DIMENSION(:),   ALLOCATABLE :: ims_size_p_buffer
#endif

  is=1
  i=0
  l_q=C_0_R
  ALLOCATE(l_buffer(isize_particle,2))

  CALL RANDOM_SEED(SIZE=size_seed) !Alberto new
  ALLOCATE(x_seed(size_seed))

#ifdef USE_MPI
  !Particle number per processor
  particle_number_each = INT( particle_number /INT(ims_npro, KIND=8) )
  IF ( ims_pro .LT. INT( MOD(particle_number, INT(ims_npro, KIND=8)) ) ) THEN
     particle_number_each = particle_number_each +1
  ENDIF
  
  !Generate seed - different seed for each processor
  !x_seed=(/1+ims_pro/)
  x_seed = (/ (j,j=1+ims_pro, size_seed+ims_pro) /) !Alberto
  CALL RANDOM_SEED(PUT=x_seed)

!#######################################################################
! GENERATE RANDOM NUMBERS
!#######################################################################
  IF (particle_rnd_mode .eq. 1) THEN

    !Each processor writes in its chunk of the array
    !DO i=ims_pro*particle_number_each+1,(ims_pro+1)*particle_number_each
    DO i=1,particle_number_each
        !z-y-x direction
        DO j=1,3
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO

        !Division of two integers needs real buffer
        real_buffer_i=ims_pro_i
        real_nbuffer_i=ims_npro_i
        real_buffer_k=ims_pro_k
        real_nbuffer_k=ims_npro_k
        !Set up random numbers
        !Each processor creates random numbers in its territory
        l_q(i,1)=((rnd_number(1)/ims_npro_i)+(real_buffer_i/real_nbuffer_i))*g(1)%scale
        l_q(i,2)= rnd_number(2) *y_particle_width +y_particle_pos -C_05_R*y_particle_width
!        l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*g(2)%scale
!        l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part))*g(2)%size
        IF ( g(3)%size .GT. 1) THEN
           l_q(i,3)=((rnd_number(3)/ims_npro_k)+(real_buffer_k/real_nbuffer_k))*g(3)%scale
        ELSE
           l_q(i,3)=1
        END IF

    END DO
    !Set up the vector containing the number of particles for every processor
    ims_size_p(1:ims_npro)=0
    ims_size_p(ims_pro+1)=particle_number_each 

!#######################################################################
  ELSEIF (particle_rnd_mode .eq. 2) THEN

!Read the scalar field
     CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d)
     
     IF (jmin_part/sbg(is)%ymean .GT. g(2)%size) THEN
        CALL IO_WRITE_ASCII(efile,'DNS_INIPART. JMIN_PART exceeds YCorrScalar value')
        CALL DNS_STOP(DNS_ERROR_PARTICLE)
     END IF
     
     DO WHILE (i .le. particle_number_each) 
!z-y-x direction
        DO j=1,3
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO
        
      
        CALL RANDOM_NUMBER(rnd_number_second)
        
        rnd_scal(1)=floor(rnd_number(1)*imax)+1
        rnd_scal(2)=floor(rnd_number(2)*(jmax_part-jmin_part+1)+jmin_part)
        real_buffer_frac = rnd_number(2)*(jmax_part-jmin_part+1) - floor(rnd_number(2)*(jmax_part-jmin_part+1))
        rnd_scal(3)=floor(rnd_number(3)*kmax)+1
        
        real_buffer_mean=sbg(is)%mean
        real_buffer_delta=sbg(is)%delta
        real_buffer_calc=((real_buffer_mean/real_buffer_delta)-0.5)

!Use the scalar field to create the particle distribution
        IF (rnd_number_second .le. abs(txc(rnd_scal(1),rnd_scal(2),rnd_scal(3), 1)/real_buffer_delta-real_buffer_calc)) THEN
           
!Division of two integers needs real buffer
           real_buffer_i=ims_pro_i
           real_nbuffer_i=ims_npro_i
           real_buffer_k=ims_pro_k
           real_nbuffer_k=ims_npro_k
           
!Set up random numbers
!Each processor creates random numbers in its territory
           l_q(i,1)=((rnd_number(1)/ims_npro_i)+(real_buffer_i/real_nbuffer_i))*g(1)%scale
!        l_q(i,2)=g(2)%nodes(rnd_scal(2)) + real_buffer_frac*g(2)%scale/jmax 
           l_q(i,2)=g(2)%nodes(rnd_scal(2)) + real_buffer_frac*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part) )
           IF ( g(3)%size .GT. 1) THEN
              l_q(i,3)=((rnd_number(3)/ims_npro_k)+(real_buffer_k/real_nbuffer_k))*g(3)%scale
           ELSE
              l_q(i,3)=1
           END IF
           i=i+1
        ELSE
           i=i
      END IF
      
   END DO
!Set up the vector containing the number of particles for every processor
   ims_size_p(1:ims_npro)=0
   ims_size_p(ims_pro+1)=particle_number_each 
   
END IF

  ALLOCATE(ims_size_p_buffer(ims_npro))

  CALL MPI_ALLGATHER(ims_size_p(ims_pro+1),1,MPI_INTEGER4,ims_size_p_buffer,1,MPI_INTEGER4,MPI_COMM_WORLD,ims_err)
  ims_size_p(:)=ims_size_p_buffer(:)

  DEALLOCATE(ims_size_p_buffer)

#else

  IF (particle_rnd_mode .eq. 1) THEN
!  x_seed=1
     x_seed = (/ (j,j=1, size_seed) /) !Alberto
     CALL RANDOM_SEED(PUT=x_seed)
      
     DO i=1,particle_number
        DO j=1,3
           CALL RANDOM_NUMBER(rnd_number(j))
        END DO
        
        l_q(i,1) = rnd_number(1) *g(1)%scale
        l_q(i,2) = rnd_number(2) *y_particle_width +y_particle_pos -C_05_R*y_particle_width
!      l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*g(2)%scale
!      l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part))*g(2)%size
        IF ( g(3)%size .GT. 1) THEN
           l_q(i,3)=((rnd_number(3)*g(3)%scale))
        ELSE
           l_q(i,3)=1
        END IF
     END DO
     
  ELSEIF (particle_rnd_mode .eq. 2) THEN
    !Read the scalar field
    CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d)
  
    IF (jmin_part/sbg(is)%ymean .GT. g(2)%size) THEN
      CALL IO_WRITE_ASCII(efile,'DNS_INIPART. JMIN_PART exceeds YCorrScalar value')
      CALL DNS_STOP(DNS_ERROR_PARTICLE)
    END IF
    DO WHILE (i .le. particle_number) 
    !z-y-x direction
      DO j=1,3
        CALL RANDOM_NUMBER(rnd_number(j))
      END DO


      CALL RANDOM_NUMBER(rnd_number_second)

      rnd_scal(1)=floor(rnd_number(1)*imax)+1
      rnd_scal(2)=floor(rnd_number(2)*(jmax_part-jmin_part+1)+jmin_part)
      real_buffer_frac = rnd_number(2)*(jmax_part-jmin_part+1) - floor(rnd_number(2)*(jmax_part-jmin_part+1))
      rnd_scal(3)=floor(rnd_number(3)*kmax)+1

      real_buffer_mean=sbg(is)%mean
      real_buffer_delta=sbg(is)%delta
      real_buffer_calc=((real_buffer_mean/real_buffer_delta)-0.5)
      
!Use the scalar field to create the particle distribution
      IF (rnd_number_second .le. abs(txc(rnd_scal(1),rnd_scal(2),rnd_scal(3), 1)/real_buffer_delta-real_buffer_calc)) THEN

       !Set up random numbers
        l_q(i,1)=(rnd_number(1))*g(1)%scale
        l_q(i,2)=g(2)%nodes(rnd_scal(2)) + real_buffer_frac*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part) )
        IF ( g(3)%size .GT. 1) THEN
           l_q(i,3)=(rnd_number(3))*g(3)%scale
        ELSE
           l_q(i,3)=1
        END IF
        i=i+1
      ELSE
        i=i
      END IF

     END DO

  END IF
#endif

! ************************************************************
! Do the rest of the scalra properties of the lagrangian field
! ************************************************************
#ifdef USE_MPI
  particle_number_local=ims_size_p(ims_pro+1)
#else
  particle_number_local=particle_number
#endif

  IF (inb_particle .GT. 3 ) THEN
     IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
         CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_field, txc, wrk3d) ! Read the scalar fields into txc
          IF ( imixture .EQ.  MIXT_TYPE_AIRWATER_LINEAR) THEN 
             CALL FIELD_TO_PARTICLE_OLD (txc(1,1,1,1),wrk1d,wrk2d,wrk3d, l_buffer(1,1), l_tags, l_hq, l_q) ! Not sure about l_q(:,4). Maybe we need a particle txc 
             CALL FIELD_TO_PARTICLE_OLD (txc(1,1,1,2),wrk1d,wrk2d,wrk3d, l_buffer(1,2), l_tags, l_hq, l_q) ! Not sure about l_q(:,4). Maybe we need a particle txc 
       
             CALL THERMO_AIRWATER_LINEAR(isize_particle,1,1,l_buffer(1,1),l_q(1,4)) !INEFFICIENT. It loops the whole array and not only until particle_number_local

            l_q(:,5) = l_q(:,4)

            ! ATTTENTION l_hq(:,6) for bil_cloud_4 is set =0 in dns_main at initialization
          ELSE
             !Give error for wrong mixture
          ENDIF
     ENDIF
  ENDIF


  RETURN
END SUBROUTINE PARTICLE_RANDOM_POSITION_OLD

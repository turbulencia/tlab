#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

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

SUBROUTINE  PARTICLE_RANDOM_POSITION(l_q,l_hq,l_tags,isize_wrk3d,wrk1d,wrk2d,wrk3d,txc)
  
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
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE    :: s  
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE        :: l_buffer
  TREAL, DIMENSION(*)                             :: wrk1d,wrk2d,wrk3d
  TREAL, DIMENSION(isize_field,*)                 :: txc
  TINTEGER  i,j, is
  TINTEGER  particle_number_local
 ! TINTEGER, DIMENSION(1)::x_seed
  TINTEGER, ALLOCATABLE                            ::x_seed(:)
  TINTEGER size_seed

  TREAL real_buffer_calc, real_buffer_mean, real_buffer_delta, real_buffer_frac
  TREAL rnd_number(4) 
  TINTEGER rnd_scal(3)
  TREAL rnd_number_second

  TINTEGER isize_wrk3d
#ifdef USE_MPI
  TINTEGER particle_number_each
  TREAL real_buffer_i, real_nbuffer_i, real_buffer_k, real_nbuffer_k
  TREAL, DIMENSION(jmax) :: dummy
#endif

  is=1
  i=0
  l_q=C_0_R
  ALLOCATE(s(imax,jmax,kmax, inb_scal_array))
  ALLOCATE(l_buffer(isize_particle,2))

  CALL RANDOM_SEED(SIZE=size_seed) !Alberto new
  ALLOCATE(x_seed(size_seed))

#ifdef USE_MPI
  !Particle number per processor
  particle_number_each=int(particle_number/INT(ims_npro, KIND=8)) 

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
        l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*g(2)%scale
!        l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part))*jmax_total
        IF ( kmax_total .GT. 1) THEN
           l_q(i,3)=((rnd_number(3)/ims_npro_k)+(real_buffer_k/real_nbuffer_k))*g(3)%scale
        ELSE
           l_q(i,3)=1
        END IF

    END DO
    !Set up the vector containing the number of particles for every processor
    particle_vector(1:ims_npro)=C_0_R
    particle_vector(ims_pro+1)=particle_number_each 

  ELSEIF (particle_rnd_mode .eq. 2) THEN

    !Read the scalar field
    CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)

    IF (jmin_part/sbg(is)%ymean .GT. jmax_total) THEN
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
      IF (rnd_number_second .le. abs(s(rnd_scal(1),rnd_scal(2),rnd_scal(3), 1)/real_buffer_delta-real_buffer_calc)) THEN

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
        IF ( kmax_total .GT. 1) THEN
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
     particle_vector(1:ims_npro)=C_0_R
     particle_vector(ims_pro+1)=particle_number_each 

  END IF

#else

  IF (particle_rnd_mode .eq. 1) THEN
  !  x_seed=1
      x_seed = (/ (j,j=1, size_seed) /) !Alberto
    CALL RANDOM_SEED(PUT=x_seed)

    DO i=1,particle_number
      DO j=1,3
        CALL RANDOM_NUMBER(rnd_number(j))
      END DO

      l_q(i,1)=((rnd_number(1)))*g(1)%scale
      l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*g(2)%scale
!      l_q(i,2)=((rnd_number(2)*y_particle_width)+(y_particle_pos-(y_particle_width/2)))*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part))*jmax_total
      IF ( kmax_total .GT. 1) THEN
        l_q(i,3)=((rnd_number(3)*g(3)%scale))
      ELSE
        l_q(i,3)=1
      END IF
    END DO

  ELSEIF (particle_rnd_mode .eq. 2) THEN
    !Read the scalar field
    CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
  
    IF (jmin_part/sbg(is)%ymean .GT. jmax_total) THEN
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
      IF (rnd_number_second .le. abs(s(rnd_scal(1),rnd_scal(2),rnd_scal(3), 1)/real_buffer_delta-real_buffer_calc)) THEN

       !Set up random numbers
        l_q(i,1)=(rnd_number(1))*g(1)%scale
        l_q(i,2)=g(2)%nodes(rnd_scal(2)) + real_buffer_frac*(g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part) )
        IF ( kmax_total .GT. 1) THEN
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
particle_number_local=particle_vector(ims_pro+1)
#else
particle_number_local=particle_number
#endif

  IF (inb_particle .GT. 3 ) THEN
     IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
         CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, txc, wrk3d) ! Read the scalar fields into txc
          IF ( imixture .EQ.  MIXT_TYPE_AIRWATER_LINEAR) THEN 
             CALL FIELD_TO_PARTICLE (txc(:,1),wrk1d,wrk2d,wrk3d, l_buffer(1,1), l_tags, l_hq, l_q) ! Not sure about l_q(:,4). Maybe we need a particle txc 
             CALL FIELD_TO_PARTICLE (txc(:,2),wrk1d,wrk2d,wrk3d, l_buffer(1,2), l_tags, l_hq, l_q) ! Not sure about l_q(:,4). Maybe we need a particle txc 
       

             CALL THERMO_AIRWATER_LINEAR(isize_particle,1,1,l_buffer(1,1),l_q(1,4)) !INEFFICIENT. It loops the whole array and not only until particle_number_local


            l_q(:,5) = l_q(:,4)

            ! ATTTENTION l_hq(:,6) for bil_cloud_4 is set =0 in dns_main at initialization
          ELSE
             !Give error for wrong mixture
          ENDIF
     ENDIF
  ENDIF


  RETURN
END SUBROUTINE PARTICLE_RANDOM_POSITION

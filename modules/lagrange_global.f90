#include "types.h"

MODULE LAGRANGE_GLOBAL
  IMPLICIT NONE
  SAVE

  TYPE particle_dt
     SEQUENCE
     TINTEGER size ! size of arrays
     TINTEGER np   ! number of particles
     LOGICAL uniform
     INTEGER(8), DIMENSION(:), ALLOCATABLE :: tags
     TINTEGER,   DIMENSION(:), ALLOCATABLE :: nodes
  END TYPE particle_dt

  TINTEGER, PARAMETER :: MAX_LAGPARAM = 10 !Maximum size of Lagrange Parameters

! ###################################################################
! Lagrange Parameter
! ###################################################################
  TYPE(particle_dt) :: l_g

  TINTEGER      :: ilagrange
  TLONGINTEGER  :: particle_number_total
  TINTEGER      :: inb_particle_interp !Total number of interpolated fields into lagrangian
  TINTEGER      :: isize_l_comm, isize_pbuffer
  
  TINTEGER      :: particle_rnd_mode !which initializing mode
  TINTEGER      :: residence_reset  !if reseidence l_q should be reset

  TREAL         :: y_particle_pos  !position where particles will be initialized
  TREAL         :: y_particle_width  !width of particle distribution
  TREAL         :: l_y_lambda !y coordinate where approx radiation begins for residence times (set in dns_main)
  TREAL         :: l_y_base   !set to be 1/3 of cloud domain between two bouyancy stratification for residence times 

  TINTEGER      :: itrajectory       ! Type of trajectories
  TINTEGER      :: isize_trajectory  ! number of saved trajectories
  TINTEGER      :: inb_trajectory    ! number of properties saved along trajectories

  TINTEGER      :: icalc_part_pdf    ! if calculation of pdf for particles
  TREAL         :: y_particle_pdf_pos
  TREAL         :: y_particle_pdf_width
  TREAL         :: x_particle_pdf_pos
  TREAL         :: x_particle_pdf_width
  TREAL         :: z_particle_pdf_pos
  TREAL         :: z_particle_pdf_width
  TREAL         :: particle_pdf_max 
  TREAL         :: particle_pdf_interval 
  TINTEGER      :: number_of_bins
  
  TREAL         :: lagrange_param(MAX_LAGPARAM)                 ! lagrange function parameters
  CHARACTER*32, DIMENSION(15) :: LAGRANGE_SPNAME             !Name of different lagrange species

END MODULE LAGRANGE_GLOBAL

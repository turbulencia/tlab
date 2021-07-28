#include "types.h"

MODULE LAGRANGE_GLOBAL
  IMPLICIT NONE
  SAVE

  TYPE particle_dt
     SEQUENCE
     TINTEGER size ! size of arrays
     TINTEGER np   ! number of particles
     LOGICAL uniform, lpadding(3)
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
  TREAL         :: y_particle_pos  !position where particles will be initialized
  TREAL         :: y_particle_width  !width of particle distribution

  TINTEGER      :: residence_reset  !if reseidence l_q should be reset
  TREAL         :: l_y_lambda !y coordinate where approx radiation begins for residence times (set in dns_main)
  TREAL         :: l_y_base   !set to be 1/3 of cloud domain between two bouyancy stratification for residence times

  TINTEGER      :: itrajectory       ! Type of trajectories
  TINTEGER      :: isize_trajectory  ! number of saved trajectories
  TINTEGER      :: inb_trajectory    ! number of properties saved along trajectories

  TINTEGER      :: icalc_part_pdf    ! if calculation of pdf for particles
  TREAL         :: particle_pdf_subdomain(6)
  TREAL         :: particle_pdf_max
  TREAL         :: particle_pdf_interval

  TREAL         :: lagrange_param(MAX_LAGPARAM)                 ! lagrange function parameters
  CHARACTER*32, DIMENSION(15) :: LAGRANGE_SPNAME             !Name of different lagrange species

END MODULE LAGRANGE_GLOBAL

MODULE LAGRANGE_ARRAYS
  IMPLICIT NONE
  SAVE
  PRIVATE

  TREAL, ALLOCATABLE, PUBLIC :: l_q(:,:)      ! Lagrangian fields, flow vartiables
  TREAL, ALLOCATABLE, PUBLIC :: l_txc(:,:)    ! Temporary space for Lagrnagian fields

END MODULE

#include "types.h"

MODULE LAGRANGE_GLOBAL
  IMPLICIT NONE
  SAVE

   TINTEGER, PARAMETER :: MAX_LAGPARAM = 10 !Maximum size of Lagrange Parameters
! ###################################################################
! Lagrange Parameter
! ###################################################################
  TINTEGER      :: ilagrange       !Type of particle
  TLONGINTEGER  :: particle_number  !particle number parameter
  TINTEGER      :: particle_rnd_mode !which initializing mode
  TINTEGER      :: icalc_trajectories  !if calculation of trajectories
  TINTEGER      :: residence_reset  !if reseidence l_q should be reset
  TINTEGER      :: nzone_max  !maximum size of vector p_buffer in particle_send_recv
  TINTEGER      :: isize_hf_1, isize_hf_2, isize_hf_3
  TINTEGER      :: isize_max_hf, isize_pbuffer
  TINTEGER      :: isize_l_comm
  TINTEGER      :: jmax_part, jmin_part
  TINTEGER      :: inb_lag_aux_field ! Number of fields needed by lagrangian (no extra memory, usually txc fields)
  TINTEGER      :: inb_scal_particle ! Number of scalar properties solved in the lagrangian
  TINTEGER      :: inb_lag_total_interp !Total number of interpolated fields into lagrangian (=inb_lag_aux_fields +3)
  TINTEGER      :: inb_particle_evolution !inb_particle  number for time runge kutta
  TINTEGER      :: inb_particle_aux   ! additional inb_particle property which is not looped in time runge kutta
  TREAL         :: y_particle_pos  !position where particles will be initialized
  TREAL         :: y_particle_width  !width of particle distribution
  TREAL         :: particle_bumper !bumper zone in particle for send/recv
  TREAL         :: l_y_lambda !y coordinate where approx radiation begins for residence times (set in dns_main)
  TREAL         :: l_y_base   !set to be 1/3 of cloud domain between two bouyancy stratification for residence times 


  TINTEGER      :: icalc_particle_pdf  !if calculation of pdf for particles
  TINTEGER      :: num_trajectories  !number of followed trajectories
  TINTEGER      :: num_dispersion  !number of pairs for dispersion comparison
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
  CHARACTER*16, DIMENSION(15) :: LAGRANGE_SPNAME             !Name of different lagrange species

END MODULE LAGRANGE_GLOBAL

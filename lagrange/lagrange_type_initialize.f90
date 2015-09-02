#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 
!# 2014/17/06 - L. Muessle  
!#
!# 
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE LAGRANGE_TYPE_INITIALIZE
  
  USE LAGRANGE_GLOBAL

  USE DNS_GLOBAL, ONLY : inb_particle, inb_particle_txc
  IMPLICIT NONE

! -------------------------------------------------------------------


    inb_particle_txc = 0

   IF   (ilagrange .EQ. LAG_TYPE_TRACER) THEN
    inb_particle_evolution = 3
    inb_particle_aux = 0          
    inb_particle_txc = 0
    inb_lag_aux_field = 0
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent

  ELSEIF  (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN
    inb_particle_evolution = 3
    inb_particle_aux = 0          
    inb_particle_txc = 0
    inb_lag_aux_field = 0
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent

  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD) THEN
    inb_particle_evolution = 5 
    inb_particle_aux = 0          
    inb_particle_txc = 1
    inb_lag_aux_field = 3 
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
    LAGRANGE_SPNAME(1) = 'droplet_diff'
    LAGRANGE_SPNAME(2) = 'droplet_nodiff'
 
  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2) THEN
    inb_particle_evolution = 5 
    inb_particle_aux = 0          
    inb_particle_txc = 1
    inb_lag_aux_field = 4 
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
    LAGRANGE_SPNAME(1) = 'droplet_diff_2'
    LAGRANGE_SPNAME(2) = 'droplet_nodiff_2'
  
  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3) THEN
    inb_particle_evolution = 5    !amount of particle properties 
    inb_particle_aux = 0          !amount of particle properties without runge kutta (only sent and sorted)
    inb_particle_txc = 1          !l_txc properties
    inb_lag_aux_field = 4         !field data on txc
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
    LAGRANGE_SPNAME(1) = 'droplet_diff_3'
    LAGRANGE_SPNAME(2) = 'droplet_nodiff_3'

  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
    inb_particle_evolution = 5    !amount of particle properties with runge kutta
    inb_particle_aux = 1          !amount of particle properties without runge kutta (only sent and sorted)
    inb_particle_txc = 1          !l_txc properties
    inb_lag_aux_field = 4         !field data on txc
    inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
    LAGRANGE_SPNAME(1) = 'droplet_diff_3'
    LAGRANGE_SPNAME(2) = 'droplet_nodiff_3'
    LAGRANGE_SPNAME(3) = 'residence_part'
    

  END IF
  
  inb_scal_particle = inb_particle_evolution - 3          !Number of scalar properties solved in the lagrangian
  inb_lag_total_interp = inb_lag_aux_field + 3  !Number of fields needed by lagrangian (no extra memory, usually txc fields)

  RETURN
END SUBROUTINE LAGRANGE_TYPE_INITIALIZE

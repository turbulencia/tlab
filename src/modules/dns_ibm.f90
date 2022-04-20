#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF MODLE
!#   IBM based on:
!#     "A simple and scalable immersed boundary method for high-ﬁdelity 
!#      simulations of ﬁxed and moving objects on a Cartesian mesh"
!#      https://doi.org/10.1016/j.apm.2021.06.026
!#                    
!#
!########################################################################

module DNS_IBM

  use TLAB_TYPES,     only: ibm_geo_dt
  use TLAB_CONSTANTS, only: MAX_NSP 

  implicit none

  save 

  ! descriptive geometry fields (saved)
  TREAL, dimension(:),allocatable, target :: eps, epsp                   ! eps indicator field
  TINTEGER, dimension(:),     allocatable :: nobi,    nobj,   nobk       ! number of objects in i/j/k 
  TINTEGER, dimension(:),     allocatable :: nobi_b,  nobj_b, nobk_b     ! beginn of objects in i/j/k 
  TINTEGER, dimension(:),     allocatable :: nobi_e,  nobj_e, nobk_e     ! end    of objects in i/j/k
  TINTEGER                                :: nobi_max, nobj_max, nobk_max, nob_max
  TINTEGER                                :: max_height_objlo, max_height_objup

  ! modified field
  TREAL, dimension(:), allocatable, target:: fld_ibm                     ! with splines in solid regions

  ! boundary values of scalar fields
  TREAL, dimension(MAX_NSP)               :: ibmscaljmin, ibmscaljmax 

  ! work array for splines
  TREAL,    dimension(:),     allocatable :: xa, xb, ya, yb

  ! flag (decides which fdm calls are with modified fields)
  logical                                 :: ibm_burgers, ibm_partial, ibm_objup

  ! read_local
  logical                                 :: ibm_restart
  TINTEGER                                :: nflu                        ! number of fluid points used for Splines 
  TINTEGER                                :: ibm_io                      ! IBM IO Type 
  
  ! array sizes
  TINTEGER                                :: isize_nobi,    isize_nobj,    isize_nobk
  TINTEGER                                :: isize_nobi_be, isize_nobj_be, isize_nobk_be
  TINTEGER                                :: nspl
  TINTEGER                                :: isize_wrk1d_ibm

  ! check IBM procs (active/idle)
  logical                                 :: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

  ! ibm_dt geometry type 
  type(ibm_geo_dt)                        :: xbars_geo                   ! create new geometry here

  ! name of io eps
  character(len=32), parameter            :: eps_name      = 'eps0.1'
  character(len=32), parameter            :: eps_name_real = eps_name(1:4)

end module DNS_IBM

!########################################################################
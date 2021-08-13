#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF MODLE
!#
!#
!#
!#                    
!#
!########################################################################

module DNS_IBM

  use TLAB_TYPES, only: ibm_geo_dt

  implicit none

  save 

  ! descriptive geometry fields (saved)
  TREAL,    dimension(:),     allocatable :: eps                         ! eps indicator field
  TINTEGER, dimension(:),     allocatable :: nobi,    nobj,   nobk       ! number of objects in i/j/k 
  TINTEGER, dimension(:),     allocatable :: nobi_b,  nobj_b, nobk_b     ! beginn of objects in i/j/k 
  TINTEGER, dimension(:),     allocatable :: nobi_e,  nobj_e, nobk_e     ! end    of objects in i/j/k
  
  ! descriptive geometry fields (deallocated after initialization of geometry)
  TREAL,    dimension(:,:,:), allocatable :: eps_aux                     ! eps_aux field 
                                                                         ! (debugging / geometry generation)
  TREAL,    dimension(:),     allocatable :: epsi, epsj, epsk            ! eps transposed in i/j/k

  ! modified field
  TREAL, dimension(:), allocatable, target:: fld_ibm                     ! with splines in solid regions

  ! work array for splines
  TREAL,    dimension(:),     allocatable :: wrk_ibm
  TINTEGER, dimension(:),     allocatable :: iwrk_ibm
  !
  TREAL,    dimension(:),     allocatable :: xa, xb, ya, yb
  logical,  dimension(:),     allocatable :: x_mask, y_mask

  ! flag - decides which fdm calls are with modiefied field or not (opr_burgers.f90)
  logical                                 :: ibm_burgers, ibm_partial

  ! read_local from dns.ini file 
  logical                                 :: ibm_spline_global, ibm_procs_idle
  TINTEGER                                :: kspl                        ! spline order kspl=[1,5] (best: 3 or 5)
  TINTEGER                                :: nflu                        ! number of fluid points used for Splines 
                                                                         ! (on one side) nflu >= kspl
  
  ! array sizes
  TINTEGER                                :: isize_nobi,    isize_nobj,    isize_nobk
  TINTEGER                                :: isize_nobi_be, isize_nobj_be, isize_nobk_be
  TINTEGER                                :: nsp, nest
  TINTEGER                                :: isize_wrk_ibm, isize_iwrk_ibm, isize_wrk1d_ibm

  ! check IBM procs (active/idle)
  logical                                 :: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

  ! ibm_dt type 
  type(ibm_geo_dt)                        :: xbars_geo                   ! create new here for new geometries

end module DNS_IBM

!########################################################################
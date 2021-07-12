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

  ! use DNS_TYPES, only: ibm_dt

  implicit none

  save 

  ! geometry fields
  TREAL, dimension(:,:,:), allocatable :: eps_aux                     ! eps_aux field (debugging / geometry generation)
  TREAL, dimension(:),     allocatable :: epsi,    epsj,   epsk,  eps ! eps    transposed in i/j/k
  TREAL, dimension(:),     allocatable :: nobi,    nobj,   nobk       ! number of objects in i/j/k 
  TREAL, dimension(:),     allocatable :: nobi_b,  nobj_b, nobk_b     ! beginn of objects in i/j/k 
  TREAL, dimension(:),     allocatable :: nobi_e,  nobj_e, nobk_e     ! end    of objects in i/j/k

  ! modified velocity fields
  TREAL, dimension(:),     allocatable :: u_ibm                       ! with splines in solid regions

  ! used in opr_burgers.f90
  logical :: ibm_burgers 

  ! informations of type of immersed objects (--> introduce ibm_type in modules/dns_types)
  TINTEGER, dimension(3)               :: xbars_geo                   ! bars in x, xbars_geo(3)=[nbars,hbar,wbar]

  ! ibm_dt type (--> introduce ibm_type in modules/dns_types)
  ! type(ibm_dt) :: xbars

end module DNS_IBM

!########################################################################
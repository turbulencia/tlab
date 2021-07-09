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

  implicit none

  save 

  TREAL, dimension(:,:,:), allocatable :: eps_aux                     ! eps_aux field (debugging / geometry generation)
  TREAL, dimension(:),     allocatable :: epsi,    epsj,   epsk,  eps ! eps    transposed in i/j/k
  TREAL, dimension(:),     allocatable :: nobi,    nobj,   nobk       ! number of objects in i/j/k 
  TREAL, dimension(:),     allocatable :: nobi_b,  nobj_b, nobk_b     ! beginn of objects in i/j/k 
  TREAL, dimension(:),     allocatable :: nobi_e,  nobj_e, nobk_e     ! end    of objects in i/j/k

end module DNS_IBM

!########################################################################
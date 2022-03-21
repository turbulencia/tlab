!########################################################################
!# HISTORY / AUTHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 in ibm_initialize()
!#    
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################

subroutine IBM_DEALLOCATE()

  use DNS_IBM, only : ibm_restart, eps_aux, epsi, epsj, epsk

  implicit none

  ! ================================================================== !

  ! deallocate here not needed arrays
  if ( .not. ibm_restart ) then
    deallocate(eps_aux)
  end if
  
  deallocate(epsi)
  
  deallocate(epsj)
  
  deallocate(epsk)

  return
end subroutine IBM_DEALLOCATE

!########################################################################
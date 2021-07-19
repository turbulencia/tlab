!########################################################################
!# HISTORY / ATHORS
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

  use DNS_IBM

  implicit none

  ! ================================================================== !

  ! deallocate here not needed arrays

  deallocate(eps_aux)
  
  deallocate(epsi)
  
  deallocate(epsj)
  
  deallocate(epsk)

  return
end subroutine IBM_DEALLOCATE

!########################################################################
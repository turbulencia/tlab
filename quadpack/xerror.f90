      subroutine xerror(messg,nmessg,nerr,level)
      
      IMPLICIT NONE

#include "types.h"
!
! abstract
!    xerror processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    (see subroutine xsetf for details.)
!
! description of parameters
!  --input--
!    messg - the hollerith message to be processed, containing
!            no more than 72 characters.
!    nmessg- the actual number of characters in messg.
!    nerr  - the error number associated with this message.
!            nerr must not be zero.
!    level - error category.
!            =2 means this is an unconditionally fatal error.
!            =1 means this is a recoverable error.  (i.e., it is
!               non-fatal if xsetf has been appropriately called.)
!            =0 means this is a warning message only.
!            =-1 means this is a warning message which is to be
!               printed at most once, regardless of how many
!               times this call is executed.
!
! examples
!    call xerror(23hsmooth -- num was zero.,23,1,2)
!    call xerror(43hinteg  -- less than full accuracy achieved.,
!                43,2,1)
!    call xerror(65hrooter -- actual zero of f found before interval
!1 fully collapsed.,65,3,0)
!    call xerror(39hexp    -- underflows being set to zero.,39,1,-1)
!
! written by ron jones, with slatec common math library subcommittee
! latest revision ---  7 feb 1979
!
      TINTEGER messg, nmessg, nerr, level

      dimension messg(nmessg)
!  call xerrwv(messg,nmessg,nerr,level,0,0,0,0,C_0_R,C_0_R)
      print*,'error reported to xerror.f'
      stop

!  return
      end

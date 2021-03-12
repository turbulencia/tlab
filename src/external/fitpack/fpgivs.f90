      subroutine fpgivs(piv,ww,cos,sin)
!  subroutine fpgivs calculates the parameters of a givens
!  transformation .
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL piv,ww,cos,sin
!  ..local scalars..
      TREAL dd,one,store
!  ..function references..
!  ..
      one = C_1_R
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end

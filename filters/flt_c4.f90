#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# 4th-Order Compact Filter from Lele [J. Comp. Phys., V103, 1992]
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################
SUBROUTINE FLT_C4(imax,jkmax, periodic, bcsimin, bcsimax, cxi, u, uf, txi)

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER,                     INTENT(IN)    :: imax,jkmax, bcsimin, bcsimax
  TREAL, DIMENSION(jkmax,imax), INTENT(IN)    :: u
  TREAL, DIMENSION(jkmax,imax), INTENT(OUT)   :: uf
  TREAL, DIMENSION(imax,6),     INTENT(IN)    :: cxi
  TREAL, DIMENSION(imax,5),     INTENT(INOUT) :: txi

! -----------------------------------------------------------------------
  TINTEGER jk, i
  TREAL vmult1, vmultmx

! ###################################################################
  vmult1 = C_1_R; vmultmx = C_1_R
  IF ( bcsimin .EQ. DNS_FILTER_BCS_ZERO ) vmult1  = C_0_R ! No filter at i=1
  IF ( bcsimax .EQ. DNS_FILTER_BCS_ZERO ) vmultmx = C_0_R ! No filter at i=imax

! #######################################################################
! Set up the left hand side and factor the matrix
! Constant alpha contained in array cxi(6)
! #######################################################################
  IF ( periodic ) THEN ! periodic
     txi(1,1) = cxi(1,6)
     txi(1,2) = C_1_R
     txi(1,3) = cxi(1,6)

     txi(imax,1) = cxi(imax,6)
     txi(imax,2) = C_1_R
     txi(imax,3) = cxi(imax,6)

  ELSE ! biased
     txi(1,1) = C_0_R
     txi(1,2) = C_1_R
     txi(1,3) = cxi(1,6)*vmult1

     txi(imax,1) = cxi(imax,6)*vmultmx
     txi(imax,2) = C_1_R
     txi(imax,3) = C_0_R

  ENDIF

  DO i = 2, imax-1
     txi(i,1) = cxi(i,6)
     txi(i,2) = C_1_R
     txi(i,3) = cxi(i,6)
  ENDDO

! #######################################################################
! Set up right hand side and DO forward/backward substitution
! #######################################################################
  IF ( periodic ) THEN ! periodic
     DO jk=1, jkmax
        uf(jk,1) = cxi(1,1)*u(jk,imax-1) + cxi(1,2)*u(jk,imax) + cxi(1,3)*u(jk,1) + &
             cxi(1,4)*u(jk,2) + cxi(1,5)*u(jk,3)

        uf(jk,2) = cxi(2,1)*u(jk,imax)   + cxi(2,2)*u(jk,1)    + cxi(2,3)*u(jk,2) + &
             cxi(2,4)*u(jk,3) + cxi(2,5)*u(jk,4)

        uf(jk,imax  ) = cxi(imax,1)  *u(jk,imax-2) + cxi(imax,2)*  u(jk,imax-1) &
             + cxi(imax,3)*u(jk,imax)     + cxi(imax,4)*u(jk,1)      + cxi(imax,5)*u(jk,2)

        uf(jk,imax-1) = cxi(imax-1,1)*u(jk,imax-3) + cxi(imax-1,2)*u(jk,imax-2) &
             + cxi(imax-1,3)*u(jk,imax-1) + cxi(imax-1,4)*u(jk,imax) + cxi(imax-1,5)*u(jk,1)
     ENDDO

  ELSE ! biased
     DO jk=1, jkmax
        uf(jk,1) = cxi(1,1)*u(jk,1) &
             + cxi(1,2)*u(jk,2) &
             + cxi(1,3)*u(jk,3) &
             + cxi(1,4)*u(jk,4) &
             + cxi(1,5)*u(jk,5)
        uf(jk,2) = cxi(2,1) *u(jk,1) &
             + cxi(2,2) *u(jk,2) &
             + cxi(2,3) *u(jk,3) &
             + cxi(2,4) *u(jk,4) &
             + cxi(2,5) *u(jk,5)

        uf(jk,imax-1) = cxi(imax-1,5)*u(jk,imax)&
             + cxi(imax-1,4)*u(jk,imax-1) &
             + cxi(imax-1,3)*u(jk,imax-2) &
             + cxi(imax-1,2)*u(jk,imax-3) &
             + cxi(imax-1,1)*u(jk,imax-4)
        uf(jk,imax) = cxi(imax,5)*u(jk,imax)&
             + cxi(imax,4)*u(jk,imax-1) &
             + cxi(imax,3)*u(jk,imax-2) &
             + cxi(imax,2)*u(jk,imax-3) &
             + cxi(imax,1)*u(jk,imax-4)
     ENDDO

     IF ( bcsimin .EQ. DNS_FILTER_BCS_ZERO ) THEN ! No filter at i=1
        DO jk=1, jkmax
           uf(jk,   1) = u(jk,   1)
        ENDDO
     ENDIF
     IF ( bcsimax .EQ. DNS_FILTER_BCS_ZERO ) THEN ! No filter at i=1
        DO jk=1, jkmax
           uf(jk,imax) = u(jk,imax) 
        ENDDO
     ENDIF

  ENDIF

  DO i=3,imax-2
     DO jk=1,jkmax
        uf(jk,i) = cxi(i,1)*u(jk,i-2)+ cxi(i,2)*u(jk,i-1) + cxi(i,3)*u(jk,i) + &
             cxi(i,4)*u(jk,i+1) + cxi(i,5)*u(jk,i+2) 
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLT_C4


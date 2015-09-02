#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2000/09/22 - J.P. Mellado
!#              Created
!# 2002/01/01 - J.P. Mellado
!#              Adding operations routines
!# 2013/01/10 - J.P. Mellado
!#              Adding general block in-place routine
!#
!########################################################################
!# DESCRIPTION
!#
!# Reduction in i of the matrices a1,a2... into the smaller b gathering
!# the i-planes given in p(np).
!#
!########################################################################
SUBROUTINE REDUCE(imax,jmax,kmax, a, np,p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np                           ! np is the number of sampled planes
  TINTEGER, DIMENSION(np)            , INTENT(IN)  :: p ! array with the i location of the planes
  TREAL,    DIMENSION(imax,jmax,kmax), INTENT(IN)  :: a ! input array (big)
  TREAL,    DIMENSION(np,  jmax,kmax), INTENT(OUT) :: b ! output array (small)

  TINTEGER i, j, k, n

  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,np
           i = p(n)
           b(n,j,k) = a(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_SUM(imax,jmax,kmax, a1,a2, np,p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np                           
  TINTEGER, DIMENSION(np)            , INTENT(IN)  :: p 
  TREAL,    DIMENSION(imax,jmax,kmax), INTENT(IN)  :: a1,a2
  TREAL,    DIMENSION(np,  jmax,kmax), INTENT(OUT) :: b

  TINTEGER i, j, k, n

  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,np
           i = p(n)
           b(n,j,k) = a1(i,j,k) + a2(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_SUM

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_SUB(imax,jmax,kmax, a1,a2, np,p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np                           
  TINTEGER, DIMENSION(np)            , INTENT(IN)  :: p 
  TREAL,    DIMENSION(imax,jmax,kmax), INTENT(IN)  :: a1,a2
  TREAL,    DIMENSION(np,  jmax,kmax), INTENT(OUT) :: b

  TINTEGER i, j, k, n

  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,np
           i = p(n)
           b(n,j,k) = a1(i,j,k) - a2(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_SUB

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_MUL(imax,jmax,kmax, a1, a2, np, p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np                           
  TINTEGER, DIMENSION(np)            , INTENT(IN)  :: p 
  TREAL,    DIMENSION(imax,jmax,kmax), INTENT(IN)  :: a1,a2
  TREAL,    DIMENSION(np,  jmax,kmax), INTENT(OUT) :: b

  TINTEGER i, j, k, n

  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,np
           i = p(n)
           b(n,j,k) = a1(i,j,k) * a2(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_MUL

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_DIV( imax, jmax, kmax, a1, a2, np, p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np                           
  TINTEGER, DIMENSION(np)            , INTENT(IN)  :: p 
  TREAL,    DIMENSION(imax,jmax,kmax), INTENT(IN)  :: a1,a2
  TREAL,    DIMENSION(np,  jmax,kmax), INTENT(OUT) :: b

  TINTEGER i, j, k, n

  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,np
           i = p(n)
           b(n,j,k) = a1(i,j,k) / a2(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_DIV

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax_dst,jmax_dst,kmax_dst, a, wrk3d)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, imax_dst,jmax_dst,kmax_dst
  TREAL, DIMENSION(imax*jmax*kmax),             INTENT(INOUT) :: a
  TREAL, DIMENSION(imax_dst*jmax_dst*kmax_dst), INTENT(INOUT) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER j,k, ijmax,ijmax_dst, ip,ip_dst

! -------------------------------------------------------------------
  ijmax     = imax*jmax
  ijmax_dst = imax_dst*jmax_dst

  DO k = 1,kmax_dst
     ip     = (k-1)*ijmax     + 1
     ip_dst = (k-1)*ijmax_dst + 1
     DO j = 1,jmax_dst
        wrk3d(ip_dst:ip_dst+imax_dst-1) = a(ip:ip+imax_dst-1)
        ip     = ip     + imax
        ip_dst = ip_dst + imax_dst
     ENDDO
  ENDDO

  ip = imax_dst*jmax_dst*kmax_dst 
  a(1:ip) = wrk3d(1:ip)

  RETURN
END SUBROUTINE REDUCE_BLOCK_INPLACE

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_Y_ALL(imax,jmax,kmax, nvar1,a1, nvar2,a2, np,p, b)

  IMPLICIT NONE

  TINTEGER imax,jmax,kmax, np, nvar1,nvar2                               ! np is the number of sampled planes
  TINTEGER, DIMENSION(np),                         INTENT(IN)  :: p      ! array with the j location of the planes
  TREAL,    DIMENSION(imax,jmax,kmax,       *   ), INTENT(IN)  :: a1, a2 ! input array (big)
  TREAL,    DIMENSION(imax,np  ,nvar1+nvar2,kmax), INTENT(OUT) :: b      ! output array (small)

  TINTEGER j, j_loc, ivar, k
  
  DO k = 1,kmax
     DO ivar = 1,nvar1
        DO j = 1,np
           j_loc = p(j)
           b(1:imax,j,ivar,k) = a1(1:imax,j_loc,k,ivar)
        ENDDO
     ENDDO
     
     DO ivar = 1+nvar1,nvar2+nvar1 ! if nvar is 0, then array a2 is not used
        DO j = 1,np
           j_loc = p(j)
           b(1:imax,j,ivar,k) = a2(1:imax,j_loc,k,ivar)
        ENDDO
     ENDDO
     
  ENDDO
  
! ip_o=1
! DO k=1,kmax
!    DO ivar = 1,inb_flow
!       DO j=1,npln_j !nsave_planes
!          j_loc = pln_j(j) !j_save(j)
!          ip_i = (k-1)*jmax*imax + (j_loc-1)*imax +1;
!          wrk3d(ip_o:ip_o+imax-1) = q(ip_i:ip_i+imax-1,ivar);  ip_o = ip_o + imax;
!       ENDDO
!    ENDDO
  
!    DO ivar = 1,inb_scal
!       DO j=1,npln_j !nsave_planes
!          j_loc = pln_j(j) !j_save(j)
!          ip_i = (k-1)*jmax*imax + (j_loc-1)*imax +1;
!          wrk3d(ip_o:ip_o+imax-1) = s(ip_i:ip_i+imax-1,ivar);  ip_o = ip_o + imax;
!       ENDDO
!    ENDDO
  
! ENDDO
  
  RETURN
END SUBROUTINE REDUCE_Y_ALL

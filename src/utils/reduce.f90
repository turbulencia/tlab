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
SUBROUTINE REDUCE_BLOCK_INPLACE(nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst, a, wrk1d)

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: a
  TREAL, DIMENSION(nx_dst),   INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER j,k, nxy,nxy_dst, ip,ip_dst

! -------------------------------------------------------------------
  nxy     = nx    *ny
  nxy_dst = nx_dst*ny_dst

  DO k = 1,nz_dst
     ip     = (k-1) *nxy    +(nz1-1) *nxy +(ny1-1)*nx +(nx1-1) +1
     ip_dst = (k-1) *nxy_dst                                   +1
     DO j = 1,ny_dst
        wrk1d(1:nx_dst)           = a(ip:ip+nx_dst-1)
        a(ip_dst:ip_dst+nx_dst-1) = wrk1d(1:nx_dst)

        ip     = ip     + nx
        ip_dst = ip_dst + nx_dst
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_BLOCK_INPLACE

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_BLOCK_INPLACE_INT1(nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst, a, wrk1d)

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst
  INTEGER(1), DIMENSION(nx*ny*nz), INTENT(INOUT) :: a
  INTEGER(1), DIMENSION(nx_dst),   INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER j,k, nxy,nxy_dst, ip,ip_dst

! -------------------------------------------------------------------
  nxy     = nx    *ny
  nxy_dst = nx_dst*ny_dst

  DO k = 1,nz_dst
     ip     = (k-1) *nxy    +(nz1-1) *nxy +(ny1-1)*nx +(nx1-1) +1
     ip_dst = (k-1) *nxy_dst                                   +1
     DO j = 1,ny_dst
        wrk1d(1:nx_dst)           = a(ip:ip+nx_dst-1)
        a(ip_dst:ip_dst+nx_dst-1) = wrk1d(1:nx_dst)

        ip     = ip     + nx
        ip_dst = ip_dst + nx_dst
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_BLOCK_INPLACE_INT1

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_BLOCK(nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst, a, wrk3d)

  IMPLICIT NONE

  TINTEGER,                               INTENT(IN)    :: nx,ny,nz, nx1,ny1,nz1, nx_dst,ny_dst,nz_dst
  TREAL, DIMENSION(nx*ny*nz),             INTENT(INOUT) :: a
  TREAL, DIMENSION(nx_dst*ny_dst*nz_dst), INTENT(INOUT) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER j,k, nxy,nxy_dst, ip,ip_dst

! -------------------------------------------------------------------
  nxy     = nx    *ny
  nxy_dst = nx_dst*ny_dst

  DO k = 1,nz_dst
     ip     = (k-1) *nxy    +(nz1-1) *nxy +(ny1-1)*nx +(nx1-1) +1
     ip_dst = (k-1) *nxy_dst                                   +1
     DO j = 1,ny_dst
        wrk3d(ip_dst:ip_dst+nx_dst-1) = a(ip:ip+nx_dst-1)

        ip     = ip     + nx
        ip_dst = ip_dst + nx_dst
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE REDUCE_BLOCK

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_Y_ALL(nx,ny,nz, nvar1,a1, nvar2,a2, aux, np,np_aux,p, b)

  IMPLICIT NONE

  TINTEGER nx,ny,nz, np,np_aux, nvar1,nvar2                           ! np is the number of sampled planes
                                                                      ! np_aux is one or zero, additional data
  TINTEGER, DIMENSION(np),                      INTENT(IN)  :: p      ! array with the j location of the planes
  TREAL,    DIMENSION(nx,ny,nz,            * ), INTENT(IN)  :: a1, a2 ! input array (big)
  TREAL,    DIMENSION(nx,   nz,nvar1+nvar2   ), INTENT(IN)  :: aux    ! additional data with different structure
  TREAL,    DIMENSION(nx,np   ,nvar1+nvar2,nz), INTENT(OUT) :: b      ! output array (small)

  TINTEGER j, j_loc, ivar, k

  DO k = 1,nz
     DO ivar = 1,nvar1
        DO j = 1,np-np_aux
           j_loc = p(j)
           b(1:nx,j,ivar,k) = a1(1:nx,j_loc,k,ivar)
        ENDDO
        IF ( np_aux .GT. 0 ) b(1:nx,j,ivar,k) = aux(1:nx,k,ivar            ) ! Additional data, if needed
     ENDDO

     DO ivar = 1,nvar2 ! if nvar is 0, then array a2 is not used
        DO j = 1,np-np_aux
           j_loc = p(j)
           b(1:nx,j,ivar+nvar1,k) = a2(1:nx,j_loc,k,ivar)
        ENDDO
        IF ( np_aux .GT. 0 ) b(1:nx,j,ivar+nvar1,k) = aux(1:nx,k,ivar+nvar1) ! Additional data, if needed
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE REDUCE_Y_ALL

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_Z_ALL(nx,ny,nz, nvar1,a1, nvar2,a2, np,p, b)

  IMPLICIT NONE

  TINTEGER nx,ny,nz, np, nvar1,nvar2                               ! np is the number of sampled planes
  TINTEGER, DIMENSION(np),                   INTENT(IN)  :: p      ! array with the j location of the planes
  TREAL,    DIMENSION(nx*ny,nz,         * ), INTENT(IN)  :: a1, a2 ! input array (big)
  TREAL,    DIMENSION(nx*ny,nvar1+nvar2,np), INTENT(OUT) :: b      ! output array (small)

  TINTEGER k, k_loc, ivar

  DO k = 1,np
     k_loc = p(k)
     DO ivar = 1,nvar1
        b(1:nx*ny,ivar,k) = a1(1:nx*ny,k_loc,ivar)
     ENDDO

     DO ivar = 1,nvar2 ! if nvar is 0, then array a2 is not used
        b(1:nx*ny,ivar+nvar1,k) = a2(1:nx*ny,k_loc,ivar)
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE REDUCE_Z_ALL

! #######################################################################
! #######################################################################
SUBROUTINE REDUCE_X_ALL(nx,ny,nz, nvar1,a1, nvar2,a2, np,p, b)

  IMPLICIT NONE

  TINTEGER nx,ny,nz, np, nvar1,nvar2                               ! np is the number of sampled planes
  TINTEGER, DIMENSION(np),                   INTENT(IN)  :: p      ! array with the j location of the planes
  TREAL,    DIMENSION(nx,ny,nz,         * ), INTENT(IN)  :: a1, a2 ! input array (big)
  TREAL,    DIMENSION(np,ny,nvar1+nvar2,nz), INTENT(OUT) :: b      ! output array (small)

  TINTEGER i, i_loc, j, ivar, k

  DO k = 1,nz
     DO ivar = 1,nvar1
        DO j = 1,ny
           DO i = 1,np
              i_loc = p(i)
              b(i,j,ivar,k) = a1(i_loc,j,k,ivar)
           ENDDO
        ENDDO
     ENDDO

     DO ivar = 1,nvar2 ! if nvar is 0, then array a2 is not used
        DO j = 1,ny
           DO i = 1,np
              i_loc = p(i)
              b(i,j,ivar+nvar1,k) = a2(i_loc,j,k,ivar)
           ENDDO
        ENDDO
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE REDUCE_X_ALL

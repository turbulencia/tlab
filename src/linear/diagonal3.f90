#include "types.h"

#ifdef SINGLE_PREC
#define AXPY_LOC SAXPY
#define SCAL_LOC SSCAL
#else
#define AXPY_LOC DAXPY
#define SCAL_LOC DSCAL
#endif

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2009/09/01 - J.P. Mellado
!#              Reviewed, OpenMP
!# 2011/11/01 - C. Ansorge
!#              OpenMP Optimization
!#
!########################################################################
!# DESCRIPTION
!#
!# Vectorized tridiagonal solver based ob Thomas algorithm. 
!#
!########################################################################

! #######################################################################
! LU factorization stage; L with diagonal unity
! #######################################################################
SUBROUTINE TRIDFS(nmax, a,b,c)

  IMPLICIT NONE
  
  TINTEGER, INTENT(IN)                  ::  nmax
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: a,b,c
  
! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
  DO n = 2,nmax
     a(n) = a(n)/b(n-1)
     b(n) = b(n) - a(n)*c(n-1)
  ENDDO
  
! Final operations
  a =-a
  b = C_1_R/b
  c =-c

  RETURN
END SUBROUTINE TRIDFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE TRIDSS(nmax,len, a,b,c, f)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif 

  IMPLICIT NONE

  TINTEGER, INTENT(IN)                      :: nmax  ! dimension of tridiagonal systems
  TINTEGER, INTENT(IN)                      :: len   ! number of systems to be solved
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c ! factored LHS
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f     ! RHS and solution

! -------------------------------------------------------------------
  TINTEGER :: n
  TINTEGER :: srt,end,siz
  TREAL    :: dummy1, dummy2

#ifdef USE_BLAS
  TREAL alpha
  INTEGER ilen
#else
  TINTEGER l
#endif

! ###################################################################
! -----------------------------------------------------------------------
! Forward sweep
! -----------------------------------------------------------------------
#ifdef USE_BLAS
  ilen = len
#endif

#ifdef USE_BLAS
!$omp parallel default(none) &
!$omp private(n,ilen,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,nmax,len)
#else
!$omp parallel default(none) &
!$omp private(n,l,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,nmax,len)
#endif

  CALL DNS_OMP_PARTITION(len,srt,end,siz)
  IF ( siz .LE. 0 ) THEN 
     GOTO 999
  END IF

#ifdef USE_BLAS
  ilen = siz
#endif 

  DO n = 2,nmax
     dummy1 = a(n)
#ifdef USE_BLAS
     CALL AXPY_LOC(ilen, dummy1, f(srt,n-1), 1, f(srt,n), 1)
#else 
     DO l=srt,end
        f(l,n) = f(l,n) + dummy1*f(l,n-1)
     ENDDO
#endif
  ENDDO
  
! -----------------------------------------------------------------------
! Backward sweep
! -----------------------------------------------------------------------
  dummy1=b(nmax) 
#ifdef USE_BLAS
  CALL SCAL_LOC(ilen, dummy1, f(srt,nmax), 1)
#else 
  DO l=srt,end 
     f(l,nmax) = f(l,nmax)*dummy1
  ENDDO
#endif

  DO n = nmax-1, 1, -1
     dummy1 = c(n)
     dummy2 = b(n)
#ifdef USE_BLAS
     CALL AXPY_LOC(ilen, dummy1, f(srt,n+1), 1, f(srt,n), 1)
     CALL SCAL_LOC(ilen, dummy2, f(srt,n), 1)
#else 
     DO l = srt, end
        f(l,n) = (f(l,n) + dummy1*f(l,n+1))*dummy2
     ENDDO
     
#endif
  ENDDO
999 CONTINUE
!$omp end parallel
 
  RETURN
END SUBROUTINE TRIDSS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE TRIDSS_ADD(nmax,len, a,b,c, f, g,h, d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif 

  IMPLICIT NONE

  TINTEGER, INTENT(IN)                      :: nmax  ! dimension of tridiagonal systems
  TINTEGER, INTENT(IN)                      :: len   ! number of systems to be solved
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c ! factored LHS
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f     ! RHS and solution
  TREAL, DIMENSION(len,nmax), INTENT(IN)    :: g,h   ! Additive term
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: d     ! Scaled Jacobian

! -------------------------------------------------------------------
  TINTEGER :: n
  TINTEGER :: srt,end,siz
  TREAL    :: dummy1, dummy2

  TINTEGER l
#ifdef USE_BLAS
  TREAL alpha
  INTEGER ilen
#endif

! ###################################################################
! -----------------------------------------------------------------------
! Forward sweep
! -----------------------------------------------------------------------
#ifdef USE_BLAS
  ilen = len
#endif


#ifdef USE_BLAS
!$omp parallel default(none) &
!$omp private(n,l,ilen,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,d,g,h,nmax,len) 
#else
!$omp parallel default(none) &
!$omp private(n,l,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,d,g,h,nmax,len) 
#endif

  CALL DNS_OMP_PARTITION(len,srt,end,siz)
  IF ( siz .LE. 0 ) THEN 
     GOTO 999
  END IF

#ifdef USE_BLAS
  ilen = siz
#endif 

  DO n = 2,nmax
     dummy1 = a(n)
#ifdef USE_BLAS
     CALL AXPY_LOC(ilen, dummy1, f(srt,n-1), 1, f(srt,n), 1)
#else 
     DO l=srt,end
        f(l,n) = f(l,n) + dummy1*f(l,n-1)
     ENDDO
#endif
  ENDDO
  
! -----------------------------------------------------------------------
! Backward sweep
! -----------------------------------------------------------------------
  dummy1=b(nmax) 
#ifdef USE_BLAS
  CALL SCAL_LOC(ilen, dummy1, f(srt,nmax), 1)
#else 
  DO l=srt,end 
     f(l,nmax) = f(l,nmax)*dummy1
  ENDDO
#endif

  DO n = nmax-1, 1, -1
     dummy1 = c(n)
     dummy2 = b(n)
#ifdef USE_BLAS
     CALL AXPY_LOC(ilen, dummy1, f(srt,n+1), 1, f(srt,n), 1)
     CALL SCAL_LOC(ilen, dummy2, f(srt,n), 1)
     DO l = srt, end
        f(l,n+1) =  f(l,n+1) - (d(n+1)+g(l,n+1)) *h(l,n+1)  ! TO BE IMPLEMENTED USING BLAS
     ENDDO
#else 
     DO l = srt, end
        f(l,n  ) = (f(l,n) + dummy1*f(l,n+1))*dummy2
        f(l,n+1) =  f(l,n+1) - (d(n+1)+g(l,n+1)) *h(l,n+1)
     ENDDO
#endif
  ENDDO
  DO l = srt, end
     f(l,1) =  f(l,1) - (d(1)+g(l,1)) *h(l,1)
  ENDDO
  
999 CONTINUE
!$omp end parallel
 
  RETURN
END SUBROUTINE TRIDSS_ADD

!########################################################################
!# DESCRIPTION
!#
!# Circulant system (periodic bcs)
!#
!########################################################################

! #######################################################################
! LU factorization stage
! #######################################################################
SUBROUTINE TRIDPFS(nmax, a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER nmax
  TREAL, DIMENSION(nmax) ::  a,b,c,d,e

! -------------------------------------------------------------------
  TINTEGER n
  TREAL sum

! ###################################################################
! Generate first elements of LU
  c(1) = c(1)/b(1)
  e(1) = a(1)/b(1)
  d(1) = c(nmax)

! Generate n=2 to n=n-2 elements of LU
  DO n = 2, nmax-2
     b(n) = b(n) - a(n)*c(n-1)
     c(n) = c(n)/b(n)               
     e(n) =-a(n)*e(n-1)/b(n)
     d(n) =-d(n-1)*c(n-1)
  ENDDO


! Generate n-1 elements
  b(nmax-1) =  b(nmax-1) - a(nmax-1)*c(nmax-2)
  e(nmax-1) = (c(nmax-1) - a(nmax-1)*e(nmax-2))/b(nmax-1)
  d(nmax-1) =  a(nmax  ) - d(nmax-2)*c(nmax-2)

! Generate the n-th element
  sum = C_0_R
  DO n = 1,nmax-1
     sum = sum + d(n)*e(n)
  ENDDO
  b(nmax) = b(nmax) - sum

! Final operations
  DO n = 1,nmax
     b(n) = C_1_R/b(n)
     a(n) =-a(n)*b(n)
     c(n) =-c(n)
     e(n) =-e(n)
  ENDDO  

  RETURN
END SUBROUTINE TRIDPFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE TRIDPSS(nmax,len, a,b,c,d,e, f, wrk)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  
  IMPLICIT NONE
  
  TINTEGER,                   INTENT(IN)    :: nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c,d,e
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f
  TREAL, DIMENSION(len)                     :: wrk

! -------------------------------------------------------------------
  TREAL                                     :: dummy1, dummy2
  TINTEGER                                  :: srt, end, siz

  TINTEGER l, n
#ifdef USE_BLAS 
  INTEGER :: ilen
#endif 

! -------------------------------------------------------------------
! Forward sweep
! -------------------------------------------------------------------

#ifdef USE_BLAS
!$omp parallel default( none ) &
!$omp private(n, l,ilen, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,wrk,nmax,a,b,c,d,e,len)
#else
!$omp parallel default( none ) &
!$omp private(n, l, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,wrk,nmax,a,b,c,d,e,len)
#endif

  CALL DNS_OMP_PARTITION(len,srt,end,siz) 
  IF ( siz .LE. 0 ) THEN 
     GOTO 999
  ENDIF

#ifdef USE_BLAS 
  ilen = siz
#endif 

  dummy1 = b(1)
#ifdef USE_BLAS 
  CALL SCAL_LOC(ilen,dummy1,f(srt,1),1)
#else
  DO l=srt,end
     f(l,1) = f(l,1)*dummy1
  ENDDO
#endif 

  DO n = 2, nmax-1
     dummy1 = a(n)
     dummy2 = b(n)
#ifdef USE_BLAS 
     CALL SCAL_LOC(ilen,dummy2,f(srt,n),1) 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n-1),1, f(srt,n), 1)
#else
     DO l=srt,end
        f(l,  n) = f(l,  n)*dummy2 + dummy1*f(l,  n-1)
     ENDDO
#endif
  ENDDO
      
  wrk(srt:end) = C_0_R
  
  DO n = 1, nmax-1
     dummy1 = d(n)
#ifdef USE_BLAS 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n),1,wrk(srt),1)
#else
     DO l=srt,end
        wrk(l) = wrk(l)   + dummy1*f(l,n)
     ENDDO
#endif
  ENDDO

  dummy1 = b(nmax)
#ifdef USE_BLAS 
  CALL SCAL_LOC(ilen, dummy1,f(srt,nmax),1)
  CALL AXPY_LOC(ilen,-dummy1,wrk(srt),1,f(srt,nmax),1)
#else
  DO l=srt,end
     f(l,nmax) = (f(l,nmax) - wrk(l))*dummy1
  ENDDO
#endif

! -------------------------------------------------------------------
! Backward sweep
! -------------------------------------------------------------------
  dummy1 = e(nmax-1)

#ifdef USE_BLAS
  CALL AXPY_LOC(ilen,dummy1,f(srt,nmax),1,f(srt,nmax-1),1)
#else
  DO l=srt,end
     f(l,nmax-1) = dummy1*f(l,nmax) + f(l,nmax-1)
  ENDDO
#endif 


  DO n = nmax-2, 1, -1
     dummy1 = c(n)
     dummy2 = e(n)
#ifdef USE_BLAS 
     CALL AXPY_LOC(ilen,dummy2,f(srt,nmax),1,f(srt,n),1) 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n+1), 1,f(srt,n),1)
#else
     DO l = srt,end
        f(l,n) = f(l,n) + dummy1*f(l,n+1) + dummy2*f(l,nmax)
     ENDDO
#endif
  ENDDO
999 CONTINUE
!$omp end parallel
  
  RETURN
END SUBROUTINE TRIDPSS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE TRIDPSS_ADD(nmax,len, a,b,c,d,e, f, g,h, wrk)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  
  IMPLICIT NONE
  
  TINTEGER,                   INTENT(IN)    :: nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c,d,e
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f
  TREAL, DIMENSION(len,nmax), INTENT(IN)    :: g,h
  TREAL, DIMENSION(len)                     :: wrk

! -------------------------------------------------------------------
  TREAL                                     :: dummy1, dummy2
  TINTEGER                                  :: srt, end, siz

  TINTEGER l, n
#ifdef USE_BLAS 
  INTEGER :: ilen
#endif 

! -------------------------------------------------------------------
! Forward sweep
! -------------------------------------------------------------------

#ifdef USE_BLAS
!$omp parallel default( none ) &
!$omp private(n, l,ilen, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,g,h,wrk,nmax,a,b,c,d,e,len)
#else
!$omp parallel default( none ) &
!$omp private(n, l, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,g,h,wrk,nmax,a,b,c,d,e,len)
#endif

  CALL DNS_OMP_PARTITION(len,srt,end,siz) 
  IF ( siz .LE. 0 ) THEN 
     GOTO 999
  ENDIF

#ifdef USE_BLAS 
  ilen = siz
#endif 

  dummy1 = b(1)
#ifdef USE_BLAS 
  CALL SCAL_LOC(ilen,dummy1,f(srt,1),1)
#else
  DO l=srt,end
     f(l,1) = f(l,1)*dummy1
  ENDDO
#endif 

  DO n = 2, nmax-1
     dummy1 = a(n)
     dummy2 = b(n)
#ifdef USE_BLAS 
     CALL SCAL_LOC(ilen,dummy2,f(srt,n),1) 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n-1),1, f(srt,n), 1)
#else
     DO l=srt,end
        f(l,  n) = f(l,  n)*dummy2 + dummy1*f(l,  n-1)
     ENDDO
#endif
  ENDDO
      
  wrk(srt:end) = C_0_R
  
  DO n = 1, nmax-1
     dummy1 = d(n)
#ifdef USE_BLAS 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n),1,wrk(srt),1)
#else
     DO l=srt,end
        wrk(l) = wrk(l)   + dummy1*f(l,n)
     ENDDO
#endif
  ENDDO

  dummy1 = b(nmax)
#ifdef USE_BLAS 
  CALL SCAL_LOC(ilen, dummy1,f(srt,nmax),1)
  CALL AXPY_LOC(ilen,-dummy1,wrk(srt),1,f(srt,nmax),1)
#else
  DO l=srt,end
     f(l,nmax) = (f(l,nmax) - wrk(l))*dummy1
  ENDDO
#endif

! -------------------------------------------------------------------
! Backward sweep
! -------------------------------------------------------------------
  dummy1 = e(nmax-1)

#ifdef USE_BLAS
  CALL AXPY_LOC(ilen,dummy1,f(srt,nmax),1,f(srt,nmax-1),1)
#else
  DO l=srt,end
     f(l,nmax-1) = dummy1*f(l,nmax) + f(l,nmax-1)
  ENDDO
#endif 


  DO n = nmax-2, 1, -1
     dummy1 = c(n)
     dummy2 = e(n)
#ifdef USE_BLAS 
     CALL AXPY_LOC(ilen,dummy2,f(srt,nmax),1,f(srt,n),1) 
     CALL AXPY_LOC(ilen,dummy1,f(srt,n+1), 1,f(srt,n),1)
     DO l = srt,end
        f(l,n+1) = f(l,n+1) - g(l,n+1)*h(l,n+1)  ! TO BE IMPLEMENTED USING BLAS
     ENDDO
#else
     DO l = srt,end
        f(l,n  ) = f(l,n  ) + dummy1*f(l,n+1) + dummy2*f(l,nmax)
        f(l,n+1) = f(l,n+1) - g(l,n+1)*h(l,n+1)
     ENDDO
#endif
  ENDDO
  DO l = srt, end
     f(l,1   ) =  f(l,1   ) - g(l,1   )*h(l,1   )
     f(l,nmax) =  f(l,nmax) - g(l,nmax)*h(l,nmax)
  ENDDO
999 CONTINUE
!$omp end parallel
  
  RETURN
END SUBROUTINE TRIDPSS_ADD

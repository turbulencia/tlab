#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2011/11/01 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Transposition using Cache-Blocking and OpenMP
!# transpose of the matrix a and place the transposed matrix in b
!# routine trans below is faster than TRANSPOSE routine from f90
!#
!########################################################################
!# ARGUMENTS 
!#
!# nra   In    Number of rows in a
!# nca   In    Number of columns in b
!# ma    In    Leading dimension on the input matrix a
!# mb    In    Leading dimension on the output matrix b
!#
!########################################################################
SUBROUTINE DNS_TRANSPOSE(a, nra, nca, ma, b, mb)
  
  IMPLICIT NONE
  
  TINTEGER jb,kb
#ifdef HLRS_HAWK
  PARAMETER(jb=16,kb=8)
#else
  PARAMETER(jb=64, kb=64)
#endif

  TINTEGER nra, nca, ma, mb
  TREAL a(ma,*),b(mb,*)

  TINTEGER :: srt,end,siz

  TINTEGER k,j,jj,kk
  TINTEGER last_k, last_j

#ifdef USE_MKL 
  CALL MKL_DOMATCOPY('c','t',nra,nca,C_1_R,a,ma,b,mb)
<<<<<<< HEAD
#else 
=======
#else
>>>>>>> upstream_CA/master
   ! use own implementation
!$omp parallel default(none) &
!$omp private(k,j,jj,kk,srt,end,siz,last_k,last_j) &
!$omp shared(a,b,nca,nra)

  CALL DNS_OMP_PARTITION(nca,srt,end,siz)

  kk=1; jj=1

  DO k=srt,end-kb+1,kb; 
     DO j=1,nra-jb+1,jb;
        DO jj=j,j+jb-1
           DO kk=k,k+kb-1
              b(kk,jj) = a(jj,kk)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  last_k = kk 
  last_j = jj 

  DO k=last_k,end
     DO j=1,nra
        b(k,j) = a(j,k)
     ENDDO
  ENDDO

  DO k=srt,end
     DO j=last_j,nra
        b(k,j) = a(j,k) 
     ENDDO
  ENDDO

!$omp end parallel 

#endif 
  RETURN
END SUBROUTINE DNS_TRANSPOSE

!########################################################################
!########################################################################
SUBROUTINE DNS_TRANSPOSE_INT1(a, nra, nca, ma, b, mb)
  
  IMPLICIT NONE
  
  TINTEGER jb,kb
  PARAMETER(jb=32, kb=32)
 
  TINTEGER nra, nca, ma, mb
  INTEGER(1) a(ma,*),b(mb,*)

  TINTEGER :: srt,end,siz

  TINTEGER k,j,jj,kk
  TINTEGER last_k, last_j

!$omp parallel default(none) &
!$omp private(k,j,jj,kk,srt,end,siz,last_k,last_j) &
!$omp shared(a,b,nca,nra)

  CALL DNS_OMP_PARTITION(nca,srt,end,siz)

  kk=1; jj=1

  DO k=srt,end-kb+1,kb; 
     DO j=1,nra-jb+1,jb;
        DO jj=j,j+jb-1
           DO kk=k,k+kb-1
              b(kk,jj) = a(jj,kk)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  last_k = kk 
  last_j = jj 

  DO k=last_k,end
     DO j=1,nra
        b(k,j) = a(j,k)
     ENDDO
  ENDDO

  DO k=srt,end
     DO j=last_j,nra
        b(k,j) = a(j,k) 
     ENDDO
  ENDDO

!$omp end parallel 

  RETURN
END SUBROUTINE DNS_TRANSPOSE_INT1

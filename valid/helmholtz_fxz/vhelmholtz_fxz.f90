#include "types.h"  
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI 
#include "dns_const_mpi.h" 
#endif 
PROGRAM VHELMHOLTZ_FXZ
  
  USE DNS_TYPES, ONLY : pointers_structure
  USE DNS_GLOBAL,  ONLY : iunifx, iunify, iunifz, imode_fdm, imax_total, jmax_total, kmax_total, &
       i1bc,j1bc,k1bc,scalex,scaley,scalez,area,volume,imax,jmax,kmax,&
       inb_grid,inb_wrk1d,inb_wrk2d,isize_wrk1d,isize_wrk2d,gfile,isize_txc_field
#ifdef USE_MPI
  USE DNS_MPI,     ONLY : ims_pro, ims_err 
#endif 

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI 
#include "mpif.h"
#endif 

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: b, c, d,g
  TREAL, DIMENSION(:,:,:,:),ALLOCATABLE:: a,f,e
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: txc
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: bcs_hb, bcs_ht 
  TREAL, DIMENSION(:),     ALLOCATABLE :: cx, cy, cz, wrk3d
  
  TINTEGER i, j, k,  ibc_x(4), ibc_y(4), ibc_z(4), ifield,imeasure
  TINTEGER nfield, nmeasure
  PARAMETER(nfield=4,nmeasure=1)
  TINTEGER opt
  TINTEGER isize_wrk3d
  TREAL dummy, error, rms, max_error, beta ,t_new,t_old 
  CHARACTER*8 :: date
  CHARACTER*10:: time1,time2 
  TINTEGER :: h1,m1,s1,n1,h2,m2,s2,n2 
  TYPE(pointers_structure),DIMENSION(nfield) ::  data 

  TARGET a


#ifdef USE_MPI  
#else 
  TINTEGER :: ims_pro 
  PARAMETER(ims_pro=0) 
#endif

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz

! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')
#ifdef USE_MPI 
  CALL DNS_MPI_INITIALIZE  
#endif 
  
  isize_wrk3d = imax*jmax*kmax
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)

! --------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(imax_total,inb_grid))
  ALLOCATE(y(jmax_total,inb_grid))
  ALLOCATE(z(kmax_total,inb_grid))
  ! ALLOCATE(dx(imax_total,inb_grid))
  ! ALLOCATE(dy(jmax_total,inb_grid))
  ! ALLOCATE(dz(kmax_total,inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,4*nfield+7))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(bcs_ht(imax,kmax,nfield),bcs_hb(imax,kmax,nfield))
  ALLOCATE(a(imax,jmax,kmax,nfield), b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(d(imax,jmax,kmax),e(imax,jmax,kmax,nfield),f(imax,jmax,kmax,nfield),g(imax,jmax,kmax))
  ALLOCATE(txc(isize_txc_field,nfield+1),wrk3d(isize_wrk3d))
  ALLOCATE(cx(6*imax),cy(6*jmax),cz(6*kmax_total))

#include "dns_read_grid.h"

  CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  
! ###################################################################
! Define forcing term
! ###################################################################

  beta = -C_1_R/0.05   
  ! Velocity BCs
  bcs_ht(:,:,1)   = C_1_R
  bcs_ht(:,:,2:3) = C_0_R  
  bcs_hb(:,:,:)   = C_0_R 
  
  ! Scalar BCs  
  IF ( nfield .GT. 3 ) THEN 
     bcs_hb(:,:,4:nfield)   = C_1_R 
     bcs_ht(:,:,4:nfield)   = C_0_R 
  ENDIF

  t_new=0.0
  t_old=0.0
  DO opt=1,2
     CALL DNS_READ_FIELDS('field.inp', i2, imax,jmax,kmax, nfield,i0, isize_wrk3d, a, wrk3d)
     f = a 
     IF ( opt .EQ. 1 ) THEN   
        DO imeasure=1,nmeasure  
           a=f  
           
           DO ifield=1,nfield
              data(ifield)%field => a(:,1,1,ifield) 
           ENDDO

           CALL DATE_AND_TIME(date,time1) 
           CALL OPR_HELMHOLTZ_FXZ_2_N(imax,jmax,kmax, nfield, i0, beta, &
                data, txc(1,1),txc(1,nfield+1), &
                bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,wrk1d(1,4*nfield+1),wrk3d)    
           CALL DATE_AND_TIME(date,time2) 
           READ(time1(1:10),'(i2,i2,i2,1x,i3)') h1,m1,s1,n1
           READ(time2(1:10),'(i2,i2,i2,1x,i3)') h2,m2,s2,n2
           t_new = t_new + (3600*(h2-h1)+60*(m2-m1)+(s2-s1) ) + (n2-n1)/1.0E+03
        ENDDO 
     ELSE IF (opt .EQ. 2) THEN 
        DO imeasure=1,nmeasure  
           a=f 
           CALL DATE_AND_TIME(date,time1)
           DO ifield=1,nfield
              CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, i0, beta, &
                   a(1,1,1,ifield), txc(1,1), txc(1,2),&
                   bcs_hb(1,1,ifield),bcs_ht(1,1,ifield), wrk1d,wrk1d(1,5),wrk3d) 
           ENDDO
           CALL DATE_AND_TIME(date,time2) 
           READ(time1(1:10),'(i2,i2,i2,1x,i3)') h1,m1,s1,n1
           READ(time2(1:10),'(i2,i2,i2,1x,i3)') h2,m2,s2,n2
           t_old = t_old+(3600*(h2-h1)+60*(m2-m1)+(s2-s1) ) + (n2-n1)/1.0E+03
        ENDDO
     ENDIF
#ifdef USE_MPI  
     dummy=t_new
     CALL MPI_Reduce(dummy, t_new, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)  
     dummy=t_old
     CALL MPI_Reduce(dummy, t_old, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)  
     IF ( ims_pro .EQ. 0 ) THEN 
#endif
        IF ( opt .EQ. 1 ) THEN 
           WRITE(*,*) 'NEW SOLVER', t_new/nmeasure, 'ms'  
        ELSE IF ( opt .EQ. 2 ) THEN 
           WRITE(*,*) 'OLD SOLVER', t_old/nmeasure, 'ms'  
        ENDIF
#ifdef USE_MPI 
     ENDIF
#endif 
     ! normalize solution  
     a(:,:,:,:) =a(:,:,:,:)*beta  

     ! ###################################################################
     ! Error
     ! ###################################################################
     DO ifield=1,nfield
        ! ------------------------------------------------------------------- 
        CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax,jmax,kmax, i1bc, &
             dx, a(1,1,1,ifield), b, i0,i0,i0,i0, g,wrk1d,wrk2d,wrk3d) 
        CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax,jmax,kmax, k1bc, &
             dz, a(1,1,1,ifield), c, i0,i0,i0,i0, g,wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_YY(i0, iunify, imode_fdm, imax,jmax,kmax, j1bc, &
             dy, a(1,1,1,ifield), d, i0,i0,i0,i0, g,wrk1d,wrk2d,wrk3d)

        error     = C_0_R
        max_error = C_0_R 
        rms       = C_0_R

        DO k=1,kmax
           DO j = 2,jmax-1
              DO i = 1,imax
                 b(i,j,k)        = a(i,j,k,ifield)+ (b(i,j,k)+c(i,j,k)+d(i,j,k))/beta
                 e(i,j,k,ifield) = b(i,j,k)-f(i,j,k,ifield)
                 error = error + e(i,j,k,ifield)*e(i,j,k,ifield)
                 rms = rms + a(i,j,k,ifield)*a(i,j,k,ifield)

                 IF ( error .GT. max_error ) max_error = error
              ENDDO
           ENDDO
        ENDDO
#ifdef USE_MPI 
         dummy=error
         CALL MPI_Reduce(dummy, error,     i1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err) 
         dummy=rms
         CALL MPI_Reduce(dummy, rms,       i1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)  
         dummy=max_error 
         CALL MPI_Reduce(dummy, max_error, i1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err) 
         IF ( ims_pro .EQ. 0 ) THEN 
#endif 
           WRITE(*,*) ims_pro, ifield, 'Relative error .............: ', sqrt(error)/sqrt(rms)
           WRITE(*,*) ims_pro, ifield, 'Maximum error ..............: ', max_error 
#ifdef USE_MPI 
         ENDIF
#endif 
     ENDDO
  ENDDO
  
  CALL DNS_END(0)
  
  STOP
END PROGRAM VHELMHOLTZ_FXZ

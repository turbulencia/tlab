#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif
PROGRAM VHELMHOLTZ_FXZ

  USE TLAB_TYPES, ONLY : pointers_dt
  USE TLAB_VARS, ONLY : imax,jmax,kmax, inb_wrk1d,inb_wrk2d,isize_wrk1d,isize_wrk2d,gfile,isize_txc_field
  USE TLAB_PROCS
  USE IO_FIELDS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: b, c, d, h
  TREAL, DIMENSION(:,:,:,:),ALLOCATABLE:: a,f,e
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: txc
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: bcs_hb, bcs_ht
  TREAL, DIMENSION(:),     ALLOCATABLE :: cx, cy, cz, wrk3d

  TINTEGER i, j, k,  ibc_x(4), ibc_y(4), ibc_z(4), ifield,imeasure
  TINTEGER nfield, nmeasure
  PARAMETER(nfield=4,nmeasure=1)
  TINTEGER opt
  TREAL dummy, error, rms, max_error, beta ,t_new,t_old
  CHARACTER*8 :: date
  CHARACTER*10:: time1,time2
  TINTEGER :: h1,m1,s1,n1,h2,m2,s2,n2
  TYPE(pointers_dt),DIMENSION(nfield) ::  data

  TARGET a


#ifdef USE_MPI
#else
  TINTEGER :: ims_pro
  PARAMETER(ims_pro=0)
#endif

! ###################################################################
  CALL TLAB_START()

  CALL IO_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  isize_wrk3d = imax*jmax*kmax
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)

! --------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,4*nfield+7))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(bcs_ht(imax,kmax,nfield),bcs_hb(imax,kmax,nfield))
  ALLOCATE(a(imax,jmax,kmax,nfield), b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(d(imax,jmax,kmax),e(imax,jmax,kmax,nfield),f(imax,jmax,kmax,nfield),h(imax,jmax,kmax))
  ALLOCATE(txc(isize_txc_field,nfield+1),wrk3d(isize_wrk3d))
  ALLOCATE(cx(6*imax),cy(6*jmax),cz(6*g(3)%size))

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

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
     CALL IO_READ_FIELDS('field.inp', IO_FLOW, imax,jmax,kmax, nfield,i0, a, wrk3d)
     f = a
     IF ( opt .EQ. 1 ) THEN
        DO imeasure=1,nmeasure
           a=f

           DO ifield=1,nfield
              data(ifield)%field => a(:,1,1,ifield)
           ENDDO

           CALL DATE_AND_TIME(date,time1)
           CALL OPR_HELMHOLTZ_FXZ_2_N(imax,jmax,kmax, h, nfield, i0, beta, &
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
              CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, h, i0, beta, &
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
        CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), a(1,1,1,ifield), b, h, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), a(1,1,1,ifield), c, h, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), a(1,1,1,ifield), d, h, wrk2d,wrk3d)

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

  CALL TLAB_STOP(0)
END PROGRAM VHELMHOLTZ_FXZ

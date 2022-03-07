#include "types.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE DNS_TOWER - save columns at every iterations to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RELATED INPUT FROM dns.ini
!
! [SaveTowers]
! Stride=<stride_x>,<stride_y>,<stride_z>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINES
!
! DNS_TOWER_INITIALIZE(x,y,z,stride)
!     - find out which towers to process (according to stride)
!     - allocate corresponding space for arrays
!     - organize buffers where tower data are saved until write
!
! DNS_TOWER_ACCUMULATE(v,index,wrk1d)
!     - PARAMTERS: v        -- data
!                  index    -- what to write (1 - flow, 2 - scalars, 4 - pressure)
!                  wrk1d    -- wrkspace for averaging
!     - accumulates data from flow arrays into tower buffers
!       until towers are written to disk
!
! DNS_TOWER_WRITE(wrk3d)
!     - writes tower data from buffers to disk
!
! DNS_TOWER_FINALIZE()
!     - nothing so far
!
! TOWER_AVG_IK_V(imax, jmax, kmax, a, avg, wrk)
!     - average 1 Variable horizontally over all planes in the vertical.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE DNS_TOWER

  USE TLAB_VARS, ONLY : inb_flow, inb_scal
  USE TLAB_PROCS

  TINTEGER tower_imax, tower_jmax,tower_kmax, tower_isize_field
  TINTEGER, TARGET :: tower_isize_plane,tower_isize_plane_total
  TINTEGER tower_imax_total, tower_jmax_total,tower_kmax_total
  TINTEGER tower_offset_i, tower_offset_j, tower_offset_k
  TINTEGER tower_isize_acc_field, tower_isize_acc_mean, tower_isize_acc_write
  TINTEGER tower_accumulation,tower_varcount
  TINTEGER tower_ncid, tower_ncmid, tower_stat
  TINTEGER tower_istride, tower_jstride, tower_kstride
  TINTEGER tower_master
  TINTEGER, DIMENSION(:),ALLOCATABLE :: tower_ipos, tower_kpos, tower_jpos
  TINTEGER ::  tower_mode_check

  INTEGER(KIND=8) :: tower_bufsize
  CHARACTER(LEN=128) :: tower_nc_name

  TREAL, DIMENSION(:), ALLOCATABLE, TARGET :: tower_buf
  TREAL, DIMENSION(:), POINTER :: tower_u, tower_um, tower_v, tower_vm, tower_w, tower_wm
  TREAL, DIMENSION(:), POINTER :: tower_p, tower_pm, tower_s, tower_sm, tower_t, tower_it

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE DNS_TOWER_INITIALIZE(stride)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE TLAB_VARS,ONLY : imax,jmax,kmax
    USE DNS_LOCAL, ONLY : nitera_save

#ifdef USE_MPI
    USE MPI
    USE TLAB_VARS,ONLY : g
    USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k,ims_pro
#endif

    IMPLICIT NONE

#ifdef USE_MPI
#else
    TINTEGER :: ims_offset_i,ims_offset_j, ims_offset_k,ims_pro
    PARAMETER(ims_offset_i=0,ims_offset_j=0,ims_offset_k=0,ims_pro=-1)
#endif

    TINTEGER, DIMENSION(3), INTENT(IN) :: stride

    TINTEGER :: istart, iend, jstart,jend,kstart,kend, ibuf,i,j,k,ii
    TINTEGER, POINTER :: tip, tip_total
    !
    IF ( ims_offset_i .EQ. 0 .AND. ims_offset_k .EQ. 0 ) THEN
       tower_master = 1
    ELSE
       tower_master = 0
    ENDIF

    tip => tower_isize_plane
    tip_total => tower_isize_plane_total

    tower_isize_acc_write = nitera_save
    tower_varcount=inb_flow+inb_scal+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DISTRIBUTE TOWERS:

    tower_istride = stride(1)
    tower_jstride = stride(2)
    tower_kstride = stride(3)

    istart=ims_offset_i; iend=ims_offset_i + imax -1
    jstart=ims_offset_j; jend=ims_offset_j + jmax -1
    kstart=ims_offset_k; kend=ims_offset_k + kmax -1

    tower_imax=0;    tower_jmax=0;    tower_kmax=0

    tower_offset_i = CEILING(DBLE(istart)/tower_istride)
    tower_offset_j = CEILING(DBLE(jstart)/tower_jstride)
    tower_offset_k = CEILING(DBLE(kstart)/tower_kstride)

    DO i = istart,iend
       IF ( MOD(i-1, tower_istride ) .EQ. 0 ) THEN
          tower_imax = tower_imax + 1
       ENDIF
    ENDDO

    DO j = jstart,jend
       IF ( MOD(j-1, tower_jstride ) .EQ. 0 ) THEN
          tower_jmax = tower_jmax + 1
       ENDIF
    ENDDO

    DO k=kstart,kend
       IF ( MOD(k-1, tower_kstride) .EQ. 0 ) THEN
          tower_kmax = tower_kmax+1
       ENDIF
    ENDDO

#ifdef USE_MPI
    CALL MPI_Allreduce(tower_imax, tower_imax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD,i)
    CALL MPI_Allreduce(tower_jmax, tower_jmax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD,i)
    CALL MPI_Allreduce(tower_kmax, tower_kmax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD,i)

    tower_imax_total = tower_imax_total / (g(3)%size/kmax * g(2)%size/jmax)
    tower_jmax_total = tower_jmax_total / (g(1)%size/imax * g(3)%size/kmax)
    tower_kmax_total = tower_kmax_total / (g(1)%size/imax * g(2)%size/jmax)
#else
    tower_imax_total = tower_imax
    tower_jmax_total = tower_jmax
    tower_kmax_total = tower_kmax
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALLOCATE SPACE FOR TOWERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tower_isize_plane = tower_imax * tower_kmax
    tower_isize_field = tower_imax * tower_jmax * tower_kmax
    tower_isize_plane_total = tower_imax_total * tower_kmax_total
    tower_isize_acc_field = tower_isize_acc_write*tower_jmax*tip
    tower_isize_acc_mean  = tower_isize_acc_write*tower_jmax
    tower_bufsize = 5*tower_isize_acc_field + 5*tower_isize_acc_mean + 2* tower_isize_acc_write

    ALLOCATE(tower_ipos(tower_imax), tower_kpos(tower_kmax), tower_jpos(tower_jmax))
    ALLOCATE(tower_buf( tower_bufsize))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! REMEMBER WHICH TOWERS TO SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ii=1
    DO i=istart,iend
       IF ( MOD(i, tower_istride ) .EQ. 0 ) THEN
          tower_ipos(ii) = i+1 - ims_offset_i
          ii=ii+1
       ENDIF
    ENDDO

    ii=1
    DO j=jstart,jend
       IF ( MOD(j, tower_jstride ) .EQ. 0 ) THEN
          tower_jpos(ii) = j+1 - ims_offset_j
          ii=ii+1
       ENDIF
    ENDDO

    ii=1
    DO k=kstart,kend
       IF ( MOD(k,tower_kstride ) .EQ. 0 ) THEN
          tower_kpos(ii) = k+1 - ims_offset_k
          ii=ii+1
       ENDIF
    ENDDO

    ibuf = 1;
    tower_u(1:)=>tower_buf(ibuf:); ibuf = ibuf+tower_isize_acc_field;
    tower_v(1:)=>tower_buf(ibuf:); ibuf = ibuf+tower_isize_acc_field;
    tower_w(1:)=>tower_buf(ibuf:); ibuf = ibuf+tower_isize_acc_field;
    tower_p(1:)=>tower_buf(ibuf:); ibuf = ibuf+tower_isize_acc_field;
    tower_s(1:)=>tower_buf(ibuf:); ibuf = ibuf+tower_isize_acc_field;
    tower_um(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_mean;
    tower_vm(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_mean;
    tower_wm(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_mean;
    tower_pm(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_mean;
    tower_sm(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_mean;
    tower_t(1:) =>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_write;
    tower_it(1:)=>tower_buf(ibuf:);ibuf = ibuf+tower_isize_acc_write;

    tower_accumulation = 1

  END SUBROUTINE DNS_TOWER_INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE DNS_TOWER_ACCUMULATE(v,index,wrk1d)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
    USE MPI
    USE TLAB_MPI_VARS, ONLY : ims_err
#endif

    USE TLAB_VARS, ONLY : imax,jmax,kmax
    USE TLAB_VARS, ONLY : itime,rtime

    IMPLICIT NONE

    ! 1 -- flow;  2 -- scalar;  4 -- pressure
    TINTEGER,                           INTENT(IN)   :: index
    TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)   :: v
    TREAL, DIMENSION(*),                INTENT(INOUT):: wrk1d

    TINTEGER  :: ii,kk,ip,ipm
    TINTEGER, POINTER :: tip


    tip => tower_isize_plane
    ip = 1+(tower_accumulation-1)*tower_jmax*tip; ipm=ip+tower_jmax-1

    ! HANDLE PRESSURE
    ! The Counter tower_accumulation is only increased when the pressure is handled
    ! - This needs to be done first!
    ! - The iteration is not increased yet when the pressure is calculated -> use it
    IF ( index .EQ. 4 ) THEN
       tower_t(tower_accumulation) = rtime
       tower_it(tower_accumulation)= itime
       DO kk=1,tower_kmax
          DO ii=1,tower_imax
             tower_p(ip:ipm) = v(tower_ipos(ii),1:tower_jmax:tower_jstride, tower_kpos(kk),1)
             ip = ip + tower_jmax; ipm= ipm+tower_jmax
          ENDDO
       ENDDO
       ip = 1+(tower_accumulation-1)*tower_jmax; ipm=ip+tower_jmax-1

       CALL TOWER_AVG_IK_V(imax,jmax,kmax,v(1,1,1,1),tower_pm(ip:ipm),wrk1d(6*jmax))
    ! HANDLY FLOW FIELDS
    ELSE IF ( index .EQ. 1 ) THEN
       DO kk=1,tower_kmax
          DO ii=1,tower_imax
             tower_u(ip:ipm) = &
                  v(tower_ipos(ii),1:tower_jmax:tower_jstride, tower_kpos(kk),1)
             tower_v(ip:ipm) = &
                  v(tower_ipos(ii),1:tower_jmax:tower_jstride, tower_kpos(kk),2)
             tower_w(ip:ipm) = &
                  v(tower_ipos(ii),1:tower_jmax:tower_jstride, tower_kpos(kk),3)
             ip = ip + tower_jmax; ipm= ipm+tower_jmax
          ENDDO
       ENDDO
       ip = 1+(tower_accumulation-1)*tower_jmax; ipm=ip+tower_jmax-1
       CALL TOWER_AVG_IK_V(imax,jmax,kmax,v(1,1,1,1),tower_um(ip:ipm),wrk1d(6*jmax))
       CALL TOWER_AVG_IK_V(imax,jmax,kmax,v(1,1,1,2),tower_vm(ip:ipm),wrk1d(6*jmax))
       CALL TOWER_AVG_IK_V(imax,jmax,kmax,v(1,1,1,3),tower_wm(ip:ipm),wrk1d(6*jmax))

    ! HANDLE SCALARS
    ELSE IF (index .EQ. 2 ) THEN
       DO kk=1,tower_kmax
          DO ii=1,tower_imax
             tower_s(ip:ipm) = &
                  v(tower_ipos(ii),1:tower_jmax:tower_jstride, tower_kpos(kk),1)
             ip = ip + tower_jmax; ipm= ipm+tower_jmax
          ENDDO
       ENDDO
       ip = 1+(tower_accumulation-1)*tower_jmax; ipm=ip+tower_jmax-1
       CALL TOWER_AVG_IK_V(imax,jmax,kmax,v(1,1,1,1),tower_sm(ip:ipm),wrk1d(6*jmax))

    ENDIF


    tower_mode_check = tower_mode_check + index

    IF ( tower_mode_check .EQ. 7 ) THEN
       tower_mode_check = 0
       tower_accumulation = tower_accumulation+1
    ELSE IF ( tower_mode_check .GT. 7 ) THEN
       WRITE(*,*) 'ERROR - tower_mode_check GREATER THAN 7'
       STOP 'INTERNAL ERROR'
    ENDIF

  END SUBROUTINE DNS_TOWER_ACCUMULATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE DNS_TOWER_WRITE(wrk3d)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE DNS_LOCAL,     ONLY : nitera_save
    USE TLAB_VARS,    ONLY : itime
    USE TLAB_CONSTANTS, ONLY : wfile
#ifdef USE_MPI
    USE MPI
    USE TLAB_MPI_VARS,   ONLY : ims_offset_i, ims_offset_j, ims_offset_k,ims_pro,ims_err
#endif
    IMPLICIT NONE

#ifdef USE_MPI
#else
    TINTEGER :: ims_offset_i, ims_offset_j, ims_offset_k, ims_pro
    PARAMETER(ims_offset_i=0, ims_offset_j=0, ims_offset_k=0,ims_pro=-1)
#endif

    TREAL, DIMENSION(*), INTENT(INOUT) :: wrk3d

    TINTEGER :: it,ivar
    CHARACTER(LEN=64) :: cdummy

    TINTEGER, POINTER :: tip
    TINTEGER :: itower,ktower,ip_skp,ip_srt,ip_end,tower_count,op_srt,op_end
    tip => tower_isize_plane

    IF ( tip .LT. 1 ) THEN
       ! DO NOTHING
    ELSE

       IF  ( itime .NE. INT(tower_it(nitera_save)+1)  )  THEN
          CALL TLAB_WRITE_ASCII(wfile,'tools/dns/dns_tower.f90 (DNS_TOWER_WRITE)')
          CALL TLAB_WRITE_ASCII(wfile,'nitera_save for towers does not match current iteration')
          !                          (But it should if the code is set-up properly)
       ENDIF

       DO ivar=1,tower_varcount
          tower_count = 0
#ifdef USE_MPI
          if ( ims_pro .EQ. 0 ) THEN
#endif
             op_srt=1; op_end=tower_jmax+2;
             ip_skp = tower_jmax*1;
             ip_srt = 1;
             ip_end = ip_srt + tower_jmax - 1;

             DO it=1,nitera_save
                wrk3d(op_srt) =  tower_t(it)
                wrk3d(op_srt+1)= tower_it(it)
                SELECT CASE(ivar)
                CASE(1)
                   wrk3d(op_srt+2:op_end) = tower_um(ip_srt:ip_end)
                CASE(2)
                   wrk3d(op_srt+2:op_end) = tower_vm(ip_srt:ip_end)
                CASE(3)
                   wrk3d(op_srt+2:op_end) = tower_wm(ip_srt:ip_end)
                CASE(4)
                   wrk3d(op_srt+2:op_end) = tower_pm(ip_srt:ip_end)
                CASE(5)
                   wrk3d(op_srt+2:op_end) = tower_sm(ip_srt:ip_end)
                CASE DEFAULT
                   ! ISSUE WARNING - NO MORE THAN ONE SCALAR
                END SELECT
                ip_srt = ip_srt + ip_skp; op_srt = op_srt + tower_jmax+2
                ip_end = ip_end + ip_skp; op_end = op_end + tower_jmax+2
             ENDDO

             WRITE(cdummy,995) &
                  INT(tower_it(1))+1,itime,ivar
995          FORMAT('tower.mean','.',I6.6,'-',I6.6,'.',I1)
             OPEN(73,FILE=TRIM(ADJUSTL(cdummy)),ACCESS='STREAM', FORM='UNFORMATTED')
             WRITE(73,POS=1) wrk3d(1:nitera_save*(tower_jmax+2))
             CLOSE(73)

#ifdef USE_MPI
          ENDIF
#endif

          DO itower=1,tower_imax
             DO ktower=1,tower_kmax
                ip_skp = tower_jmax*tip;
                ip_srt = tower_count*tower_jmax + 1;
                ip_end = ip_srt + tower_jmax - 1;
                op_srt=1; op_end=tower_jmax+2;

                DO it=1,nitera_save
                   wrk3d(op_srt) =  tower_t(it)
                   wrk3d(op_srt+1)= tower_it(it)
                   SELECT CASE(ivar)
                   CASE(1)
                      wrk3d(op_srt+2:op_end) = tower_u(ip_srt:ip_end)
                   CASE(2)
                      wrk3d(op_srt+2:op_end) = tower_v(ip_srt:ip_end)
                   CASE(3)
                      wrk3d(op_srt+2:op_end) = tower_w(ip_srt:ip_end)
                   CASE(4)
                      wrk3d(op_srt+2:op_end) = tower_p(ip_srt:ip_end)
                   CASE(5)
                      wrk3d(op_srt+2:op_end) = tower_s(ip_srt:ip_end)
                   CASE DEFAULT
                      ! ISSUE WARNING - NO MORE THAN ONE SCALAR
                   END SELECT
                   ip_srt = ip_srt + ip_skp; op_srt = op_srt + tower_jmax+2
                   ip_end = ip_end + ip_skp; op_end = op_end + tower_jmax+2
                ENDDO
                tower_count = tower_count + 1

                WRITE(cdummy,997) &
                     ims_offset_i+tower_ipos(itower), &
                     ims_offset_k+tower_kpos(ktower),&
                     INT(tower_it(1))+1,itime,ivar
997             FORMAT('tower.',I6.6,'x',I6.6,'.',I6.6,'-',I6.6,'.',I1)
                OPEN(73,FILE=TRIM(ADJUSTL(cdummy)),ACCESS='STREAM', FORM='UNFORMATTED')
                WRITE(73,POS=1) wrk3d(1:nitera_save*(tower_jmax+2))
                CLOSE(73)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    tower_accumulation = 1

  END SUBROUTINE DNS_TOWER_WRITE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE DNS_TOWER_FINALIZE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

  END SUBROUTINE DNS_TOWER_FINALIZE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE TOWER_AVG_IK_V(imax, jmax, kmax, a, avg, wrk)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "types.h"

#ifdef USE_MPI
    USE MPI
#endif
    USE TLAB_VARS, ONLY : area, g

    IMPLICIT NONE

    TINTEGER imax, jmax, kmax
    TREAL a(imax, jmax, kmax)
    TREAL avg(tower_jmax), wrk(tower_jmax)
#ifdef USE_MPI
    INTEGER ims_err, len
#endif

    TINTEGER i, j, k

    DO j=1,tower_jmax
       avg(j) = C_0_R
    ENDDO

    DO k = 1, kmax
       DO j=1,tower_jmax
          DO i = 1, imax
             avg(j) = avg(j) + a(i,tower_jpos(j),k) *g(1)%jac(i,1) *g(3)%jac(k,1)
          ENDDO
       ENDDO
    ENDDO

#ifdef USE_MPI
    CALL MPI_REDUCE(avg, wrk, tower_jmax, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    DO j=1,tower_jmax
       avg(j) = wrk(j)/area
    ENDDO
#else
    DO j=1,tower_jmax
       avg(j) = avg(j)/area
    ENDDO
#endif

    RETURN
  END SUBROUTINE TOWER_AVG_IK_V

END MODULE DNS_TOWER

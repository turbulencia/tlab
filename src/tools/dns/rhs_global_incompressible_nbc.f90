#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "info_vars.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#ifdef USE_PSFFT
#include "nb3dfft_defines.inc"
#endif

SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC(u,v,w,s,&
     tmpu,tmpw,tmp11,tmp12,tmp21,tmp22,tmp31,tmp32,tmp41,tmp42,&
     bt1,bt2,bt3,bt4,&
     h1,h2,h3,hs,&
     wrk1d,wrk2d,wrk3d)
  USE, INTRINSIC :: iso_c_binding, ONLY : c_int,c_loc,c_ptr,c_f_pointer

  USE OMP_LIB,    ONLY : omp_get_thread_num

  USE TLAB_CONSTANTS, ONLY : lfile,wfile,efile,tfile
  !
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imode_eqns
  USE TLAB_VARS, ONLY : inb_flow,inb_scal,inb_scal_array
  USE TLAB_VARS, ONLY : isize_field, isize_wrk1d, imax,jmax,kmax
  USE TLAB_VARS, ONLY : rbackground, ribackground
  !
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS
  USE TIME,  ONLY : rkm_substep,rkm_endstep, dte
  USE DNS_LOCAL,  ONLY : use_tower

  USE DNS_TOWER

#ifdef USE_PSFFT
  USE DNS_LOCAL,  ONLY : nbcsetup
#endif

  USE TLAB_MPI_VARS,    ONLY : ims_npro, ims_pro, ims_err,ims_size_i,ims_size_k

  USE NB3DFFT,    ONLY : nb3dfft_nbc_prepare,nb3dfft_nbc_finish,nb3dfft_infoType
  USE NB3DFFT,    ONLY : nb3dfft_nbc_schedl_start, nb3dfft_nbc_worker_start
  USE NB3DFFT,    ONLY : nb3dfft_nbc_schedl_stop,  nb3dfft_schedlType
  USE NB3DFFT,    ONLY : nb3dfft_r2r_yxcomm,       nb3dfft_r2r_yzcomm
  USE NB3DFFT,    ONLY : nb3dfft_r2r_xycomm,       nb3dfft_r2r_zycomm
  USE NB3DFFT,    ONLY : nb3dfft_r2r_xunpack,      nb3dfft_r2r_zunpack
  USE NB3DFFT,    ONLY : nb3dfft_r2r_y1unpack,     nb3dfft_r2r_y2unpack
  USE NB3DFFT,    ONLY : nb3dfft_r2r_ready,        mytype

  USE MPI

  IMPLICIT NONE

#include "integers.h"
  !
  ! PARAMETERS
  !
  TINTEGER, PARAMETER :: nmeasure=3

  TREAL, DIMENSION(isize_field),                INTENT(IN)   :: u,v,w
  TREAL, DIMENSION(isize_field,inb_scal_array), INTENT(IN)   :: s

  TREAL, DIMENSION(isize_field),                INTENT(INOUT):: h1,h2,h3
  TREAL, DIMENSION(isize_field,inb_scal),TARGET,INTENT(OUT)  :: hs
  TREAL, DIMENSION(isize_field),                INTENT(INOUT):: tmpu,tmpw,tmp11,tmp12,tmp21,tmp22,tmp31,tmp32,tmp41,tmp42
  TREAL, DIMENSION(isize_field)  :: bt1,bt2,bt3,bt4
  TREAL, DIMENSION(isize_wrk1d,*):: wrk1d
  TREAL, DIMENSION(*)            :: wrk2d,wrk3d
  TREAL, DIMENSION(:),POINTER :: p_h
  !
  ! LOCAL VARIABLES
  !
  TINTEGER :: nxy_trans,nyz_trans,nxy,id,imeasure,ij,k,is,commID,iq
  TINTEGER :: finished,ip_b,ip_t,ibc, bcs(2,2)
  TREAL tdummy
  TREAL, DIMENSION(:), POINTER :: p_bcs
  !
  TYPE(nb3dfft_infoType), DIMENSION(24) :: info
  INTEGER(KIND=4),DIMENSION(2) :: cur_time
  REAL(KIND=mytype) :: t_comp,t_test,t_ser,t_init,t_wait,t_tmp
  REAL(KIND=8),DIMENSION(6)      :: t_snd, t_rcv
  TYPE(nb3dfft_schedlType) :: nbcsetup_
  INTEGER :: pkg_cnt
  TREAL rtime, rtime_loc,t_run,ptime,ctime_loc

  TARGET h1,h2,h3


#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'ENTERING SUBROUTINE, RHS_GLOBAL_INCOMPRESSIBLE_NBC')
#endif

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero
  BcsScalJmin%ref(:,:,:) = C_0_R
  BcsScalJmax%ref(:,:,:) = C_0_R

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     DO is =1,inb_scal
        IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
             p_bcs => hs(ip_b:,is); BcsScalJmin%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
        IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
             p_bcs => hs(ip_t:,is); BcsScalJmax%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     ENDDO
     ip_b = ip_b + nxy ! bottom BC address
     ip_t = ip_t + nxy ! top BC address
  ENDDO


  nbcsetup_ = nbcsetup
  pkg_cnt=24*ims_npro
  t_comp=0
  t_test=0
  t_ser =0
  ptime =0

  nxy=imax*jmax

  rtime = -MPI_WTime()
  t_run = rtime
  t_init= rtime

! #######################################################################
! Advection-diffusion terms
! #######################################################################
  CALL NB3DFFT_NBC_PREPARE(24,.FALSE.)

!$omp parallel num_threads(2)  &
!$omp default(shared)
  IF (omp_get_thread_num() == 0) THEN
     CALL NB3DFFT_NBC_SCHEDL_START(nbcsetup_)
  ELSE IF (omp_get_thread_num() == 1) THEN
     CALL NB3DFFT_NBC_WORKER_START()
     ;                      rtime_loc=MPI_WTime();         commID=1;
     info(FUYX)%timer_start=rtime_loc;info(FUYX)%id=commID;commID=commID+1     ! STEP1
     info(BUXY)%timer_start=rtime_loc;info(BUXY)%id=commID;commID=commID+1     ! STEP1
     info(FWYZ)%timer_start=rtime_loc;info(FWYZ)%id=commID;commID=commID+1     ! STEP1
     info(BWZY)%timer_Start=rtime_loc;info(BWZY)%id=commID;commID=commID+1     ! STEP1
     !
     IF ( inb_scal .GT. 0 ) THEN
        info(F1YX)%timer_start=rtime_loc;info(F1YX)%id=commID;commID=commID+1  ! STEP2
        info(B1XY)%timer_start=rtime_loc;info(B1XY)%id=commID;commID=commID+1  ! STEP2
        info(F1YZ)%timer_start=rtime_loc;info(F1YZ)%id=commID;commID=commID+1  ! STEP2
        info(B1ZY)%timer_start=rtime_loc;info(B1ZY)%id=commID;commID=commID+1  ! STEP2
     ENDIF

     info(FUYZ)%timer_start=rtime_loc; info(FUYZ)%id=commID;commID=commID+1    ! STEP3
     info(BUZY)%timer_start=rtime_loc; info(BUZY)%id=commID;commID=commID+1    ! STEP3
     info(FVYZ)%timer_start=rtime_loc; info(FVYZ)%id=commID;commID=commID+1    ! STEP3
     info(BVZY)%timer_start=rtime_loc; info(BVZY)%id=commID;commID=commID+1    ! STEP3

     IF ( inb_scal .EQ. 2 ) THEN
        info(F2YX)%timer_start=rtime_loc; info(F2YX)%id=commID;commID=commID+1 ! STEP4
        info(B2XY)%timer_start=rtime_loc; info(B2XY)%id=commID;commID=commID+1 ! STEP4
        info(F2YZ)%timer_start=rtime_loc; info(F2YZ)%id=commID;commID=commID+1 ! STEP4
        info(B2ZY)%timer_start=rtime_loc; info(B2ZY)%id=commID;commID=commID+1 ! STEP4
     ENDIF

     info(FVYX)%timer_start=rtime_loc; info(FVYX)%id=commID;commID=commID+1    ! STEP5
     info(BVXY)%timer_start=rtime_loc; info(BVXY)%id=commID;commID=commID+1    ! STEP5
     info(FWYX)%timer_start=rtime_loc; info(FWYX)%id=commID;commID=commID+1    ! STEP5
     info(BWXY)%timer_start=rtime_loc; info(BWXY)%id=commID;commID=commID+1    ! STEP5

     info(FPYX)%timer_start=rtime_loc; info(FPYX)%id=commID;commID=commID+1    ! STEP6
     info(BPXY)%timer_start=rtime_loc; info(BPXY)%id=commID;commID=commID+1    ! STEP6
     info(FPYZ)%timer_start=rtime_loc; info(FPYZ)%id=commID;commID=commID+1    ! STEP6
     info(BPZY)%timer_start=rtime_loc; info(BPZY)%id=commID;                   ! STEP6

     t_init = t_init + MPI_WTime()

     id = TLAB_MPI_I_PARTIAL;   nyz_trans = ims_size_i(id)
     id = TLAB_MPI_K_PARTIAL;   nxy_trans = ims_size_k(id)
     !
     ! kick off transpose U y->x and W y->z
     CALL NB3DFFT_R2R_YXCOMM(u,bt1,  bt1,  tmp11,info(FUYX),t_tmp);t_comp=t_comp+t_tmp
     CALL NB3DFFT_R2R_YZCOMM(w,tmpw, tmpw, bt2,  info(FWYZ),t_tmp);t_comp=t_comp+t_tmp
     CALL NB3DFFT_R2R_YXCOMM(w,bt3,  bt3,  tmp31,info(FWYX),t_tmp);t_comp=t_comp+t_tmp
     !
     ! Vertical derivatives, and Vertical advection
     !
     t_tmp = -MPI_WTime()
     CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), v,v,v,     tmp21, tmp22, wrk2d,wrk3d) ! store v transposed in tmp22
     h2 = h2 + tmp21
     CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), u,v,tmp22, tmp21, tmpu,  wrk2d,wrk3d) ! using tmp22
     h1 = h1 + tmp21
     CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), w,v,tmp22, tmp21, tmpu,  wrk2d,wrk3d) ! using tmp22
     h3 = h3 + tmp21
     t_ser = t_ser + (t_tmp +MPI_WTime())

     CALL NB3DFFT_R2R_YZCOMM(u,tmp41,tmp41,bt4,  info(FUYZ),t_tmp);t_comp=t_comp+t_tmp

     t_tmp = -MPI_WTime()
     DO is = 1,inb_scal
        CALL OPR_BURGERS_Y(i1,is, imax,jmax,kmax, bcs, g(2), s(1,is),v,tmp22, tmp21, tmpu, wrk2d,wrk3d) ! using tmp22
        hs(:,is) = hs(:,is) + tmp21
     ENDDO
     t_ser = t_ser + (t_tmp +MPI_WTime())

     finished = 0
     DO WHILE ( finished /= 2 )
        IF ( NB3DFFT_R2R_READY(info(FUYX), t_tmp ) ) THEN
           t_test = t_test + t_tmp
           ! u du/dx + 1/Re d2u/dx2
           CALL NB3DFFT_R2R_XUNPACK(bt1,tmp11, info(FUYX), t_tmp);t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt1,g(1)%size,nyz_trans,g(1)%size,tmpu,nyz_trans)
           CALL OPR_BURGERS(0, nyz_trans, bcs, g(1), tmpu,tmpu,tmp11, wrk2d,wrk3d)
           CALL DNS_TRANSPOSE(tmp11,nyz_trans,g(1)%size,nyz_trans,bt1,g(1)%size)
           t_ser = t_ser + (t_tmp + MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt1,bt1,tmp12,tmp11,info(BUXY), t_tmp);t_comp=t_comp+t_tmp
           finished = finished + 1
        ENDIF

        IF ( NB3DFFT_R2R_READY(info(FWYZ), t_tmp ) ) THEN
           t_test = t_test +t_tmp
           ! w dw/dz + 1/Re d2w/dz2
           CALL NB3DFFT_R2R_ZUNPACK(tmpw,bt2,info(FWYZ),t_tmp);t_comp = t_comp +t_tmp
           !
           t_tmp = -MPI_WTime()
           CALL OPR_BURGERS(0, nxy_trans, bcs, g(3), tmpw,tmpw,bt2, wrk2d,wrk3d)
           t_ser = t_ser + (t_tmp + MPI_WTime())
           !
           CALL NB3DFFT_R2R_ZYCOMM(bt2,bt2,tmp22,tmp21,info(BWZY),t_tmp); t_comp=t_comp+t_tmp
           finished = finished + 1
        ENDIF
     ENDDO
     !
     DO WHILE ( finished /= 4 )
        IF ( NB3DFFT_R2R_READY(info(FWYX), t_tmp) )  THEN
           t_test=t_test+t_tmp
           ! u dw/dx + 1/Re d2w/dx2
           CALL NB3DFFT_R2R_XUNPACK(bt3,tmp31,info(FWYX),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt3,g(1)%size,nyz_trans,g(1)%size,tmp31,nyz_trans)
           CALL OPR_BURGERS(0, nyz_trans, bcs, g(1), tmp31,tmpu,tmp32, wrk2d,wrk3d)
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,g(1)%size,nyz_trans,bt3,g(1)%size)
           t_ser = t_ser + (t_tmp + MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt3,bt3,tmp32,tmp31,info(BWXY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(FUYZ), t_tmp) ) THEN
           t_test=t_test+t_tmp
           ! w du/dz + 1/Re d2u/dz2
           CALL NB3DFFT_R2R_ZUNPACK(tmp41,bt4,info(FUYZ),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL OPR_BURGERS(0, nxy_trans, bcs, g(3), tmp41,tmpw,bt4, wrk2d,wrk3d)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_ZYCOMM(bt4,bt4,tmp42,tmp41,info(BUZY),t_tmp);t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
     ENDDO
     !
     !
     DO WHILE ( finished /= 8 )
        IF ( NB3DFFT_R2R_READY(info(BUXY),t_tmp) ) THEN
           t_test=t_test+t_tmp;
           CALL NB3DFFT_R2R_Y1UNPACK(bt1,tmp11,info(BUXY),t_tmp);     t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h1=h1+bt1
           t_ser = t_ser + (t_tmp + MPI_WTime())
           !
           CALL NB3DFFT_R2R_YXCOMM(v,bt1,bt1,tmp11,info(FVYX),t_tmp); t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(BWZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt2,tmp21,info(BWZY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h3=h3+bt2
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_YZCOMM(v,tmp21,tmp21,bt2,info(FVYZ),t_tmp);t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(BWXY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y1UNPACK(bt3,tmp31,info(BWXY),t_tmp);         t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h3=h3+bt3
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_YXCOMM(s(:,1),bt3,bt3,tmp31,  info(F1YX),t_tmp);t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(BUZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt4,tmp41,info(BUZY),t_tmp);         t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h1=h1+bt4
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_YZCOMM(s(:,1),tmp41,tmp41,bt4,info(F1YZ),t_tmp);t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
     ENDDO
     ! u and w are finished and v and s1 are initiated
     !
     DO WHILE ( finished /= 12 )
        IF ( NB3DFFT_R2R_READY(info(FVYX), t_tmp) )  THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_XUNPACK(bt1,tmp11,info(FVYX),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt1,g(1)%size,nyz_trans,g(1)%size,tmp11,nyz_trans)
           CALL OPR_BURGERS(0, nyz_trans, bcs, g(1), tmp11,tmpu,tmp12, wrk2d,wrk3d)
           CALL DNS_TRANSPOSE(tmp12,nyz_trans,g(1)%size,nyz_trans,bt1,g(1)%size)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt1,bt1,tmp12,tmp11,info(BVXY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(FVYZ), t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_ZUNPACK(tmp21,bt2,info(FVYZ),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL OPR_BURGERS(0, nxy_trans, bcs, g(3), tmp21,tmpw,bt2, wrk2d,wrk3d)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_ZYCOMM(bt2,bt2,tmp22,tmp21,info(BVZY),t_tmp);t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(F1YX), t_tmp) )  THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_XUNPACK(bt3,tmp31,info(F1YX),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt3,g(1)%size,nyz_trans,g(1)%size,tmp31,nyz_trans)
           CALL OPR_BURGERS(0, nyz_trans, bcs, g(1), tmp31,tmpu,tmp32, wrk2d,wrk3d)
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,g(1)%size,nyz_trans,bt3,g(1)%size)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt3,bt3,tmp32,tmp31,info(B1XY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(F1YZ), t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_ZUNPACK(tmp41,bt4,info(F1YZ),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL OPR_BURGERS(0, nxy_trans, bcs, g(3), tmp41,tmpw,bt4, wrk2d,wrk3d)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_ZYCOMM(bt4,bt4,tmp42,tmp41,info(B1ZY),t_tmp);t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
     ENDDO
     !
     DO WHILE ( finished /= 16 )
        IF ( NB3DFFT_R2R_READY(info(BVXY),t_tmp) ) THEN
           t_test=t_test+t_tmp;
           CALL NB3DFFT_R2R_Y1UNPACK(bt1,tmp11,info(BVXY),t_tmp);     t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h2=h2+bt1
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           IF ( inb_scal .GT. 1 ) &
                CALL NB3DFFT_R2R_YXCOMM(s(:,2),bt1,bt1,tmp11,info(F2YX),t_tmp);
           t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(BVZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt2,tmp21,info(BVZY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           h2=h2+bt2
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           IF ( inb_scal .GT. 1 ) &
                CALL NB3DFFT_R2R_YZCOMM(s(:,2),tmp21,tmp21,bt2,info(F2YZ),t_tmp);
           t_comp=t_comp+t_tmp
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(B1XY),t_tmp) ) THEN
           t_test=t_test+t_tmp;
           CALL NB3DFFT_R2R_Y1UNPACK(bt3,tmp31,info(B1XY),t_tmp);     t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           hs(:,1)=hs(:,1)+bt3
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(B1ZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt4,tmp41,info(B1ZY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           hs(:,1)=hs(:,1)+bt4
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           finished = finished+1
        ENDIF
     ENDDO
     ! u,v,w,s1 are finished and s2 initiated
     ! we can prepare the pressure solver, and update tendencies

     t_tmp = -MPI_WTime()
     !
     ! Source terms
     !
     CALL FI_SOURCES_FLOW(u,s, h1, tmp31,       wrk1d,wrk2d,wrk3d)
     CALL FI_SOURCES_SCAL(  s, hs, tmp31,tmp32, wrk1d,wrk2d,wrk3d)
     !
     ! Impose buffer zone as relaxation terms
     !
     IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
        CALL BOUNDARY_BUFFER_RELAX_FLOW(u,h1)
     ENDIF
     !
     ! Calculate divergence for pressure solver
     !
     tdummy = C_1_R / dte
     tmp32 = h1 + u*tdummy
     tmp42 = h3 + w*tdummy
     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp32)
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp42)
     ENDIF
     t_ser = t_ser + (t_tmp+MPI_WTime())

     CALL NB3DFFT_R2R_YXCOMM(tmp32,bt3,  bt3,  tmp31,info(FPYX),t_tmp);t_comp=t_comp+t_tmp
     CALL NB3DFFT_R2R_YZCOMM(tmp42,tmp41,tmp41,bt4,  info(FPYZ),t_tmp);t_comp=t_comp+t_tmp

     ! Oy source term for pressure solver below
     !
     IF ( inb_scal .LT. 2 ) &
          finished = finished + 2

     DO WHILE ( finished /= 20 )
        IF ( inb_scal .GT. 1 ) THEN
           IF ( NB3DFFT_R2R_READY(info(F2YX), t_tmp) )  THEN
              t_test=t_test+t_tmp
              CALL NB3DFFT_R2R_XUNPACK(bt1,tmp11,info(F2YX),t_tmp); t_comp=t_comp+t_tmp;
              !
              t_tmp = -MPI_WTime()
              CALL DNS_TRANSPOSE(bt1,g(1)%size,nyz_trans,g(1)%size,tmp11,nyz_trans)
              CALL OPR_BURGERS(0, nyz_trans, bcs, g(1), tmp11,tmpu,tmp12, wrk2d,wrk3d)
              CALL DNS_TRANSPOSE(tmp12,nyz_trans,g(1)%size,nyz_trans,bt1,g(1)%size)
              t_ser = t_ser + (t_tmp+MPI_WTime())
              !
              CALL NB3DFFT_R2R_XYCOMM(bt1,bt1,tmp12,tmp11,info(B2XY),t_tmp); t_comp=t_comp+t_tmp;
              finished = finished+1
           ENDIF
           IF ( NB3DFFT_R2R_READY(info(F2YZ), t_tmp) ) THEN
              t_test=t_test+t_tmp
              CALL NB3DFFT_R2R_ZUNPACK(tmp21,bt2,info(F2YZ),t_tmp); t_comp=t_comp+t_tmp;
              !
              t_tmp = -MPI_WTime()
              CALL OPR_BURGERS(0, nxy_trans, bcs, g(3), tmp21,tmpw,bt2, wrk2d,wrk3d)
              t_ser = t_ser + (t_tmp+MPI_WTime())
              !
              CALL NB3DFFT_R2R_ZYCOMM(bt2,bt2,tmp22,tmp21,info(B2ZY),t_tmp);t_comp=t_comp+t_tmp;
              finished = finished+1
           ENDIF
        ENDIF
        !
        IF ( NB3DFFT_R2R_READY(info(FPYX), t_tmp) )  THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_XUNPACK(bt3,tmp31,info(FPYX),t_tmp); t_comp=t_comp+t_tmp;
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt3,g(1)%size,nyz_trans,g(1)%size,tmp31,nyz_trans)
           CALL OPR_PARTIAL1(nyz_trans, bcs, g(1), tmp31,tmp32, wrk2d)
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,g(1)%size,nyz_trans,bt3,g(1)%size)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt3,bt3,tmp32,tmp31,info(BPXY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(FPYZ), t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_ZUNPACK(tmp41,bt4,info(FPYZ),t_tmp); t_comp=t_comp+t_tmp;
           t_tmp = -MPI_WTime()
           CALL OPR_PARTIAL1(nxy_trans, bcs, g(3), tmp41,bt4, wrk2d)
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_ZYCOMM(bt4,bt4,tmp42,tmp41,info(BPZY),t_tmp);t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
     ENDDO
     !
     DO WHILE ( finished /= 22 )
        IF ( NB3DFFT_R2R_READY(info(B2XY),t_tmp) ) THEN
           t_test=t_test+t_tmp;
           CALL NB3DFFT_R2R_Y1UNPACK(bt1,tmp11,info(B2XY),t_tmp);     t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           hs(:,2)=hs(:,2)+bt1
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(B2ZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt2,tmp21,info(B2ZY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           hs(:,2)=hs(:,2)+bt2
           t_ser = t_ser + t_tmp + MPI_WTime()
           !
           finished = finished+1
        ENDIF
     ENDDO
     !
     t_tmp = -MPI_WTime()
     tdummy=C_1_R/dte
     tmp11=h2+v*tdummy
     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp11)
     ENDIF
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp11,tmp12, wrk3d, wrk2d,wrk3d)
     t_ser = t_ser + (t_tmp+MPI_WTime())
     !
     DO WHILE ( finished /= 24 )
        IF ( NB3DFFT_R2R_READY(info(BPXY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y1UNPACK(bt3,tmp31,info(BPXY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           tmp12=tmp12+bt3
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(BPZY),t_tmp) ) THEN
           t_test=t_test+t_tmp
           CALL NB3DFFT_R2R_Y2UNPACK(bt4,tmp41,info(BPZY),t_tmp);      t_comp=t_comp+t_tmp
           !
           t_tmp = -MPI_WTime()
           tmp12=tmp12+bt4
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           finished = finished+1
        ENDIF
     ENDDO
     CALL NB3DFFT_NBC_SCHEDL_STOP(nbcsetup_)
  ELSE
     PRTERR1("there MUST NOT be more threads than 2 on first level")
  ENDIF
!$omp end parallel
  t_run = t_run + MPI_WTime()
  !               rhs-stuff    testing    packing     init
  t_wait= t_run - t_ser        - t_test   - t_comp  - t_init
  CALL NB3DFFT_NBC_FINISH(nbcsetup_,t_run,pkg_cnt)
  ! CALL NB3DFFT_NBC_FINISH()
  nbcsetup = nbcsetup_

  t_snd=(/t_run,t_ser,t_comp,t_test,t_init,t_wait/)
  CALL MPI_Reduce(t_snd,t_rcv,6,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ims_err)

#ifdef USE_PROFILE
  IF ( ims_pro .EQ. 0 ) THEN
     t_rcv=t_rcv/ims_npro
     t_run=t_rcv(1);  t_ser= t_rcv(2); t_comp=t_rcv(3);
     t_test=t_rcv(4); t_init=t_rcv(5); t_wait=t_rcv(6);
     WRITE(*,903)  ims_npro,t_run, &
          t_ser/t_run, t_comp/t_run, t_test/t_run, t_init/t_run, t_wait/t_run
903  FORMAT(&
          'NBC-time nproc', i6, &
          '(total [s], serial[1], comp[1], test[1], init[1], wait [1]);', &
          F10.3, 5(F6.3,';'))
  ENDIF
#endif

  ptime = -MPI_WTime()

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     p_bcs => h2(ip_b:); BcsFlowJmin%ref(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => h2(ip_t:); BcsFlowJmax%ref(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! Adding density in BCs
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     BcsFlowJmin%ref(:,:,2) = BcsFlowJmin%ref(:,:,2) *rbackground(1)
     BcsFlowJmax%ref(:,:,2) = BcsFlowJmax%ref(:,:,2) *rbackground(g(2)%size)
  ENDIF

! pressure in tmp12, Oy derivative in tmp11
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
       tmp12,tmp11, tmp41,tmp42, BcsFlowJmin%ref(1,1,2),BcsFlowJmax%ref(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  IF ( use_tower .AND. rkm_substep .EQ. rkm_endstep ) THEN
     CALL DNS_TOWER_ACCUMULATE(tmp12,i4,wrk1d)
  ENDIF

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp12,tmp41, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp12,tmp42, wrk3d, wrk2d,wrk3d)

  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp41, h1)
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp11, h2)
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp42, h3)

  ELSE
     h1 = h1 - tmp41
     h2 = h2 - tmp11
     h3 = h3 - tmp42

  ENDIF

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  BcsFlowJmin%ref = C_0_R ! default is no-slip (dirichlet)
  BcsFlowJmax%ref = C_0_R

  DO iq = 1,inb_flow
     ibc = 0
     IF ( BcsFlowJmin%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsFlowJmax%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( iq .EQ. 1 ) p_h => h1(1:)
     IF ( iq .EQ. 2 ) p_h => h2(1:)
     IF ( iq .EQ. 3 ) p_h => h3(1:)
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), p_h, &
             BcsFlowJmin%ref(1,1,iq),BcsFlowJmax%ref(1,1,iq), wrk1d,tmp11,wrk3d)
     ENDIF
  ENDDO

  DO is = 1,inb_scal
     ibc = 0
     IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs(1,is), &
             BcsScalJmin%ref(1,1,is),BcsScalJmax%ref(1,1,is), wrk1d,tmp11,wrk3d)
     ENDIF

     IF ( BcsScalJmin%type(is) .NE. DNS_SFC_STATIC .OR. &
          BcsScalJmax%type(is) .NE. DNS_SFC_STATIC ) THEN
        CALL BOUNDARY_SURFACE_J(is,bcs,s,hs,tmp11,tmp12,tmp21,wrk1d,wrk2d,wrk3d)
     ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     h1(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,2)
     h3(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,3)
     DO is = 1,inb_scal
        hs(ip_b:ip_b+imax-1,is) = BcsScalJmin%ref(1:imax,k,is)
     ENDDO
     ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     h1(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,2)
     h3(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,3)
     DO is = 1,inb_scal
        hs(ip_t:ip_t+imax-1,is) = BcsScalJmax%ref(1:imax,k,is)
     ENDDO
     ip_t = ip_t + nxy
  ENDDO

  ptime = ptime + MPI_WTime()
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC')
#endif
  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC

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

SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC(dte,etime,x,y,z,dx,dy,dz,&
     u,v,w,s,&
     tmpu,tmpw,tmp11,tmp12,tmp21,tmp22,tmp31,tmp32,tmp41,tmp42,&
     bt1,bt2,bt3,bt4,&
     h1,h2,h3,hs,&
     bcs_hb,bcs_ht,b_ref,vaux,&
     wrk1d,wrk2d,wrk3d) 
  USE, INTRINSIC :: iso_c_binding, ONLY : c_int,c_loc,c_ptr,c_f_pointer 

  USE OMP_LIB,    ONLY : omp_get_thread_num 

  USE DNS_CONSTANTS, ONLY : lfile,wfile,efile
  !
  USE DNS_GLOBAL, ONLY : g, i1bc,iunifx,j1bc,iunify,k1bc
  USE DNS_GLOBAL, ONLY : imax_total,jmax_total,kmax_total,imode_fdm
  USE DNS_GLOBAL, ONLY : iunifz,inb_flow,inb_vars,inb_scal,inb_scal_array,visc,schmidt,prandtl 
  USE DNS_GLOBAL, ONLY : isize_field, isize_wrk1d, imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : rotn_param,rotn_vector,body_param,body_vector 
  USE DNS_GLOBAL, ONLY : icoriolis 
  ! 
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE DNS_LOCAL,  ONLY : bcs_scal_jmin, bcs_scal_jmax
  USE DNS_LOCAL,  ONLY : idivergence
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type 
  USE DNS_LOCAL,  ONLY : rkm_substep,rkm_endstep,tower_mode 

  USE DNS_TOWER

#ifdef USE_PSFFT 
  USE DNS_LOCAL,  ONLY : nbcsetup 
#endif 

  USE DNS_MPI,    ONLY : ims_npro, ims_pro, ims_err,ims_size_i,ims_size_k  

  USE NB3DFFT,    ONLY : nb3dfft_nbc_prepare,nb3dfft_nbc_finish,nb3dfft_infoType
  USE NB3DFFT,    ONLY : nb3dfft_nbc_schedl_start, nb3dfft_nbc_worker_start
  USE NB3DFFT,    ONLY : nb3dfft_nbc_schedl_stop,  nb3dfft_schedlType
  USE NB3DFFT,    ONLY : nb3dfft_r2r_yxcomm,       nb3dfft_r2r_yzcomm
  USE NB3DFFT,    ONLY : nb3dfft_r2r_xycomm,       nb3dfft_r2r_zycomm
  USE NB3DFFT,    ONLY : nb3dfft_r2r_xunpack,      nb3dfft_r2r_zunpack
  USE NB3DFFT,    ONLY : nb3dfft_r2r_y1unpack,     nb3dfft_r2r_y2unpack
  USE NB3DFFT,    ONLY : nb3dfft_r2r_ready,        mytype
  
  IMPLICIT NONE

#include "integers.h"
#include "mpif.h" 
  !
  ! PARAMETERS 
  ! 
  TINTEGER, PARAMETER :: nmeasure=3
  
  TREAL,                                 INTENT(IN) :: dte,etime
  
  TREAL, DIMENSION(*),                   INTENT(IN) :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field),         INTENT(IN) :: u,v,w
  TREAL, DIMENSION(isize_field,inb_scal_array),INTENT(IN) :: s

  TREAL, DIMENSION(isize_field),         INTENT(INOUT):: h1,h2,h3 
  TREAL, DIMENSION(isize_field,inb_scal),INTENT(OUT):: hs 
  TREAL, DIMENSION(imax,kmax,inb_vars)              :: bcs_hb, bcs_ht 
  TREAL, DIMENSION(jmax),                INTENT(IN) :: b_ref
  TREAL, DIMENSION(isize_field),         INTENT(INOUT):: tmpu,tmpw,tmp11,tmp12,tmp21,tmp22,tmp31,tmp32,tmp41,tmp42
  TREAL, DIMENSION(isize_field)  :: bt1,bt2,bt3,bt4
  TREAL, DIMENSION(isize_wrk1d,*):: wrk1d
  TREAL, DIMENSION(*)            :: wrk2d,wrk3d,vaux
  !
  ! LOCAL VARIABLES 
  ! 
  TINTEGER :: nxy_trans,nyz_trans,nxy,id,imeasure,ij,k,is,commID
  TINTEGER :: finished,ip_b,ip_t,ibc
  TREAL tdummy!,bdummy,fdummy,u_geo,w_geo 
  TREAL, DIMENSION(:), POINTER :: p_bcs 
  !
  TREAL, DIMENSION(inb_scal)      :: err_s, diff  
  ! 
  TYPE(nb3dfft_infoType), DIMENSION(24) :: info
  INTEGER(KIND=4),DIMENSION(2) :: cur_time
  REAL(KIND=mytype) :: t_comp,t_test,t_ser,t_init,t_wait,t_tmp 
  REAL(KIND=8),DIMENSION(6)      :: t_snd, t_rcv
  TYPE(nb3dfft_schedlType) :: nbcsetup_
  INTEGER :: pkg_cnt
  TREAL err_x,err_y,err_z,kx,ky,kz,dum 
  TREAL rtime, rtime_loc,t_run,ptime,ctime_loc

  TARGET h2, tmp42 

  IF ( inb_scal .GT. 2 ) THEN   
     CALL IO_WRITE_ASCII(efile,&
          'Nonblocking Communication not implemented >2 scalars' )
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
  ELSE IF ( inb_scal .LT. 1 ) THEN 
     CALL IO_WRITE_ASCII(efile,&
          'Nonblocking Communication require at least 1 scalar')  
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  ! IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN .OR.&
  !      bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN 
  !    CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
  ! ENDIF

  nbcsetup_ = nbcsetup
  pkg_cnt=24*ims_npro
  t_comp=0
  t_test=0
  t_ser =0
  ptime =0

  nxy=imax*jmax 
  ! u_geo = COS(rotn_param(1))
  ! w_geo =-SIN(rotn_param(1)) 
  ! IF ( icoriolis .EQ. EQNS_COR_NORMALIZED ) THEN 
  !    fdummy= rotn_vector(2) 
  ! ELSE 
  !    fdummy = C_0_R 
  ! ENDIF

  DO is=1,inb_scal
     diff(is) = visc*prandtl*schmidt(is) 
  ENDDO

  rtime = -MPI_WTime()  
  t_run = rtime
  t_init= rtime  

! #######################################################################
! Source terms
! #######################################################################
  ! CALL FI_SOURCES_FLOW(u,s, h1, b_ref, wrk1d,wrk3d)
  ! CALL FI_SOURCES_SCAL(y,dy, s, hs, tmp11,tmp12, wrk1d,wrk2d,wrk3d)

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

     id = DNS_MPI_I_PARTIAL;   nyz_trans = ims_size_i(id)  
     id = DNS_MPI_K_PARTIAL;   nxy_trans = ims_size_k(id)   
     !
     ! kick off transpose U y->x and W y->z
     CALL NB3DFFT_R2R_YXCOMM(u,bt1, bt1, tmp11,info(FUYX),t_tmp);  t_comp=t_comp+t_tmp 
     CALL NB3DFFT_R2R_YZCOMM(w,tmpw,tmpw,bt2,  info(FWYZ),t_tmp);  t_comp=t_comp+t_tmp  
     CALL NB3DFFT_R2R_YXCOMM(w,bt3,  bt3,  tmp31,info(FWYX),t_tmp);t_comp=t_comp+t_tmp   
     !
     ! Coriolis Force, Vertical derivatives, and Vertical advection  
     !   
     t_tmp = -MPI_WTime()
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc, & 
          dy, u, tmp41, i0,i0, i0,i0, tmp42, wrk1d,wrk2d,wrk3d)  
!     h1=h1+visc*tmp41 - v*tmp42 + fdummy*(w_geo-w)
     h1=h1+visc*tmp41 - v*tmp42 !+ fdummy*(w_geo-w)

     ! -----------------------------------------------------------------------
     ! BCs s.t. derivative d/dy(u) is zero
     ! -----------------------------------------------------------------------
     IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
        ip_b =                 1
        DO k = 1,kmax
           p_bcs => tmp42(ip_b:); bcs_hb(1:imax,k,1) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        ENDDO
     ENDIF
     IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
        ip_t = imax*(jmax-1) + 1
        DO k = 1,kmax
           p_bcs => tmp42(ip_t:); bcs_ht(1:imax,k,1) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
        ENDDO
     ENDIF

     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc, & 
          dy, v, tmp41, i0,i0, i0,i0, tmp42, wrk1d,wrk2d,wrk3d)  
     h2=h2+visc*tmp41 - v*tmp42
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc, & 
          dy, w, tmp41, i0,i0, i0,i0, tmp42, wrk1d,wrk2d,wrk3d)  
!     h3=h3+visc*tmp41 - v*tmp42 + fdummy*(u-u_geo)
     h3=h3+visc*tmp41 - v*tmp42 !+ fdummy*(u-u_geo)

     ! -----------------------------------------------------------------------
     ! BCs s.t. derivative d/dy(w) is zero
     ! -----------------------------------------------------------------------
     IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
        ip_b =                 1
        DO k = 1,kmax
           p_bcs => tmp42(ip_b:); bcs_hb(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        ENDDO
     ENDIF
     IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
        ip_t = imax*(jmax-1) + 1
        DO k = 1,kmax
           p_bcs => tmp42(ip_t:); bcs_ht(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
        ENDDO
     ENDIF

     CALL FI_SOURCES_FLOW(u,s, h1, b_ref, wrk1d,wrk3d)
     CALL FI_SOURCES_SCAL(y,dy, s, hs, tmp41,tmp42, wrk1d,wrk2d,wrk3d)

     t_ser    = t_ser + (t_tmp +MPI_WTime())
     !
     CALL NB3DFFT_R2R_YZCOMM(u,tmp41,tmp41,bt4,  info(FUYZ),t_tmp);t_comp=t_comp+t_tmp  

     finished = 0
     DO WHILE ( finished /= 2 )   
        IF ( NB3DFFT_R2R_READY(info(FUYX), t_tmp ) ) THEN     
           t_test = t_test + t_tmp
           ! u du/dx + 1/Re d2u/dx2
           CALL NB3DFFT_R2R_XUNPACK(bt1,tmp11, info(FUYX), t_tmp);t_comp=t_comp+t_tmp    
           !
           t_tmp = -MPI_WTime() 
           CALL DNS_TRANSPOSE(bt1,imax_total,nyz_trans,imax_total,tmpu,nyz_trans)
!           CALL OPR_BURGERS(iunifx,imode_fdm,0,nyz_trans,imax_total,i1bc, dx, &
           CALL OPR_BURGERS(imode_fdm, 0, nyz_trans, g(1), &
                tmpu,tmpu,0,0,0,0,tmp11,wrk2d,wrk3d)  
           CALL DNS_TRANSPOSE(tmp11,nyz_trans,imax_total,nyz_trans,bt1,imax_total) 
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
!           CALL OPR_BURGERS(iunifz,imode_fdm,0,nxy_trans,kmax_total,k1bc,dz, & 
           CALL OPR_BURGERS(imode_fdm, 0, nxy_trans, g(3), & 
                tmpw,tmpw,0,0,0,0,bt2,wrk2d,wrk3d)     
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
           CALL NB3DFFT_R2R_XUNPACK(bt3,tmp31,info(FWYX),t_tmp); t_comp=t_comp+t_tmp;  
           !
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt3,imax_total,nyz_trans,imax_total,tmp31,nyz_trans) 
!           CALL OPR_BURGERS(iunifx,imode_fdm,0,nyz_trans,imax_total,i1bc,dx, &
           CALL OPR_BURGERS(imode_fdm, 0, nyz_trans, g(1), &
                tmp31,tmpu,0,0,0,0,tmp32,wrk2d,wrk3d) 
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,imax_total,nyz_trans,bt3,imax_total)   
           t_ser = t_ser + (t_tmp + MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt3,bt3,tmp32,tmp31,info(BWXY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(FUYZ), t_tmp) ) THEN  
           t_test=t_test+t_tmp  
           CALL NB3DFFT_R2R_ZUNPACK(tmp41,bt4,info(FUYZ),t_tmp); t_comp=t_comp+t_tmp;  
           ! 
           t_tmp = -MPI_WTime() 
!           CALL OPR_BURGERS(iunifz,imode_fdm,0,nxy_trans,kmax_total,k1bc,dz, &
           CALL OPR_BURGERS(imode_fdm, 0, nxy_trans, g(3), &
                tmp41,tmpw,0,0,0,0,bt4,wrk2d,wrk3d)   
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
     ! u and w are finished 
     !
     DO WHILE ( finished /= 12 ) 
        IF ( NB3DFFT_R2R_READY(info(FVYX), t_tmp) )  THEN
           t_test=t_test+t_tmp 
           CALL NB3DFFT_R2R_XUNPACK(bt1,tmp11,info(FVYX),t_tmp); t_comp=t_comp+t_tmp;  
           ! 
           t_tmp = -MPI_WTime()
           CALL DNS_TRANSPOSE(bt1,imax_total,nyz_trans,imax_total,tmp11,nyz_trans) 
!           CALL OPR_BURGERS(iunifx,imode_fdm,0,nyz_trans,imax_total,i1bc,dx, &
           CALL OPR_BURGERS(imode_fdm, 0, nyz_trans, g(1), &
                tmp11,tmpu,0,0,0,0,tmp12,wrk2d,wrk3d) 
           CALL DNS_TRANSPOSE(tmp12,nyz_trans,imax_total,nyz_trans,bt1,imax_total)  
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
!           CALL OPR_BURGERS(iunifz,imode_fdm,0,nxy_trans,kmax_total,k1bc,dz, &
           CALL OPR_BURGERS(imode_fdm, 0, nxy_trans, g(3), &
                tmp21,tmpw,0,0,0,0,bt2,wrk2d,wrk3d)   
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
           CALL DNS_TRANSPOSE(bt3,imax_total,nyz_trans,imax_total,tmp31,nyz_trans) 
!           CALL OPR_BURGERS(iunifx,imode_fdm,0,nyz_trans,imax_total,i1bc,dx, &
           CALL OPR_BURGERS(imode_fdm, 0, nyz_trans, g(1), &
                tmp31,tmpu,0,0,0,0,tmp32,wrk2d,wrk3d) 
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,imax_total,nyz_trans,bt3,imax_total) 
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
!           CALL OPR_BURGERS(iunifz,imode_fdm,0,nxy_trans,kmax_total,k1bc,dz, &
           CALL OPR_BURGERS(imode_fdm, 0, nxy_trans, g(3), &
                tmp41,tmpw,0,0,0,0,bt4,wrk2d,wrk3d)   
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
     ! u,v,w,s1  finished 
     ! we can prepare the pressure solver, and update tendencies 
     !
     t_tmp = -MPI_WTime() 
     CALL PARTIAL_YY(i1,iunify,imode_fdm,imax,jmax,kmax, j1bc, & 
          dy, s(1,1), tmp31, i0,i0,i0,i0,tmp32, wrk1d,wrk2d,wrk3d) 
     hs(:,1) = hs(:,1) + diff(1)*tmp31-v*tmp32 
     !
     IF ( inb_scal .GT. 1 ) THEN
     CALL PARTIAL_YY(i1,iunify,imode_fdm,imax,jmax,kmax, j1bc, & 
          dy, s(1,2), tmp31, i0,i0,i0,i0,tmp32, wrk1d,wrk2d,wrk3d) 
     hs(:,2) = hs(:,2) + diff(2)*tmp31-v*tmp32 
     ENDIF

     ! #######################################################################
     ! Impose buffer zone as relaxation terms
     ! #######################################################################
     IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
        CALL BOUNDARY_BUFFER_RELAXATION_FLOW(&
             vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
             vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), x,y, u,h1)
     ENDIF

     tdummy = C_1_R / dte   
     ! #######################################################################
     ! Calculate divergence for pressure solver 
     ! #######################################################################
     tmp32 = h1 + u*tdummy 
     tmp42 = h3 + w*tdummy 
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
              CALL DNS_TRANSPOSE(bt1,imax_total,nyz_trans,imax_total,tmp11,nyz_trans) 
!              CALL OPR_BURGERS(iunifx,imode_fdm,0,nyz_trans,imax_total,i1bc,dx, &
              CALL OPR_BURGERS(imode_fdm, 0, nyz_trans, g(1), &
                   tmp11,tmpu,0,0,0,0,tmp12,wrk2d,wrk3d) 
              CALL DNS_TRANSPOSE(tmp12,nyz_trans,imax_total,nyz_trans,bt1,imax_total) 
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
!              CALL OPR_BURGERS(iunifz,imode_fdm,0,nxy_trans,kmax_total,k1bc,dz, &
              CALL OPR_BURGERS(imode_fdm, 0, nxy_trans, g(3), &
                   tmp21,tmpw,0,0,0,0,bt2,wrk2d,wrk3d)  
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
           CALL DNS_TRANSPOSE(bt3,imax_total,nyz_trans,imax_total,tmp31,nyz_trans) 
           CALL PARTIAL(imode_fdm,nyz_trans,imax_total,i1bc,dx, &
                tmp31,tmp32,0,0,wrk1d,wrk2d,wrk3d) 
           CALL DNS_TRANSPOSE(tmp32,nyz_trans,imax_total,nyz_trans,bt3,imax_total) 
           t_ser = t_ser + (t_tmp+MPI_WTime())
           !
           CALL NB3DFFT_R2R_XYCOMM(bt3,bt3,tmp32,tmp31,info(BPXY),t_tmp); t_comp=t_comp+t_tmp;
           finished = finished+1
        ENDIF
        IF ( NB3DFFT_R2R_READY(info(FPYZ), t_tmp) ) THEN  
           t_test=t_test+t_tmp  
           CALL NB3DFFT_R2R_ZUNPACK(tmp41,bt4,info(FPYZ),t_tmp); t_comp=t_comp+t_tmp;  
           t_tmp = -MPI_WTime()
           CALL PARTIAL(imode_fdm,nxy_trans,kmax_total,k1bc,dz, &
                tmp41,bt4,0,0,wrk1d,wrk2d,wrk3d)  
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
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, &
          tmp11, tmp12, i0,i0, wrk1d,wrk2d,wrk3d) 
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
     p_bcs => h2(ip_b:); bcs_hb(1:imax,k,3) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => h2(ip_t:); bcs_ht(1:imax,k,3) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp12, Oy derivative in tmp11
  CALL OPR_POISSON_FXZ(imode_fdm,i2,i3, imax,jmax,kmax,  &
       y,dx,dy,dz, tmp12,tmp11, tmp41,tmp42, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

  IF ( tower_mode .EQ. 1 .AND. rkm_substep .EQ. rkm_endstep ) THEN 
     CALL DNS_TOWER_ACCUMULATE(tmp12,i4,dx,dy,dz,wrk1d) 
  ENDIF


  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i0,dx,tmp12,tmp41,i0,i0,wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, i0,dz,tmp12,tmp42,i0,i0,wrk1d,wrk2d,wrk3d) 
  
  h1 = h1 - tmp41 
  h2 = h2 - tmp11
  h3 = h3 - tmp42 

  bcs_hb(:,:,1:inb_vars) = C_0_R  ! default is no-slip 
  bcs_ht(:,:,1:inb_vars) = C_0_R
  

! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  ibc = 0
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h1, &
          bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp11,wrk3d)
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h3, &
          bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp11,wrk3d)
  ENDIF

  DO is = 1,inb_scal
  ibc = 0
  IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, hs(1,is), &
          bcs_hb(1,1,is+inb_flow),bcs_ht(1,1,is+inb_flow), wrk1d,tmp11,wrk3d)
  ENDIF
  ENDDO


  ip_b = 1 
  DO k=1,kmax 
     DO is=1,inb_scal 
        hs(ip_b:ip_b+imax-1,is) = bcs_hb(1:imax,k,is+inb_flow)
     ENDDO
     h1(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = C_0_R ! no penetration
     h3(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2);  
     ip_b=ip_b+nxy 
  ENDDO

  ip_t = imax*(jmax-1)+1
  DO k=1,kmax  
     DO is=1,inb_scal 
        hs(ip_t:ip_t+imax-1,is) = bcs_ht(1:imax,k,is+inb_flow)
     ENDDO
     h1(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = C_0_R ! no penetration
     h3(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2);  
     ip_t = ip_t + nxy 
  ENDDO

  ptime = ptime + MPI_WTime()
!
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_NBC

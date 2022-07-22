#ifdef IBM_DEBUG
#include "types.h"
#ifdef USE_MPI 
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   debug IO for all relevant geometry fields for IBM routines
!#     nobi,    nobj,   nobk   : number of objects in i/j/k 
!#     nobi_b,  nobj_b, nobk_b : beginn of objects in i/j/k 
!#     nobi_e,  nobj_e, nobk_e : end    of objects in i/j/k
!#    
!# 
!########################################################################
!# ARGUMENTS                         
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!#
!########################################################################

subroutine IBM_GEOMETRY_DEBUG_IO(epsi, epsj, epsk, tmp1, tmp2, tmp3, wrk3d)
  
  use IBM_VARS
  use IO_FIELDS
  use TLAB_VARS,     only : g, imax, jmax, kmax, isize_field
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_PROCS
  use TLAB_MPI_VARS, only : ims_pro
  use TLAB_MPI_VARS, only : ims_size_i, ims_size_k    
  use TLAB_MPI_VARS, only : ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i   
  use TLAB_MPI_VARS, only : ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k   
  use TLAB_MPI_VARS, only : ims_npro_i, ims_npro_k , ims_err      
#endif  
  
  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field), intent(in   ) :: epsi, epsj, epsk
  TREAL, dimension(isize_field), intent(inout) :: tmp1, tmp2, tmp3, wrk3d
  
#ifdef USE_MPI 
  TINTEGER, parameter                          :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter                          :: idk = TLAB_MPI_K_PARTIAL 
#endif
  TINTEGER                                     :: i, j, k, ij, ik, jk, ip, inum
  TINTEGER                                     :: nyz, nxz, nxy
#ifdef USE_MPI 
#else
  TINTEGER, parameter                          :: ims_pro = 0         
#endif

  ! ================================================================== !
  ! npages
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    nyz = ims_size_i(idi)
  else
#endif
  nyz = jmax * kmax 
#ifdef USE_MPI
  end if
#endif

  nxz = imax * kmax     

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    nxy = ims_size_k(idk)
  else
#endif
  nxy = imax * jmax
#ifdef USE_MPI
  end if
#endif

  tmp1(:) = C_0_R; tmp2(:) = C_0_R; tmp3(:) = C_0_R

  ! ================================================================== !
  ! number of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((ip == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for objects
        tmp1(jk) = C_1_R
      end if 
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! check for interface 
        tmp1(jk) = tmp1(jk) + C_1_R
      end if
    end do
    ip = ip + nyz
  end do

  call DNS_TRANSPOSE(tmp1, nyz, g(1)%size, nyz,        tmp2, g(1)%size)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(tmp2, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  endif
  call IO_WRITE_FIELDS('nobi3d',  IO_FLOW, imax,jmax,kmax, i1, tmp1, wrk3d)
#else
  call IO_WRITE_FIELDS('nobi3d',  IO_FLOW, imax,jmax,kmax, i1, tmp2, wrk3d)
#endif
  tmp1(:) = C_0_R; tmp2(:) = C_0_R

  ! ================================================================== !
  ! number of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((ip == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for objects
        tmp1(ik) = C_1_R
      end if 
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! check for interface 
        tmp1(ik) = tmp1(ik) + C_1_R
      end if
    end do
    ip = ip + nxz
  end do

  call DNS_TRANSPOSE(tmp1, kmax, imax * jmax, kmax, tmp2, imax * jmax)
  call IO_WRITE_FIELDS('nobj3d',  IO_FLOW, imax,jmax,kmax, i1, tmp2, wrk3d)
  tmp1(:) = C_0_R; tmp2(:) = C_0_R

  ! ================================================================== !
  ! number of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
        tmp1(ij) = C_1_R
      end if 
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
        tmp1(ij) = tmp1(ij) + C_1_R
      end if
    end do
    ip = ip + nxy
  end do

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(tmp1, tmp2, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  endif
  call IO_WRITE_FIELDS('nobk3d',  IO_FLOW, imax,jmax,kmax, i1, tmp2, wrk3d)
#else
  call IO_WRITE_FIELDS('nobk3d',  IO_FLOW, imax,jmax,kmax, i1, tmp1, wrk3d)
#endif
  if (ims_pro == 0) write(*,*) '============= Writing all geometry fields ==============='
  if (ims_pro == 0) write(*,*) 'done writing files: nobi3d,   nobj3d,   nobk3d'
  tmp1(:) = C_0_R; tmp2(:) = C_0_R
  
  ! ================================================================== !
  ! begin and end of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((i == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for interface
        tmp1(jk) = dble(i) ! nobi_b
      end if
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then     ! nobi_b check for interface 
        inum = i0
        do while (tmp1(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp1(inum+jk) = dble(i + 1)
      elseif((epsi(ip+jk-1) == C_1_R) .and. (epsi(ip+jk-1+nyz) == C_0_R)) then ! nobi_e check for interface 
        inum = i0
        do while (tmp2(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp2(inum+jk) = dble(i)        
      end if
      if((i == (g(1)%size - 1)) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp2(inum+jk) = dble(g(1)%size)    
      end if
    end do
    ip = ip + nyz
  end do

  call DNS_TRANSPOSE(tmp1, nyz, g(1)%size, nyz,        tmp3, g(1)%size)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(tmp3, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  endif
  call IO_WRITE_FIELDS('nobi3d_b',  IO_FLOW, imax,jmax,kmax, i1, tmp1, wrk3d)
#else
  call IO_WRITE_FIELDS('nobi3d_b',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
#endif

  call DNS_TRANSPOSE(tmp2, nyz, g(1)%size, nyz,        tmp3, g(1)%size)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(tmp3, tmp2, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  endif
  call IO_WRITE_FIELDS('nobi3d_e',  IO_FLOW, imax,jmax,kmax, i1, tmp2, wrk3d)
#else
  call IO_WRITE_FIELDS('nobi3d_e',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
#endif
  if (ims_pro == 0) write(*,*) 'done writing files: nobi3d_b, nobi3d_e'
  tmp1(:) = C_0_R; tmp2(:) = C_0_R; tmp3(:) = C_0_R

  ! ================================================================== !
  ! begin and end of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((j == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for interface
        tmp1(ik) = dble(j) ! nobj_b
      end if
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then     ! nobj_b check for interface 
        inum = i0
        do while (tmp1(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        tmp1(inum+ik) = dble(j + 1)
      elseif((epsj(ip+ik-1) == C_1_R) .and. (epsj(ip+ik-1+nxz) == C_0_R)) then ! nobj_e check for interface 
        inum = i0
        do while (tmp2(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        tmp2(inum+ik) = dble(j)        
      end if
      if((j == (g(2)%size - 1)) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+ik) /= C_0_R)
          inum = inum + nxz           
        end do 
        tmp2(inum+ik) = dble(g(2)%size)    
      end if
    end do
    ip = ip + nxz
  end do

  call DNS_TRANSPOSE(tmp1, kmax, imax * jmax, kmax, tmp3, imax * jmax)
  call IO_WRITE_FIELDS('nobj3d_b',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
  call DNS_TRANSPOSE(tmp2, kmax, imax * jmax, kmax, tmp3, imax * jmax)
  call IO_WRITE_FIELDS('nobj3d_e',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
  if (ims_pro == 0) write(*,*) 'done writing files: nobj3d_b, nobj3d_e'
  tmp1(:) = C_0_R; tmp2(:) = C_0_R; tmp3(:) = C_0_R

  ! ================================================================== !
  ! begin and end of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
        tmp1(ij) = dble(k)    ! nobk_b
      end if
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
        inum = i0
        do while (tmp1(ij+inum) /= C_0_R)
          inum = inum + nxy
        end do
        tmp1(ij+inum) = dble(k + 1)
      elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
        inum = i0
        do while (tmp2(ij+inum) /= C_0_R)
          inum = inum + nxy
        end do
        tmp2(ij+inum) = dble(k)
      end if
      if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+ij) /= C_0_R)
          inum = inum + nxy           
        end do 
        tmp2(inum+ij) = dble(g(3)%size)    
      end if
    end do
    ip = ip + nxy
  end do

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(tmp1, tmp3, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  endif
  call IO_WRITE_FIELDS('nobk3d_b',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(tmp2, tmp3, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  endif
  call IO_WRITE_FIELDS('nobk3d_e',  IO_FLOW, imax,jmax,kmax, i1, tmp3, wrk3d)
#else
  call IO_WRITE_FIELDS('nobk3d_b',  IO_FLOW, imax,jmax,kmax, i1, tmp1, wrk3d)
  call IO_WRITE_FIELDS('nobk3d_e',  IO_FLOW, imax,jmax,kmax, i1, tmp2, wrk3d)
#endif
  if (ims_pro == 0) write(*,*) 'done writing files: nobk3d_b, nobk3d_e'
  if (ims_pro == 0) write(*,*) '========================================================='
  tmp1(:) = C_0_R; tmp2(:) = C_0_R; tmp3(:) = C_0_R

  return
end subroutine IBM_GEOMETRY_DEBUG_IO

!########################################################################

#endif
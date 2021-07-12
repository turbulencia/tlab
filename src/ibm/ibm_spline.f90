#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  
!#  
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################

subroutine IBM_SPLINZ(txc, wrk3d) 
  
  use DNS_IBM
  use DNS_GLOBAL, only: g   
  use DNS_GLOBAL, only: imax, jmax, kmax 
  use DNS_GLOBAL, only: isize_field, inb_txc !,isize_txc_field
  
#ifdef USE_MPI
  use DNS_MPI,    only: ims_pro, ims_npro
  use DNS_MPI,    only: ims_size_i, ims_size_j, ims_size_k    
  use DNS_MPI,    only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i 
  use DNS_MPI,    only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use DNS_MPI,    only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
#endif    
  
  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
  TINTEGER, parameter                            :: idi        = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter                            :: idj        = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter                            :: idk        = DNS_MPI_K_PARTIAL 
  TINTEGER, parameter                            :: idi_nob    = DNS_MPI_I_IBM_NOB 
  TINTEGER, parameter                            :: idk_nob    = DNS_MPI_K_IBM_NOB 
  TINTEGER, parameter                            :: idi_nob_be = DNS_MPI_I_IBM_NOB_BE 
  TINTEGER, parameter                            :: idk_nob_be = DNS_MPI_K_IBM_NOB_BE
  TINTEGER                                       :: ims_err, max_nob
  TREAL, dimension(ims_size_i(idi))              :: nobi_out    ! DEBUG
  TREAL, dimension(ims_size_j(idj))              :: nobj_out    ! DEBUG
  TREAL, dimension(ims_size_k(idk))              :: nobk_out    ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo(1)) :: nobi_b_out  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo(1)) :: nobj_b_out  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo(1)) :: nobk_b_out  ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo(1)) :: nobi_e_out  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo(1)) :: nobj_e_out  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo(1)) :: nobk_e_out  ! DEBUG
#else
  TINTEGER, parameter                            :: ims_pro=0, ims_npro=1
  TREAL, dimension(jmax * kmax)                  :: nobi_out    ! DEBUG
  TREAL, dimension(imax * kmax)                  :: nobj_out    ! DEBUG
  TREAL, dimension(imax * jmax)                  :: nobk_out    ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo(1))   :: nobi_b_out  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo(1))   :: nobj_b_out  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo(1))   :: nobk_b_out  ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo(1))   :: nobi_e_out  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo(1))   :: nobj_e_out  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo(1))   :: nobk_e_out  ! DEBUG
#endif

  TREAL, dimension(isize_field), intent(inout)   :: wrk3d 

  TINTEGER                                       :: nobi_max, nobj_max, nobk_max
  TINTEGER                                       :: i, j, k, ij, ik, jk, ip, inum
  TINTEGER                                       :: nyz, nxz, nxy
  TINTEGER                                       :: nob_max

  CHARACTER(32)                                  :: fname

  ! DEBUG
  TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc 
  TREAL, dimension(isize_field)                        :: tmp1, tmp2, tmp3, tmp4 

  ! ================================================================== !

  ! npages (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI
  nyz = ims_size_i(idi) ! local
  nxz = ims_size_j(idj) 
  nxy = ims_size_k(idk) 
#else
  nyz = jmax * kmax     ! global
  nxz = imax * kmax     
  nxy = imax * jmax     
#endif
  
!   ! ================================================================== !
!   ! number of objects in z-direction
!   ip = i1
!   do k = 1, g(3)%size - 1     ! contiguous k-lines
!     do ij = 1, nxy            ! pages of   k-lines
!       if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
!         nobk(ij) = C_1_R
!       end if 
!       if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
!         nobk(ij) = nobk(ij) + C_1_R
!       end if
!     end do
!     ip = ip + nxy
!   end do

!   nobk_max = int(maxval(nobk))
! #ifdef USE_MPI
!   max_nob = nobk_max
!   call MPI_ALLREDUCE(max_nob, nobk_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
! #endif

!   ! ------------------------------------------------------------------ !
!   ! DEBUG

! #ifdef USE_MPI
!   if ( ims_npro_k > 1 ) then
!     call DNS_MPI_TRPB_K(nobk, nobk_out, ims_ds_k(1,idk_nob), ims_dr_k(1,idk_nob), ims_ts_k(1,idk_nob), ims_tr_k(1,idk_nob))
!   endif
! #else 
!   nobk_out = nobk
! #endif

!   if (ims_pro == 0) then
!     write(*,*) 'nobk_b     =', size(nobk_b)
!     write(*,*) 'nobk_b_out =', size(nobk_b_out)
!     write(*,*) 'nobk_e     =', size(nobk_e)
!     write(*,*) 'nobk_e_out =', size(nobk_e_out)
!   end if

!   ! ! ================================================================== !
!   ! ! begin and end of objects in z-direction
!   ip = i1
!   do k = 1, g(3)%size - 1     ! contiguous k-lines
!     do ij = 1, nxy            ! pages of   k-lines
!       if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
!         nobk_b(ij) = dble(k) ! nobj_b
!       end if
!       if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
!         inum = i0
!         do while (nobk_b(inum+ij) /= C_0_R)
!           inum = inum + nxy            
!         end do 
!         nobk_b(inum+ij) = dble(k)
!       elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
!         inum = i0
!         do while (nobk_e(inum+ij) /= C_0_R)
!           inum = inum + nxy            
!         end do 
!         nobk_e(inum+ij) = dble(k)        
!       end if
!       if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
!         inum = i0
!         do while (nobk_e(inum+ij) /= C_0_R)
!           inum = inum + nxy           
!         end do 
!         nobk_e(inum+ij) = dble(g(3)%size)    
!       end if
!     end do
!     ip = ip + nxy
!   end do

!   ! ------------------------------------------------------------------ !
!   ! DEBUG

! #ifdef USE_MPI
!   if ( ims_npro_k > 1 ) then
!     call DNS_MPI_TRPB_K(nobk_b, nobk_b_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
!     call DNS_MPI_TRPB_K(nobk_e, nobk_e_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
!   endif
! #else 
!   nobk_b_out = nobk_b
!   nobk_e_out = nobk_e
! #endif

 

  return
end subroutine IBM_SPLINZ

!########################################################################

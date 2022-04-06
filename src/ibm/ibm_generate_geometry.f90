#include "types.h"
#ifdef USE_MPI 
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   generates relevant geometry fields for IBM routines
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

subroutine IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)
  
  use TLAB_VARS,     only : g, imax, jmax, kmax, isize_field
  use DNS_IBM
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS, only : ims_size_i, ims_size_j, ims_size_k    
  use TLAB_MPI_PROCS
#endif    
#ifdef IBM_DEBUG
#ifdef USE_MPI
  use TLAB_MPI_VARS, only : ims_pro, ims_err
#endif  
#endif
  
  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: epsi, epsj, epsk
  
#ifdef USE_MPI 
  TINTEGER, parameter                          :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter                          :: idj = TLAB_MPI_J_PARTIAL 
  TINTEGER, parameter                          :: idk = TLAB_MPI_K_PARTIAL 
  TINTEGER                                     :: dummy
#endif
  TINTEGER                                     :: i, j, k, ij, ik, jk, ip, inum
  TINTEGER                                     :: nyz, nxz, nxy
  character(len=32)                            :: fname
#ifdef IBM_DEBUG
#ifdef USE_MPI 
#else
  TINTEGER, parameter                          :: ims_pro = 0         
#endif
#endif

  ! ================================================================== !
  ! npages
#ifdef USE_MPI
  nyz = ims_size_i(idi) ! local
  nxz = ims_size_j(idj) 
  nxy = ims_size_k(idk) 
#else
  nyz = jmax * kmax     ! global
  nxz = imax * kmax     
  nxy = imax * jmax     
#endif

  ! ================================================================== !
  ! initialize 
  nobi(:)   = i0; nobj(:)   = i0; nobk(:)   = i0
  nobi_b(:) = i0; nobj_b(:) = i0; nobk_b(:) = i0
  nobi_e(:) = i0; nobj_e(:) = i0; nobk_e(:) = i0
  nobi_max  = i0; nobj_max  = i0; nobk_max  = i0

  ! ================================================================== !
  ! number of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((ip == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for objects
        nobi(jk) = i1
      end if 
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! check for interface 
        nobi(jk) = nobi(jk) + i1
      end if
    end do
    ip = ip + nyz
  end do

#ifdef IBM_DEBUG
  nobi_max = maxval(nobi)
#ifdef USE_MPI
  dummy = nobi_max
  call MPI_ALLREDUCE(dummy, nobi_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif
#endif

  ! ================================================================== !
  ! begin and end of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((i == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for interface
        nobi_b(jk) = i ! nobi_b
      end if
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then     ! nobi_b check for interface 
        inum = i0
        do while (nobi_b(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_b(inum+jk) = i + i1
      elseif((epsi(ip+jk-1) == C_1_R) .and. (epsi(ip+jk-1+nyz) == C_0_R)) then ! nobi_e check for interface 
        inum = i0
        do while (nobi_e(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = i        
      end if
      if((i == (g(1)%size - 1)) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobi_e(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = g(1)%size    
      end if
    end do
    ip = ip + nyz
  end do

  ! ================================================================== !
  ! number of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((ip == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for objects
        nobj(ik) = i1
      end if 
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! check for interface 
        nobj(ik) = nobj(ik) + i1
      end if
    end do
    ip = ip + nxz
  end do

#ifdef IBM_DEBUG
  nobj_max = maxval(nobj)
#ifdef USE_MPI
  dummy = nobj_max
  call MPI_ALLREDUCE(dummy, nobj_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif
#endif

  ! ================================================================== !
  ! begin and end of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((j == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for interface
        nobj_b(ik) = j ! nobj_b
      end if
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then     ! nobj_b check for interface 
        inum = i0
        do while (nobj_b(inum+ik) /= i0)
          inum = inum + nxz            
        end do 
        nobj_b(inum+ik) = j + i1
      elseif((epsj(ip+ik-1) == C_1_R) .and. (epsj(ip+ik-1+nxz) == C_0_R)) then ! nobj_e check for interface 
        inum = i0
        do while (nobj_e(inum+ik) /= i0)
          inum = inum + nxz            
        end do 
        nobj_e(inum+ik) = j        
      end if
      if((j == (g(2)%size - 1)) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobj_e(inum+ik) /= i0)
          inum = inum + nxz           
        end do 
        nobj_e(inum+ik) = g(2)%size    
      end if
    end do
    ip = ip + nxz
  end do

  ! ================================================================== !
  ! number of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
        nobk(ij) = i1
      end if 
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
        nobk(ij) = nobk(ij) + i1
      end if
    end do
    ip = ip + nxy
  end do

#ifdef IBM_DEBUG
  nobk_max = maxval(nobk)
#ifdef USE_MPI
  dummy = nobk_max
  call MPI_ALLREDUCE(dummy, nobk_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif
#endif

  ! ================================================================== !
  ! begin and end of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
        nobk_b(ij) = k ! nobj_b
      end if
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
        inum = i0
        do while (nobk_b(inum+ij) /= i0)
          inum = inum + nxy            
        end do 
        nobk_b(inum+ij) = k + i1
      elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
        inum = i0
        do while (nobk_e(inum+ij) /= i0)
          inum = inum + nxy            
        end do 
        nobk_e(inum+ij) = k        
      end if
      if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobk_e(inum+ij) /= i0)
          inum = inum + nxy           
        end do 
        nobk_e(inum+ij) = g(3)%size    
      end if
    end do
    ip = ip + nxy
  end do

  ! ================================================================== !
#ifdef IBM_DEBUG
  if (ims_pro == 0) then
    write(*,*) '======== Max number of objects in each direction ========'
    write(*,*) 'max number of objects in x = ', nobi_max
    write(*,*) 'max number of objects in y = ', nobj_max
    write(*,*) 'max number of objects in z = ', nobk_max
  end if
#endif

  return
end subroutine IBM_GENERATE_GEOMETRY

!########################################################################
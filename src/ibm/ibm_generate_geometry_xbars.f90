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
!#   square bars in x-/streamwise direction,
!#   equally distributed + spaced in z-/spanwise direction
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#  epsilon field 'eps' is an indicator field:
!#    eps(i,j,k) = 0   for  flow
!#    eps(i,j,k) = 1   for  solid
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#   interfaces of bars have to be on gridpoints (center doesn't matter), if:
!#     kmax_total/(2*nbars) = integer :: wbar = even   (center is on gridpoints)                       
!#     kmax_total/(2*nbars) = x * 1/2 :: wbar = uneven (center is between grid points)                    
!#
!#
!########################################################################

subroutine IBM_GENERATE_GEOMETRY_XBARS(wrk3d)

  use DNS_IBM
  use DNS_GLOBAL, only: g
  use DNS_GLOBAL, only: imax, jmax, kmax 
  use DNS_GLOBAL, only: isize_field

#ifdef USE_MPI 
  use DNS_MPI,    only: ims_offset_i, ims_offset_j, ims_offset_k
  use DNS_MPI,    only: ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k  ! each number of each proc
  use DNS_MPI,    only: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k ! total numbers of proc
#endif 

  implicit none

#include "integers.h"

#ifdef USE_MPI
#include "mpif.h"
#else
  TINTEGER, parameter :: ims_offset_i=0, ims_offset_j=0, ims_offset_k=0
  TINTEGER, parameter :: ims_pro_i=0,    ims_pro_j=0,    ims_pro_k=0,    ims_pro=0  
  TINTEGER, parameter :: ims_npro_i=1,   ims_npro_j=1,   ims_npro_k=1,   ims_npro=0 
#endif

  TINTEGER                                     :: nbars, hbar, wbar
  TREAL                                        :: zcenter_bar
  TINTEGER, dimension(xbars_geo%number)        :: zstart_bar, zend_bar
  TINTEGER                                     :: istart, iend, jstart, jend, kstart, kend
  TINTEGER                                     :: i,j,k,l

  ! DEBUG 
  character(32)                                :: fname
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

  ! ================================================================== !

  ! global array indicies for each mpi task (indices start with 0)
  istart = ims_offset_i; iend = ims_offset_i + imax - 1    
  jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
  kstart = ims_offset_k; kend = ims_offset_k + kmax - 1

  ! DEBUG
  if ( ims_pro == 0 ) then 
    write(*,*) '======== Initialization of grid and decomposition ======='
    write(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax, ' x ', kmax*ims_npro_k
    ! write(*,*) 'GRID g(...): ', g(1)%size,      ' x ', g(2)%size,' x ', g(3)%size 
    write(*,*) 'DECOMP: ranks', ims_npro_i,     ' x ', 1,    ' x ', ims_npro_k 
    write(*,*) '        grid ', imax,           ' x ', jmax, ' x ', kmax 
  end if
  do l = 0, ims_npro - 1
    if (ims_pro == l) then
      write(*,'(1X,A,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A,I3)') & 
      'Task:', ims_pro, &
      ' |i-part:', istart, '-', iend,&
      ' |j-part:', jstart, '-', jend,&
      ' |k-part:', kstart, '-', kend 
    end if 
  end do 
    
  ! geometry (from dns.ini)
  nbars=xbars_geo%number; hbar=xbars_geo%height; wbar=xbars_geo%width  
  
  ! global z-positions of bars, equally distributed on gridpoints with equal spacing
  do l = 1, nbars
    zcenter_bar   = g(3)%size / nbars * (l - 0.5)
    zstart_bar(l) = int(zcenter_bar - 0.5 * wbar)
    zend_bar(l)   = int(zcenter_bar + 0.5 * wbar)
  end do

  ! DEBUG
  if ( ims_pro == 0 ) then
    write(*,*) '======== Z - Positions of streamwise aligned bars ======='
    do l = 1, nbars
      write(*,*)'bar nr.', l, ' start:', zstart_bar(l) + i1, ' end:', zend_bar(l)
    end do
  end if
  
  ! ini eps_aux
  eps_aux(:,:,:) = C_0_R
  
  ! streamwise aligned square bars, equaly spaced, interfaces on gridpoints (define own geometry)
  !! ============================= Comment =========================== !!
  !! (k+kstart)>zstart_bar(l) ==> objects start at  zstart_bar(l) + 1  !!
  !! ================================================================= !!
  do j = 1, hbar   
    do k = 1, kmax 
      do l = 1, nbars 
        if( ((k+kstart)>zstart_bar(l)) .and. ((k+kstart)<=zend_bar(l)) ) then 
          do i = 1, imax
            eps_aux(i,j,k) = C_1_R
          end do
        end if
      end do 
    end do
  end do 

  ! reshape 3D-eps_aux field into 1D-eps
  eps = reshape(eps_aux,(/isize_field/))

  ! write eps field
  write(fname,*) i0; 
  fname = trim(adjustl('eps'))//trim(adjustl(fname))
  if (ims_pro == 0) then
    write(*,*) '======== Write eps field ================================'
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, i1, imax*jmax*kmax, eps, wrk3d)

  return
end subroutine IBM_GENERATE_GEOMETRY_XBARS

!########################################################################
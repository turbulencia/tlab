#include "dns_error.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   square bars in x-/streamwise direction,
!#   equally distributed + spaced in z-/spanwise direction
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#  epsilon field 'eps' is an indicator field:
!#    eps(i,j,k) = 0 for fluid domain
!#    eps(i,j,k) = 1 for solid domain
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

  use IBM_VARS
  use TLAB_VARS,      only : g, imax, jmax, kmax, isize_field
  use IO_FIELDS
  use TLAB_CONSTANTS, only : wi, wp
#ifdef USE_MPI 
  use MPI
  use TLAB_MPI_VARS,  only : ims_offset_i, ims_offset_j, ims_offset_k
#ifdef IBM_DEBUG
  use TLAB_MPI_VARS,  only : ims_pro, ims_npro, ims_npro_i, ims_npro_k 
#endif
#endif 

  implicit none

  real(wp), dimension(imax, jmax, kmax), intent(inout) ::  wrk3d

#ifdef USE_MPI
#else
  integer(wi), parameter                               :: ims_offset_i=0, ims_offset_j=0, ims_offset_k=0
#ifdef IBM_DEBUG
  integer(wi), parameter                               :: ims_pro=0, ims_npro_i=1, ims_npro_k=1, ims_npro=0 
#endif
#endif
  integer(wi)                                          :: nbars, hbar, wbar
  real(wp)                                             :: zcenter_bar
  integer(wi), dimension(xbars_geo%number)             :: zstart_bar, zend_bar
  integer(wi)                                          :: istart, iend, jstart, jend, kstart, kend
  integer(wi)                                          :: i,j,k,l

  ! ================================================================== !
  ! global array indicies for each mpi task (indices start with 0)
  istart = ims_offset_i; iend = ims_offset_i + imax - 1    
  jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
  kstart = ims_offset_k; kend = ims_offset_k + kmax - 1

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) then 
    write(*,*) '======== Initialization of grid and decomposition ======='
    write(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax, ' x ', kmax*ims_npro_k
    write(*,*) 'DECOMP: ranks', ims_npro_i,     ' x ', 1,    ' x ', ims_npro_k 
    write(*,*) '        grid ', imax,           ' x ', jmax, ' x ', kmax 
  end if
  ! do l = 0, ims_npro - 1
  !   if (ims_pro == l) then
  !     write(*,'(1X,A,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A,I3)') & 
  !     'Task:', ims_pro, &
  !     ' |i-part:', istart, '-', iend,&
  !     ' |j-part:', jstart, '-', jend,&
  !     ' |k-part:', kstart, '-', kend 
  !   end if 
  ! end do 
#endif

  ! geometry
  nbars=xbars_geo%number; hbar=xbars_geo%height; wbar=xbars_geo%width  
  wrk3d(:,:,:) = 0.0_wp

  ! global z-positions of bars, equally distributed on gridpoints with equal spacing
  do l = 1, nbars
    zcenter_bar   = g(3)%size / nbars * (l - 0.5)
    zstart_bar(l) = int(zcenter_bar - 0.5 * wbar)
    zend_bar(l)   = int(zcenter_bar + 0.5 * wbar)
  end do

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) then
    write(*,*) '======== Z - Positions of streamwise aligned bars ======='
    if (xbars_geo%mirrored) write(*,*) 'Bars are mirrored on upper wall!'
    do l = 1, nbars
      write(*,*)'bar nr.', l, ' start:', zstart_bar(l) + 1, ' end:', zend_bar(l)
    end do
  end if
#endif
  
  ! streamwise aligned square bars, equaly spaced, interfaces on gridpoints
  ! (k+kstart)>zstart_bar(l) ==> objects start at  zstart_bar(l) + 1 
  
  ! geometries on lower boundary 
  do j = 1, hbar   
    do k = 1, kmax 
      do l = 1, nbars 
        if( ((k+kstart)>zstart_bar(l)) .and. ((k+kstart)<=zend_bar(l)) ) then 
          do i = 1, imax
            wrk3d(i,j,k) = 1.0_wp
          end do
        end if
      end do 
    end do
  end do
  
  ! geometries on upper boundary
  if (xbars_geo%mirrored) then
    do j = jmax-hbar+1, jmax   
      do k = 1, kmax 
        do l = 1, nbars 
          if( ((k+kstart)>zstart_bar(l)) .and. ((k+kstart)<=zend_bar(l)) ) then 
            do i = 1, imax
              wrk3d(i,j,k) = 1.0_wp
            end do
          end if
        end do 
      end do
    end do
  end if

  ! reshape 3D-field into 1D-field
  eps          = reshape(wrk3d,(/isize_field/))
  wrk3d(:,:,:) = 0.0_wp

  ! io of eps
  select case( ibm_io )
  case ( IBM_IO_REAL )
    call IO_WRITE_FIELDS(eps_name_real, IO_FLOW, imax,jmax,kmax, 1, eps, wrk3d)
  case ( IBM_IO_INT  )
    call IBM_IO_WRITE_INT_GEOMETRY(wrk3d)
  case ( IBM_IO_BIT  )
    call IBM_IO_WRITE_BIT_GEOMETRY(wrk3d)
  end select 
  
  return
end subroutine IBM_GENERATE_GEOMETRY_XBARS

!########################################################################
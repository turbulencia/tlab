#include "dns_error.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2023/09/08 - S. Deshpande
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   3D SINE wave,
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
!#
!########################################################################

subroutine IBM_GENERATE_GEOMETRY_VALLEY(wrk3d)

  use IBM_VARS
  use FDM,      only : g
  use TLab_Memory, only: imax, jmax, kmax, isize_field
  use TLab_Constants, only : wi, wp, pi_wp
  use TLab_WorkFlow,      only : stagger_on
#ifdef USE_MPI 
  use TLabMPI_VARS,  only : ims_offset_i, ims_offset_j, ims_offset_k
#ifdef IBM_DEBUG
  use TLabMPI_VARS,  only : ims_pro, ims_npro, ims_npro_i, ims_npro_k 
#endif
#endif 

implicit none

real(wp), dimension(imax, jmax, kmax), intent(inout) ::  wrk3d

#ifdef USE_MPI
#else
  integer(wi), parameter                             :: ims_offset_i=0, ims_offset_j=0, ims_offset_k=0
#ifdef IBM_DEBUG
  integer(wi), parameter                             :: ims_pro=0, ims_npro_i=1, ims_npro_k=1, ims_npro=0 
#endif
#endif
  integer(wi)                                        :: ncrest, hhill, whill, hill_slope
  integer(wi)                                        :: istart, iend, jstart, jend, kstart, kend
  integer(wi)                                        :: i, j, k
  real(wp)                                           :: dx
  ! ================================================================== !
  ! Decleration of variables
  dx = (2*pi_wp/(g(1)%size))                ! Unit cell width of hill base

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
  
#endif

  ! geometry
  ncrest=ibm_geo%number; hhill=ibm_geo%height; whill=ibm_geo%width; hill_slope=ibm_geo%hill_slope 
  wrk3d(:,:,:) = 0.0_wp
  do i = 1,imax
    do k = 1,kmax
      do j = 1,hhill
        if ((j+jstart) < ((hhill/(2**hill_slope))*(1 + cos(dx*(i + istart)))**hill_slope)) then
          wrk3d(i,j,k) = 1.0_wp
        else
          wrk3d(i,j,k) = 0.0_wp
        end if 
      end do
    end do
  end do
  
  ! reshape 3D-field into 1D-field
  eps = reshape(wrk3d,(/isize_field/))
  call IBM_IO_WRITE(wrk3d, .false.)
  wrk3d(:,:,:) = 0.0_wp

  ! ================================================================== !
  if (stagger_on) then
    ! create epsp field
    do i = 1,imax
      do k = 1,kmax
        do j = 1,hhill
          if ((j+jstart) < ((hhill/(2**hill_slope))*(1 + cos(dx*(i + istart)))**hill_slope)) then
            wrk3d(i,j,k) = 1.0_wp
          else
            wrk3d(i,j,k) = 0.0_wp
          end if 
        end do
      end do
    end do

    ! reshape 3D-field into 1D-field
    epsp = reshape(wrk3d,(/isize_field/))
    
    ! export geometry to confirm the shape and size
    call IBM_IO_WRITE(wrk3d, .true.)
    wrk3d(:,:,:) = 0.0_wp
  end if

  return
end subroutine IBM_GENERATE_GEOMETRY_VALLEY
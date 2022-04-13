#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
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
!#   verify that the geometry used satisfies all restrictions
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

subroutine IBM_VERIFY_GEOMETRY()

  use DNS_IBM 
  use TLAB_VARS,      only : g, imode_ibm_scal
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_size_i, ims_size_j, ims_size_k, ims_err
#ifdef IBM_DEBUG
  use TLAB_MPI_VARS,  only : ims_pro
#endif
#else
  use TLAB_VARS,      only : imax, jmax, kmax
#endif    
   
  implicit none

#include "integers.h"

#ifdef USE_MPI 
  TINTEGER, parameter :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter :: idj = TLAB_MPI_J_PARTIAL 
  TINTEGER, parameter :: idk = TLAB_MPI_K_PARTIAL 
  TREAL               :: dummy
#else
#ifdef IBM_DEBUG
  TINTEGER, parameter :: ims_pro = 0 
#endif
#endif
  TINTEGER            :: nyz, nxz, nxy    
  TREAL               :: ob_min
  ! ================================================================== !

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) write(*,*) '================ Verifying the geometry ================='
#endif

  ! nlines
#ifdef USE_MPI 
  nyz = ims_size_i(idi)
  nxz = ims_size_j(idj) 
  nxy = ims_size_k(idk) 
#else
  nyz = jmax * kmax  
  nxz = imax * kmax     
  nxy = imax * jmax
#endif

  ! check if any objects are present
  ob_min = sum(eps)
#ifdef USE_MPI
  dummy = ob_min
  call MPI_ALLREDUCE(dummy, ob_min, i0, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  if ( ob_min == 0 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY no objects in flow.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  endif

  ! check all objects in each direction
  call IBM_VERIFY(g(1), nyz, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e) ! x
  call IBM_VERIFY(g(2), nxz, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e) ! y
  call IBM_VERIFY(g(3), nxy, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e) ! z

  ! check if objects on upper boundary are present
  if ( imode_ibm_scal == 1 ) then
    call IBM_VERIFY_SCAL(eps)
  end if

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) write(*,*) 'Verification successfull!'
#endif

  return
end subroutine IBM_VERIFY_GEOMETRY

!########################################################################

subroutine IBM_VERIFY(g, nlines, isize_nob, isize_nob_be, nob, nob_b, nob_e)

  use DNS_IBM,        only : nflu
  use TLAB_CONSTANTS, only : efile
  use TLAB_TYPES,     only : grid_dt
  use TLAB_PROCS

  implicit none
  
#include "integers.h"

  type(grid_dt),                     intent(in) :: g
  TINTEGER,                          intent(in) :: nlines, isize_nob, isize_nob_be
  TINTEGER, dimension(isize_nob),    intent(in) :: nob
  TINTEGER, dimension(isize_nob_be), intent(in) :: nob_b, nob_e

  TINTEGER                                      :: ii, ip, ipp, iob
  TINTEGER                                      :: fp_min, fp_l, fp_r, fp_inter
  TINTEGER                                      :: sp_min, sp_ob
  
  ! ================================================================== !
  
  ! min vals solid/fluid
  fp_min = nflu ! 3rd point can be next interface/boundary
  sp_min = i4   ! 3 solid + 2 interface points - 1

  ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
  do ii = 1, nlines        ! index of ii-plane, loop over plane and check for objects in each line
    if ( nob(ii) /= 0 ) then ! if line contains immersed object(s) --yes-->  spline interpolation
      ip = i0        
      do iob = 1, nob(ii)  ! loop over immersed object(s)
        ! ================================================================== !
        ! check number of fluid points (to boarders and between objects)
        if ( iob == 1 ) then                    ! left
          fp_l = nob_b(ip+ii) - i1
          if ( fp_l < fp_min .and. fp_l /= 0 ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY not enough fluid points between left boarder and first object.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        end if 
        if ( iob > 1 .and. iob < nob(ii) ) then ! in between objects
          ipp = ip + nlines ! next obj
          fp_inter = nob_b(ipp+ii) - nob_e(ip+ii)
          if ( fp_inter < fp_min ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY not enough fluid points between objects.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        end if
        if ( iob == nob(ii) ) then              ! right
          fp_r = g%size - nob_e(ip+ii)
          if ( fp_r < fp_min .and. fp_r /= 0 ) then
            call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY not enough fluid points between right boarder and first object.')
            call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
          end if
        end if 
        ! ================================================================== !
        ! check number of solid points of object
        sp_ob = nob_e(ip+ii) - nob_b(ip+ii)
        if ( sp_ob < sp_min ) then
          call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY not enough solid points in object(s).')
          call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
        end if
        ! ================================================================== !        
        ip = ip + nlines
      end do
    end if
  end do

  return
end subroutine IBM_VERIFY

!########################################################################

subroutine IBM_VERIFY_SCAL(eps)

  use TLAB_VARS,      only : isize_field, imax,jmax,kmax
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_err
#endif    

  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field), intent(in) :: eps
  
  TINTEGER                                     :: ip_t, k, nxy
  TREAL                                        :: dummy, top

  ! ================================================================== !

  ! check that no objects are present on upper 
  nxy  = imax*jmax
  ip_t = imax*(jmax-1) + 1
  do k = 1, kmax
    top = sum(eps(ip_t:ip_t+imax-1)) 
    ip_t = ip_t + nxy
  end do

#ifdef USE_MPI
  dummy = top
  call MPI_ALLREDUCE(dummy, top, i0, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

  if ( top == 0 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY no objects on upper domain allowed if IBM is turned on for scalars.')
    call TLAB_STOP(DNS_ERROR_IBM_GEOMETRY)
  endif

end subroutine IBM_VERIFY_SCAL
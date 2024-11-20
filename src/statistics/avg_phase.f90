#include "types.h"
#include "dns_const_mpi.h"
#include "dns_error.h"
#include "dns_const.h"

module AVG_PHASE
  
  use TLab_WorkFlow
  use TLab_Constants, only : wp, wi, longi, efile
  use TLAB_CONSTANTS, only : sizeofint, sizeofreal
  use FDM, only : g
  use TLAB_VARS, only : imax, jmax, kmax, isize_field
  use TLAB_VARS, only : rtime
  use TLAB_VARS, only : visc, froude, rossby, prandtl, mach
  use TLAB_VARS, only : imode_eqns
  use TLAB_VARS, only : inb_flow, inb_scal
  use TLAB_ARRAYS, only: q, s
  use TLab_Arrays, only: wrk2d, wrk3d
  use Thermodynamics, only: gama0
  use IO_FIELDS
  use TLab_Memory, only: Tlab_Allocate_Real_LONG
  use Avg_Vars
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

  interface AvgPhaseSpace 
    module procedure AvgPhaseSpaceFieldPtr, AvgPhaseSpaceIndex
  end interface AvgPhaseSpace

  ! type(phaseavg_dt) :: PhAvg
  ! real(wp), dimension(:), allocatable, target :: avg_flow, avg_stress, avg_p, avg_scal
  ! integer(wi) :: nxy, nxz, nyz, nz_total
  ! integer(wi) :: avg_planes
  ! character(len=32), parameter :: avgu_name = 'avg_flow'
  ! character(len=32), parameter :: avgstr_name = 'avg_stress'
  ! character(len=32), parameter :: avgp_name = 'avg_p'
  ! character(len=32), parameter :: avgs_name = 'avg_scal'

  ! public :: AvgPhaseSpace
  ! public :: avg_flow, avg_p, avg_scal, avg_stress, avg_planes
  ! public :: PhAvg
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AvgPhaseInitializeMemory(C_FILE_LOC, restart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only : ims_size_i, ims_size_k, ims_pro, ims_npro_i, ims_pro_k
#endif

    implicit none
    character(len=*), intent(in) :: C_FILE_LOC
    integer(wi), intent(in)    :: restart
    integer(longi) :: alloc_size
      ! ================================================================== !

    nxy = imax*jmax
    nyz = jmax*kmax
    nxz = imax*kmax

    if ( restart == -1 ) then ! used for calls from outside dns_main / dns.x 
      avg_planes = 0 
    elseif (mod(restart,PhAvg%stride) == 0) then
      avg_planes = restart/PhAvg%stride
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
      call TLAB_STOP(DNS_ERROR_AVG_PHASE)
    end if
    
#ifdef USE_MPI
    if (ims_pro_k == 0) then
#endif
      alloc_size = imax*jmax*(avg_planes+1)
      call Tlab_Allocate_Real_LONG(C_FILE_LOC,   avg_flow,      [alloc_size*inb_flow], 'avgflow.')
      call Tlab_Allocate_Real_LONG(C_FILE_LOC,   avg_stress,    [alloc_size*       6], 'avgstr.' ) ! allocated not yet coded
      call Tlab_Allocate_Real_LONG(C_FILE_LOC,   avg_p,         [alloc_size*       1], 'avgp.'   )
      call Tlab_Allocate_Real_LONG(C_FILE_LOC,   avg_scal,      [alloc_size*inb_scal], 'avgscal.')
      
      avg_flow(:)   = 0.0_wp
      avg_stress(:) = 0.0_wp
      avg_p(:)      = 0.0_wp
      avg_scal(:)   = 0.0_wp

#ifdef USE_MPI
    end if
#endif
    return
  end subroutine AvgPhaseInitializeMemory

  subroutine AvgPhaseSpaceFieldPtr(localsum, nfield, itr, it_first, it_save, field)
    implicit none
    real(wp), dimension(imax, jmax), intent(inout) :: localsum
    integer(wi), intent(in) :: nfield
    integer(wi), intent(in) :: itr, it_first, it_save
    real(wp), pointer, intent(in) :: field(:)
  
    integer :: index_loc = 4 ! Index needs to be set to appropriate value for pressure. Needed later for if statement

    call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index_loc, field)
  end subroutine AvgPhaseSpaceFieldPtr


  subroutine AvgPhaseSpaceIndex(localsum, nfield, itr, it_first, it_save, index)
    implicit none
    real(wp), dimension(imax,jmax), intent(inout) :: localsum
    integer(wi), intent(in) :: nfield
    integer(wi), intent(in) :: itr, it_first, it_save, index
    real(wp), pointer :: field_loc(:) => null()

    call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field_loc)    
  end subroutine AvgPhaseSpaceIndex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI 
    use MPI
    use TLabMPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

    implicit none
    real(wp), dimension(imax*jmax), intent(inout) :: localsum
    integer(wi) , intent(in) :: nfield
    integer(wi) , intent(in) :: itr, it_first, it_save, index
    real(wp), dimension(imax, jmax, kmax, 1), target, intent(in) :: field

    integer(wi) :: k, ifld, plane_id
    real(wp), dimension(:), pointer :: avg_ptr
    real(wp), dimension(:), pointer :: loc_field
    integer(wi) :: ipl_srt, ipl_end, iavg_srt, iavg_end, lpl_srt, lpl_end
      ! ================================================================== !
      ! Calculation of the plane id to write the spatial average 
    plane_id = 1
    if (it_save /= 0) plane_id = mod((itr-1) - (it_first), it_save) + 1

      ! Determing the tendency to be written
    if (index == 1) then
      avg_ptr => avg_flow
      call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 2) then
      avg_ptr => avg_scal
      call c_f_pointer(c_loc(s), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 4) then
      avg_ptr => avg_p
      call c_f_pointer(c_loc(field), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 5) then
      avg_ptr => avg_stress
      call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
      ! Not yet coded
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in AvgPhaseSpaceExec')
      call TLAB_STOP(DNS_ERROR_AVG_PHASE)
    end if

    if ((index == 1) .or. (index == 2) .or. (index == 4)) then
      do ifld = 1, nfield
        localsum = 0.0_wp
        ! Computing the space average
        do k = 1, kmax
          ! Computing the start and end of the plane in field
          ipl_srt = (ifld - 1)*isize_field + nxy*(k-1) + 1
          ipl_end = (ifld - 1)*isize_field + nxy*k

          localsum = localsum + loc_field(ipl_srt:ipl_end) / g(3)%size !loc_field(:,:,k,ifld)/g(3)%size
        end do
        
        ! Computing the local sum from start and end of the field for accumulating the space averages
        iavg_srt = (ifld - 1)*nxy*(avg_planes+1) + nxy*(plane_id-1) + 1
        iavg_end = (ifld - 1)*nxy*(avg_planes+1) + nxy*plane_id  
#ifdef USE_MPI
        if (ims_pro_k == 0) then
          call MPI_Reduce(localsum, avg_ptr(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_ptr(imax*jmax*restarts*fld)
        else
          call MPI_Reduce(localsum, MPI_IN_PLACE,                nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        end if
#else
        avg_ptr(iavg_srt:iavg_end) = localsum
#endif
      
        lpl_srt = (ifld - 1)*nxy*(avg_planes+1) + nxy*avg_planes + 1
        lpl_end = ifld*nxy*(avg_planes+1)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
          avg_ptr(lpl_srt:lpl_end) = avg_ptr(lpl_srt:lpl_end) + avg_ptr(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI 
        endif 
#endif 
      end do
    end if
  return
  end subroutine AvgPhaseSpaceExec

  subroutine AvgPhaseStress(q, itr, it_first, it_save)
    real(wp), dimension(:,:), intent(in) :: q
    integer(wi), intent(in) :: itr
    integer(wi) , intent(in) :: it_first
    integer(wi) , intent(in) :: it_save

    real(wp), dimension(:), pointer :: u,v,w
    integer(wi) :: plane_id

    target q

    u => q(:, 1)
    v => q(:, 2)
    w => q(:, 3)

    ! uu 1
    ! uv 2
    ! uw 5
    ! vv 3
    ! vw 4
    ! ww 6

    plane_id = 1
    if (it_save /= 0) plane_id = mod((itr-1) - (it_first), it_save) + 1

    call AvgPhaseCalcStress(u, u, 1, plane_id)
    call AvgPhaseCalcStress(u, v, 2, plane_id)
    call AvgPhaseCalcStress(v, v, 4, plane_id)
    call AvgPhaseCalcStress(v, w, 5, plane_id)
    call AvgPhaseCalcStress(u, w, 3, plane_id)
    call AvgPhaseCalcStress(w, w, 6, plane_id)

  end subroutine AvgPhaseStress

  subroutine AvgPhaseCalcStress(field1, field2, stress_id, plane_id)
#ifdef USE_MPI 
    use MPI
    use TLabMPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif
    real(wp), pointer, intent(in) :: field1(:)
    real(wp), pointer, intent(in) :: field2(:)
    integer(wi), intent(in) :: stress_id
    integer(wi), intent(in) :: plane_id

    integer(wi) :: k, ipl_srt, ipl_end, iavg_srt, iavg_end, lpl_srt, lpl_end

    wrk3d(:) = 0.0_wp
    wrk2d(:,:) = 0.0_wp 

    wrk3d(:) = field1(:)*field2(:)

    do k = 1, kmax
      ipl_srt = nxy*(k-1) + 1
      ipl_end = nxy*k

      wrk2d(:,1) = wrk2d(:,1) + (wrk3d(ipl_srt:ipl_end))/g(3)%size
    end do

    iavg_srt = (stress_id -1)*nxy*(avg_planes+1) + nxy*(plane_id - 1) + 1
    iavg_end = (stress_id -1)*nxy*(avg_planes+1) + nxy*(plane_id)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
          call MPI_Reduce(wrk2d, avg_stress(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_ptr(imax*jmax*restarts*fld)
        else
          call MPI_Reduce(wrk2d, MPI_IN_PLACE,                nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        end if
#else
        avg_stress(iavg_srt:iavg_end) = wrk2d(:,1)
#endif

    lpl_srt = (stress_id -1)*nxy*(avg_planes+1) + nxy*avg_planes + 1
    lpl_end = (stress_id)*nxy*(avg_planes+1)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
          avg_stress(lpl_srt:lpl_end) = avg_stress(lpl_srt:lpl_end) + avg_stress(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI 
        endif 
#endif 
    
  end subroutine AvgPhaseCalcStress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AvgPhaseResetVariable()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI 
    use MPI
    use TLabMPI_VARS,  only : ims_comm_x, ims_err, ims_pro, ims_pro_k
#endif

#ifdef USE_MPI
    if (ims_pro_k == 0) then
#endif
      avg_flow(:)   = 0.0_wp
      avg_stress(:) = 0.0_wp
      avg_p(:)      = 0.0_wp
      avg_scal(:)   = 0.0_wp
#ifdef USE_MPI
    end if
#endif
  end subroutine AvgPhaseResetVariable
end module

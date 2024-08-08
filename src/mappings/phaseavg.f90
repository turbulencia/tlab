#include "types.h"
#include "dns_const_mpi.h"
#include "dns_error.h"
#include "dns_const.h"

module PHASEAVG
  
  USE TLAB_PROCS
  USE TLAB_TYPES
  use TLAB_CONSTANTS, only : wp, wi
  use TLAB_VARS,      only : imax, jmax, kmax, g
  use TLAB_VARS,      only : rtime
  use TLAB_VARS,      only : visc, froude, rossby, prandtl, mach, gama0
  use TLAB_VARS,      only : imode_eqns
  use TLAB_VARS,      only : inb_flow, inb_scal
  use TLAB_VARS,      only : phAvg
  use TLAB_CONSTANTS, only : sizeofint, sizeofreal

  implicit none

  real(wp),   dimension(:,:,:,:), allocatable   :: avg_flow, avg_stress, avg_p, avg_scal
  integer(wi)                                   :: nx_total, ny_total, nz_total
  integer(wi)                                   :: nxy, nxz, nyz
  integer(wi)                                   :: num_strides
  character(len=32), parameter                  :: avgu_name       = 'avg_flow'
  character(len=32), parameter                  :: avgstr_name     = 'avg_stress'
  character(len=32), parameter                  :: avgp_name       = 'avg_p'
  character(len=32), parameter                  :: avgs_name       = 'avg_scal'
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Allocate(C_FILE_LOC, restart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS,  only : ims_size_i, ims_size_k, ims_pro, ims_npro_i, ims_pro_k
#endif

    implicit none
    character(len=*), intent(in)    :: C_FILE_LOC
    integer(wi)     , intent(in)    :: restart

    integer(wi)                     :: isize_planej, isize_planek
      ! ================================================================== !

#ifdef USE_MPI    
    nx_total = imax*ims_npro_i ! allocate to the 0th process
#else
    nx_total = imax
#endif
    nxy = imax*jmax
    nyz = jmax*kmax
    nxz = imax*kmax

    if ( restart == -1 ) then ! used for calls from outside dns_main / dns.x 
      num_strides = 0 
    elseif (mod(restart,phAvg%stride) == 0) then
      num_strides = restart/phAvg%stride
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
      call TLAB_STOP(DNS_ERROR_PHASEAVG)
    end if
    
#ifdef USE_MPI
    if (ims_pro_k == 0) then
#endif
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_flow,      [imax,jmax,inb_flow,num_strides+1], 'avgflow.')
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_stress,    [imax,jmax,6       ,num_strides+1], 'avgstr.' ) ! allocated not yet coded
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_p,         [imax,jmax,1       ,num_strides+1], 'avgp.'   )
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_scal,      [imax,jmax,inb_scal,num_strides+1], 'avgscal.')

      avg_flow   = 0.0_wp
      avg_stress = 0.0_wp
      avg_p      = 0.0_wp
      avg_scal   = 0.0_wp

#ifdef USE_MPI
    end if
#endif
    return
  end subroutine PhaseAvg_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  subroutine PhaseAvg_Initialize()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use IO_FIELDS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only : ims_comm_x, ims_pro_k, ims_pro
#endif
    use TLAB_VARS, only : rtime
    use TLAB_VARS, only : visc, froude, rossby, prandtl, mach, gama0
    use TLAB_VARS, only : imode_eqns

    implicit none

    integer(wi)                                     :: isize, iheader, header_offset, id
    integer(wi),       parameter                    :: isize_max = 20
    real(wp)                                        :: params(isize_max)
      ! ================================================================== !

    isize = 0
    isize = isize + 1; params(isize) = rtime
    isize = isize + 1; params(isize) = visc ! inverse of reynolds
    if (iheader == IO_SCAL) then
        isize = isize + 1 + 1                     ! prepare space for schmidt and damkohler

    else if (iheader == IO_FLOW) then
        isize = isize + 1; params(isize) = froude
        isize = isize + 1; params(isize) = rossby
        if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
            isize = isize + 1; params(isize) = gama0
            isize = isize + 1; params(isize) = prandtl
            isize = isize + 1; params(isize) = mach
        end if
    end if

    ! INITIALIZATION OF MPI TYPES SHOULD ONLY BE CARRIED OUT ONCE *AND* DURING INITIALIZATION
    header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

    id = IO_AVERAGE_PLANE
    io_aux(id)%offset = header_offset
    io_aux(id)%precision = IO_TYPE_DOUBLE
#ifdef USE_MPI
    if (ims_pro_k == 0) THEN 
        io_aux(id)%active = .true.
        io_aux(id)%communicator = ims_comm_x
        io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOY(imax, jmax, MPI_REAL8)
    else
        io_aux(id)%active = .false.  
        io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOY(imax, jmax, MPI_REAL8)

    end if
#endif
  end subroutine PhaseAvg_Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Space_Accumulate(localsum, field, nfield, itime, it_first, it_save, index)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_VARS, only : imax, jmax, kmax
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

    implicit none
    real(wp), dimension(imax, jmax)                           , intent(inout)        :: localsum
    real(wp), dimension(imax, jmax, kmax, *)                  , intent(in)           :: field
    integer(wi)                                               , intent(in)           :: nfield
    integer(wi)                                               , intent(in)           :: itime, it_first, it_save, index
    

    integer(wi)                                               :: k, fld, ifld, jfld, plane_id
      ! ================================================================== !

    plane_id = 1
    if (it_save /= 0) plane_id = mod((itime-1) - (it_first), it_save) + 1 
    if ((index == 1) .or. (index == 2) .or. (index == 4)) then
      do ifld = 1, nfield
          localsum = 0.0_wp
          do k = 1, kmax
            localsum = localsum + field(:,:,k,ifld)/g(3)%size
          end do
#ifdef USE_MPI
          if (ims_pro_k == 0) then
            if (index == 1) then
              call MPI_Reduce(localsum, avg_flow(:,:,ifld,plane_id), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
              call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_flow)
            else if (index == 2) then
              call MPI_Reduce(localsum, avg_scal(:,:,ifld,plane_id), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
              call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_scal)
            else if (index == 4) then
              call MPI_Reduce(localsum, avg_p(:,:,ifld,plane_id), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
              call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_p)
            else if (index == 5) then ! allocated not yet coded
              call MPI_Reduce(localsum, avg_stress(:,:,ifld,plane_id), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
              call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_stress)
            end if
          else
            call MPI_Reduce(localsum, MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
          end if
#else
          if (index == 1) then
            avg_flow(:,:,ifld,plane_id) = localsum
            call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_scal)
          else if (index == 2) then
            avg_scal(:,:,ifld,plane_id) = localsum
            call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_scal)
          else if (index == 4) then
            avg_p(:,:,ifld,plane_id) = localsum
            call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_p)
          else if (index == 5) then ! allocated not yet coded
            avg_stress(:,:,ifld,plane_id) = localsum
            call PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg_stress)
          end if
#endif
      end do
    end if
  return
  end subroutine PhaseAvg_Space_Accumulate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_PhasePlane(nfield, it_save, ifld, plane_id, localsum, avg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

    implicit none
    integer(wi)                                               , intent(in)           :: nfield, it_save
    integer(wi)                                               , intent(in)           :: ifld
    integer(wi)                                               , intent(in)           :: plane_id 
    real(wp), dimension(imax, jmax)                           , intent(in)           :: localsum
    real(wp), dimension(imax, jmax, num_strides+1, nfield)    , intent(inout)        :: avg

#ifdef USE_MPI
    if ( ims_pro_k == 0) then
#endif
      if (it_save /= 0) then
        avg(:,:,ifld,it_save + 1) = avg(:,:,ifld,it_save + 1) + avg(:,:,ifld,plane_id)/it_save
      end if
#ifdef USE_MPI
    end if
#endif
  return
  end subroutine PhaseAvg_PhasePlane

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Write(nfield, avg, iheader, it_save, basename, avg_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use IO_FIELDS
    use TLAB_VARS, only : imax, jmax, kmax
    use TLAB_VARS, only : rtime, itime
    use TLAB_VARS, only : visc, froude, rossby, prandtl, mach, gama0
    use TLAB_VARS, only : imode_eqns
    use TLAB_CONSTANTS, only : sizeofint, sizeofreal
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_npro_i, ims_pro_i, ims_pro, ims_comm_x, ims_pro_k
#endif
    implicit none
    integer(wi),                                          intent(in)          :: nfield
    real(wp),     dimension(imax,jmax,it_save+1, nfield), intent(in)          :: avg
    integer(wi),                                          intent(in)          :: it_save
    character(len=*),                                     intent(in)          :: basename
    integer(wi),                                          intent(in),optional :: avg_start

    character(len=128)                              :: name
    character(len=32)                               :: varname(1)
    integer(wi)                                     :: sizes(5)
    integer(wi),       parameter                    :: isize_max = 20
    real(wp)                                        :: params(isize_max)
    integer(wi)                                     :: isize, iheader, header_offset, id, ifld
    character(len=10)                               :: start, end, fld_id
  
#ifdef USE_MPI
    nx_total = imax*ims_npro_i
#else
    nx_total = imax
#endif
    ny_total = jmax
    nz_total = it_save + 1

    isize = 0
    isize = isize + 1; params(isize) = rtime
    isize = isize + 1; params(isize) = visc ! inverse of reynolds
    if (iheader == IO_SCAL) then
        isize = isize + 1 + 1                     ! prepare space for schmidt and damkohler

    else if (iheader == IO_FLOW) then
        isize = isize + 1; params(isize) = froude
        isize = isize + 1; params(isize) = rossby
        if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
            isize = isize + 1; params(isize) = gama0
            isize = isize + 1; params(isize) = prandtl
            isize = isize + 1; params(isize) = mach
        end if
    end if
    ! INITIALIZATION OF MPI TYPES SHOULD ONLY BE CARRIED OUT ONCE *AND* DURING INITIALIZATION
    header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

    id = IO_AVERAGE_PLANE
    io_aux(id)%offset = header_offset
    io_aux(id)%precision = IO_TYPE_DOUBLE

    ! Define the sizes array
    sizes(1) = imax * jmax * (it_save+1) ! total size of the 3D field
    sizes(2) = 1                         ! lower bound (starting index in the data array)
    sizes(3) = sizes(1)                  ! upper bound (ending index in the data array)
    sizes(4) = 1                         ! stride (1 means every element)
    sizes(5) = 1                         ! number of variables
    
    if (present(avg_start)) then
      write(start, '(I10)') (avg_start)
      write(end,   '(I10)') itime
    else
      write(start, '(I10)') (itime - it_save*phAvg%stride)
      write(end,   '(I10)') itime
    end if

    do ifld = 1, nfield
        write(fld_id,   '(I10)') ifld
        varname(1) = ''
        if (start == end) then ! write single iteration
          name =  trim(adjustl(basename)) // trim(adjustl(start)) // '.' // trim(adjustl(fld_id))
        else ! write multiple iteration including phase average
          name =  trim(adjustl(basename)) // trim(adjustl(start)) &
          // '_' // trim(adjustl(end)) // '.' // trim(adjustl(fld_id))
        end if

#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
#define LOC_STATUS "unknown"
#define LOC_UNIT_ID 75
#include "dns_open_file.h"
            call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, itime, params)
            close(LOC_UNIT_ID)
#ifdef USE_MPI
        end if
#endif
      call IO_WRITE_SUBARRAY(io_aux(id), name, varname, avg(:,:,ifld,:), sizes)
    end do
    return
  end subroutine PhaseAvg_Write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_ResetVariable()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_CONSTANTS, only : wp
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

#ifdef USE_MPI
    if (ims_pro_k == 0) then
#endif
    avg_flow   = 0.0_wp
    avg_stress = 0.0_wp
    avg_p      = 0.0_wp
    avg_scal   = 0.0_wp
#ifdef USE_MPI
    end if
#endif
  end subroutine PhaseAvg_ResetVariable
end module
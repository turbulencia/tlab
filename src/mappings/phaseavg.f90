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

  interface PhaseAvg_Space 
    module procedure PhaseAvg_Space_FieldPtr, PhaseAvg_Space_index
  end interface PhaseAvg_Space


  real(wp),   dimension(:),       allocatable, target   :: avg_flow, avg_stress, avg_p, avg_scal
  integer(wi)                                           :: nxy, nxz, nyz, nz_total
  integer(wi)                                           :: avg_planes
  character(len=32), parameter                          :: avgu_name       = 'avg_flow'
  character(len=32), parameter                          :: avgstr_name     = 'avg_stress'
  character(len=32), parameter                          :: avgp_name       = 'avg_p'
  character(len=32), parameter                          :: avgs_name       = 'avg_scal'

  public :: PhaseAvg_Space
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Allocate(C_FILE_LOC, restart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_VARS,      only : g
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS,  only : ims_size_i, ims_size_k, ims_pro, ims_npro_i, ims_pro_k
#endif

    implicit none
    character(len=*), intent(in)    :: C_FILE_LOC
    integer(wi)     , intent(in)    :: restart
      ! ================================================================== !

    nxy = imax*jmax
    nyz = jmax*kmax
    nxz = imax*kmax

    if ( restart == -1 ) then ! used for calls from outside dns_main / dns.x 
      avg_planes = 0 
    elseif (mod(restart,phAvg%stride) == 0) then
      avg_planes = restart/phAvg%stride
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
      call TLAB_STOP(DNS_ERROR_PHASEAVG)
    end if
    
#ifdef USE_MPI
    if (ims_pro_k == 0) then
#endif
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_flow,      [imax*jmax*(avg_planes+1)*inb_flow], 'avgflow.')
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_stress,    [imax*jmax*(avg_planes+1)*       6], 'avgstr.' ) ! allocated not yet coded
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_p,         [imax*jmax*(avg_planes+1)*       1], 'avgp.'   )
      call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_scal,      [imax*jmax*(avg_planes+1)*inb_scal], 'avgscal.')

      avg_flow(:)   = 0.0_wp
      avg_stress(:) = 0.0_wp
      avg_p(:)      = 0.0_wp
      avg_scal(:)   = 0.0_wp

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

  subroutine PhaseAvg_Space_FieldPtr(localsum, nfield, itime, it_first, it_save, field)
    implicit none
    real(wp), dimension(imax, jmax)                           , intent(inout)        :: localsum
    integer(wi)                                               , intent(in)           :: nfield
    integer(wi)                                               , intent(in)           :: itime, it_first, it_save
    real(wp), pointer                                         , intent(in)           :: field(:)
  
    integer :: index_loc = 4 ! Index needs to be set to appropriate value for pressure. Needed later for if statement

    call PhaseAvg_Space_Exec(localsum, nfield, itime, it_first, it_save, index_loc, field)
  end subroutine PhaseAvg_Space_FieldPtr


  subroutine PhaseAvg_Space_index(localsum, nfield, itime, it_first, it_save, index)
    implicit none
    real(wp), dimension(imax,jmax)                            , intent(inout)        :: localsum
    integer(wi)                                               , intent(in)           :: nfield
    integer(wi)                                               , intent(in)           :: itime, it_first, it_save, index
    real(wp), pointer                                                                :: field_loc(:) => null()

    call PhaseAvg_Space_Exec(localsum, nfield, itime, it_first, it_save, index, field_loc)    
  end subroutine PhaseAvg_Space_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Space_Exec(localsum, nfield, itime, it_first, it_save, index, field)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_VARS, only : imax, jmax, kmax, isize_field
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif
    use TLAB_ARRAYS, only: q, s
    use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

    implicit none
    real(wp), dimension(imax*jmax)                            , intent(inout)        :: localsum
    integer(wi)                                               , intent(in)           :: nfield
    integer(wi)                                               , intent(in)           :: itime, it_first, it_save, index
    real(wp), dimension(imax, jmax, kmax, 1), target          , intent(in)           :: field

    integer(wi)                                               :: k, ifld, plane_id
    real(wp), dimension(:), pointer                           :: avg_type
    real(wp), dimension(:), pointer                           :: loc_field
    integer(wi)                                               :: ipl_srt, ipl_end, iavg_srt, iavg_end, lpl_srt, lpl_end
      ! ================================================================== !
      ! Calculation of hte plane id to write the spatial average 
    plane_id = 1
    if (it_save /= 0) plane_id = mod((itime-1) - (it_first), it_save) + 1

      ! Determing the tendency to be written
    if (index == 1) then
      avg_type => avg_flow
      call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 2) then
      avg_type => avg_scal
      call c_f_pointer(c_loc(s), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 4) then
      avg_type => avg_p
      call c_f_pointer(c_loc(field), loc_field, shape=[imax*jmax*kmax*nfield])
    elseif (index == 5) then
      avg_type => avg_stress
      call c_f_pointer(c_loc(q), loc_field, shape=[imax*jmax*kmax*nfield])
      ! Not yet coded
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in PhaseAvg_Space_Exec')
      call TLAB_STOP(DNS_ERROR_PHASEAVG)
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
        
        ! Computing the start and end of the field for accumulating the space averages
        iavg_srt = (ifld - 1)*nxy*(avg_planes+1) + nxy*(plane_id-1) + 1
        iavg_end = (ifld - 1)*nxy*(avg_planes+1) + nxy*plane_id  
#ifdef USE_MPI
        if (ims_pro_k == 0) then
          call MPI_Reduce(localsum, avg_type(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_type(imax*jmax*restarts*fld)
        else
          call MPI_Reduce(localsum, MPI_IN_PLACE,                nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        end if
#else
        avg_type(iavg_srt:iavg_end) = localsum
#endif
      
        lpl_srt = (ifld - 1)*nxy*(avg_planes+1) + nxy*avg_planes + 1
        lpl_end = ifld*nxy*(avg_planes+1)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
          avg_type(lpl_srt:lpl_end) = avg_type(lpl_srt:lpl_end) + avg_type(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI 
        endif 
#endif 
      end do
    end if
  return
  end subroutine PhaseAvg_Space_Exec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_Write(nfield, iheader, it_save, basename, index, avg_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use IO_FIELDS
    use TLAB_VARS, only : imax, jmax
    use TLAB_VARS, only : rtime, itime
    use TLAB_VARS, only : visc, froude, rossby, prandtl, mach, gama0
    use TLAB_VARS, only : imode_eqns
    use TLAB_CONSTANTS, only : sizeofint, sizeofreal
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_x, ims_err, ims_npro_i, ims_pro_i, ims_pro, ims_comm_z, ims_pro_k
#endif
    implicit none
    integer(wi),                                          intent(in)          :: nfield
    integer(wi),                                          intent(in)          :: it_save
    character(len=*),                                     intent(in)          :: basename
    integer(wi),                                          intent(in)          :: index
    integer(wi),                                          intent(in),optional :: avg_start

    real(wp), dimension(:), pointer                 :: avg_type
    character(len=128)                              :: name
    character(len=32)                               :: varname(1)
    integer(wi)                                     :: sizes(5)
    integer(wi),       parameter                    :: isize_max = 20
    real(wp)                                        :: params(isize_max)
    integer(wi)                                     :: isize, iheader, header_offset, id, ifld, ifld_srt, ifld_end
    character(len=10)                               :: start, end, fld_id

    if (index == 1) then
      avg_type => avg_flow
    elseif (index == 2) then
      avg_type => avg_scal
    elseif (index == 4) then
      avg_type => avg_p
    elseif (index == 5) then
      avg_type => avg_stress
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in PhaseAvg_Write')
      call TLAB_STOP(DNS_ERROR_PHASEAVG)
    end if

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
    sizes(1) = imax * jmax * (avg_planes+1) ! size of the 3D field
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
        ifld_srt = (ifld-1)*nxy*(avg_planes+1) + 1
        ifld_end = ifld*nxy*(avg_planes+1)
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
            call IO_WRITE_HEADER(LOC_UNIT_ID, isize, g(1)%size, g(2)%size, nz_total, itime, params)
            close(LOC_UNIT_ID)
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            call MPI_BARRIER(ims_comm_x,ims_err)
            call IO_WRITE_SUBARRAY(io_aux(id), name, varname, avg_type(ifld_srt:ifld_end), sizes)
#ifdef USE_MPI
        end if
#endif
    end do
    return
  end subroutine PhaseAvg_Write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PhaseAvg_ResetVariable()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_CONSTANTS, only : wp
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_x, ims_err, ims_pro, ims_pro_k
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
  end subroutine PhaseAvg_ResetVariable
end module

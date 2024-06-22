#include "types.h"
#include "dns_const_mpi.h"
#include "dns_error.h"
#include "dns_const.h"
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

module DNS_PHASEAVG
  
  USE TLAB_PROCS
  USE TLAB_TYPES
  use TLAB_CONSTANTS, ONLY : wp, wi
  use TLAB_VARS,      ONLY : imax, jmax, kmax, g
  use TLAB_VARS,      ONLY : rtime, itime
  use TLAB_VARS,      ONLY : visc, froude, rossby, prandtl, mach, gama0
  use TLAB_VARS,      ONLY : imode_eqns
  use TLAB_CONSTANTS, ONLY : sizeofint, sizeofreal
  use DNS_LOCAL,      ONLY : nitera_first, nitera_save

  implicit none

  real(wp),   dimension(:,:,:), allocatable     :: avg_u, avg_v, avg_w, avg_p, avg_uu, avg_uv, avg_uw, avg_vv, avg_vw, avg_ww
  integer(wi)                                   :: nx_total, ny_total, nz_total
  integer(wi)                                   :: nxy, nxz, nyz
  character(len=32), parameter                  :: avgu_name       = 'avgu'
  character(len=32), parameter                  :: avgv_name       = 'avgv'
  character(len=32), parameter                  :: avgw_name       = 'avgw'
  character(len=32), parameter                  :: avguu_name      = 'avguu'
  character(len=32), parameter                  :: avguv_name      = 'avguv'
  character(len=32), parameter                  :: avguw_name      = 'avguw'
  character(len=32), parameter                  :: avgvv_name      = 'avgvv'
  character(len=32), parameter                  :: avgvw_name      = 'avgvw'
  character(len=32), parameter                  :: avgww_name      = 'avgww'
  character(len=32), parameter                  :: avgp_name       = 'avgp'
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine DNS_AVG_ALLOCATE(C_FILE_LOC, restart)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS,  only : ims_size_i, ims_size_k, ims_pro, ims_npro_i
#endif

    implicit none

    character(len=*), intent(in)    :: C_FILE_LOC
    integer(wi)     , intent(in)    :: restart

    integer(wi)                     :: isize_planej, isize_planek

#ifdef USE_MPI    
    nx_total = imax*ims_npro_i
#else
    nx_total = imax
#endif
    nxy = imax*jmax
    nyz = jmax*kmax
    nxz = imax*kmax

    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_u,     [imax,jmax,restart+1], 'avg_u'    )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_v,     [imax,jmax,restart+1], 'avg_v'    )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_w,     [imax,jmax,restart+1], 'avg_w'    )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_uu,    [imax,jmax,restart+1], 'avg_uu'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_uv,    [imax,jmax,restart+1], 'avg_uv'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_uw,    [imax,jmax,restart+1], 'avg_uw'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_vv,    [imax,jmax,restart+1], 'avg_vv'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_vw,    [imax,jmax,restart+1], 'avg_vw'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_ww,    [imax,jmax,restart+1], 'avg_ww'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_p,     [imax,jmax,restart+1], 'avg_p'    )

    avg_u = 0.0_wp
    avg_v = 0.0_wp
    avg_w = 0.0_wp
    avg_uu = 0.0_wp
    avg_uv = 0.0_wp
    avg_uw = 0.0_wp
    avg_vv = 0.0_wp
    avg_vw = 0.0_wp
    avg_ww = 0.0_wp
    avg_p = 0.0_wp

    return
  end subroutine DNS_AVG_ALLOCATE
    
  subroutine DNS_PHASEAVG_INITILIZE()
    use IO_FIELDS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, ONLY : ims_comm_x, ims_pro_k, ims_pro
#endif
    use TLAB_VARS, ONLY : rtime
    use TLAB_VARS, ONLY : visc, froude, rossby, prandtl, mach, gama0
    use TLAB_VARS, ONLY : imode_eqns

    implicit none

    integer(wi)                                     :: isize, iheader, header_offset, id
    integer(wi),       parameter                    :: isize_max = 20
    real(wp)                                        :: params(isize_max)
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
  end subroutine DNS_PHASEAVG_INITILIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DNS_SPACE_AVG(field, avg, localsum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_VARS, ONLY : imax, jmax, kmax
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

    implicit none
    real(wp), dimension(imax, jmax, kmax),           intent(in)           :: field
    real(wp), dimension(imax, jmax, nitera_save+1),  intent(inout)        :: avg
    real(wp), dimension(imax, jmax)                                       :: localsum
    integer(wi)                                                           :: plane_id, k
      ! ================================================================== !
    localsum = 0.0_wp
    plane_id = mod(itime - nitera_first-1, nitera_save) + 1
    do k = 1, kmax
        localsum = localsum + field(:,:,k)/g(3)%size
    end do

#ifdef USE_MPI
    if (ims_pro_k == 0) then
      call MPI_Reduce(localsum, avg(:,:,plane_id), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)

    else
      ! using 2nd plane of local sum (usually wrk2d) as dummy array on processes 
      call MPI_Reduce(localsum,MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
    end if
#else
    avg(:,:,plane_id) = localsum
#endif
#ifdef USE_MPI
  if ( ims_pro_k == 0) &
    avg(:,:,nitera_save+ 1) = avg(:,:,nitera_save+1) + avg(:,:,plane_id)/nitera_save
#else
  avg(:,:,nitera_save+ 1) = avg(:,:,nitera_save+1) + avg(:,:,plane_id)/nitera_save
#endif
  return
  end subroutine DNS_SPACE_AVG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DNS_WRITE_AVG(avg, iheader, nfield, name, index)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use IO_FIELDS
    use TLAB_VARS, ONLY : imax, jmax, kmax
    use TLAB_VARS, ONLY : rtime
    use TLAB_VARS, ONLY : visc, froude, rossby, prandtl, mach, gama0
    use TLAB_VARS, ONLY : imode_eqns
    use TLAB_CONSTANTS, ONLY : sizeofint, sizeofreal
    use PLANES, ONLY : kplanes
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_npro_i, ims_pro_i, ims_pro, ims_comm_x, ims_pro_k
#endif
    implicit none
    real(wp),     dimension(imax,jmax,nitera_save+1), intent(in)      :: avg
    integer(wi),                                      intent(in)      :: nfield
    character(len=*),                                 intent(in)      :: name

    character(len=128)                              :: fname
    character(len=32)                               :: varname(1)
    integer(wi)                                     :: sizes(5)
    integer(wi),       parameter                    :: isize_max = 20
    real(wp)                                        :: params(isize_max)
    integer(wi)                                     :: isize, iheader, header_offset, id, index
    character(len=10)                               :: start, end, ind
  
#ifdef USE_MPI
    nx_total = imax*ims_npro_i
#else
    nx_total = imax
#endif
    ny_total = jmax
    nz_total = nitera_save

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
    sizes(1) = imax * jmax * nitera_save ! total size of the 3D field
    sizes(2) = 1                         ! lower bound (starting index in the data array)
    sizes(3) = sizes(1)                  ! upper bound (ending index in the data array)
    sizes(4) = 1                         ! stride (1 means every element)
    sizes(5) = 1                         ! number of variables
    
    write(start, '(I10)') (itime - nitera_save)
    write(end,   '(I10)') itime

    if (index == 1) then
      varname(1) = ''
      fname =  trim(adjustl(name)) // trim(adjustl(start)) // '_' // trim(adjustl(end))

#ifdef USE_MPI
      if (ims_pro == 0) then
#endif
        call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, itime, params)
#ifdef USE_MPI
      end if
#endif

      call IO_WRITE_SUBARRAY(io_aux(id), fname, varname, avg, sizes)
    else if ((index == 4) .or. (index ==  2)) then
      varname(1) = ''
      fname =  trim(adjustl(name)) // trim(adjustl(start)) // '_' // trim(adjustl(end))
#ifdef USE_MPI
      if (ims_pro == 0) then
#endif
        call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, itime, params)
#ifdef USE_MPI
      end if
#endif          
      call IO_WRITE_SUBARRAY(io_aux(id), fname, varname, avg, sizes)
    end if
    return
  end subroutine DNS_WRITE_AVG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DNS_RESET_VARIABLE()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_CONSTANTS, only : wp
    avg_u = 0.0_wp
    avg_v = 0.0_wp
    avg_w = 0.0_wp
    avg_uu = 0.0_wp
    avg_uv = 0.0_wp
    avg_uw = 0.0_wp
    avg_vv = 0.0_wp
    avg_vw = 0.0_wp
    avg_ww = 0.0_wp
    avg_p = 0.0_wp
  end subroutine DNS_RESET_VARIABLE   
end module
#include "types.h"
#include "dns_const_mpi.h"
#include "dns_error.h"
#include "dns_const.h"
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

module PHASEAVG
  
  USE TLAB_PROCS
  USE TLAB_TYPES
  use TLAB_CONSTANTS, only : wp, wi
  use TLAB_VARS,      only : imax, jmax, kmax, g
  use TLAB_VARS,      only : rtime
  use TLAB_VARS,      only : visc, froude, rossby, prandtl, mach, gama0
  use TLAB_VARS,      only : imode_eqns
  use TLAB_VARS,      only : inb_flow, inb_scal
  use TLAB_VARS,      only : phaseaverage
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
  subroutine AVG_ALLOCATE(C_FILE_LOC, restart)
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
    if (mod(restart,phaseaverage%stride) == 0) then
      num_strides = restart/phaseaverage%stride
    else
      call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
      call TLAB_STOP(DNS_ERROR_PHASEAVG)
    end if
    
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_flow,      [imax,jmax,num_strides+1,inb_flow], 'avgflow.')
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_stress,    [imax,jmax,num_strides+1,6],        'avgstr.' )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_p,         [imax,jmax,num_strides+1,1],        'avgp.'   )
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   avg_scal,      [imax,jmax,num_strides+1,inb_scal], 'avgscal.')

    avg_flow   = 0.0_wp
    avg_stress = 0.0_wp
    avg_p      = 0.0_wp
    avg_scal   = 0.0_wp
    return
  end subroutine AVG_ALLOCATE
    
  subroutine PHASEAVG_INITIALIZE()
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
  end subroutine PHASEAVG_INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SPACE_AVG(field, avg, nfield, localsum, itime, it_first, it_save, index)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_VARS, only : imax, jmax, kmax
#ifdef USE_MPI 
    use MPI
    use TLAB_MPI_VARS,  only : ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

    implicit none
    real(wp), dimension(imax, jmax, kmax, *)                  , intent(in)           :: field
    real(wp), dimension(imax, jmax, num_strides+1, nfield)    , intent(inout)        :: avg
    integer(wi)                                               , intent(in)           :: nfield
    integer(wi)                                               , intent(in)           :: itime, it_first, it_save
    real(wp), dimension(imax, jmax)                                                  :: localsum

    integer(wi)                                               :: k, fld, ifld, jfld, plane_id, index
      ! ================================================================== !
    if ((index == 1) .or. (index == 2) .or. (index == 4)) then
      do ifld = 1, nfield
          localsum = 0.0_wp
          plane_id = ((mod((itime-1) - (it_first), it_save) + 1)/phaseaverage%stride)
          do k = 1, kmax
            localsum = localsum + field(:,:,k,ifld)/g(3)%size
          end do

#ifdef USE_MPI
          if (ims_pro_k == 0) then
            call MPI_Reduce(localsum, avg(:,:,plane_id, ifld), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
          else
            call MPI_Reduce(localsum, MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
          end if
#else
            avg(:,:,plane_id, ifld) = localsum
#endif
#ifdef USE_MPI
          if ( ims_pro_k == 0) &
              avg(:,:,it_save + 1, ifld) = avg(:,:,it_save+1, ifld) + avg(:,:,plane_id, ifld)/it_save
#else
          avg(:,:,it_save+ 1, ifld) = avg(:,:,it_save+1, ifld) + avg(:,:,plane_id, ifld)/it_save
#endif
      end do
    else if (index == 5) then
      do ifld = 1, inb_flow
        do jfld = 1, inb_flow
          if (ifld <= jfld) then
            fld = ifld + jfld
            localsum = 0.0_wp
            plane_id = mod(itime - (it_first), it_save+1)
            do k = 1, kmax
              localsum = localsum + (field(:,:,k,ifld)*field(:,:,k,jfld))/(g(3)%size)
            end do
#ifdef USE_MPI
            if (ims_pro_k == 0) then
              call MPI_Reduce(localsum, avg(:,:,plane_id, fld), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
            else
              call MPI_Reduce(localsum, MPI_IN_PLACE, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
            end if
#else
            avg(:,:,plane_id, fld) = localsum
#endif
#ifdef USE_MPI
            if ( ims_pro_k == 0) &
                avg(:,:,it_save + 1, fld) = avg(:,:,it_save + 1, fld) + avg(:,:,plane_id, fld)/it_save
#else
            avg(:,:,it_save + 1, fld) = avg(:,:,it_save + 1, fld) + avg(:,:,plane_id, fld)/it_save
#endif
          end if
        end do
      end do
    end if
  return
  end subroutine SPACE_AVG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine WRITE_AVG(nfield, avg, iheader, it_save, name)
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
    integer(wi),                                          intent(in)      :: nfield
    real(wp),     dimension(imax,jmax,it_save+1, nfield), intent(in)      :: avg
    integer(wi),                                          intent(in)      :: it_save
    character(len=*),                                     intent(in)      :: name

    character(len=128)                              :: fname
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
    nz_total = it_save

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
    
    write(start, '(I10)') (itime - it_save)
    write(end,   '(I10)') itime

    do ifld = 1, nfield
        write(fld_id,   '(I10)') ifld
        varname(1) = ''
        fname =  trim(adjustl(name)) // '_' // trim(adjustl(start)) &
        // '_' // trim(adjustl(end)) // '.' // trim(adjustl(fld_id))

#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
            call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, itime, params)
#ifdef USE_MPI
        end if
#endif

        call IO_WRITE_SUBARRAY(io_aux(id), fname, varname, avg(:,:,:,ifld), sizes)
    end do
    return
  end subroutine WRITE_AVG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine RESET_VARIABLE()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use TLAB_CONSTANTS, only : wp
    avg_flow   = 0.0_wp
    avg_stress = 0.0_wp
    avg_p      = 0.0_wp
    avg_scal   = 0.0_wp
  end subroutine RESET_VARIABLE  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine STRIDE_CHECK(itime_size_max, itime_vec, stride)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(wi)     , intent(in)      :: itime_size_max
    integer(wi)     , intent(in)      :: itime_vec(itime_size_max)
    integer(wi)     , intent(in)      :: stride

    integer(wi)                       :: i

    do i = 1, itime_size_max - 1
      if (itime_vec(i) /= stride) then
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Stride does not match the plane integer list')
        call TLAB_STOP(DNS_ERROR_PHASEAVG)
      end if
    end do
  end subroutine STRIDE_CHECK
end module
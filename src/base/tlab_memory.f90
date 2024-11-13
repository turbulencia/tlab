! ###################################################################
! ###################################################################
module TLab_Arrays
    use TLab_Constants, only: wp
    implicit none
    save

    real(wp), allocatable :: x(:, :), y(:, :), z(:, :)      ! Grid and associated arrays
    real(wp), allocatable :: q(:, :)                        ! Eulerian fields, flow vartiables
    real(wp), allocatable :: s(:, :)                        ! Eulerian fields, scalar variables
    real(wp), allocatable :: txc(:, :)                      ! Temporary space for Eulerian fields
    real(wp), allocatable :: wrk1d(:, :)                    ! Work arrays (scratch space)
    real(wp), allocatable :: wrk2d(:, :)                    ! Work arrays (scratch space)
    real(wp), allocatable :: wrk3d(:)                       ! Work arrays (scratch space)
    real(wp), allocatable :: wrkdea(:, :)                   ! Work arrays for dealiasing (scratch space)

    target q, s, txc, wrk1d, wrk2d, wrk3d, wrkdea

end module TLab_Arrays

! ###################################################################
! ###################################################################
module TLab_Pointers
    use TLab_Constants, only: wp
    implicit none

    type :: pointers_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:)
    end type pointers_dt

    real(wp), pointer :: u(:) => null()
    real(wp), pointer :: v(:) => null()
    real(wp), pointer :: w(:) => null()
    real(wp), pointer :: e(:) => null()
    real(wp), pointer :: rho(:) => null()
    real(wp), pointer :: p(:) => null()
    real(wp), pointer :: T(:) => null()
    real(wp), pointer :: vis(:) => null()

    real(wp), pointer :: tmp1(:) => null()
    real(wp), pointer :: tmp2(:) => null()
    real(wp), pointer :: tmp3(:) => null()
    real(wp), pointer :: tmp4(:) => null()
    real(wp), pointer :: tmp5(:) => null()
    real(wp), pointer :: tmp6(:) => null()
    real(wp), pointer :: tmp7(:) => null()
    real(wp), pointer :: tmp8(:) => null()
    real(wp), pointer :: tmp9(:) => null()

end module TLab_Pointers

! ###################################################################
! ###################################################################
module TLab_Pointers_3D
    use TLab_Constants, only: wp
    implicit none

    type :: pointers3d_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:, :, :)
    end type pointers3d_dt

    real(wp), pointer :: u(:, :, :) => null()
    real(wp), pointer :: v(:, :, :) => null()
    real(wp), pointer :: w(:, :, :) => null()
    real(wp), pointer :: e(:, :, :) => null()
    real(wp), pointer :: rho(:, :, :) => null()
    real(wp), pointer :: p(:, :, :) => null()
    real(wp), pointer :: T(:, :, :) => null()
    real(wp), pointer :: vis(:, :, :) => null()

    real(wp), pointer :: p_q(:, :, :, :) => null()
    real(wp), pointer :: p_s(:, :, :, :) => null()
    real(wp), pointer :: p_wrk1d(:, :) => null()
    real(wp), pointer :: p_wrk2d(:, :, :) => null()
    real(wp), pointer :: p_wrk3d(:, :, :) => null()

    real(wp), pointer :: tmp1(:, :, :) => null()
    real(wp), pointer :: tmp2(:, :, :) => null()
    real(wp), pointer :: tmp3(:, :, :) => null()
    real(wp), pointer :: tmp4(:, :, :) => null()
    real(wp), pointer :: tmp5(:, :, :) => null()
    real(wp), pointer :: tmp6(:, :, :) => null()
    real(wp), pointer :: tmp7(:, :, :) => null()
    real(wp), pointer :: tmp8(:, :, :) => null()
    real(wp), pointer :: tmp9(:, :, :) => null()

end module TLab_Pointers_3D

! ###################################################################
! ###################################################################
module TLab_Pointers_C
    use TLab_Constants, only: wp
    implicit none

    complex(wp), pointer :: c_wrk1d(:, :) => null()
    complex(wp), pointer :: c_wrk3d(:, :) => null()

end module TLab_Pointers_C

! ###################################################################
! ###################################################################
#include "dns_const.h"
#include "dns_error.h"

module TLab_Memory
    use TLab_Constants, only: sp, wp, wi, longi, lfile, efile
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: isize_field, inb_flow_array, inb_scal_array
    use TLAB_VARS, only: isize_txc_field, inb_txc, isize_txc_dimx, isize_txc_dimz
    use TLAB_VARS, only: isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d
    use TLAB_VARS, only: Dealiasing
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private
    save

    character*128 :: str, line
    integer :: ierr

    interface TLAB_ALLOCATE_LOG
        module procedure TLAB_ALLOCATE_LOG_SHORT, TLAB_ALLOCATE_LOG_LONG
    end interface TLAB_ALLOCATE_LOG

#ifdef NO_ASSUMED_RANKS
    interface TLab_Allocate_SINGLE
        module procedure TLab_Allocate_SINGLE1, TLab_Allocate_SINGLE2, TLab_Allocate_SINGLE3, TLab_Allocate_SINGLE4
    end interface TLab_Allocate_SINGLE

    interface TLab_Allocate_Real
        module procedure TLab_Allocate_Real1, TLab_Allocate_Real2, TLab_Allocate_Real3, TLab_Allocate_Real4
    end interface TLab_Allocate_Real

    interface TLab_Allocate_INT
        module procedure TLab_Allocate_INT1, TLab_Allocate_INT2, TLab_Allocate_INT3, TLab_Allocate_INT4
    end interface TLab_Allocate_INT

    interface TLab_Allocate_LONG_INT
        module procedure TLab_Allocate_LONG_INT1, TLab_Allocate_LONG_INT2, TLab_Allocate_LONG_INT3, TLab_Allocate_LONG_INT4
    end interface TLab_Allocate_LONG_INT
#endif

    public :: TLab_Initialize_Memory
    public :: TLab_Allocate_SINGLE
    public :: TLab_Allocate_Real
    public :: TLab_Allocate_INT
    public :: TLab_Allocate_LONG_INT
    public :: Tlab_Allocate_Real_Long
contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Initialize_Memory(C_FILE_LOC)
        use TLab_Arrays

        character(len=*), intent(in) :: C_FILE_LOC

! loop counters over the whole domain are integer*4
        if (isize_txc_field > huge(imax)) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Integer model of 4 bytes is not big enough.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call TLab_Allocate_Real(C_FILE_LOC, q, [isize_field, inb_flow_array], 'flow')
        call TLab_Allocate_Real(C_FILE_LOC, s, [isize_field, inb_scal_array], 'scal')

        call TLab_Allocate_Real(C_FILE_LOC, txc, [isize_txc_field, inb_txc], 'txc')
        call TLab_Allocate_Real(C_FILE_LOC, wrk1d, [isize_wrk1d, inb_wrk1d], 'wrk1d')
        call TLab_Allocate_Real(C_FILE_LOC, wrk2d, [isize_wrk2d, inb_wrk2d], 'wrk2d')
        call TLab_Allocate_Real(C_FILE_LOC, wrk3d, [isize_wrk3d], 'wrk3d')

        if (any(Dealiasing(:)%type /= DNS_FILTER_NONE)) then
            call TLab_Allocate_Real(C_FILE_LOC, wrkdea, [isize_field, 2], 'wrk-dealiasing')
        end if

        call TLab_Set_Pointers()

        call TLab_Set_Pointers_3D()

        call TLab_Set_Pointers_C()

        return
    end subroutine TLab_Initialize_Memory

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Set_Pointers()
        use TLab_Arrays
        use TLab_Pointers

        integer(wi) idummy(2)

        idummy = shape(q)
        if (idummy(2) >= 1) u(1:isize_field) => q(1:isize_field, 1)
        if (idummy(2) >= 2) v(1:isize_field) => q(1:isize_field, 2)
        if (idummy(2) >= 3) w(1:isize_field) => q(1:isize_field, 3)
        ! compressible flows variables
        if (idummy(2) >= 4) e(1:isize_field) => q(1:isize_field, 4)
        if (idummy(2) >= 5) rho(1:isize_field) => q(1:isize_field, 5)
        if (idummy(2) >= 6) p(1:isize_field) => q(1:isize_field, 6)
        if (idummy(2) >= 7) T(1:isize_field) => q(1:isize_field, 7)
        if (idummy(2) >= 8) vis(1:isize_field) => q(1:isize_field, 8)

        idummy = shape(txc)
        if (idummy(2) >= 1) tmp1(1:isize_field) => txc(1:isize_field, 1)
        if (idummy(2) >= 2) tmp2(1:isize_field) => txc(1:isize_field, 2)
        if (idummy(2) >= 3) tmp3(1:isize_field) => txc(1:isize_field, 3)
        if (idummy(2) >= 4) tmp4(1:isize_field) => txc(1:isize_field, 4)
        if (idummy(2) >= 5) tmp5(1:isize_field) => txc(1:isize_field, 5)
        if (idummy(2) >= 6) tmp6(1:isize_field) => txc(1:isize_field, 6)
        if (idummy(2) >= 7) tmp7(1:isize_field) => txc(1:isize_field, 7)
        if (idummy(2) >= 8) tmp8(1:isize_field) => txc(1:isize_field, 8)
        if (idummy(2) >= 9) tmp9(1:isize_field) => txc(1:isize_field, 9)

        return
    end subroutine TLab_Set_Pointers

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Set_Pointers_3D()
        use TLab_Arrays
        use TLab_Pointers_3D

        integer(wi) idummy(2)

        idummy = shape(q)
        if (idummy(2) >= 1) u(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 1)
        if (idummy(2) >= 2) v(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 2)
        if (idummy(2) >= 3) w(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 3)
        ! compressible flows variables
        if (idummy(2) >= 4) e(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 4)
        if (idummy(2) >= 5) rho(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 5)
        if (idummy(2) >= 6) p(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 6)
        if (idummy(2) >= 7) T(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 7)
        if (idummy(2) >= 8) vis(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 8)

        if (allocated(q)) p_q(1:imax, 1:jmax, 1:kmax, 1:inb_flow_array) => q(1:isize_field*inb_flow_array, 1)
        if (allocated(s)) p_s(1:imax, 1:jmax, 1:kmax, 1:inb_scal_array) => s(1:isize_field*inb_scal_array, 1)
        if (allocated(wrk3d)) p_wrk3d(1:imax, 1:jmax, 1:kmax) => wrk3d(1:isize_field)
        if (allocated(wrk2d)) p_wrk2d(1:imax, 1:kmax, 1:inb_wrk2d) => wrk2d(1:imax*kmax*inb_wrk2d, 1)    ! this is the most common wrk2d dimensions
        if (allocated(wrk1d)) p_wrk1d(1:jmax, 1:inb_wrk1d) => wrk1d(1:jmax*inb_wrk1d, 1)                 ! this is the most common wrk1d dimensions

        idummy = shape(txc)
        if (idummy(2) >= 1) tmp1(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 1)
        if (idummy(2) >= 2) tmp2(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 2)
        if (idummy(2) >= 3) tmp3(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 3)
        if (idummy(2) >= 4) tmp4(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 4)
        if (idummy(2) >= 5) tmp5(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 5)
        if (idummy(2) >= 6) tmp6(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 6)
        if (idummy(2) >= 7) tmp7(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 7)
        if (idummy(2) >= 8) tmp8(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 8)
        if (idummy(2) >= 9) tmp9(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 9)

        return
    end subroutine TLab_Set_Pointers_3D

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Set_Pointers_C()
        use TLab_Arrays
        use TLab_Pointers_C
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

        if (allocated(wrk1d)) call c_f_pointer(c_loc(wrk1d), c_wrk1d, shape=[jmax, inb_wrk1d/2])
        if (allocated(wrk3d)) call c_f_pointer(c_loc(wrk3d), c_wrk3d, shape=[isize_txc_dimz/2, kmax])

        return
    end subroutine TLab_Set_Pointers_C

    ! ######################################################################
    ! ######################################################################
#ifndef NO_ASSUMED_RANKS
    subroutine TLab_Allocate_Real(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(wp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        ! if (any(dims <= 0)) return; better allocate to zero than not allocate; error in supermuc
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        rank (4)
            allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        rank default
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Rank too large while allocating memory space for '//trim(adjustl(s))//'.')
            call TLab_Stop(DNS_ERROR_ALLOC)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLab_Allocate_Real

! ### DOUBLE ALLOCATION ROUTINES FOR LARGE 1D ARRAYS
    subroutine Tlab_Allocate_Real_Long(C_FILE_LOC, a, dims, s)
      character(len=*), intent(in) :: C_FILE_LOC,s
      real(8), allocatable, intent(inout) :: a(:)
      integer(8), intent(in) :: dims(1)
      integer id
      !#####################################################################
      call TLAB_ALLOCATE_LOG_LONG(lfile,dims,s)
      allocate (a(dims(1)), stat=ierr)
      call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine Tlab_Allocate_Real_Long

   ! ### DOUBLE ALLOCATION ROUTINES FOR LARGE 1D ARRAYS
    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1_LONG(C_FILE_LOC, a, dims, s)
      character(len=*), intent(in) :: C_FILE_LOC,s
      real(8), allocatable, intent(inout) :: a(:)
      integer(8), intent(in) :: dims(1)
      integer id
      !#####################################################################
      call TLAB_ALLOCATE_LOG_LONG(lfile,dims,s)
      allocate (a(dims(1)), stat=ierr)
      call TLAB_ALLOCATE_ERR(C_FILE_LOC,efile,s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1_LONG
    

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Allocate_SINGLE(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(sp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLab_Allocate_SINGLE

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Allocate_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLab_Allocate_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLab_Allocate_LONG_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(longi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLab_Allocate_LONG_INT
#endif

#ifdef NO_ASSUMED_RANKS
    ! #######################################################
    ! ### INSTANCES FOR INTERFACES TO ALLOCATION ROUTINES
    ! ###
    ! ### SINGLE ALLOCATION ROUTINES
    subroutine TLab_Allocate_SINGLE1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(sp), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_SINGLE1

    subroutine TLab_Allocate_SINGLE2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(sp), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_SINGLE2

    subroutine TLab_Allocate_SINGLE3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(sp), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_SINGLE3

    subroutine TLab_Allocate_SINGLE4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(sp), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_SINGLE4

    ! ### DOUBLE ALLOCATION ROUTINES
    subroutine TLab_Allocate_Real1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(wp), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_Real1

    subroutine TLab_Allocate_Real2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(wp), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_Real2

    subroutine TLab_Allocate_Real3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(wp), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_Real3

    subroutine TLab_Allocate_Real4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(wp), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_Real4

    ! # INTEGER ALLOCATION ROUTINES
    subroutine TLab_Allocate_INT1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_INT1

    subroutine TLab_Allocate_INT2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        !   if (any(dims <= 0)) return; better allocate to zero than not allocate; error in supermuc
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_INT2

    subroutine TLab_Allocate_INT3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_INT3

    subroutine TLab_Allocate_INT4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_INT4

    ! # LONG INTEGER ALLOCATION ROUTINES
    subroutine TLab_Allocate_LONG_INT1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_LONG_INT1

    subroutine TLab_Allocate_LONG_INT2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_LONG_INT2

    subroutine TLab_Allocate_LONG_INT3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_LONG_INT3

    subroutine TLab_Allocate_LONG_INT4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLab_Allocate_LONG_INT4

#endif

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_ALLOCATE_LOG_SHORT(log_file, dims, s)
        integer(wi), intent(IN) :: dims(:)
        character(len=*), intent(IN) :: log_file, s
        integer id
        !#####################################################################

        if (any(dims < 0)) then
            ierr = DNS_ERROR_ALLOC
            call TLAB_ALLOCATE_ERR('TLAB_ALLOCATE_LOG', efile, s)
        end if

        DO id=1,size(dims) 
            dims_long(id) = dims(id)
        ENDDO
        CALL TLAB_ALLOCATE_LOG_LONG(log_file,dims_long,s)
      
        if (any(dims == 0)) return      ! do not print out lines when allocation a zero-space array

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//' x '//trim(adjustl(str))
        end do
        call TLab_Write_ASCII(log_file, line)
    end subroutine TLAB_ALLOCATE_LOG_SHORT

    subroutine TLAB_ALLOCATE_LOG_LONG(log_file, dims, s)
        integer(longi), intent(IN) :: dims(:)
        character(len=*), intent(IN) :: log_file, s
        integer id
        !#####################################################################

        if (any(dims < 0)) then
            ierr = DNS_ERROR_ALLOC
            call TLAB_ALLOCATE_ERR('TLAB_ALLOCATE_LOG', efile, s)
        end if

        if (any(dims == 0)) return      ! do not print out lines when allocation a zero-space array

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//' x '//trim(adjustl(str))
        end do
        call TLab_Write_ASCII(log_file, line)
    end subroutine TLAB_ALLOCATE_LOG_LONG

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_ALLOCATE_ERR(C_FILE_LOC, log_file, s)
        character(len=*) :: C_FILE_LOC, log_file, s

        if (ierr /= 0) then
            call TLab_Write_ASCII(log_file, C_FILE_LOC//'. Error while allocating memory space for'//trim(adjustl(s))//'.')
            call TLab_Stop(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ERR

end module TLab_Memory

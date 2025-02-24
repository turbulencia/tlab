#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   IO of geometry field eps.
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

subroutine IBM_IO_READ_INT_GEOMETRY(wrk3d, stag)

    use  TLab_Constants, only: wp, wi
    use  TLab_Memory, only: imax, jmax, kmax, isize_field
    use  IBM_VARS, only: epsp, eps, eps_name, epsp_name
    use  IO_Fields, only: IO_Read_Field_INT1
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none

    real(wp), dimension(isize_field), target, intent(inout) :: wrk3d
    logical, intent(in) :: stag

    integer(1), pointer :: int_wrk(:) => null()
    real(wp) :: params(0)
    character(len=32) :: name

    ! ================================================================== !

    if (stag) then
        name = epsp_name
    else
        name = eps_name
    end if

    ! pass memory address from double precision array to int1 array
    call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax*jmax*kmax])
    int_wrk(:) = int(wrk3d(:), 1)

    ! header without params
    ! read eps field as int(1)
    call IO_Read_Field_INT1(name, imax, jmax, kmax, 0, int_wrk, params)

    ! type casting
    if (stag) then
        epsp(:) = real(int_wrk(:), wp)
    else
        eps(:) = real(int_wrk(:), wp)
    end if

    nullify (int_wrk)

    return
end subroutine IBM_IO_READ_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_INT_GEOMETRY(wrk3d, stag)

    use TLab_Constants, only: wp, wi
    use IBM_VARS,       only: epsp, eps, eps_name, epsp_name
    use TLab_Memory, only: isize_field, imax, jmax, kmax
    use IO_Fields,      only: IO_Write_Field_INT1

    implicit none

    integer(1), dimension(isize_field), intent(inout) :: wrk3d
    logical, intent(in) :: stag

    character(len=32) :: name

    ! ================================================================== !

    ! wp to int1
    if (stag) then
        name = epsp_name
        wrk3d(:) = int(epsp(:), 1)
    else
        name = eps_name
        wrk3d(:) = int(eps(:), 1)
    end if

    ! header (offset, nx, ny, nz, nt == 20 byte)
    ! write eps field as int(1)
    call IO_Write_Field_INT1(name, imax, jmax, kmax, 0, wrk3d)

    return
end subroutine IBM_IO_WRITE_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_BIT_GEOMETRY(wrk3d, stag)
    use TLab_Constants, only: wp, wi
    use IBM_VARS, only: epsp, eps, eps_name, epsp_name
    use TLab_Memory, only: isize_field, imax, jmax, kmax
    use IO_Fields,      only: IO_Write_Field_INT1

    implicit none

    integer(1), dimension(isize_field), target, intent(inout) :: wrk3d
    logical, intent(in) :: stag

    integer(1), dimension(:), pointer :: eps_bit

    integer(wi) :: bsize_field, imax_bit
    character(len=32) :: name

    ! ================================================================== !

    ! size of bit-array
    bsize_field = isize_field/8
    imax_bit = imax/8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

    ! assign to scratch
    eps_bit => wrk3d(1:bsize_field)

    ! wp real to bitwise int1
    if (stag) then
        name = epsp_name
        call IBM_IO_R2B(isize_field, bsize_field, epsp, eps_bit)
    else
        name = eps_name
        call IBM_IO_R2B(isize_field, bsize_field, eps, eps_bit)
    end if

    ! header (offset, nx, ny, nz, nt == 20 byte)
    ! write bitwise eps_bit field as int(1)
    call IO_Write_Field_INT1(name, imax_bit, jmax, kmax, 0, eps_bit)

    nullify (eps_bit)

    return
end subroutine IBM_IO_WRITE_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_READ_BIT_GEOMETRY(wrk3d, stag)
    use TLab_Constants, only: wp, wi
    use IBM_VARS, only: epsp, eps, eps_name, epsp_name
    use TLab_Memory, only: isize_field, imax, jmax, kmax
    use IO_Fields, only: IO_Read_Field_INT1
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none

    real(wp), dimension(isize_field), target, intent(inout) :: wrk3d
    logical, intent(in) :: stag

    integer(1), pointer :: int_wrk(:) => null()
    integer(wi) :: bsize_field, imax_bit
    real(wp) :: params(0)
    character(len=32) :: name
    ! ================================================================== !

    if (stag) then
        name = epsp_name
    else
        name = eps_name
    end if

    ! size of bit-array
    bsize_field = isize_field/8
    imax_bit = imax/8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

    ! pass memory address from double precision array to int1 array
    call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax_bit*jmax*kmax])
    int_wrk(:) = int(wrk3d(:), 1)

    ! header without params
    ! read eps field as int(1)
    call IO_Read_Field_INT1(name, imax_bit, jmax, kmax, 0, int_wrk, params)

    ! wp real to bitwise int1
    if (stag) then
        call IBM_IO_B2R(bsize_field, isize_field, int_wrk, epsp)
    else
        call IBM_IO_B2R(bsize_field, isize_field, int_wrk, eps)
    end if

    ! type casting
    ! eps(:) = real(int_wrk(:), dp)

    nullify (int_wrk)

    return
end subroutine IBM_IO_READ_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_R2B(rsize, bsize, r, b)

    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(in) :: rsize
    integer(wi), intent(in) :: bsize
    real(wp), dimension(rsize), intent(in) :: r
    integer(1), dimension(bsize), intent(out) :: b

    integer(wi) :: i, ip, ib
    ! ================================================================== !

    ! initialize bit-array
    b(1:bsize) = int(0, 1)

    ! convert
    do i = 1, bsize
        ip = (i - 1)*8
        do ib = 1, 8
            if (r(ip + ib) > 0.0_wp) then
                b(i) = ibset(b(i), ib - 1)
            end if
        end do
    end do

    return
end subroutine IBM_IO_R2B

!########################################################################

subroutine IBM_IO_B2R(bsize, rsize, b, r)

    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(in) :: bsize
    integer(wi), intent(in) :: rsize
    integer(1), dimension(bsize), intent(in) :: b
    real(wp), dimension(rsize), intent(out) :: r

    integer(wi) :: i, ip, ib
    ! ================================================================== !

    ! initialize real-array
    r(:) = 0.0_wp

    ! convert
    do i = 1, bsize
        ip = (i - 1)*8
        do ib = 1, 8
            if (btest(b(i), ib - 1)) then
                r(ip + ib) = 1.0_wp
            end if
        end do
    end do

    return
end subroutine IBM_IO_B2R

!########################################################################

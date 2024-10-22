#include "dns_error.h"

!########################################################################
!########################################################################
subroutine IO_READ_GRID(name, imax, jmax, kmax, scalex, scaley, scalez, x, y, z, area)
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow

    implicit none

    character*(*) name
    integer(wi), intent(in) :: imax, jmax, kmax
    real(wp) scalex, scaley, scalez
    real(wp) x(imax), y(jmax), z(kmax)
    real(wp), optional :: area

    ! -----------------------------------------------------------------------
    integer(wi) imax_loc, jmax_loc, kmax_loc
    real(wp) scale(3)
    character*(32) line

    ! #######################################################################
    open (50, file=name, status='old', form='unformatted')
    rewind (50)

    ! -----------------------------------------------------------------------
    read (50) imax_loc, jmax_loc, kmax_loc
    read (50) scale
    scalex = scale(1)
    scaley = scale(2)
    scalez = scale(3)

    ! -----------------------------------------------------------------------
    if (imax_loc /= imax .or. jmax_loc /= jmax .or. kmax_loc /= kmax) then
        close (50)
        write (line, 100) imax_loc, jmax_loc, kmax_loc
        call TLAB_WRITE_ASCII(efile, 'IO_READ_GRID. Dimensions ('//trim(line)//') unmatched.')
        call TLAB_STOP(DNS_ERROR_DIMGRID)
    end if

    ! -----------------------------------------------------------------------
    read (50) x
    read (50) y
    read (50) z
    close (50)

    if (present(area)) then
        area = scalex
        if (kmax > 1) area = area*scalez ! 3D case
    end if

    return

100 format(I5, ',', I5, ',', I5)

end subroutine IO_READ_GRID

!########################################################################
!########################################################################
subroutine IO_WRITE_GRID(name, imax, jmax, kmax, scalex, scaley, scalez, x, y, z)
    use TLab_Constants, only: wp, wi
    implicit none

    character*(*) name
    integer(wi) imax, jmax, kmax
    real(wp) scalex, scaley, scalez
    real(wp) x(imax), y(jmax), z(kmax)

    ! -----------------------------------------------------------------------
    real(wp) scale(3)

    !########################################################################
    open (unit=51, file=name, form='unformatted', status='unknown')

    ! -----------------------------------------------------------------------
    scale(1) = scalex
    scale(2) = scaley
    scale(3) = scalez
    write (51) imax, jmax, kmax
    write (51) scale

    ! -----------------------------------------------------------------------
    write (51) x
    write (51) y
    write (51) z
    close (51)

    return
end subroutine IO_WRITE_GRID

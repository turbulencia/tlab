#include "dns_error.h"

module TLab_Grid
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    real(wp), allocatable, target, public :: x(:), y(:), z(:)

    public :: TLab_Grid_Read
    public :: TLab_Grid_Write

contains
!########################################################################
!########################################################################
    subroutine TLab_Grid_Read(name, x, y, z, sizes, scales)
        character*(*) name
        real(wp), allocatable, intent(out) :: x(:), y(:), z(:)
        integer(wi), intent(in), optional :: sizes(3)
        real(wp), intent(out), optional :: scales(3)

        ! -----------------------------------------------------------------------
        integer(wi) locSizes(3)
        real(wp) locScales(3)
        character*(32) line

        ! #######################################################################
        open (50, file=name, status='old', form='unformatted')
        rewind (50)

        ! -----------------------------------------------------------------------
        read (50) locSizes(1:3)

        ! -----------------------------------------------------------------------
        if (present(sizes) .and. any(sizes /= locSizes)) then
            close (50)
            write (line, 100) locSizes
            call TLab_Write_ASCII(efile, __FILE__//'. Dimensions ('//trim(line)//') unmatched.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        read (50) locScaleS(1:3)
        if (present(scales)) scales = locScales
        
        ! -----------------------------------------------------------------------
        allocate (x(locSizes(1)), y(locSizes(2)), z(locSizes(3)))
        read (50) x
        read (50) y
        read (50) z
        close (50)

        return

100     format(I5, ',', I5, ',', I5)

    end subroutine TLab_Grid_Read

!########################################################################
!########################################################################
    subroutine TLab_Grid_Write(name, imax, jmax, kmax, scalex, scaley, scalez, x, y, z)
        character*(*) name
        integer(wi) imax, jmax, kmax
        real(wp) scalex, scaley, scalez
        real(wp), intent(in) :: x(imax), y(jmax), z(kmax)

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
    end subroutine TLab_Grid_Write

end module TLab_Grid

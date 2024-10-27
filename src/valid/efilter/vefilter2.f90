program VEFILTER2

    use TLAB_VARS
    use IO_FIELDS
    use OPR_FILTERS

    implicit none

#include "types.h"

    TREAL, dimension(:, :), allocatable :: x, y, z
    TREAL, dimension(:), allocatable :: a, cx, cy, cz
    TREAL, dimension(:, :), allocatable :: wrk3d
    TINTEGER :: i

! ###################################################################
    call DNS_START

    call TLab_Initialize_Parameters('tlab.ini')

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

    allocate (wrk3d(imax*jmax*kmax, 2), a(imax*jmax*kmax))
    allocate (cx(imax*5), cy(jmax*5), cz(kmax_total*5))

! ###################################################################
    call IO_READ_GRID(gfile, imax, jmax, kmax_total, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z)

    ! CALL FLT4E_INI(g(1)%scale, x, cx)
    ! CALL FLT4E_INI(g(2)%scale, y, cy)
    ! CALL FLT4E_INI(g(3)%scale, z, cz)

    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, 1, a)

!  CALL OPR_FILTER(i4, imax, jmax, kmax,  i1bc, j1bc, k1bc, i1, i1, i1, i1, a, &
!       cx, cy, cz, wrk3d)

    ! to be rewritten with new arrays
    call OPR_FILTER(i3, imax, jmax, kmax, i1bc, j1bc, k1bc, i1, i1, i1, i1, a, &
                    cx, cy, cz, wrk3d(1, 1))
    call OPR_FILTER(i3, imax, jmax, kmax, i1bc, j1bc, k1bc, i1, i1, i1, i1, wrk3d(1, 1), &
                    cx, cy, cz, wrk3d(1, 2))
    do i = 1, imax*jmax*kmax
        wrk3d(i, 2) = wrk3d(i, 2) + C_3_R*(a(i) - wrk3d(i, 1))
    end do
    call OPR_FILTER(i3, imax, jmax, kmax, i1bc, j1bc, k1bc, i1, i1, i1, i1, wrk3d(1, 2), &
                    cx, cy, cz, a)

    call IO_WRITE_FIELDS('field.out', IO_SCAL, imax, jmax, kmax, 1, 1, a)

    stop
end program VEFILTER2

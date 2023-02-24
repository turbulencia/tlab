#include "types.h"
#include "dns_const.h"

subroutine SL_NORMAL_GRADIENT(isl, nmax, istep, kstep, ibuffer_npy, &
                              u, v, w, z1, a, sl, profiles, txc, wrk1d, wrk2d, wrk3d)

    use TLAB_VARS
    use OPR_PARTIAL

    implicit none

#define L_NFIELDS_MAX 1

    TINTEGER isl, nmax, istep, kstep, ibuffer_npy
    TREAL u(*), v(*), w(*), z1(*), a(*), sl(*)
    TREAL profiles(L_NFIELDS_MAX, nmax, imax/istep, kmax/kstep)
    TREAL txc(imax*jmax*kmax, *)
    TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
    TREAL vmin, vmax
    TINTEGER ij, i, k, n, ifield, nfield, jmin_loc, jmax_loc
    character*32 fname

! ###################################################################
    nfield = L_NFIELDS_MAX
    jmin_loc = max(1, 2*ibuffer_npy)
    jmax_loc = min(jmax, jmax - 2*ibuffer_npy + 1)

! Calculate scalar gradient field sqrt(G_iG_i), put it in array z1
    call FI_GRADIENT(imax, jmax, kmax, z1, a, txc(1, 1))
    do ij = 1, imax*jmax*kmax
        a(ij) = sqrt(a(ij))
    end do

! Calculate boundaries, upper or lower depending on flag isl
    call MINMAX(imax, jmax, kmax, a, vmin, vmax)
    vmin = vmin + C_1EM2_R*(vmax - vmin)
    if (isl == 1) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1, 1), sl, wrk2d)
    else if (isl == 2) then
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1, 1), sl, wrk2d)
    end if

! -------------------------------------------------------------------
! Normal analysis
! -------------------------------------------------------------------
! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, u)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, v)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, w)

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
    write (fname, *) itime; fname = 'slg'//trim(adjustl(fname))

    open (unit=21, file=fname)
    write (21, '(A8,E14.7E3)') 'RTIME = ', rtime
    write (21, *) 'I J N G'

    do k = 1, kmax/kstep
        do i = 1, imax/istep
            do n = 1, nmax
                write (21, 1020) i, k, M_REAL(n - 1 - nmax/2)*(g(1)%nodes(2) - g(1)%nodes(1)), &
                    (profiles(ifield, n, i, k), ifield=1, nfield)
            end do
        end do
    end do
1020 format(I3, 1x, I3, 1x, E10.3E3, L_NFIELDS_MAX(1x, E10.3E3))

    close (21)

    return
end subroutine SL_NORMAL_GRADIENT

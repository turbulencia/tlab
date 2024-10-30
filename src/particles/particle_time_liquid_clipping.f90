#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2015/03 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Sets negative particle liquid to zero
!# Sets particle liquid with no eulerian liquid surrounded to zero
!#
!########################################################################
subroutine PARTICLE_TIME_LIQUID_CLIPPING(s, l_q, l_txc)
    use TLab_Constants, only: wp, wi
    use TLab_Types, only: pointers_dt, pointers3d_dt
    use TLAB_VARS, only: imax, jmax, kmax
    use PARTICLE_VARS, only: isize_part, inb_part_array
    use TLAB_VARS, only: isize_field, inb_scal_array
    use PARTICLE_ARRAYS, only: l_g
    use PARTICLE_INTERPOLATE
    
    implicit none

    real(wp), dimension(isize_field, *), target :: s
    real(wp), dimension(isize_part, *) :: l_q
    real(wp), dimension(isize_part), target :: l_txc

! -------------------------------------------------------------------
    integer(wi) is, i, nvar
    type(pointers3d_dt), dimension(1) :: data
    type(pointers_dt), dimension(1) :: data_out

! ###################################################################
! If negative liquid set lagrange liquid 0
! ###################################################################
    do is = 4, inb_part_array
        do i = 1, l_g%np
            if (l_q(i, is) < 0.0_wp) then
                l_q(i, is) = 0.0_wp
            end if
        end do
    end do

! ###################################################################
! If no liquid around in Eulerian, set liquid droplet to zero
! ###################################################################
    nvar = 0
    nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => s(:, inb_scal_array); data_out(nvar)%field => l_txc(:)
    l_txc = 0.0_wp
    call FIELD_TO_PARTICLE(data(1:nvar), data_out(1:nvar), l_g, l_q)

    do i = 1, l_g%np
        if (l_txc(i) < 0.00001) then
            do is = 4, inb_part_array
                l_q(i, is) = 0.0_wp
            end do
        end if
    end do

    return
end subroutine PARTICLE_TIME_LIQUID_CLIPPING

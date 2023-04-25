! Locate the closest node from below to a list of vertical coordinates
subroutine LOCATE_Y(pmax, y_part, j_part, jmax, y_grid)
    use TLAB_CONSTANTS, only: wp, wi
    implicit none

    integer(wi), intent(in)  :: pmax, jmax
    real(wp),    intent(in)  :: y_part(pmax)
    integer(wi), intent(out) :: j_part(pmax)
    real(wp),    intent(in)  :: y_grid(jmax)

    integer(wi) ip, jm, jp, jc

    do ip = 1, pmax
        jp = jmax
        jm = 1
        jc = (jm + jp)/2
        do while ((y_part(ip) - y_grid(jc))*(y_part(ip) - y_grid(jc + 1)) > 0.0_wp .and. jc > jm)
            if (y_part(ip) < y_grid(jc)) then; jp = jc; 
            else; jm = jc; end if
            jc = (jm + jp)/2
        end do
        j_part(ip) = jc
!     WRITE(*,'(i,3f)') ip, y_grid(jc), y_part(ip), y_grid(jc+1)
    end do

    return
end subroutine LOCATE_Y

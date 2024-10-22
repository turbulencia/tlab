#include "dns_error.h"

!########################################################################
!# Chops string into list of strings
!########################################################################
subroutine LIST_STRING(line, n, a)
    use TLab_Constants, only: wp, wi
    implicit none

    character*(*), intent(IN) :: line
    integer(wi), intent(INOUT) :: n
    character*(*), intent(OUT) :: a(n)

    ! -------------------------------------------------------------------
    integer(wi) i, l1, l2, lmax

    ! ###################################################################
    i = 0                                             ! number of items

    lmax = len_trim(line)
    if (lmax > 0) then
        l1 = lmax - len_trim(adjustl(line(1:lmax))) + 1   ! absolute position of first nonblank in remaining string
        do
            l2 = index(line(l1:lmax), ' ')                 ! relative position of first blank in remaining string

            i = i + 1
            if (l2 == 0) then                           ! we found the last element
                a(i) = line(l1:lmax)
                exit
            else
                a(i) = line(l1:l1 + l2 - 1)
                l1 = lmax - len_trim(adjustl(line(l1 + l2:lmax))) + 1
            end if
            if (i == n) exit

        end do
    end if

    n = i                                             ! return the number of items

    return
end subroutine LIST_STRING

!########################################################################
!# Chops string into list of integers
!########################################################################
subroutine LIST_INTEGER(line, n, a)
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLAB_STOP
    implicit none

    character*(*), intent(IN) :: line
    integer(wi), intent(INOUT) :: n
    integer(wi), dimension(n), intent(OUT) :: a

    ! -------------------------------------------------------------------
    integer(wi) i, lloc
    logical iread
    integer(wi) incr, itmax
    integer(wi) lfirst, ilast
    integer(wi) l1, l2

    ! ###################################################################
    l2 = len_trim(line)
    if (l2 == 0) then ! empty string
        n = 0
        return
    else
        l1 = l2 - len_trim(adjustl(line)) + 1
    end if
    lloc = index(line(l1:l2), ':')

    ! -------------------------------------------------------------------
    ! List separated by commas
    ! -------------------------------------------------------------------
    if (lloc == 0) then
        i = 0
        lfirst = l1 - 1
        do while (.true.)
            lfirst = lfirst + 1
            if ((line(lfirst:lfirst) /= ' ' .and. line(lfirst:lfirst) /= ',')) then
                ! beggining of an item
                do lloc = lfirst, l2
                    iread = .false.
                    if (line(lloc:lloc) == ' ' .or. line(lloc:lloc) == ',') then
                        iread = .true.
                        ilast = lloc - 1
                    end if
                    if (lloc == l2) then
                        iread = .true.
                        ilast = lloc
                    end if
                    if (iread) then
                        i = i + 1
                        ! check the array is big enough
                        if (i > n) then
                            call TLAB_STOP(DNS_ERROR_PARAMETER)
                        end if
                        read (line(lfirst:ilast), *) a(i)
                        lfirst = lloc
                        goto 111
                    end if
                end do
            end if
111         if (lfirst == l2) goto 222
        end do

        ! assign the correct size
222     n = i

        ! -------------------------------------------------------------------
        ! Matlab notation (first:step:last)
        ! -------------------------------------------------------------------
    else
        read (line(l1:lloc - 1), *) a(1)
        l1 = lloc + 1
        lloc = index(line(l1:l2), ':')
        if (lloc == 0) then
            call TLAB_STOP(DNS_ERROR_PARAMETER)
        end if
        lloc = l1 + lloc - 1
        read (line(l1:lloc - 1), *) incr
        l1 = lloc + 1
        read (line(l1:l2), *) itmax

        ! check the array is big enough
        if ((itmax - a(1))/incr + 1 > n) then
            call TLAB_STOP(DNS_ERROR_PARAMETER)
        else
            n = (itmax - a(1))/incr + 1
        end if

        do i = 2, n
            a(i) = a(1) + (i - 1)*incr
        end do

    end if

    return
end subroutine LIST_INTEGER

!########################################################################
!# Chops string into list of real numbers
!########################################################################
subroutine LIST_REAL(line, n, a)
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLAB_STOP
    implicit none

    character*(*), intent(IN) :: line
    integer(wi), intent(INOUT) :: n
    real(wp), dimension(n), intent(OUT) :: a

    ! -------------------------------------------------------------------
    integer(wi) i, lloc
    logical iread
    real(wp) aincr, amax
    integer(wi) lfirst, ilast
    integer(wi) l1, l2

    ! ###################################################################
    l2 = len_trim(line)
    if (l2 == 0) then ! empty string
        n = 0
        return
    else
        l1 = l2 - len_trim(adjustl(line)) + 1
    end if
    lloc = index(line(l1:l2), ':')

    ! -------------------------------------------------------------------
    ! List separated by commas
    ! -------------------------------------------------------------------
    if (lloc == 0) then
        i = 0
        lfirst = l1 - 1
        do while (.true.)
            lfirst = lfirst + 1
            if ((line(lfirst:lfirst) /= ' ' .and. line(lfirst:lfirst) /= ',')) then
                ! beggining of an item
                do lloc = lfirst, l2
                    iread = .false.
                    if (line(lloc:lloc) == ' ' .or. line(lloc:lloc) == ',') then
                        iread = .true.
                        ilast = lloc - 1
                    end if
                    if (lloc == l2) then
                        iread = .true.
                        ilast = lloc
                    end if
                    if (iread) then
                        i = i + 1
                        ! check the array is big enough
                        if (i > n) then
                            call TLAB_STOP(DNS_ERROR_PARAMETER)
                        end if
                        read (line(lfirst:ilast), *) a(i)
                        lfirst = lloc
                        goto 111
                    end if
                end do
            end if
111         if (lfirst == l2) goto 222
        end do

        ! assign the correct size
222     n = i

        ! -------------------------------------------------------------------
        ! Matlab notation (first:step:last)
        ! -------------------------------------------------------------------
    else
        read (line(l1:lloc - 1), *) a(1)
        l1 = lloc + 1
        lloc = index(line(l1:l2), ':')
        if (lloc == 0) then
            call TLAB_STOP(DNS_ERROR_PARAMETER)
        end if
        lloc = l1 + lloc - 1
        read (line(l1:lloc - 1), *) aincr
        l1 = lloc + 1
        read (line(l1:l2), *) amax

        ! check the array is big enough
        if (int((amax - a(1))/aincr) + 1 > n) then
            call TLAB_STOP(DNS_ERROR_PARAMETER)
        else
            n = int((amax - a(1))/aincr) + 1
        end if

        do i = 2, n
            a(i) = a(1) + (i - 1)*aincr
        end do

    end if

    return
end subroutine LIST_REAL

!########################################################################
!########################################################################
subroutine SORT_INTEGER(n, a) ! Sorting elements in array from min to max
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: n
    integer(wi), dimension(n), intent(INOUT) :: a

    ! -------------------------------------------------------------------
    integer(wi) i, j
    integer(wi) dummy

    ! ###################################################################
    do i = 1, n - 1
        do j = i + 1, n
            if (a(j) < a(i)) then
                dummy = a(i)
                a(i) = a(j)
                a(j) = dummy
            end if
        end do
    end do

    return
end subroutine SORT_INTEGER

!########################################################################
!########################################################################
subroutine SORT_REAL(n, a)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: n
    real(wp), dimension(n), intent(INOUT) :: a

    ! -------------------------------------------------------------------
    integer(wi) i, j
    real(wp) dummy

    ! ###################################################################
    do i = 1, n - 1
        do j = i + 1, n
            if (a(j) < a(i)) then
                dummy = a(i)
                a(i) = a(j)
                a(j) = dummy
            end if
        end do
    end do

    return
end subroutine SORT_REAL

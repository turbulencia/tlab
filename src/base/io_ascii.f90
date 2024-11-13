!########################################################################
!#
!# Scan ASCII file for inputs
!# Strings are always forced into lower-case to make the input case-insensitive
!#
!########################################################################

! #######################################################################
! Scan file for an integer value
! #######################################################################
subroutine ScanFile_Int(ofile, ifile, title, name, default, value)
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII
    implicit none

    character*(*), intent(IN) :: ofile, ifile, title, name, default
    integer(wi), intent(OUT) :: value

    character*(128) StrValue

    call TLab_Read_ASCII(ifile, title, name, StrValue, default)
    read (StrValue, *) value

    call TLab_Write_ASCII(ofile, trim(adjustl(name))//'='//trim(adjustl(StrValue)))

    return
end subroutine ScanFile_Int

! #######################################################################
! Scan file for an integer value
! #######################################################################
subroutine ScanFile_LongInt(ofile, ifile, title, name, default, value)
    use TLab_Constants, only: wp, wi, longi
    use TLab_WorkFlow, only: TLab_Write_ASCII
    implicit none

    character*(*), intent(IN) :: ofile, ifile, title, name, default
    integer(longi), intent(OUT) :: value

    character*(128) StrValue

    call TLab_Read_ASCII(ifile, title, name, StrValue, default)
    read (StrValue, *) value

    call TLab_Write_ASCII(ofile, trim(adjustl(name))//'='//trim(adjustl(StrValue)))

    return
end subroutine ScanFile_LongInt

! #######################################################################
! Scan file for an real value
! #######################################################################
subroutine ScanFile_Real(ofile, ifile, title, name, default, value)
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII
    implicit none

    character*(*), intent(IN) :: ofile, ifile, title, name, default
    real(wp), intent(OUT) :: value

    character*(128) StrValue

    call TLab_Read_ASCII(ifile, title, name, StrValue, default)
    read (StrValue, *) value

    call TLab_Write_ASCII(ofile, trim(adjustl(name))//'='//trim(adjustl(StrValue)))

    return
end subroutine ScanFile_Real

! #######################################################################
! Scan file for an char value
! #######################################################################
subroutine ScanFile_Char(ofile, ifile, title, name, default, value)
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII
    implicit none

    character*(*), intent(IN) :: ofile, ifile, title, name, default
    character*(*), intent(OUT) :: value

    call TLab_Read_ASCII(ifile, title, name, value, default)

    call TLab_Write_ASCII(ofile, trim(adjustl(name))//'='//trim(adjustl(value)))

    return
end subroutine ScanFile_Char

! #######################################################################
! Scan file for a string
! #######################################################################
subroutine TLab_Read_ASCII(fname, title, name, value, default)
    use TLab_Constants, only: wp, wi

#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_pro, ims_err
#endif
    implicit none

    character*(*), intent(IN) :: fname, title, name, default
    character*(*), intent(OUT) :: value

! -----------------------------------------------------------------------
    character*512 line
    character*128 str, tag1, tag2, tag3
    integer(wi) equal, code, n

! #######################################################################
! make case-insensitive: title, name and default tags.
    tag1 = title; tag2 = name; tag3 = default
    do n = 1, len(tag1)
        code = iachar(tag1(n:n)); if (code >= 65 .and. code <= 90) tag1(n:n) = achar(code + 32)
    end do
    do n = 1, len(tag2)
        code = iachar(tag2(n:n)); if (code >= 65 .and. code <= 90) tag2(n:n) = achar(code + 32)
    end do
    do n = 1, len(tag3)
        code = iachar(tag3(n:n)); if (code >= 65 .and. code <= 90) tag3(n:n) = achar(code + 32)
    end do

    value = tag3 !default

! -----------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

        open (unit=45, file=fname, status='OLD')

20      continue
        read (45, '(A512)', end=50) line; line = trim(adjustl(line))
        do n = 1, len(line)
            code = iachar(line(n:n)); if (code >= 65 .and. code <= 90) line(n:n) = achar(code + 32)
        end do

        if (trim(adjustl(line)) == '['//trim(adjustl(tag1))//']') then ! Scan within block

30          continue
            read (45, '(A512)', end=50) line; line = trim(adjustl(line))
            if (line(1:1) == '[') goto 20 ! Scape sequence
            if (line(1:1) == '#') goto 30 ! Scape sequence
            equal = index(line, '='); str = trim(adjustl(line(1:equal - 1)))
            do n = 1, len(str)
                code = iachar(str(n:n)); if (code >= 65 .and. code <= 90) str(n:n) = achar(code + 32)
            end do

            if (str == tag2) then
                value = trim(adjustl(line(equal + 1:)))
                do n = 1, len(value)
                    code = iachar(value(n:n)); if (code >= 65 .and. code <= 90) value(n:n) = achar(code + 32)
                end do
                goto 50
            end if
            goto 30

        end if
        goto 20

50      close (unit=45)

#ifdef USE_MPI
    end if
    n = len(value)
    call MPI_BCast(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCast(value, n, MPI_CHAR, 0, MPI_COMM_WORLD, ims_err)
#endif

    return
end subroutine TLab_Read_ASCII

! #######################################################################
! Write ASCII data; complete fields
! #######################################################################
subroutine TLab_Write_ASCII_FIELD(fname, imax, jmax, kmax, u)
    use TLab_Constants, only: wp, wi

#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro, ims_offset_i, ims_offset_k
#endif

    implicit none

    character*(*), intent(IN) :: fname
    integer(wi), intent(IN) :: imax, jmax, kmax
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: u

! -----------------------------------------------------------------------
    integer(wi) idsp, kdsp, i, j, k
    character*32 name_loc

! #######################################################################
#ifdef USE_MPI
    write (name_loc, *) ims_pro; name_loc = trim(adjustl(fname))//'-'//trim(adjustl(name_loc))
    idsp = ims_offset_i; kdsp = ims_offset_k
#else
    name_loc = trim(adjustl(fname))
    idsp = 0; kdsp = 0
#endif

    open (unit=31, file=name_loc, status='unknown')

    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                write (31, *) i + idsp, j, k + kdsp, u(i, j, k)
            end do
        end do
    end do

    close (31)

    return
end subroutine TLab_Write_ASCII_FIELD

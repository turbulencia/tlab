#include "dns_error.h"

#define NMOMS_MAX 10

!########################################################################
!#
!# Calculating the first nm moments of the nv fields defined by the pointers array vars
!#
!########################################################################
subroutine AVG_N_XZ(fname, itime, rtime, nx, ny, nz, nv, nm, vars, igate, gate, y, avg)
    use TLAB_CONSTANTS, only: efile, lfile, wp, wi
    use TLAB_TYPES, only: pointers_dt
    use TLAB_PROCS
    use AVGS, only: AVG1V2d, AVG1V2D1G

    implicit none

    character*(*) fname
    integer(wi) itime
    real(wp) rtime
    integer(wi), intent(IN) :: nx, ny, nz, nv, nm ! Number of moments to consider in the analysis
    type(pointers_dt), intent(IN) :: vars(nv)         ! Array of pointer to the fields to be processed
    integer(1), intent(IN) :: gate(*), igate   ! discrete conditioning criteria
    real(wp), intent(IN) :: y(ny)            ! heights of each plane
    real(wp), intent(OUT) :: avg(ny, nm, nv)

    ! -------------------------------------------------------------------
    integer(wi) j, iv, im
    real(wp) moments(nm)

    character*32 str
    character*32 groupname(1) ! 3 groups for consistency with old TkStat format
    character*250 varname(1)  ! to be reduce to just 1 group in the future

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do j = 1, ny
        do iv = 1, nv
            do im = 1, nm
                if (igate > 0) then
                    moments(im) = AVG1V2D1G(nx, ny, nz, j, igate, im, vars(iv)%field, gate)
                else
                    moments(im) = AVG1V2D(nx, ny, nz, j, im, vars(iv)%field)
                end if
            end do
            if (nm > 1) call RAW_TO_CENTRAL(nm, moments)
            avg(j, 1:nm, iv) = moments(1:nm)
        end do
    end do

    ! ###################################################################
    ! Output
    ! ###################################################################
    varname(:) = ''
    groupname(:) = ''
    do iv = 1, nv
        do im = 1, nm
            varname(1) = trim(adjustl(varname(1)))//' '//trim(adjustl(vars(iv)%tag))
            if (im > 1) then
                write (str, *) im; varname(1) = trim(adjustl(varname(1)))//'.'//trim(adjustl(str))
            end if
        end do
    end do
    ! for consistency with old format I pass 2 groups; no longer (Cedrick Jul 3rd 2021)
    call IO_WRITE_AVERAGES(fname, itime, rtime, ny, nv*nm, 1, y, varname, groupname, avg)

    return
end subroutine AVG_N_XZ

! ###################################################################
subroutine RAW_TO_CENTRAL(nm, moments)
    use TLAB_CONSTANTS, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: nm
    real(wp), intent(INOUT) :: moments(nm)

    integer(wi) im, i, k
    real(wp) aux(NMOMS_MAX), num, den

    ! -------------------------------------------------------------------
    aux(1:nm) = moments(1:nm)

    do im = 2, nm
        moments(im) = 0.0_wp
        do k = 0, im - 1
            num = 1.0_wp
            do i = im, im - k + 1, -1
                num = num*real(i, wp)
            end do
            den = 1.0_wp
            do i = k, 1, -1
                den = den*real(i, wp)
            end do
            moments(im) = moments(im) + num/den*aux(im - k)*((-aux(1))**k)
        end do
        moments(im) = moments(im) + (-aux(1))**im
    end do

    return
end subroutine RAW_TO_CENTRAL

!########################################################################
!#
!# Intermittency factors, i.e., area fraction occupied by each gate level
!#
!########################################################################
subroutine INTER_N_XZ(fname, itime, rtime, nx, ny, nz, np, parname, gate, y, inter)
    use TLAB_CONSTANTS, only: efile, lfile, wp, wi
    use TLAB_PROCS
    use AVGS, only: INTER1V2D

    implicit none

    character*(*) fname, parname(np)
    integer(wi) itime
    real(wp) rtime
    integer(wi), intent(IN) :: nx, ny, nz, np   ! npar is the number of partitions in gate field
    real(wp), intent(IN) :: y(ny)          ! heights of each plane
    integer(1), intent(IN) :: gate(*)        ! field with partitions
    real(wp), intent(OUT) :: inter(ny, np)   ! intermittency factor

    ! -------------------------------------------------------------------
    integer(wi) ip, j
    integer(1) gate_level

    character*32 groupname(1) ! Two groups for consistency with old TkStat format
    character*250 varname(1)  ! to be reduce to just 1 group in the future

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do j = 1, ny
        do ip = 1, np
            gate_level = int(ip, KIND=1)
            inter(j, ip) = INTER1V2D(nx, ny, nz, j, gate_level, gate)
        end do
    end do

    ! ###################################################################
    ! Output
    ! ###################################################################
    varname(:) = ''
    groupname(:) = ''
    do ip = 1, np
        varname(1) = trim(adjustl(varname(1)))//' '//trim(adjustl(parname(ip)))
    end do
    ! for consistency with old format I pass 2 groups
    call IO_WRITE_AVERAGES(fname, itime, rtime, ny, np, 1, y, varname, groupname, inter)

    return
end subroutine INTER_N_XZ

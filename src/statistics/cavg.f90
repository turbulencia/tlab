!########################################################################
!#
!# Calcualte {<a|u_i>, i=1,...,nv] in ny planes.
!# Volume calculation in ny+1, if needed.
!#
!########################################################################
subroutine CAVG1V_N(fname, time, nx, ny, nz, nv, nbins, ibc, umin, umax, u, igate, gate, a, y, avg)
    use TLab_Constants, only: lfile, wp, wi
    use TLab_Pointers, only: pointers_dt
    use TLab_Arrays, only: wrk1d
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use PDFS
#ifdef USE_MPI
    use MPI
#endif

    implicit none

    character*(*), intent(IN) :: fname
    real(wp), intent(IN) :: time
    integer(wi), intent(IN) :: nx, ny, nz, nv, nbins, ibc(nv) ! ibc=0 for external interval, 1 for local
    real(wp), intent(IN) :: umin(nv), umax(nv)            ! Random variables
    type(pointers_dt), intent(IN) :: u(nv)
    integer(1), intent(IN) :: gate(*), igate               ! discrete conditioning criteria
    real(wp), intent(IN) :: a(nx*ny*nz)
    real(wp), intent(IN) :: y(ny)                        ! heights of each plane
    real(wp), intent(OUT) :: avg(nbins + 2, ny + 1, nv)         ! last 2 bins contain the interval bounds

    ! -------------------------------------------------------------------
    integer(wi) j, iv
    character*64 name

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do iv = 1, nv

        do j = 1, ny             ! calculation in planes
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx, ny, nz, j, umin(iv), umax(iv), u(iv)%field, nbins, avg(1, j, iv), wrk1d, a, wrk1d(1, 2))
            else
   call PDF1V2D1G(ibc(iv), nx, ny, nz, j, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, avg(1, j, iv), wrk1d, a, wrk1d(1, 2))
            end if
            avg(1:nbins, j, iv) = wrk1d(1:nbins, 2)
        end do

        if (ny > 1) then   ! calculation in whole volume, saved as plane j=ny+1
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx*ny, 1, nz, 1, umin(iv), umax(iv), u(iv)%field, nbins, avg(1, j, iv), wrk1d, a, wrk1d(1, 2))
            else
 call PDF1V2D1G(ibc(iv), nx*ny, 1, nz, 1, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, avg(1, j, iv), wrk1d, a, wrk1d(1, 2))
            end if
            avg(1:nbins, j, iv) = wrk1d(1:nbins, 2)
        end if

    end do

    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
        do iv = 1, nv
            name = trim(adjustl(fname))
            if (u(iv)%tag /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(u(iv)%tag))
            call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
#include "dns_open_file.h"
            if (ny > 1) then
                write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(avg(:, :, iv))
            else
                write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(avg(:, 1, iv))
            end if
            close (LOC_UNIT_ID)
        end do

#ifdef USE_MPI
    end if
#endif

    return

end subroutine CAVG1V_N

!########################################################################
!########################################################################
subroutine CAVG2V(fname, time, nx, ny, nz, nbins, u, v, a, y, avg)
    use TLab_Constants, only: lfile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Arrays, only: wrk2d
    use PDFS
#ifdef USE_MPI
    use MPI
#endif

    implicit none

    character*(*), intent(IN) :: fname
    real(wp), intent(IN) :: time
    integer(wi), intent(IN) :: nx, ny, nz, nbins(2)
    real(wp), intent(IN) :: u(nx*ny*nz), v(nx*ny*nz), a(nx*ny*nz)
    real(wp), intent(IN) :: y(ny)
    real(wp), intent(OUT) :: avg(nbins(1)*nbins(2) + 2 + 2*nbins(1), ny + 1)

    ! -------------------------------------------------------------------
    integer(wi) j
    character*64 name

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do j = 1, ny             ! calculation in planes
        call PDF2V2D(nx, ny, nz, j, u, v, nbins, avg(1, j), wrk2d, a, wrk2d(1, 2))
        avg(1:nbins(1)*nbins(2), j) = wrk2d(1:nbins(1)*nbins(2), 2)
    end do

    if (ny > 1) then      ! calculation in whole volume, saved as plane j=ny+1
        call PDF2V2D(nx*ny, 1, nz, 1, u, v, nbins, avg(1, j), wrk2d, a, wrk2d(1, 2))
        avg(1:nbins(1)*nbins(2), j) = wrk2d(1:nbins(1)*nbins(2), 2)
    end if

    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
        name = trim(adjustl(fname))
        call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
#include "dns_open_file.h"
        if (ny > 1) then
            write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(avg(:, :))
        else
            write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(avg(:, 1))
        end if
        close (LOC_UNIT_ID)
#ifdef USE_MPI
    end if
#endif

    return

end subroutine CAVG2V

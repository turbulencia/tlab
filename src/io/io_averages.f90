#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Write N profiles in netCDF or ASCII format
!#
!########################################################################
subroutine IO_WRITE_AVERAGES(fname, itime, rtime, ny, nv, ng, y, varnames, groupnames, avg)
    use TLab_Constants, only: lfile, efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII
#ifdef USE_MPI
    use MPI
#endif
#ifdef USE_NETCDF
    use NETCDF
#endif
    implicit none

    character(LEN=*), intent(IN) :: fname, varnames(ng), groupnames(ng)
    integer(wi), intent(IN) :: itime
    real(wp), intent(IN) :: rtime
    integer(wi), intent(IN) :: ny, nv, ng
    real(wp), intent(IN) :: y(ny)
    real(wp), intent(IN) :: avg(ny, nv)

    ! -------------------------------------------------------------------
    integer iv, ig
#ifdef USE_NETCDF
    integer nvg
    integer fid, dtid, dyid, tid, yid, itid, vid(nv)
    character(LEN=32) vname(nv), gname(nv)
#else
    integer j
    character*1300 line
#endif

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Writing '//trim(adjustl(fname))//'...')

#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

        ! -----------------------------------------------------------------------
        ! Using NetCDF format
#ifdef USE_NETCDF
        iv = 1
        do ig = 1, ng
            nvg = nv - iv + 1
            call LIST_STRING(varnames(ig), nvg, vname(iv:nv))
            gname(iv:iv + nvg - 1) = groupnames(ig)
            iv = iv + nvg
            if (iv > nv) exit
        end do

        call NC_CHECK(NF90_CREATE(trim(adjustl(fname))//'.nc', NF90_NETCDF4, fid))

        call NC_CHECK(NF90_DEF_DIM(fid, "t", NF90_UNLIMITED, dtid))
        call NC_CHECK(NF90_DEF_DIM(fid, "y", ny, dyid))

        call NC_CHECK(Nf90_DEF_VAR(fid, "t", NF90_FLOAT, (/dtid/), tid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "y", NF90_FLOAT, (/dyid/), yid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "it", NF90_INT, (/dtid/), itid))
        do iv = 1, nv
            call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(vname(iv))), NF90_FLOAT, (/dyid, dtid/), vid(iv)))
            call NC_CHECK(NF90_PUT_ATT(fid, vid(iv), 'group', gname(iv)))
        end do

        call NC_CHECK(NF90_ENDDEF(fid))

        call NC_CHECK(NF90_PUT_VAR(fid, tid, SNGL(rtime)))
        call NC_CHECK(NF90_PUT_VAR(fid, itid, itime))
        call NC_CHECK(NF90_PUT_VAR(fid, yid, SNGL(y)))
        do iv = 1, nv
            call NC_CHECK(NF90_PUT_VAR(fid, vid(iv), SNGL(avg(1:ny, iv))))
        end do

        call NC_CHECK(NF90_CLOSE(fid))

        ! -----------------------------------------------------------------------
        ! Using ASCII format
#else
#define L_AVGMAX 250

        if (L_AVGMAX < nv) then
            call TLab_Write_ASCII(efile, 'IO_WRITE_AVERAGES. Not enough space in format definition.')
            call TLab_Stop(DNS_ERROR_AVGTMP)
        end if

#define LOC_UNIT_ID 23
#define LOC_STATUS 'unknown'
#ifdef USE_RECLEN
        open (UNIT=LOC_UNIT_ID, RECL=1050, FILE=fname, STATUS=LOC_STATUS) ! this is probably outdated
#else
        open (UNIT=LOC_UNIT_ID, FILE=fname, STATUS=LOC_STATUS)
#endif

        write (LOC_UNIT_ID, '(A8,E14.7E3)') 'RTIME = ', rtime
        line = ''
        do ig = 1, ng
            line = trim(adjustl(line))//' '//trim(adjustl(varnames(ig)))
            write (LOC_UNIT_ID, 1010) 'GROUP = '//trim(adjustl(groupnames(ig)))//' '//trim(adjustl(varnames(ig)))
        end do
        write (LOC_UNIT_ID, 1010) 'I J Y '//trim(adjustl(line)) ! Independent, dependent variables

        do j = 1, ny
            write (LOC_UNIT_ID, 1020) 1, j, y(j), (avg(j, iv), iv=1, nv)
        end do

        close (LOC_UNIT_ID)

1010    format(A)
1020    format(I5, (1x, I5), L_AVGMAX(1x, G_FORMAT_R))

#endif

#ifdef USE_MPI
    end if
#endif

    return
end subroutine IO_WRITE_AVERAGES

! ###################################################################
#ifdef USE_NETCDF
subroutine NC_CHECK(status)
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop 
    use NETCDF

    implicit none

    integer, intent(IN) :: status

    if (status /= nf90_noerr) then
        call TLab_Write_ASCII(efile, 'NETCDF error signal '//trim(adjustl(NF90_STRERROR(status))))
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

    return
end subroutine NC_CHECK
#endif

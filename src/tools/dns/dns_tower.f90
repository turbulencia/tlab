!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE DNS_TOWER - save columns at every iterations to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RELATED INPUT FROM tlab.ini
!
! [SaveTowers]
! Stride=<stride_x>,<stride_y>,<stride_z>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINES
!
! DNS_TOWER_INITIALIZE(x,y,z,stride)
!     - find out which towers to process (according to stride)
!     - allocate corresponding space for arrays
!     - organize buffers where tower data are saved until write
!
! DNS_TOWER_ACCUMULATE(v,index,wrk1d)
!     - PARAMTERS: v        -- data
!                  index    -- what to write (1 - flow, 2 - scalars, 4 - pressure)
!                  wrk1d    -- wrkspace for averaging
!     - accumulates data from flow arrays into tower buffers
!       until towers are written to disk
!
! DNS_TOWER_WRITE(wrk3d)
!     - writes tower data from buffers to disk
!
! DNS_TOWER_FINALIZE()
!     - nothing so far
!
! TOWER_AVG_IK_V(imax, jmax, kmax, a, avg, wrk)
!     - average 1 Variable horizontally over all planes in the vertical.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module DNS_TOWER
    use TLab_Constants, only: wp, wi
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII

    integer(wi) tower_imax, tower_jmax, tower_kmax, tower_isize_field
    integer(wi), target :: tower_isize_plane, tower_isize_plane_total
    integer(wi) tower_imax_total, tower_jmax_total, tower_kmax_total
    integer(wi) tower_offset_i, tower_offset_j, tower_offset_k
    integer(wi) tower_isize_acc_field, tower_isize_acc_mean, tower_isize_acc_write
    integer(wi) tower_accumulation, tower_varcount
    integer(wi) tower_ncid, tower_ncmid, tower_stat
    integer(wi) tower_istride, tower_jstride, tower_kstride
    integer(wi) tower_master
    integer(wi), dimension(:), allocatable :: tower_ipos, tower_kpos, tower_jpos
    integer(wi) :: tower_mode_check

    integer(KIND=8) :: tower_bufsize
    character(LEN=128) :: tower_nc_name

    real(wp), dimension(:), allocatable, target :: tower_buf
    real(wp), dimension(:), pointer :: tower_u, tower_um, tower_v, tower_vm, tower_w, tower_wm
    real(wp), dimension(:), pointer :: tower_p, tower_pm, tower_s, tower_sm, tower_t, tower_it

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine DNS_TOWER_INITIALIZE(stride)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use TLab_Memory, only: imax, jmax, kmax
        use DNS_LOCAL, only: nitera_save

#ifdef USE_MPI
        use MPI
        use FDM, only: g
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_j, ims_offset_k, ims_pro
#endif

        implicit none

#ifdef USE_MPI
#else
        integer(wi) :: ims_offset_i, ims_offset_j, ims_offset_k, ims_pro
        parameter(ims_offset_i=0, ims_offset_j=0, ims_offset_k=0, ims_pro=-1)
#endif

        integer(wi), dimension(3), intent(IN) :: stride

        integer(wi) :: istart, iend, jstart, jend, kstart, kend, ibuf, i, j, k, ii
        integer(wi), pointer :: tip, tip_total
        !
        if (ims_offset_i == 0 .and. ims_offset_k == 0) then
            tower_master = 1
        else
            tower_master = 0
        end if

        tip => tower_isize_plane
        tip_total => tower_isize_plane_total

        tower_isize_acc_write = nitera_save
        tower_varcount = inb_flow + inb_scal + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DISTRIBUTE TOWERS:

        tower_istride = stride(1)
        tower_jstride = stride(2)
        tower_kstride = stride(3)

        istart = ims_offset_i; iend = ims_offset_i + imax - 1
        jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
        kstart = ims_offset_k; kend = ims_offset_k + kmax - 1

        tower_imax = 0; tower_jmax = 0; tower_kmax = 0

        tower_offset_i = ceiling(dble(istart)/tower_istride)
        tower_offset_j = ceiling(dble(jstart)/tower_jstride)
        tower_offset_k = ceiling(dble(kstart)/tower_kstride)

        do i = istart, iend
            if (mod(i - 1, tower_istride) == 0) then
                tower_imax = tower_imax + 1
            end if
        end do

        do j = jstart, jend
            if (mod(j - 1, tower_jstride) == 0) then
                tower_jmax = tower_jmax + 1
            end if
        end do

        do k = kstart, kend
            if (mod(k - 1, tower_kstride) == 0) then
                tower_kmax = tower_kmax + 1
            end if
        end do

#ifdef USE_MPI
        call MPI_Allreduce(tower_imax, tower_imax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, i)
        call MPI_Allreduce(tower_jmax, tower_jmax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, i)
        call MPI_Allreduce(tower_kmax, tower_kmax_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, i)

        tower_imax_total = tower_imax_total/(g(3)%size/kmax*g(2)%size/jmax)
        tower_jmax_total = tower_jmax_total/(g(1)%size/imax*g(3)%size/kmax)
        tower_kmax_total = tower_kmax_total/(g(1)%size/imax*g(2)%size/jmax)
#else
        tower_imax_total = tower_imax
        tower_jmax_total = tower_jmax
        tower_kmax_total = tower_kmax
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALLOCATE SPACE FOR TOWERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tower_isize_plane = tower_imax*tower_kmax
        tower_isize_field = tower_imax*tower_jmax*tower_kmax
        tower_isize_plane_total = tower_imax_total*tower_kmax_total
        tower_isize_acc_field = tower_isize_acc_write*tower_jmax*tip
        tower_isize_acc_mean = tower_isize_acc_write*tower_jmax
        tower_bufsize = 5*tower_isize_acc_field + 5*tower_isize_acc_mean + 2*tower_isize_acc_write

        allocate (tower_ipos(tower_imax), tower_kpos(tower_kmax), tower_jpos(tower_jmax))
        allocate (tower_buf(tower_bufsize))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! REMEMBER WHICH TOWERS TO SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ii = 1
        do i = istart, iend
            if (mod(i, tower_istride) == 0) then
                tower_ipos(ii) = i + 1 - ims_offset_i
                ii = ii + 1
            end if
        end do

        ii = 1
        do j = jstart, jend
            if (mod(j, tower_jstride) == 0) then
                tower_jpos(ii) = j + 1 - ims_offset_j
                ii = ii + 1
            end if
        end do

        ii = 1
        do k = kstart, kend
            if (mod(k, tower_kstride) == 0) then
                tower_kpos(ii) = k + 1 - ims_offset_k
                ii = ii + 1
            end if
        end do

        ibuf = 1; 
        tower_u(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_field; 
        tower_v(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_field; 
        tower_w(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_field; 
        tower_p(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_field; 
        tower_s(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_field; 
        tower_um(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_mean; 
        tower_vm(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_mean; 
        tower_wm(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_mean; 
        tower_pm(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_mean; 
        tower_sm(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_mean; 
        tower_t(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_write; 
        tower_it(1:) => tower_buf(ibuf:); ibuf = ibuf + tower_isize_acc_write; 
        tower_accumulation = 1

    end subroutine DNS_TOWER_INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine DNS_TOWER_ACCUMULATE(v, index, wrk1d)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use MPI
        use TLabMPI_VARS, only: ims_err
#endif

        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Time, only: itime, rtime

        implicit none

        ! 1 -- flow;  2 -- scalar;  4 -- pressure
        integer(wi), intent(IN) :: index
        real(wp), dimension(imax, jmax, kmax, *), intent(IN) :: v
        real(wp), dimension(*), intent(INOUT) :: wrk1d

        integer(wi) :: ii, kk, ip, ipm
        integer(wi), pointer :: tip

        tip => tower_isize_plane
        ip = 1 + (tower_accumulation - 1)*tower_jmax*tip; ipm = ip + tower_jmax - 1

        ! HANDLE PRESSURE
        ! The Counter tower_accumulation is only increased when the pressure is handled
        ! - This needs to be done first!
        ! - The iteration is not increased yet when the pressure is calculated -> use it
        if (index == 4) then
            tower_t(tower_accumulation) = rtime
            tower_it(tower_accumulation) = itime
            do kk = 1, tower_kmax
                do ii = 1, tower_imax
                    tower_p(ip:ipm) = v(tower_ipos(ii), 1:tower_jmax:tower_jstride, tower_kpos(kk), 1)
                    ip = ip + tower_jmax; ipm = ipm + tower_jmax
                end do
            end do
            ip = 1 + (tower_accumulation - 1)*tower_jmax; ipm = ip + tower_jmax - 1

            call TOWER_AVG_IK_V(imax, jmax, kmax, v(1, 1, 1, 1), tower_pm(ip:ipm), wrk1d(6*jmax))
            ! HANDLE FLOW FIELDS
        else if (index == 1) then
            do kk = 1, tower_kmax
                do ii = 1, tower_imax
                    tower_u(ip:ipm) = &
                        v(tower_ipos(ii), 1:tower_jmax:tower_jstride, tower_kpos(kk), 1)
                    tower_v(ip:ipm) = &
                        v(tower_ipos(ii), 1:tower_jmax:tower_jstride, tower_kpos(kk), 2)
                    tower_w(ip:ipm) = &
                        v(tower_ipos(ii), 1:tower_jmax:tower_jstride, tower_kpos(kk), 3)
                    ip = ip + tower_jmax; ipm = ipm + tower_jmax
                end do
            end do
            ip = 1 + (tower_accumulation - 1)*tower_jmax; ipm = ip + tower_jmax - 1
            call TOWER_AVG_IK_V(imax, jmax, kmax, v(1, 1, 1, 1), tower_um(ip:ipm), wrk1d(6*jmax))
            call TOWER_AVG_IK_V(imax, jmax, kmax, v(1, 1, 1, 2), tower_vm(ip:ipm), wrk1d(6*jmax))
            call TOWER_AVG_IK_V(imax, jmax, kmax, v(1, 1, 1, 3), tower_wm(ip:ipm), wrk1d(6*jmax))

            ! HANDLE SCALARS
        else if (index == 2) then
            do kk = 1, tower_kmax
                do ii = 1, tower_imax
                    tower_s(ip:ipm) = &
                        v(tower_ipos(ii), 1:tower_jmax:tower_jstride, tower_kpos(kk), 1)
                    ip = ip + tower_jmax; ipm = ipm + tower_jmax
                end do
            end do
            ip = 1 + (tower_accumulation - 1)*tower_jmax; ipm = ip + tower_jmax - 1
            call TOWER_AVG_IK_V(imax, jmax, kmax, v(1, 1, 1, 1), tower_sm(ip:ipm), wrk1d(6*jmax))

        end if

        tower_mode_check = tower_mode_check + index

        if (tower_mode_check == 7) then
            tower_mode_check = 0
            tower_accumulation = tower_accumulation + 1
        else if (tower_mode_check > 7) then
            write (*, *) 'ERROR - tower_mode_check GREATER THAN 7'
            stop 'INTERNAL ERROR'
        end if

    end subroutine DNS_TOWER_ACCUMULATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine DNS_TOWER_WRITE(wrk3d)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use DNS_LOCAL, only: nitera_save
        use TLab_Time, only: itime
        use TLab_Constants, only: wfile
#ifdef USE_MPI
        use MPI
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_j, ims_offset_k, ims_pro, ims_err
#endif

#ifdef USE_H5
        use HDF5
#endif
        implicit none

#ifdef USE_MPI
#else
        integer(wi) :: ims_offset_i, ims_offset_j, ims_offset_k, ims_pro
        parameter(ims_offset_i=0, ims_offset_j=0, ims_offset_k=0, ims_pro=-1)
#endif

        real(wp), dimension(*), intent(INOUT) :: wrk3d

        integer(wi) :: it, ivar
        character(LEN=64) :: cdummy

        integer(wi) :: itower, ktower, ip_skp, ip_srt, ip_end, tower_count, op_srt, op_end
#ifdef USE_H5
        integer :: h5_err
        integer(HID_T) :: h5_avgFileID, h5_Space1ID, h5_Space2ID, h5_avgDsetID, h5_avgTsetID, h5_avgIsetID
        integer(HID_T) :: h5_dID_loc
        integer(HID_T), dimension(tower_imax, tower_kmax) :: h5_fileID
        integer(HSIZE_T), dimension(2) :: h5_varDim
        character(LEN=64), dimension(5) :: vname
        character(LEN=64), dimension(2) :: vname1d
        character(LEN=128) :: vname_loc
#endif
        integer(wi) :: include_global
        integer(wi), pointer :: tip        
        tip => tower_isize_plane

        if (tip < 1) then
            ! NOTHING TO DO FOR THIS TASK
        else

            if (itime /= int(tower_it(nitera_save) + 1)) then
                call TLab_Write_ASCII(wfile, 'tools/dns/dns_tower.f90 (DNS_TOWER_WRITE)')
                call TLab_Write_ASCII(wfile, 'nitera_save for towers does not match current iteration')
                !                          (But it should if the code is set-up properly)
            end if

#ifdef USE_H5
            h5_varDim = (/tower_jmax, nitera_save/)
            vname = (/'u', 'v', 'w', 'p', 's'/)
            vname1d = (/'time', 'iter'/)

            call h5open_f(h5_err)
            call h5screate_simple_f(2_4, h5_vardim(1:2), h5_Space2ID, h5_err)
            call h5screate_simple_f(1_4, h5_vardim(2:), h5_Space1ID, h5_err)

            write (cdummy, 993) int(tower_it(1)) + 1, itime
993         format('tower.mean', '.', I6.6, '-', I6.6, '.h5')
            call h5fcreate_f(cdummy, H5F_ACC_TRUNC_F, h5_avgFileID, h5_err)
            call h5dcreate_f(h5_avgFileID, vname1D(1), H5T_NATIVE_DOUBLE, h5_Space1ID, h5_dID_loc, h5_err)
            call h5dwrite_f(h5_dID_loc, H5T_NATIVE_DOUBLE, tower_t(1:nitera_save), h5_vardim(2:), h5_err)
            call h5dcreate_f(h5_avgFileID, vname1D(2), H5T_NATIVE_INTEGER, h5_Space1ID, h5_dID_loc, h5_err)
            call h5dwrite_f(h5_dID_loc, H5T_NATIVE_INTEGER, int(tower_it(1:nitera_save)), h5_vardim(2:), h5_err)

            do itower = 1, tower_imax
                do ktower = 1, tower_kmax
                    write (cdummy, 994) &
                        ims_offset_i + tower_ipos(itower), &
                        ims_offset_k + tower_kpos(ktower), &
                        int(tower_it(1)) + 1, itime
994                 format('tower.', I6.6, 'x', I6.6, '.', I6.6, '-', I6.6, '.h5')
                    call h5fcreate_f(cdummy, H5F_ACC_TRUNC_F, h5_fileID(itower, ktower), h5_err)
                end do
            end do
#endif
            do ivar = 1, tower_varcount
                tower_count = 0
#ifdef USE_MPI
                if (ims_pro == 0) then
#endif
#ifdef USE_H5
                    include_global = 0
#else
                    include_global = 2
#endif
                    op_srt = 1; op_end = tower_jmax + include_global; 
                    ip_skp = tower_jmax*1; 
                    ip_srt = 1; 
                    ip_end = ip_srt + tower_jmax - 1; 
                    do it = 1, nitera_save
#ifndef USE_H5
                        wrk3d(op_srt) = tower_t(it)
                        wrk3d(op_srt + 1) = tower_it(it)
#endif
                        select case (ivar)
                        case (1)
                            wrk3d(op_srt + include_global:op_end) = tower_um(ip_srt:ip_end)
                        case (2)
                            wrk3d(op_srt + include_global:op_end) = tower_vm(ip_srt:ip_end)
                        case (3)
                            wrk3d(op_srt + include_global:op_end) = tower_wm(ip_srt:ip_end)
                        case (4)
                            wrk3d(op_srt + include_global:op_end) = tower_pm(ip_srt:ip_end)
                        case (5)
                            wrk3d(op_srt + include_global:op_end) = tower_sm(ip_srt:ip_end)
                        case DEFAULT
                            ! ISSUE WARNING - NO MORE THAN ONE SCALAR
                        end select
                        ip_srt = ip_srt + ip_skp; op_srt = op_srt + tower_jmax + include_global
                        ip_end = ip_end + ip_skp; op_end = op_end + tower_jmax + include_global
                    end do
#ifdef USE_H5
                    ! IMPLICIT Type casting -- saving the variable as real, but passing double values (see below)
                    call h5dcreate_f(h5_avgFileID, vname(ivar), H5T_NATIVE_REAL, h5_space2ID, h5_avgDsetID, h5_err)
                    call h5dwrite_f(h5_avgDsetID, H5T_NATIVE_DOUBLE, wrk3d(1:tower_jmax*nitera_save), h5_vardim, h5_err)
#else
                    write (cdummy, 995) &
                        int(tower_it(1)) + 1, itime, ivar
995                 format('tower.mean', '.', I6.6, '-', I6.6, '.', I1)
                    open (73, FILE=trim(adjustl(cdummy)), ACCESS='STREAM', FORM='UNFORMATTED')
                    write (73, POS=1) wrk3d(1:nitera_save*(tower_jmax + 2))
                    close (73)
#endif
#ifdef USE_MPI
                end if
#endif

                do itower = 1, tower_imax
                    do ktower = 1, tower_kmax
#ifdef USE_H5
                        include_global = 0
#else
                        include_global = 2
#endif
                        ip_skp = tower_jmax*tip; 
                        ip_srt = tower_count*tower_jmax + 1; 
                        ip_end = ip_srt + tower_jmax - 1; 
                        op_srt = 1; op_end = tower_jmax + include_global; 
                        do it = 1, nitera_save
#ifndef USE_H5
                            wrk3d(op_srt) = tower_t(it)
                            wrk3d(op_srt + 1) = tower_it(it)
#endif
                            select case (ivar)
                            case (1)
                                wrk3d(op_srt + include_global:op_end) = tower_u(ip_srt:ip_end)
                            case (2)
                                wrk3d(op_srt + include_global:op_end) = tower_v(ip_srt:ip_end)
                            case (3)
                                wrk3d(op_srt + include_global:op_end) = tower_w(ip_srt:ip_end)
                            case (4)
                                wrk3d(op_srt + include_global:op_end) = tower_p(ip_srt:ip_end)
                            case (5)
                                wrk3d(op_srt + include_global:op_end) = tower_s(ip_srt:ip_end)
                            case DEFAULT
                                ! COULD ISSUE WARNING SOMEWHERE - NO MORE THAN ONE SCALAR FOR TOWERS
                            end select
                            ip_srt = ip_srt + ip_skp; op_srt = op_srt + tower_jmax + include_global
                            ip_end = ip_end + ip_skp; op_end = op_end + tower_jmax + include_global
                        end do
#ifdef USE_H5
                        call h5dcreate_f(h5_fileID(itower, ktower), vname(ivar), H5T_NATIVE_REAL, h5_Space2ID, h5_dID_loc, h5_err)
                        call h5dwrite_f(h5_dID_loc, H5T_NATIVE_DOUBLE, wrk3d(1:tower_jmax*nitera_save), h5_vardim, h5_err)
#else
                        write (cdummy, 997) &
                            ims_offset_i + tower_ipos(itower), &
                            ims_offset_k + tower_kpos(ktower), &
                            int(tower_it(1)) + 1, itime, ivar
997                     format('tower.', I6.6, 'x', I6.6, '.', I6.6, '-', I6.6, '.', I1)
                        open (73, FILE=trim(adjustl(cdummy)), ACCESS='STREAM', FORM='UNFORMATTED')
                        write (73, POS=1) wrk3d(1:nitera_save*(tower_jmax + 2))
                        close (73)
#endif
                        tower_count = tower_count + 1
                    end do
                end do
            end do

#ifdef USE_H5
            call h5fclose_f(h5_avgFileID, h5_err)
            do itower = 1, tower_imax
                do ktower = 1, tower_kmax
                    call h5fclose_f(h5_fileID(itower, ktower), h5_err)
                end do
            end do
            call h5close_f(h5_err)
#endif
        end if
        tower_accumulation = 1

    end subroutine DNS_TOWER_WRITE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine DNS_TOWER_FINALIZE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        implicit none

    end subroutine DNS_TOWER_FINALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine TOWER_AVG_IK_V(imax, jmax, kmax, a, avg, wrk)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_MPI
        use MPI
        use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
#endif

        implicit none

        integer(wi) imax, jmax, kmax
        real(wp) a(imax, jmax, kmax)
        real(wp) avg(tower_jmax), wrk(tower_jmax)
#ifdef USE_MPI
        integer ims_err, len
#endif

        integer(wi) i, j, k

        do j = 1, tower_jmax
            avg(j) = 0.0_wp
        end do

        do k = 1, kmax
            do j = 1, tower_jmax
                do i = 1, imax
                    avg(j) = avg(j) + a(i, tower_jpos(j), k)
                end do
            end do
        end do

        do j = 1, tower_jmax
            avg(j) = avg(j)/real(imax*kmax)
        end do
#ifdef USE_MPI
        call MPI_REDUCE(avg, wrk, tower_jmax, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
        do j = 1, tower_jmax
            avg(j) = wrk(j)/real(ims_npro_i*ims_npro_k)
        end do
#else
#endif

        return
    end subroutine TOWER_AVG_IK_V

end module DNS_TOWER

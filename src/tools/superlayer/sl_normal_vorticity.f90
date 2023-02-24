#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# The calculation of the normal direction is reapeted for different
!# groupd of variables to save memory
!#
!########################################################################
subroutine SL_NORMAL_VORTICITY(isl, ith, iavg, nmax, istep, kstep, nfield, itxc_size, &
                               threshold, ibuffer_npy, u, v, w, p, z1, a, sl, profiles, txc, mean, wrk1d, wrk2d, wrk3d)

    use TLAB_VARS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS
#endif
    use OPR_PARTIAL
    use FI_STRAIN_EQN

    implicit none

#define L_NFIELDS_MAX 13

    TINTEGER isl, ith, nmax, istep, kstep, nfield, itxc_size, iavg, ibuffer_npy
    TREAL threshold
    TREAL u(*), v(*), w(*), p(*), z1(*), a(*), sl(imax, kmax)
    TREAL profiles(L_NFIELDS_MAX, nmax, imax/istep, kmax/kstep)
    TREAL mean(L_NFIELDS_MAX, nmax, 2)
    TREAL txc(imax*jmax*kmax, 7)
    TREAL wrk1d(*), wrk2d(3, imax, kmax), wrk3d(*)

! -------------------------------------------------------------------
    TREAL vmin, vmax, vmean, AVG_IK, diff, normal_factor, dn_u, norm
    TINTEGER ij, i, k, n, ifield, nfield_loc, ipfield, jmin_loc, jmax_loc
    character*32 fname
#ifdef USE_MPI
    TINTEGER ioffset, ip
    integer mpio_ip, mpio_locsize
    integer status(MPI_STATUS_SIZE)
#endif
! ###################################################################
    jmin_loc = max(1, 2*ibuffer_npy)
    jmax_loc = min(jmax, jmax - 2*ibuffer_npy + 1)

    if (nfield < L_NFIELDS_MAX) then
        call TLAB_WRITE_ASCII(efile, 'SL_NORMAL_VORTICITY. Profiles array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    else
        nfield = L_NFIELDS_MAX
    end if
    if (itxc_size < imax*jmax*kmax*7) then
        call TLAB_WRITE_ASCII(efile, 'SL_NORMAL_VORTICITY. Txc array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    end if

! Calculate vorticiy field w_iw_i
    call FI_VORTICITY(imax, jmax, kmax, u, v, w, a, txc(1, 1), txc(1, 2))

! -------------------------------------------------------------------
! Calculate boundaries
! -------------------------------------------------------------------
! threshold w.r.t w_max, therefore threshold^2 w.r.t. w^2_max
    if (ith == 1) then
        call MINMAX(imax, jmax, kmax, a, vmin, vmax)
        vmin = threshold*threshold*vmax
! threshold w.r.t w_mean, therefore threshold^2 w.r.t. w^2_mean
    else if (ith == 2) then
        ij = jmax/2
        vmean = AVG_IK(imax, jmax, kmax, ij, a, dx, dz, area)
        vmin = threshold*threshold*vmean
    end if
! upper or lower depending on flag isl
    if (isl == 1) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, y, a, txc(1, 1), sl, wrk2d)
    else if (isl == 2) then
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, y, a, txc(1, 1), sl, wrk2d)
    end if

! grid size factor for the normal direction
    normal_factor = C_05_R

! ###################################################################
! Normal analysis:
! txc1 ....: vorticity w_i w_i
! txc2 ....: scalar gradient G_i G_i
! txc3 ....: rate-of-strain 2 s_ij s_ij
! ###################################################################
    ipfield = 1
    nfield_loc = 3

    call FI_STRAIN(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 1), txc(1, 2))
    do ij = 1, imax*jmax*kmax
        txc(ij, 3) = C_2_R*txc(ij, 3)
    end do
    call FI_GRADIENT(imax, jmax, kmax, z1, txc(1, 2), txc(1, 1))
    do ij = 1, imax*jmax*kmax
        txc(ij, 1) = a(ij)
    end do

! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, txc(1, 4), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, txc(1, 5), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, txc(1, 6), wrk3d, wrk2d, wrk3d)

    call SL_NORMAL_SAMPLE &
        (imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
         g(1)%scale, g(3)%scale, normal_factor, &
         g(1)%nodes, g(2)%nodes, g(3)%nodes, sl, txc, profiles(ipfield, 1, 1, 1), txc(1, 4), txc(1, 5), txc(1, 6))

! ###################################################################
! Normal analysis:
! txc1 ....: first invariant P
! txc2 ....: second invariant Q
! txc3 ....: third invariant R
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 3

    call FI_INVARIANT_R(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 1), txc(1, 2))
    call FI_INVARIANT_Q(imax, jmax, kmax, u, v, w, txc(1, 2), txc(1, 4), txc(1, 5), txc(1, 6))
    call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 4))

! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, txc(1, 4), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, txc(1, 5), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, txc(1, 6), wrk3d, wrk2d, wrk3d)

    call SL_NORMAL_SAMPLE &
        (imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
         g(1)%scale, g(3)%scale, normal_factor, &
         g(1)%nodes, g(2)%nodes, g(3)%nodes, sl, txc, profiles(ipfield, 1, 1, 1), txc(1, 4), txc(1, 5), txc(1, 6))

! ###################################################################
! Normal analysis:
! txc1 ....: vorticity production
! txc2 ....: vorticity diffusion
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 2

    call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                 txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
    call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, u, v, w, txc(1, 2), &
                                txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
    do ij = 1, imax*jmax*kmax
        txc(ij, 2) = txc(ij, 2)*visc
    end do

! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, txc(1, 4), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, txc(1, 5), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, txc(1, 6), wrk3d, wrk2d, wrk3d)

    call SL_NORMAL_SAMPLE &
        (imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
         g(1)%scale, g(3)%scale, normal_factor, &
         g(1)%nodes, g(2)%nodes, g(3)%nodes, sl, txc, profiles(ipfield, 1, 1, 1), txc(1, 4), txc(1, 5), txc(1, 6))

! ###################################################################
! Normal analysis:
! txc1 ....: scalar gradient production
! txc2 ....: scalar gradient diffusion
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 2

    call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, z1, u, v, w, &
                                txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
    call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, z1, &
                               txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
    diff = visc/schmidt(inb_scal)
    do ij = 1, imax*jmax*kmax
        txc(ij, 2) = txc(ij, 2)*diff
    end do

! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, txc(1, 4), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, txc(1, 5), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, txc(1, 6), wrk3d, wrk2d, wrk3d)

    call SL_NORMAL_SAMPLE &
        (imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
         g(1)%scale, g(3)%scale, normal_factor, &
         g(1)%nodes, g(2)%nodes, g(3)%nodes, sl, txc, profiles(ipfield, 1, 1, 1), txc(1, 4), txc(1, 5), txc(1, 6))

! ###################################################################
! Normal analysis:
! txc1 ....: strain production
! txc2 ....: strain diffusion
! txc3 ....: strain-pressure
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 3

    call FI_STRAIN_PRODUCTION(imax, jmax, kmax, u, v, w, &
                              txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
    do ij = 1, imax*jmax*kmax
        txc(ij, 1) = txc(ij, 1)*C_2_R
    end do
    call FI_STRAIN_DIFFUSION(imax, jmax, kmax, u, v, w, &
                             txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
    do ij = 1, imax*jmax*kmax
        txc(ij, 2) = txc(ij, 2)*visc*C_2_R
    end do
    call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, p, &
                            txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
    do ij = 1, imax*jmax*kmax
        txc(ij, 3) = txc(ij, 3)*C_2_R
    end do

! Calculate gradient of conditioning field; normal stored in u,v,w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, txc(1, 4), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, txc(1, 5), wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, txc(1, 6), wrk3d, wrk2d, wrk3d)

    call SL_NORMAL_SAMPLE &
        (imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
         g(1)%scale, g(3)%scale, normal_factor, &
         g(1)%nodes, g(2)%nodes, g(3)%nodes, sl, txc, profiles(ipfield, 1, 1, 1), txc(1, 4), txc(1, 5), txc(1, 6))

! ###################################################################
! Output averages
! This option is note developed in parallel
! ###################################################################
    if (iavg == 1) then

        do n = 1, nfield*nmax
! mean
            mean(n, 1, 1) = C_0_R
            do k = 1, kmax/kstep
                do i = 1, imax/istep
                    mean(n, 1, 1) = mean(n, 1, 1) + profiles(n, 1, i, k)
                end do
            end do
            mean(n, 1, 1) = mean(n, 1, 1)/M_REAL(imax/istep*kmax/kstep)

! standard deviation
            mean(n, 1, 2) = C_0_R
            do k = 1, kmax/kstep
                do i = 1, imax/istep
                    mean(n, 1, 2) = mean(n, 1, 2) + profiles(n, 1, i, k)*profiles(n, 1, i, k)
                end do
            end do
            mean(n, 1, 2) = mean(n, 1, 2)/M_REAL(imax/istep*kmax/kstep) - mean(n, 1, 1)*mean(n, 1, 1)
            mean(n, 1, 2) = sqrt(mean(n, 1, 2))

        end do

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
! mean value of grid spacing
        dn_u = (x(2) - x(1) + y(2) - y(1) + z(2) - z(1))/C_3_R*normal_factor

        write (fname, *) itime; fname = 'avgSl'//trim(adjustl(fname))

        open (unit=21, file=fname)
        if (ith == 1) then
            write (21, '(A12,E14.7E3)') '# w^2_max = ', vmax
        else if (ith == 2) then
            write (21, '(A12,E14.7E3)') '# w^2_avg = ', vmean
        end if
        write (21, '(A14,E14.7E3)') '# Threshold = ', vmin
        if (isl == 1) then
            write (21, '(A)') '# Upper envelope surface '
        else if (isl == 2) then
            write (21, '(A)') '# Lower envelope surface '
        end if
        write (21, '(A8,E14.7E3)') 'RTIME = ', rtime
        write (21, '(A)') 'GROUP = Mean ' &
            //'rW2 rG2 r2S2 rP rQ rR rP_W rD_W rP_G rD_G r2P_S r2D_S r2SijPij'
        write (21, '(A)') 'GROUP = Sigma ' &
            //'sW2 sG2 s2S2 sP sQ sR sP_W sD_W sP_G sD_G s2P_S s2D_S s2SijPij'
        write (21, '(A)') 'I J N ' &
            //'rW2 rG2 r2S2 rP rQ rR rP_W rD_W rP_G rD_G r2P_S r2D_S r2SijPij ' &
            //'sW2 sG2 s2S2 sP sQ sR sP_W sD_W sP_G sD_G s2P_S s2D_S s2SijPij'

        do n = 1, nmax
            write (21, 1040) i1, i1, &
                M_REAL(n - 1 - nmax/2)*dn_u, &
                (mean(ifield, n, 1), ifield=1, nfield), &
                (mean(ifield, n, 2), ifield=1, nfield)
        end do

! ###################################################################
! Output profiles
! ###################################################################
    else
! sample the normal vector into array wrk2d for output file
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i3, i3, y, sl, txc(1, 4), wrk2d)

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
#ifdef USE_MPI
        call TLAB_MPI_TAGUPDT

        if (ims_pro == 0) then
#endif

! mean value of grid spacing
            dn_u = (x(2) - x(1) + y(2) - y(1) + z(2) - z(1))/C_3_R*normal_factor

            write (str, *) itime; fname = 'slw'//trim(adjustl(str))

            open (unit=21, file=fname)
            if (ith == 1) then
                write (21, '(A12,E14.7E3)') '# w^2_max = ', vmax
            else if (ith == 2) then
                write (21, '(A12,E14.7E3)') '# w^2_avg = ', vmean
            end if
            write (21, '(A14,E14.7E3)') '# Threshold = ', vmin
            if (isl == 1) then
                write (21, '(A)') '# Upper envelope surface '
            else if (isl == 2) then
                write (21, '(A)') '# Lower envelope surface '
            end if
            write (21, '(A8,E14.7E3)') 'RTIME = ', rtime
            write (21, '(A)') 'I J N W2 G2 2S2 P Q R ' &
                //'P_W D_W P_G D_G 2P_S 2D_S 2SijPij Px Py Pz Nx Ny Nz'

            do k = 1, kmax/kstep
                do i = 1, imax/istep
                    norm = sqrt(wrk2d(1, i, k)**2 + wrk2d(2, i, k)**2 + wrk2d(3, i, k)**2)
                    do n = 1, nmax - 1
                        write (21, 1020) istep*i, kstep*k, &
                            M_REAL(n - 1 - nmax/2)*dn_u, &
                            (profiles(ifield, n, i, k), ifield=1, nfield)
                    end do
                    n = nmax
                    write (21, 1030) istep*i, kstep*k, &
                        M_REAL(n - 1 - nmax/2)*dn_u, &
                        (profiles(ifield, n, i, k), ifield=1, nfield), &
                        x(istep*i), sl(i, k), z(kstep*k), &
                        -wrk2d(1, i, k)/norm, -wrk2d(2, i, k)/norm, -wrk2d(3, i, k)/norm
                end do
            end do

#ifdef USE_MPI
            ioffset = 0
            do ip = 2, ims_npro
                mpio_ip = ip - 1
                mpio_locsize = nfield*(imax/istep)*nmax*(kmax/kstep)
                call MPI_RECV(profiles, mpio_locsize, MPI_REAL8, mpio_ip, &
                              ims_tag, MPI_COMM_WORLD, status, ims_err)

                ioffset = ioffset + kmax
                do k = 1, kmax/kstep
                    do i = 1, imax/istep
                        norm = sqrt(wrk2d(1, i, k)**2 + wrk2d(2, i, k)**2 + wrk2d(3, i, k)**2)
                        do n = 1, nmax
                            write (21, 1020) istep*i, kstep*k + ioffset, &
                                M_REAL(n - 1 - nmax/2)*dn_u, &
                                (profiles(ifield, n, i, k), ifield=1, nfield)
                        end do
                        n = nmax
                        write (21, 1030) istep*i, kstep*k + ioffset, &
                            M_REAL(n - 1 - nmax/2)*dn_u, &
                            (profiles(ifield, n, i, k), ifield=1, nfield), &
                            x(istep*i), sl(i, k), z(kstep*k + ioffset), &
                            -wrk2d(1, i, k)/norm, -wrk2d(2, i, k)/norm, -wrk2d(3, i, k)/norm
                    end do
                end do
            end do
#endif

            close (21)

#ifdef USE_MPI
        else
            mpio_locsize = nfield*(imax/istep)*nmax*(kmax/kstep)
            call MPI_SEND &
                (profiles, mpio_locsize, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
        end if
#endif

    end if

1020 format(I3, 1x, I3, 1x, E10.3E3, L_NFIELDS_MAX(1x, E10.3E3))
1030 format(I3, 1x, I3, 1x, E10.3E3, L_NFIELDS_MAX(1x, E10.3E3), 6(1x, E10.3E3))
1040 format(I3, 1x, I3, 1x, E10.3E3, L_NFIELDS_MAX(1x, E10.3E3), L_NFIELDS_MAX(1x, E10.3E3))

    return
end subroutine SL_NORMAL_VORTICITY

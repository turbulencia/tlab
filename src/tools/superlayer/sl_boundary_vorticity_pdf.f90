#include "types.h"
#include "dns_error.h"

subroutine SL_BOUNDARY_VORTICITY_PDF(isl, ith, np, nfield, itxc_size, threshold, ibuffer_npy, &
                                     u, v, w, z1, a, sl, samples, pdf, txc, wrk1d, wrk2d, wrk3d)

    use TLAB_VARS
    use TLAB_AVGS
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    
    implicit none

#define L_NFIELDS_MAX 5

    TREAL threshold
    TINTEGER isl, ith, nfield, itxc_size, np, ibuffer_npy
    TREAL u(*), v(*), w(*), z1(*), a(*), sl(imax*kmax, 2)
    TREAL samples(L_NFIELDS_MAX*imax*kmax*2), pdf(np, L_NFIELDS_MAX)
    TREAL txc(imax*jmax*kmax, 6)
    TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
    TREAL vmin, vmax, vmean
    TINTEGER ij, ikmax, ipfield, nfield_loc, ioffset, isize, iv, ip, jmin_loc, jmax_loc
    character*32 fname
    character*32 varname(L_NFIELDS_MAX)

! ###################################################################
    jmin_loc = max(1, 2*ibuffer_npy)
    jmax_loc = min(jmax, jmax - 2*ibuffer_npy + 1)

    if (nfield < L_NFIELDS_MAX) then
        call TLAB_WRITE_ASCII(efile, 'SL_VORTICITY_PDF. Samples array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    else
        nfield = L_NFIELDS_MAX
    end if
    if (itxc_size < imax*jmax*kmax*6) then
        call TLAB_WRITE_ASCII(efile, 'SL_VORTICITY_PDF. Txc array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    end if

! Offset to be used in SL_BOUNDARY_SAMPLE
    if (isl == 1 .or. isl == 2) then
        ioffset = nfield
    else
        ioffset = 2*nfield
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
        vmean = AVG_IK(imax, jmax, kmax, ij, a, g(1)%jac, g(3)%jac, area)
        vmin = threshold*threshold*vmean
    end if
! upper/lower/both depending on flag isl
    if (isl == 1) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1, 1), sl, wrk2d)
    else if (isl == 2) then
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1, 1), sl, wrk2d)
    else if (isl == 3) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1, 1), sl(1, 1), wrk2d)
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1, 1), sl(1, 2), wrk2d)
    end if

! ###################################################################
! Sample along the surface
! txc1 ....: log vorticity W_iW_i
! txc2 ....: log gradient G_i Gi
! txc3 ....: log strain 2 s_ij s_ij
! ###################################################################
    ipfield = 1
    nfield_loc = 3

    do ij = 1, imax*jmax*kmax
        txc(ij, 1) = log(a(ij))
    end do
    varname(1) = 'log(W2)'

    call FI_GRADIENT(imax, jmax, kmax, z1, txc(1, 2), txc(1, 3))
    do ij = 1, imax*jmax*kmax
        txc(ij, 2) = log(txc(ij, 2))
    end do
    varname(2) = 'log(G2)'

    call FI_STRAIN(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5))
    do ij = 1, imax*jmax*kmax
        txc(ij, 3) = log(C_2_R*txc(ij, 3))
    end do
    varname(3) = 'log(2S2)'

    if (isl == 1 .or. isl == 2) then
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                g(2)%nodes, sl, txc, samples(ipfield))
    else
! txc? in upper and lower layer consecutive in samples array
        do iv = 1, nfield_loc
            call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                    g(2)%nodes, sl(1, 1), txc(1, iv), samples(ipfield + 2*(iv - 1)))
            call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                    g(2)%nodes, sl(1, 2), txc(1, iv), samples(ipfield + 2*(iv - 1) + 1))
        end do
! correct the number of fields from this section
        nfield_loc = nfield_loc*2
    end if

! ###################################################################
! Sample along the surface
! txc1 ....: cos(grad w, grad G)
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 1

    call FI_GRADIENT(imax, jmax, kmax, z1, txc(1, 2), txc(1, 3))

    call FI_ISOSURFACE_ANGLE(imax, jmax, kmax, a, txc(1, 2), txc(1, 1), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

    varname(4) = 'cos(gradG,gradW)'

    if (isl == 1 .or. isl == 2) then
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                g(2)%nodes, sl, txc, samples(ipfield))
    else
! txc? in upper and lower layer consecutive in samples array
        do iv = 1, nfield_loc
            call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                    g(2)%nodes, sl(1, 1), txc(1, iv), samples(ipfield + 2*(iv - 1)))
            call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
                                    g(2)%nodes, sl(1, 2), txc(1, iv), samples(ipfield + 2*(iv - 1) + 1))
        end do
! correct the number of fields from this section
        nfield_loc = nfield_loc*2
    end if

! ###################################################################
! Vertical distance to the centerplane
! Data already in arrays sl
! ###################################################################
    ipfield = ipfield + nfield_loc
    nfield_loc = 1

    varname(5) = 'height'

    if (isl == 1) then
        do ij = 1, imax*kmax
            ip = ipfield + (ij - 1)*ioffset
            samples(ip) = sl(ij, 1) - qbg(1)%ymean
        end do
    else if (isl == 2) then
        do ij = 1, imax*kmax
            ip = ipfield + (ij - 1)*ioffset
            samples(ip) = qbg(1)%ymean - sl(ij, 1)
        end do
    else
        do ij = 1, imax*kmax
            ip = ipfield + (ij - 1)*ioffset
            samples(ip) = sl(ij, 1) - qbg(1)%ymean
            samples(ip + 1) = qbg(1)%ymean - sl(ij, 2)
        end do
! correct the number of fields from this section
        nfield_loc = nfield_loc*2
    end if

! ###################################################################
! Calculate PDFs
! ###################################################################
! transpose data from sample into txc space to have nfield as last index
! ioffset is equal to the number of fields
    ikmax = imax*kmax
    call DNS_TRANSPOSE(samples, ioffset, ikmax, ioffset, txc, ikmax)

    if (isl == 1 .or. isl == 2) then
        isize = 1
    else
        isize = 2
    end if
    write (fname, *) itime; fname = 'pdfSl'//trim(adjustl(fname))
    ! CALL PDF3D_N(fname, rtime, i1, rtime, &
    !      imax, isize, kmax, nfield, np, txc, pdf, wrk1d)
    ! Need to adapt to the new formulation of PDF1V_N
    ! CALL PDF1V_N(fname, rtime, imax*isize, 1, kmax, &
    !    nfield, opt_bins(1), ibc, amin,amax,data, gate_level,gate, y_aux, pdf)

    return
end subroutine SL_BOUNDARY_VORTICITY_PDF

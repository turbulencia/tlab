#include "types.h"
#include "dns_error.h"

subroutine SL_BOUNDARY_VORTICITY_JPDF(iopt, isl, ith, np, nfield, itxc_size, &
                                      threshold, ibuffer_npy, u, v, w, sl, samples, txc, wrk1d, wrk2d, wrk3d)

    use TLAB_VARS
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_VORTICITY_EQN

    implicit none

#define L_NFIELDS_MAX 4

    TREAL threshold
    TINTEGER iopt, isl, ith, nfield, itxc_size, np, ibuffer_npy
    TREAL u(*), v(*), w(*), sl(imax*kmax, *)
    TREAL samples(L_NFIELDS_MAX*imax*kmax)
    TREAL txc(imax*jmax*kmax, 6)
    TREAL wrk1d(*), wrk2d(imax*kmax, *), wrk3d(*)

! -------------------------------------------------------------------
    TREAL vmin, vmax, vmean, AVG_IK
    TINTEGER ij, ikmax, nfield_loc, isize, jmin_loc, jmax_loc
    integer(1) igate
    character*32 fname
    character*16 suffix

! ###################################################################
    jmin_loc = max(1, 2*ibuffer_npy)
    jmax_loc = min(jmax, jmax - 2*ibuffer_npy + 1)

    if (nfield < L_NFIELDS_MAX) then
        call TLAB_WRITE_ASCII(efile, 'SL_VORTICITY_JPDF. Samples array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    else
        nfield = L_NFIELDS_MAX
    end if
    if (itxc_size < imax*jmax*kmax*6) then
        call TLAB_WRITE_ASCII(efile, 'SL_VORTICITY_JPDF. Txc array size.')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    end if

! ###################################################################
! Calculate fields
! ###################################################################
! -------------------------------------------------------------------
! RQ PDF
! txc1 ....: third invariant R
! txc2 ....: second invariant Q
! -------------------------------------------------------------------
    if (iopt == 3) then
        call TLAB_WRITE_ASCII(lfile, 'Computing invariant R...')
        call FI_INVARIANT_R(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
        call TLAB_WRITE_ASCII(lfile, 'Computing invariant Q...')
        call FI_INVARIANT_Q(imax, jmax, kmax, u, v, w, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        suffix = 'RQ '

! -------------------------------------------------------------------
! WS PDF !! this one does not make sense !!
! txc1 ....: vorticity w_i w_i
! txc2 ....: strain 2 s_ij s_ij
! -------------------------------------------------------------------
    else if (iopt == 4) then
        call TLAB_WRITE_ASCII(lfile, 'Computing vorticity...')
        call FI_VORTICITY(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3))
        call TLAB_WRITE_ASCII(lfile, 'Computing strain...')
        call FI_STRAIN(imax, jmax, kmax, u, v, w, txc(1, 2), txc(1, 3), txc(1, 4))
        do ij = 1, imax*jmax*kmax
            txc(ij, 2) = C_2_R*txc(ij, 2)
        end do
        suffix = 'WS '

    end if

! ###################################################################
! Calculate vorticiy w_iw_i as conditioning field and boundaries
! Array txc3, and sl
! ###################################################################
    call FI_VORTICITY(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5))

! -------------------------------------------------------------------
! Calculate boundaries
! -------------------------------------------------------------------
! threshold w.r.t w_max, therefore threshold^2 w.r.t. w^2_max
    if (ith == 1) then
        call MINMAX(imax, jmax, kmax, txc(1, 3), vmin, vmax)
        vmin = threshold*threshold*vmax
! threshold w.r.t w_mean, therefore threshold^2 w.r.t. w^2_mean
    else if (ith == 2) then
        ij = jmax/2
        vmean = AVG_IK(imax, jmax, kmax, ij, txc(1, 3), g(1)%jac, g(3)%jac, area)
        vmin = threshold*threshold*vmean
    end if
! upper/lower/both depending on flag isl
    if (isl == 1) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, txc(1, 3), txc(1, 4), sl, wrk2d)
    else if (isl == 2) then
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, txc(1, 3), txc(1, 4), sl, wrk2d)
    else if (isl == 3) then
        call SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, txc(1, 3), txc(1, 4), sl(1, 1), wrk2d)
        call SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, txc(1, 3), txc(1, 4), sl(1, 2), wrk2d)
    end if

! ###################################################################
! Sample along the surface in sl
! ###################################################################
    if (isl == 1 .or. isl == 2) then
        nfield_loc = 2
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i2, nfield_loc, g(2)%nodes, sl, txc, samples)

    else
        nfield_loc = 4
! txc1 in upper and lower layer consecutive in samples array
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i1, nfield_loc, g(2)%nodes, sl(1, 1), txc(1, 1), samples(1))
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i1, nfield_loc, g(2)%nodes, sl(1, 2), txc(1, 1), samples(2))
! txc2 in upper and lower layer consecutive in samples array
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i1, nfield_loc, g(2)%nodes, sl(1, 1), txc(1, 2), samples(3))
        call SL_BOUNDARY_SAMPLE(imax, jmax, kmax, i1, nfield_loc, g(2)%nodes, sl(1, 2), txc(1, 2), samples(4))

    end if

! ###################################################################
! Calculate JPDF
! ###################################################################
! make ifields the last variable, putting first the imax*kmax
    ikmax = imax*kmax
    call DNS_TRANSPOSE(samples, nfield_loc, ikmax, nfield_loc, wrk2d, ikmax)

    isize = nfield_loc/2
    write (fname, *) itime; fname = 'jpdf'//trim(adjustl(suffix))//trim(adjustl(fname))
    igate = 0
    ! CALL JPDF3D(fname, i0, igate, i0, imax, isize, kmax, i0, i0,&
    !      txc(1,3), wrk2d(1,1+isize), wrk2d(1,1), np, np, wrk2d(1,5), wrk2d(1,6), wrk2d(1,7), wrk1d)
    ! Check, need to pass gate to th new formulation PDF2V of joint pdfs
    ! We pass ny=1 and it only calculates 3D pdfs (twice, but it allows us to reuse existing routines)
    call PDF2V(fname, rtime, imax*isize, 1, kmax, opt_bins, y_aux, wrk2d(1, 1 + isize), wrk2d(1, 1), pdf)

    return
end subroutine SL_BOUNDARY_VORTICITY_JPDF

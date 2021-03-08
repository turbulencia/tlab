#include "types.h"
#include "dns_error.h"

SUBROUTINE SL_BOUNDARY_VORTICITY_JPDF(iopt, isl, ith, np, nfield, itxc_size, &
     threshold, ibuffer_npy, u,v,w, sl, samples, txc, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

#define L_NFIELDS_MAX 4

  TREAL threshold
  TINTEGER iopt, isl, ith, nfield, itxc_size, np, ibuffer_npy
  TREAL u(*), v(*), w(*), sl(imax*kmax,*)
  TREAL samples(L_NFIELDS_MAX*imax*kmax)
  TREAL txc(imax*jmax*kmax,6)
  TREAL wrk1d(*), wrk2d(imax*kmax,*), wrk3d(*)

! -------------------------------------------------------------------
  TREAL vmin, vmax, vmean, AVG_IK
  TINTEGER ij, ikmax, nfield_loc, isize, jmin_loc, jmax_loc
  INTEGER(1) igate
  CHARACTER*32 fname
  CHARACTER*16 suffix
#ifdef USE_MPI
  TINTEGER ioffset, ip
  INTEGER mpio_ip, mpio_locsize
  INTEGER status(MPI_STATUS_SIZE)
#endif

! ###################################################################
  jmin_loc = MAX(1,2*ibuffer_npy)
  jmax_loc = MIN(jmax,jmax - 2*ibuffer_npy +1)

  IF ( nfield .LT. L_NFIELDS_MAX ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_VORTICITY_JPDF. Samples array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ELSE
     nfield = L_NFIELDS_MAX
  ENDIF
  IF ( itxc_size .LT. imax*jmax*kmax*6 ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_VORTICITY_JPDF. Txc array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! ###################################################################
! Calculate fields
! ###################################################################
! -------------------------------------------------------------------
! RQ PDF
! txc1 ....: third invariant R
! txc2 ....: second invariant Q
! -------------------------------------------------------------------
  IF ( iopt .EQ. 3 ) THEN
     CALL IO_WRITE_ASCII(lfile,'Computing invariant R...')
     CALL FI_INVARIANT_R(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
     CALL IO_WRITE_ASCII(lfile,'Computing invariant Q...')
     CALL FI_INVARIANT_Q(imax,jmax,kmax, u,v,w, txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
     suffix = 'RQ '

! -------------------------------------------------------------------
! WS PDF !! this one does not make sense !!
! txc1 ....: vorticity w_i w_i
! txc2 ....: strain 2 s_ij s_ij
! -------------------------------------------------------------------
  ELSE IF ( iopt .EQ. 4 ) THEN
     CALL IO_WRITE_ASCII(lfile,'Computing vorticity...')
     CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
     CALL IO_WRITE_ASCII(lfile,'Computing strain...')
     CALL FI_STRAIN(imax,jmax,kmax, u,v,w, txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
     DO ij = 1,imax*jmax*kmax
        txc(ij,2) = C_2_R*txc(ij,2)
     ENDDO
     suffix='WS '

  ENDIF

! ###################################################################
! Calculate vorticiy w_iw_i as conditioning field and boundaries
! Array txc3, and sl
! ###################################################################
  CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4),txc(1,5), wrk2d,wrk3d)

! -------------------------------------------------------------------
! Calculate boundaries
! -------------------------------------------------------------------
! threshold w.r.t w_max, therefore threshold^2 w.r.t. w^2_max
  IF ( ith .EQ. 1 ) THEN
     CALL MINMAX(imax,jmax,kmax, txc(1,3), vmin,vmax)
     vmin = threshold*threshold*vmax
! threshold w.r.t w_mean, therefore threshold^2 w.r.t. w^2_mean
  ELSE IF ( ith .EQ. 2 ) THEN
     ij = jmax/2
     vmean = AVG_IK(imax, jmax, kmax, ij, txc(1,3), g(1)%jac,g(3)%jac, area)
     vmin = threshold*threshold*vmean
  ENDIF
! upper/lower/both depending on flag isl
  IF ( isl .EQ. 1 ) THEN
     CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, g(2)%nodes, txc(1,3), txc(1,4), sl,      wrk2d)
  ELSE IF ( isl .EQ. 2 ) THEN
     CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, g(2)%nodes, txc(1,3), txc(1,4), sl,      wrk2d)
  ELSE IF ( isl .EQ. 3 ) THEN
     CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, g(2)%nodes, txc(1,3), txc(1,4), sl(1,1), wrk2d)
     CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, g(2)%nodes, txc(1,3), txc(1,4), sl(1,2), wrk2d)
  ENDIF

! ###################################################################
! Sample along the surface in sl
! ###################################################################
  IF ( isl .EQ. 1 .OR. isl .EQ. 2 ) THEN
     nfield_loc = 2
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i2, nfield_loc, g(2)%nodes, sl, txc, samples)

  ELSE
     nfield_loc = 4
! txc1 in upper and lower layer consecutive in samples array
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i1, nfield_loc, g(2)%nodes, sl(1,1), txc(1,1), samples(1))
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i1, nfield_loc, g(2)%nodes, sl(1,2), txc(1,1), samples(2))
! txc2 in upper and lower layer consecutive in samples array
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i1, nfield_loc, g(2)%nodes, sl(1,1), txc(1,2), samples(3))
     CALL SL_BOUNDARY_SAMPLE(imax,jmax,kmax, i1, nfield_loc, g(2)%nodes, sl(1,2), txc(1,2), samples(4))

  ENDIF

! ###################################################################
! Calculate JPDF
! ###################################################################
! make ifields the last variable, putting first the imax*kmax
  ikmax = imax*kmax
  CALL DNS_TRANSPOSE(samples, nfield_loc, ikmax, nfield_loc, wrk2d, ikmax)

  isize = nfield_loc/2
  WRITE(fname,*) itime; fname = 'jpdf'//TRIM(ADJUSTL(suffix))//TRIM(ADJUSTL(fname))
  igate = 0
  ! CALL JPDF3D(fname, i0, igate, i0, imax, isize, kmax, i0, i0,&
  !      txc(1,3), wrk2d(1,1+isize), wrk2d(1,1), np, np, wrk2d(1,5), wrk2d(1,6), wrk2d(1,7), wrk1d)
  ! Check, need to pass gate to th new formulation JPDF2D of joint pdfs
  ! We pass ny=1 and it only calculates 3D pdfs (twice, but it allows us to reuse existing routines)
  CALL JPDF2D(fname, imax*isize, 1, kmax, opt_bins, y_aux, wrk2d(1,1+isize), wrk2d(1,1), pdf, wrk2d )

  RETURN
END SUBROUTINE SL_BOUNDARY_VORTICITY_JPDF

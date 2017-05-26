#include "types.h"
#include "dns_error.h"

SUBROUTINE SL_BOUNDARY_VORTICITY_PDF(isl, ith, np, nfield, itxc_size, threshold, ibuffer_npy, &
     u, v, w, z1, a, sl, samples, pdf, txc, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

#define L_NFIELDS_MAX 5

  TREAL threshold
  TINTEGER isl, ith, nfield, itxc_size, np, ibuffer_npy
  TREAL u(*), v(*), w(*), z1(*), a(*), sl(imax*kmax,2)
  TREAL samples(L_NFIELDS_MAX*imax*kmax*2), pdf(np,L_NFIELDS_MAX)
  TREAL txc(imax*jmax*kmax,6)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TREAL vmin, vmax, vmean, AVG_IK, ycenter
  TINTEGER ij, ikmax, ipfield, nfield_loc, ioffset, isize, iv, ip, jmin_loc, jmax_loc
  CHARACTER*32 fname
  CHARACTER*32 varname(L_NFIELDS_MAX)

! ###################################################################
  jmin_loc = MAX(1,2*ibuffer_npy)
  jmax_loc = MIN(jmax,jmax - 2*ibuffer_npy +1)

  IF ( nfield .LT. L_NFIELDS_MAX ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_VORTICITY_PDF. Samples array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ELSE
     nfield = L_NFIELDS_MAX
  ENDIF
  IF ( itxc_size .LT. imax*jmax*kmax*6 ) THEN
     CALL IO_WRITE_ASCII(efile, 'SL_VORTICITY_PDF. Txc array size.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! Offset to be used in SL_BOUNDARY_SAMPLE
  IF ( isl .EQ. 1 .OR. isl .EQ. 2 ) THEN
     ioffset = nfield
  ELSE
     ioffset = 2*nfield
  ENDIF

! Calculate vorticiy field w_iw_i
  CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, a, txc(1,1),txc(1,2), wrk2d,wrk3d)

! -------------------------------------------------------------------
! Calculate boundaries
! -------------------------------------------------------------------
! threshold w.r.t w_max, therefore threshold^2 w.r.t. w^2_max
  IF ( ith .EQ. 1 ) THEN
     CALL MINMAX(imax,jmax,kmax, a, vmin,vmax)
     vmin = threshold*threshold*vmax
! threshold w.r.t w_mean, therefore threshold^2 w.r.t. w^2_mean
  ELSE IF ( ith .EQ. 2 ) THEN
     ij = jmax/2
     vmean = AVG_IK(imax, jmax, kmax, ij, a, g(1)%jac,g(3)%jac, area)
     vmin = threshold*threshold*vmean
  ENDIF
! upper/lower/both depending on flag isl
  IF ( isl .EQ. 1 ) THEN
     CALL SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1,1), sl, wrk2d)
  ELSE IF ( isl .EQ. 2 ) THEN
     CALL SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1,1), sl, wrk2d)
  ELSE IF ( isl .EQ. 3 ) THEN
     CALL SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1,1), sl(1,1), wrk2d)
     CALL SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1,1), sl(1,2), wrk2d)
  ENDIF

! ###################################################################
! Sample along the surface
! txc1 ....: log vorticity W_iW_i
! txc2 ....: log gradient G_i Gi
! txc3 ....: log strain 2 s_ij s_ij
! ###################################################################
  ipfield    = 1
  nfield_loc = 3

  DO ij = 1,imax*jmax*kmax
     txc(ij,1) = log(a(ij))
  ENDDO
  varname(1) = 'log(W2)'

  CALL FI_GRADIENT(imax,jmax,kmax, z1,txc(1,2), txc(1,3), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,2) = log(txc(ij,2))
  ENDDO
  varname(2) = 'log(G2)'

  CALL FI_STRAIN(imax,jmax,kmax, u,v,w, txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     txc(ij,3) = log(C_2_R*txc(ij,3))
  ENDDO
  varname(3) = 'log(2S2)'

  IF ( isl .EQ. 1 .OR. isl .EQ. 2 ) THEN
     CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
          g(2)%nodes, sl, txc, samples(ipfield))
  ELSE
! txc? in upper and lower layer consecutive in samples array
     DO iv = 1,nfield_loc
        CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
             g(2)%nodes, sl(1,1), txc(1,iv), samples(ipfield+2*(iv-1)  ))
        CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
             g(2)%nodes, sl(1,2), txc(1,iv), samples(ipfield+2*(iv-1)+1))
     ENDDO
! correct the number of fields from this section
     nfield_loc = nfield_loc*2
  ENDIF

! ###################################################################
! Sample along the surface
! txc1 ....: cos(grad w, grad G)
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 1

  CALL FI_GRADIENT(imax,jmax,kmax, z1,txc(1,2), txc(1,3), wrk2d,wrk3d)

  CALL FI_ISOSURFACE_ANGLE(imax,jmax,kmax, a,txc(1,2), txc(1,1), &
       txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

  varname(4) = 'cos(gradG,gradW)'

  IF ( isl .EQ. 1 .OR. isl .EQ. 2 ) THEN
     CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
          g(2)%nodes, sl, txc, samples(ipfield))
  ELSE
! txc? in upper and lower layer consecutive in samples array
     DO iv = 1,nfield_loc
        CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
             g(2)%nodes, sl(1,1), txc(1,iv), samples(ipfield+2*(iv-1)  ))
        CALL SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, ioffset, &
             g(2)%nodes, sl(1,2), txc(1,iv), samples(ipfield+2*(iv-1)+1))
     ENDDO
! correct the number of fields from this section
     nfield_loc = nfield_loc*2
  ENDIF

! ###################################################################
! Vertical distance to the centerplane
! Data already in arrays sl
! ###################################################################
  ipfield    = ipfield + nfield_loc
  nfield_loc = 1

  varname(5) = 'height'

  ycenter = y(1) + qbg(1)%ymean *g(2)%scale
  IF ( isl .EQ. 1 ) THEN
     DO ij = 1,imax*kmax
        ip = ipfield + (ij-1)*ioffset
        samples(ip) = sl(ij,1) - ycenter
     ENDDO
  ELSE IF (isl .EQ. 2 ) THEN
     DO ij = 1,imax*kmax
        ip = ipfield + (ij-1)*ioffset
        samples(ip) = ycenter - sl(ij,1)
     ENDDO
  ELSE
     DO ij = 1,imax*kmax
        ip = ipfield + (ij-1)*ioffset
        samples(ip  ) = sl(ij,1) - ycenter
        samples(ip+1) = ycenter - sl(ij,2)
     ENDDO
! correct the number of fields from this section
     nfield_loc = nfield_loc*2
  ENDIF

! ###################################################################
! Calculate PDFs
! ###################################################################
! transpose data from sample into txc space to have nfield as last index
! ioffset is equal to the number of fields
  ikmax = imax*kmax
  CALL DNS_TRANSPOSE(samples, ioffset, ikmax, ioffset, txc, ikmax)

  IF ( isl .EQ. 1 .OR. isl .EQ. 2 ) THEN
     isize = 1
  ELSE
     isize = 2
  ENDIF
  WRITE(fname,*) itime; fname='pdfSl'//TRIM(ADJUSTL(fname))
  CALL PDF3D_N(fname, varname, i1, i1, rtime, &
       imax, isize, kmax, nfield, np, txc, pdf, wrk1d)

  RETURN
END SUBROUTINE SL_BOUNDARY_VORTICITY_PDF

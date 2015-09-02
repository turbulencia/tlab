#include "types.h"

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/07/11 - J.P. Mellado
!#              Cleaned. Written in terms of the conditional one.
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!#
!########################################################################
!# DESCRIPTION
!#
!# The joint PDF between U and V is computed. 
!# The calculation is done using the relation P(V,U)=P(V|U)P(U); the
!# field V is discretized between min and max and for a fixed U, then
!# PDF of V is computed.
!# Note that V goes into the horizontal axis X, U in the vertical axis Y
!#
!########################################################################
!# ARGUMENTS 
!#
!# igate     In     Gate level. If 0, no intermittency considered
!# iubc/ivbc In     BCs for variable u/v: 0 None
!#                                        1 Drop minimum
!#                                        2 Drop maximum
!#                                        3 Drop both
!########################################################################
SUBROUTINE JPDF3D(fname, inorm, igate, ianalyze, imax, jmax, kmax, &
     iubc, ivbc, gate, u, v, nxp, nyp, x, y, pdf, wrk1d)

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*32 fname
  TINTEGER inorm, ianalyze
  TINTEGER imax, jmax, kmax
  TINTEGER iubc, ivbc
  TREAL v(*)
  TREAL u(*)
  TINTEGER nxp, nyp
  TREAL x(nyp,nxp)
  TREAL y(nyp,nxp)
  TREAL pdf(nyp,nxp,2)
  TREAL wrk1d(nyp,*)

  INTEGER(1) gate(*), igate

! -------------------------------------------------------------------
  TINTEGER nyp_total, ixp, iyp, np_total, nplim, nxplot, ij
  TREAL ucenter, udelta, umin, umax, xmin, xmax, vmin, vmax
  TREAL ylim, c_ratio, c_ratio_min

  ! CHARACTER*32 str1, str2, str3, str4, str5, varname(1)
  ! TINTEGER l1, l2, iv, ll1, ll2

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  TREAL umin_p, umax_p
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

! ###################################################################
! calculate min/max of u
  IF ( igate .GT. 0 ) THEN
     umin = C_BIG_R
     umax =-C_BIG_R
     DO ij = 1,imax*jmax*kmax
        IF ( gate(ij) .EQ. igate ) THEN
           umin = MIN(umin, u(ij))
           umax = MAX(umax, u(ij))
        ENDIF
     ENDDO

#ifdef USE_MPI
     CALL MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
     umin = umin_p
     CALL MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     umax = umax_p
#endif
     
  ELSE
     CALL MINMAX(imax,jmax,kmax, u, umin,umax)

  ENDIF
  udelta = (umax-umin)/M_REAL(nxp)

! BCs in variable U: eliminate the outer bins.
  IF ( iubc .EQ. 1 .OR. iubc .EQ. 3 ) THEN
     umin = umin + udelta
  ENDIF
  IF ( iubc .EQ. 2 .OR. iubc .EQ. 3 ) THEN
     umax = umax - udelta
  ENDIF
  udelta = (umax-umin)/M_REAL(nxp)

! ###################################################################
! Joint PDf calculation based on conditional PDF
! ###################################################################
  ylim = C_1EM3_R
  c_ratio_min = C_05_R
!  c_ratio_min = C_1EM2_R
  np_total = 0 
  nxplot = 0

! loop on U variable (vertical axis)
  DO ixp = 1, nxp
     ucenter = umin + udelta*C_05_R + M_REAL(ixp-1)*udelta

     IF ( igate .GT. 0 ) THEN
        CALL CPDF2V3D1G(inorm, i1, imax, jmax, kmax, igate, vmin, vmax, gate, u, v,&
             ucenter, udelta, nyp, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), nyp_total)
        xmin = wrk1d(1,  1)
        xmax = wrk1d(nyp,1)
     ELSE
        CALL CPDF2V3D(inorm, i1, imax, jmax, kmax, vmin, vmax, u, v,&
             ucenter, udelta, nyp, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), nyp_total)
        xmin = wrk1d(1,  1)
        xmax = wrk1d(nyp,1)
     ENDIF

     IF ( ianalyze .EQ. 1 ) THEN
        CALL PDF_ANALIZE(nyp, ivbc, xmin, xmax, wrk1d(1,2), ylim, vmin, vmax, nplim)
! depending on filling ratio, neglect pdf
        c_ratio = M_REAL(nplim)/M_REAL(nyp)
        IF ( c_ratio .GT. c_ratio_min ) THEN
           IF ( igate .GT. 0 ) THEN
              CALL CPDF2V3D1G(inorm, i2, imax, jmax, kmax, igate, vmin, vmax, gate, u, v,&
                   ucenter, udelta, nyp, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), nyp_total)
           ELSE
              CALL CPDF2V3D(inorm, i2, imax, jmax, kmax, vmin, vmax, u, v,&
                   ucenter, udelta, nyp, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), nyp_total)
           ENDIF
        ELSE
           DO iyp = 1,nyp
              wrk1d(iyp,2) = C_0_R
           ENDDO
        ENDIF
        xmin = wrk1d(1,  1)
        xmax = wrk1d(nyp,1)
     ENDIF

! save data for OpenDx file; V in horizontal axis X, U in vertical axis Y
! if, for the given ubin, there is no point, xmax-xmin is negative, as set
! in CPDF2V3D1G. In that case, do not write data for OpenDx
     IF ( xmax-xmin .GT. C_SMALL_R ) THEN
        nxplot = nxplot + 1
        DO iyp = 1,nyp
           x(iyp,nxplot) = wrk1d(iyp,1)
           y(iyp,nxplot) = ucenter
           pdf(iyp,nxplot,2) = wrk1d(iyp,2)
           pdf(iyp,nxplot,1) = wrk1d(iyp,2)*nyp_total
        ENDDO
        np_total = np_total + nyp_total
     ENDIF
  ENDDO

  IF ( inorm  .EQ. 1 ) THEN
     DO ixp = 1, nxplot
        DO iyp = 1,nyp
           pdf(iyp,ixp,1) = pdf(iyp,ixp,1)/(M_REAL(np_total)*udelta)
        ENDDO
     ENDDO
  ELSE
     DO ixp = 1, nxplot
        DO iyp = 1,nyp
           pdf(iyp,ixp,1) = pdf(iyp,ixp,2)
        ENDDO
     ENDDO
  ENDIF

! ###################################################################
! OpenDx file
! Using OpenDx ASCII format, since these PDF files are expected to
! be relatively small
! ###################################################################
! #ifdef USE_MPI
!   IF ( ims_pro .EQ. 0 ) THEN
! #endif

!      OPEN(i55,TRIM(ADJUSTL(fname))//'.dx',STATUS='unknown')

! ! comment section
!      IF ( inorm .EQ. 0 ) THEN
!         WRITE(i55,'(A)') '# Histogram (no normalization)'
!      ELSE
!         WRITE(i55,'(A)') '# PDF (normalization s.t. integral is 1)'
!      ENDIF
!      WRITE(i55,'(A33,I8)') '# Number of points in the sample ', np_total
!      IF ( igate .EQ. 0 ) THEN
!         WRITE(i55,'(A)') '# No conditioning on intermittency'
!      ELSE
!         WRITE(i55,'(A39,I8)') '# Conditioning on intermittency, level ', igate
!      ENDIF
!      IF ( ianalyze .EQ. 0 ) THEN
!         WRITE(i55,'(A)') '# No PDF analysis'
!      ELSE
!         WRITE(i55,'(A)') '# PDF analysis'
!      ENDIF

! ! preparing for general case with nvar
!      iv = 1
!      varname(1) = 'pdf'

! ! positions object
!      str1='positions '
!      CALL CONCATI(str1,iv)
!      CALL CLIP(str1,l1,l2)
!      WRITE(str5,'(I10)') nyp*nxplot
!      CALL CLIP(str5,ll1,ll2)
!      WRITE(i55,'(A)') 'object "'//str1(l1:l2)//'" '&
!           //'class array type float rank 1 shape 2 items '//str5(ll1:ll2)&
!           //' data follows'
!      DO ixp = 1,nxplot
!         DO iyp = 1,nyp
!            WRITE(i55,'(2(1X,E10.3E3))') SNGL(x(iyp,ixp)), SNGL(y(iyp,ixp))
!         ENDDO
!      ENDDO
!      WRITE(i55,'(A)') ' '

! ! connections object
!      str2='connections '
!      CALL CONCATI(str2,iv)
!      CALL CLIP(str2,l1,l2)
!      WRITE(i55,'(A,I10,I10)') 'object "'//str2(l1:l2)//'" '&
!           //'class gridconnections counts ',nxplot,nyp
!      WRITE(i55,'(A)') ' '

! ! joint pdf data object
!      str3='dataJoint '
!      CALL CONCATI(str3,iv)
!      CALL CLIP(str3,l1,l2)
!      WRITE(i55,'(A)') 'object "'//str3(l1:l2)//'" '&
!           //'class array type float rank 0 items '//str5(ll1:ll2)&
!           //' data follows'
!      DO ixp = 1,nxplot
!         DO iyp = 1,nyp
!            WRITE(i55,'(1X,E10.3E3)') SNGL(pdf(iyp,ixp,1))
!         ENDDO
!      ENDDO
!      WRITE(i55,'(A)') 'attribute "dep" string "positions"'
!      WRITE(i55,'(A)') ' '

! ! joint pdf field object
!      CALL CLIP(varname(iv),l1,l2)
!      WRITE(i55,'(A)') &
!           'object "'//varname(iv)(l1:l2)//'" class field'
!      CALL CLIP(str1,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "positions" value "'//str1(l1:l2)//'"'
!      CALL CLIP(str2,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "connections" value "'//str2(l1:l2)//'"'
!      CALL CLIP(str3,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "data" value "'//str3(l1:l2)//'"'
!      WRITE(i55,'(A)') ' '

! ! conditional pdf data object
!      str4='dataCond '
!      CALL CONCATI(str4,iv)
!      CALL CLIP(str4,l1,l2)
!      WRITE(i55,'(A)') 'object "'//str4(l1:l2)//'" '&
!           //'class array type float rank 0 items '//str5(ll1:ll2)&
!           //' data follows'
!      DO ixp = 1,nxplot
!         DO iyp = 1,nyp
!            WRITE(i55,'(1X,E10.3E3)') SNGL(pdf(iyp,ixp,2))
!         ENDDO
!      ENDDO
!      WRITE(i55,'(A)') 'attribute "dep" string "positions"'
!      WRITE(i55,'(A)') ' '

! ! conditional pdf field object
!      CALL CLIP(varname(iv),l1,l2)
!      WRITE(i55,'(A)') &
!           'object "'//varname(iv)(l1:l2)//'_c" class field'
!      CALL CLIP(str1,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "positions" value "'//str1(l1:l2)//'"'
!      CALL CLIP(str2,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "connections" value "'//str2(l1:l2)//'"'
!      CALL CLIP(str4,l1,l2)
!      WRITE(i55,'(A)') &
!           ' component "data" value "'//str4(l1:l2)//'"'
!      WRITE(i55,'(A)') ' '

! ! group object
!      WRITE(i55,'(A)') 'object "group" class group'

!      iv = 1

!      CALL CLIP(varname(iv),l1,l2)
!      WRITE(i55,'(A)') ' member "'//varname(iv)(l1:l2)//'" value "'//varname(iv)(l1:l2)//'"'
!      WRITE(i55,'(A)') &
!           ' member "'//varname(iv)(l1:l2)//'_c" value "'//varname(iv)(l1:l2)//'_c"'

!      CLOSE(i55)

! #ifdef USE_MPI
!   ENDIF
! #endif

  RETURN
END SUBROUTINE JPDF3D

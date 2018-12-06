#include "types.h"

!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/07/11 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Recalculating max/min for a given threshold. 
!# Adding nplim, the count of points above threshold
!# BCs flag ibc to drop extreme points.
!#
!########################################################################
SUBROUTINE PDF_ANALIZE(nbins, ibc, pdf, ylim, umin, umax, nplim)

  IMPLICIT NONE

  TINTEGER nbins, nplim, ibc
  TREAL umin, umax, ylim
  TREAL pdf(nbins+2) ! Space at the end for the min and max values in the sample variable

! -------------------------------------------------------------------
  TINTEGER i, ip, npmin, npmax
  TREAL pmax, pdfstep, yl

! ###################################################################
  pdfstep = (pdf(nbins+2)-pdf(nbins+1))/M_REAL(nbins-1)

! eliminate the outer bins according to BCs
  npmin = 1
  npmax = nbins
  IF ( ibc .EQ. 1 .OR. ibc .EQ. 3 ) THEN
     npmin = npmin + 1
  ENDIF
  IF ( ibc .EQ. 2 .OR. ibc .EQ. 3 ) THEN
     npmax = npmax - 1
  ENDIF
  umin = pdf(nbins+1) - C_05_R*pdfstep + pdfstep*M_REAL(npmin-1)
  umax = pdf(nbins+1) - C_05_R*pdfstep + pdfstep*M_REAL(npmax)

! get maximum
  pmax = pdf(npmin)
  DO i = npmin+1, npmax
     pmax = MAX(pmax, pdf(i))
  ENDDO

! -------------------------------------------------------------------
  yl = ylim*pmax

  nplim = 0
  DO i = npmin,npmax
     IF ( pdf(i) .GT. yl ) THEN
        nplim = nplim + 1
     ENDIF
  ENDDO

! -------------------------------------------------------------------
  ip = 0
  DO i = npmin,npmax
     IF ( pdf(i) .GT. yl ) THEN
        ip = i
        GOTO 123
     ENDIF
  ENDDO

123 IF ( ip .GT. 0 ) THEN
     umin = pdf(nbins+1) - C_05_R*pdfstep + pdfstep*M_REAL(ip-1)
  ENDIF

! -------------------------------------------------------------------
  ip = 0
  DO i = npmax,npmin,-1
     IF ( pdf(i) .GT. yl ) THEN
        ip = i
        GOTO 124
     ENDIF
  ENDDO

124 IF ( ip .GT. 0 ) THEN
     umax = pdf(nbins+1) - C_05_R*pdfstep + pdfstep*M_REAL(ip)
  ENDIF

  RETURN
END SUBROUTINE PDF_ANALIZE

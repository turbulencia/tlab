#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Initialize data of reference thermodynamic profiles
!#
!########################################################################
SUBROUTINE FI_BACKGROUND_INITIALIZE(wrk1d)

  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_VARS, ONLY : inb_scal, inb_scal_array, imax,jmax,kmax, imode_eqns
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : pbg, sbg, damkohler,froude,schmidt
  USE TLAB_VARS, ONLY : rbackground, ribackground, bbackground, pbackground, tbackground, epbackground
  USE TLAB_VARS, ONLY : buoyancy
  USE TLAB_PROCS
  USE THERMO_VARS, ONLY : imixture, GRATIO
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(g(2)%size,*), INTENT(INOUT) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER is, j, ip, nlines, offset
  TREAL ycenter, PROFILES

! #######################################################################
  ALLOCATE(pbackground(g(2)%size))
  ALLOCATE(rbackground(g(2)%size))
  ALLOCATE(ribackground(g(2)%size))
  ALLOCATE(bbackground(g(2)%size))
  ALLOCATE(tbackground(g(2)%size))
  ALLOCATE(epbackground(g(2)%size))

! #######################################################################
! Thermodynamic background profiles
! #######################################################################
  rbackground = C_1_R ! defaults
  ribackground= C_1_R
  pbackground = C_1_R
  tbackground = C_1_R
  epbackground= C_0_R

! Construct given thermodynamic profiles
  DO is = 1,inb_scal
     ycenter = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean
     DO j = 1,g(2)%size
        wrk1d(j,is) = PROFILES(sbg(is), ycenter, g(2)%nodes(j))
     ENDDO
!     wrk1d(:,is) = sbg(is)%reference
  ENDDO

  IF ( pbg%parameters(5) .GT. C_0_R ) THEN
! Calculate derived thermodynamic profiles
     epbackground = (g(2)%nodes - g(2)%nodes(1) - g(2)%scale *pbg%ymean) *GRATIO /pbg%parameters(5)

     IF ( buoyancy%active(2) ) THEN
!        CALL FI_HYDROSTATIC_H_OLD(g(2)%size, g(2)%nodes, wrk1d, epbackground, tbackground, pbackground, wrk1d(1,4))
        CALL FI_HYDROSTATIC_H(g(2), wrk1d, epbackground, tbackground, pbackground, wrk1d(1,inb_scal_array+1))
     ENDIF

  ENDIF

  IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Calculate q_l
     CALL THERMO_AIRWATER_PH(i1,g(2)%size,i1, wrk1d(1,2),wrk1d(1,1), epbackground,pbackground )
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
     CALL THERMO_AIRWATER_LINEAR(i1,g(2)%size,i1, wrk1d, wrk1d(1,inb_scal_array))
  ENDIF

  IF ( pbg%parameters(5) .GT. C_0_R ) THEN
     CALL THERMO_ANELASTIC_DENSITY(i1,g(2)%size,i1, wrk1d, epbackground,pbackground, rbackground)
     ribackground = C_1_R /rbackground
  ENDIF

! Calculate buoyancy profile
  IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
     CALL THERMO_ANELASTIC_BUOYANCY(i1,g(2)%size,i1, wrk1d, epbackground,pbackground,rbackground, bbackground)
  ELSE
     wrk1d(:,inb_scal_array+1) = C_0_R
     CALL FI_BUOYANCY(buoyancy, i1,g(2)%size,i1, wrk1d, bbackground, wrk1d(1,inb_scal_array+1))
  ENDIF

! #######################################################################
! Add diagnostic fields to reference profile data, if any
! #######################################################################
  DO is = inb_scal+1,inb_scal_array ! Add diagnostic fields, if any
     sbg(is)%mean = C_1_R; sbg(is)%delta = C_0_R; sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)
  ENDDO
! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
  sbg(is)%mean  =    (bbackground(1)+bbackground(g(2)%size)) /froude
  sbg(is)%delta = ABS(bbackground(1)-bbackground(g(2)%size)) /froude
  sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     is = is + 1
     CALL THERMO_ANELASTIC_THETA_L(i1,g(2)%size,i1, wrk1d, epbackground,pbackground, wrk1d(1,inb_scal_array+1))
     sbg(is)%mean  =    (wrk1d(1,inb_scal_array+1)+wrk1d(g(2)%size,inb_scal_array+1)) * C_05_R
     sbg(is)%delta = ABS(wrk1d(1,inb_scal_array+1)-wrk1d(g(2)%size,inb_scal_array+1))
     sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)
  ENDIF

! #######################################################################
! Anelastic density correction term in burgers operator
! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL TLAB_WRITE_ASCII(lfile,'Initialize anelastic density correction in burgers operator.')

! Density correction term in the burgers operator along X
     g(1)%anelastic = .TRUE.
#ifdef USE_MPI
     IF ( ims_npro_i .GT. 1 ) THEN
        nlines = ims_size_i(TLAB_MPI_I_PARTIAL)
        offset = nlines *ims_pro_i
     ELSE
#endif
        nlines = jmax*kmax
        offset = 0
#ifdef USE_MPI
     ENDIF
#endif
     ALLOCATE(g(1)%rhoinv(nlines))
     DO j = 1,nlines
        ip = MOD(offset +j -1, g(2)%size) +1
        g(1)%rhoinv(j) = ribackground(ip)
     ENDDO

! Density correction term in the burgers operator along Y; see fdm_initialize
! implemented directly in the tridiagonal system
     ip = 0
     DO is = 0,inb_scal ! case 0 for the velocity
        g(2)%lu2d(:,           ip+2) = g(2)%lu2d(:,ip+2) *ribackground(:)  ! matrix U; 1/diagonal
        g(2)%lu2d(:g(2)%size-1,ip+3) = g(2)%lu2d(:,ip+3) *rbackground (2:) ! matrix U; 1. superdiagonal
        ip = ip + 3
     ENDDO

! Density correction term in the burgers operator along Z
     g(3)%anelastic = .TRUE.
#ifdef USE_MPI
     IF ( ims_npro_k .GT. 1 ) THEN
        nlines = ims_size_k(TLAB_MPI_K_PARTIAL)
        offset = nlines *ims_pro_k
     ELSE
#endif
        nlines = imax*jmax
        offset = 0
#ifdef USE_MPI
     ENDIF
#endif
     ALLOCATE(g(3)%rhoinv(nlines))
     DO j = 1,nlines
        ip = (offset +j -1) /imax +1
        g(3)%rhoinv(j) = ribackground(ip)
     ENDDO

  END IF

  RETURN
END SUBROUTINE FI_BACKGROUND_INITIALIZE

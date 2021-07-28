#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

! ###################################################################
! ###################################################################
SUBROUTINE OPR_FILTER_INITIALIZE(g, f, wrk1d)

  USE DNS_TYPES, ONLY : grid_dt, filter_dt

  IMPLICIT NONE

  TYPE(grid_dt),       INTENT(IN)    :: g
  TYPE(filter_dt),     INTENT(INOUT) :: f
  TREAL, DIMENSION(*), INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------

! ###################################################################
  IF ( f%inb_filter .GT. 0 ) &
       ALLOCATE( f%coeffs( f%size, f%inb_filter ) )

  SELECT CASE( f%type )

  CASE( DNS_FILTER_4E, DNS_FILTER_ADM )
     CALL FLT_E4_INI(g%scale, g%nodes, f)

  CASE( DNS_FILTER_TOPHAT )
     CALL FLT_T1_INI(g%scale, g%nodes, f, wrk1d)

  CASE( DNS_FILTER_COMPACT )
     CALL FLT_C4_INI(g%jac, f)

  CASE( DNS_FILTER_HELMHOLTZ )
     f%parameters(2) =-C_1_R /( f%parameters(1) )**2

  END SELECT

  RETURN
END SUBROUTINE OPR_FILTER_INITIALIZE

! ###################################################################
! Filter kernel along one direction
! ###################################################################
SUBROUTINE OPR_FILTER_1D(nlines, f, u,result, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : filter_dt

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nlines     ! # of lines to be solved
  TYPE(filter_dt),                 INTENT(IN)    :: f
  TREAL, DIMENSION(nlines,f%size), INTENT(IN)    :: u          ! field to be filtered
  TREAL, DIMENSION(nlines,f%size), INTENT(OUT)   :: result     ! filtered filed
  TREAL, DIMENSION(f%size,5),      INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER delta

! ###################################################################
  delta = INT(f%parameters(1))

  SELECT CASE( f%type )

  CASE( DNS_FILTER_COMPACT )
     CALL FLT_C4(f%size,nlines, f%periodic, f%bcsmin,f%bcsmax, f%coeffs, u,result, wrk1d)
     IF ( f%periodic ) THEN
        CALL TRIDPFS(f%size,        wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(f%size,nlines, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), result, wrk2d)
     ELSE
        CALL TRIDFS(f%size,        wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(f%size,nlines, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), result)
     ENDIF

  CASE( DNS_FILTER_6E      )
     CALL FLT_E6 (f%size,nlines, f%periodic, f%bcsmin, f%bcsmax,           u,result)

  CASE( DNS_FILTER_4E      )
     CALL FLT_E4 (f%size,nlines, f%periodic,                     f%coeffs, u,result)

  CASE( DNS_FILTER_ADM     )
     CALL FLT_ADM(f%size,nlines, f%periodic,                     f%coeffs, u,result, wrk3d)

  CASE( DNS_FILTER_TOPHAT )
     IF ( f%periodic ) THEN
        IF ( f%uniform ) THEN
           IF      ( delta .EQ. 2 ) THEN; CALL FLT_T1PD2(f%size,nlines,        u,result)
           ELSE IF ( delta .EQ. 4 ) THEN; CALL FLT_T1PD4(f%size,nlines,        u,result)
           ELSE;                          CALL FLT_T1P  (f%size,nlines, delta, u,result)
           ENDIF
        ELSE
           CALL FLT_T1P_ND(f%size,nlines, delta, f%coeffs, u,result)
        ENDIF
     ELSE
        IF ( f%uniform ) THEN
           CALL FLT_T1(f%size,nlines, delta, f%coeffs, u,result)
        ELSE
           IF      ( delta .EQ. 2 ) THEN; CALL FLT_T1NDD2(f%size,nlines,        f%coeffs, u,result)
           ELSE IF ( delta .EQ. 4 ) THEN; CALL FLT_T1NDD4(f%size,nlines,        f%coeffs, u,result)
           ELSE IF ( delta .EQ. 6 ) THEN; CALL FLT_T1NDD6(f%size,nlines,        f%coeffs, u,result)
           ELSE;                          CALL FLT_T1ND  (f%size,nlines, delta, f%coeffs, u,result)
           ENDIF
        ENDIF
     ENDIF

  END SELECT

  RETURN
END SUBROUTINE OPR_FILTER_1D

! ###################################################################
! Filter in Ox direction
! ###################################################################
SUBROUTINE OPR_FILTER_X(nx,ny,nz, f, u, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : filter_dt
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)            :: nx,ny,nz
  TYPE(filter_dt),            INTENT(IN)            :: f
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(ny*nz),    INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(f%size,*), INTENT(INOUT)         :: wrk1d
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type

! -------------------------------------------------------------------
  TINTEGER nyz

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
  id = f%mpitype ! DNS_MPI_I_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, wrk3d, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => wrk3d
     p_b => u
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => u
     p_b => wrk3d
     nyz = ny*nz
#ifdef USE_MPI
  ENDIF
#endif

! -------------------------------------------------------------------
! Make  x  direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, f%size, f%size, nyz,        p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, f%size, nyz,        f%size, p_b, nyz)
#endif

! ###################################################################
  CALL OPR_FILTER_1D(nyz, f, p_b,p_a, wrk1d,wrk2d,tmp)

! ###################################################################
! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, nyz, nyz,        f%size, p_b, f%size)
#else
  CALL DNS_TRANSPOSE(p_a, nyz, f%size, nyz,        p_b, f%size)
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_b, p_a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_FILTER_X

! ###################################################################
! Filter in Oy direction
! ###################################################################
SUBROUTINE OPR_FILTER_Y(nx,ny,nz, f, u, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : filter_dt

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)            :: nx,ny,nz
  TYPE(filter_dt),            INTENT(IN)            :: f
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(*),        INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(f%size,*), INTENT(INOUT)         :: wrk1d
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type

! -----------------------------------------------------------------------
  TINTEGER nxy, nxz

  TREAL, DIMENSION(:), POINTER :: p_org, p_dst

! #######################################################################
  nxy = nx*ny
  nxz = nx*nz

! -------------------------------------------------------------------
! Make y-direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org => u
     p_dst => wrk3d
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (u, nxy, nxy, nz, wrk3d, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, wrk3d, nz)
#endif
     p_org => wrk3d
     p_dst => u
  ENDIF

! ###################################################################
  CALL OPR_FILTER_1D(nxz, f, p_org,p_dst, wrk1d,wrk2d,tmp)

! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst, nz, nz, nxy, p_org, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, p_org, nxy)
#endif
  ENDIF

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_FILTER_Y

! ###################################################################
! Filter in Oz direction
! ###################################################################
SUBROUTINE OPR_FILTER_Z(nx,ny,nz, f, u, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : filter_dt
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)            :: nx,ny,nz
  TYPE(filter_dt),            INTENT(IN)            :: f
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(nx*ny),    INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(f%size,*), INTENT(INOUT)         :: wrk1d
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type

! -------------------------------------------------------------------
  TINTEGER nxy

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
  id = f%mpitype ! DNS_MPI_K_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(u, wrk3d, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => wrk3d
     p_b => u
     nxy = ims_size_k(id)
  ELSE
#endif
     p_a => u
     p_b => wrk3d
     nxy = nx*ny
#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
  CALL OPR_FILTER_1D(nxy, f, p_a,p_b, wrk1d,wrk2d,tmp)

! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_b, p_a, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_FILTER_Z

!########################################################################
! Spectral band filter
!########################################################################
SUBROUTINE OPR_FILTER_BAND_2D(nx,ny,nz, spc_param, a)

  USE TLAB_VARS, ONLY : isize_txc_dimz
  USE TLAB_VARS, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(*)                                  :: spc_param
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(INOUT) :: a

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip
  TREAL fi,fk,f

! #######################################################################
  DO k = 1,nz
#ifdef USE_MPI
     kglobal = k + ims_offset_k
#else
     kglobal = k
#endif
     IF ( kglobal .LE. g(3)%size/2+1 ) THEN; fk = M_REAL(kglobal-1)/g(3)%scale
     ELSE;                                   fk =-M_REAL(g(3)%size+1-kglobal)/g(3)%scale; ENDIF

     DO i = 1,nx/2+1
#ifdef USE_MPI
        iglobal = i + ims_offset_i/2
#else
        iglobal = i
#endif
        fi = M_REAL(iglobal-1)/g(1)%scale

        f = SQRT(fi**2 + fk**2)

! apply spectral cutoff
        DO j = 1,ny
           ip = (j-1)*(nx/2+1) + i
           IF ( (f-spc_param(1))*(spc_param(2)-f) .LT. C_0_R ) a(ip,k) = C_0_R
        ENDDO

     ENDDO
  ENDDO

  RETURN
END SUBROUTINE OPR_FILTER_BAND_2D

!########################################################################
!# DESCRIPTION
!#
!# Spectral filter with smooth (error-function) transition.
!# The error function filter-response function is imposed
!# in logarithmic wavenumber space.
!#
!########################################################################
!# ARGUMENTS
!#
!#    spc_param(1) physical frequency for transition
!#                 > 0: High-pass
!#                 < 0: Low-pass
!#    spc_param(2) width of transition in logarithmic wavenumber space
!#
!########################################################################
SUBROUTINE OPR_FILTER_ERF_2D(nx,ny,nz, spc_param, a)

  USE TLAB_VARS, ONLY : isize_txc_dimz
  USE TLAB_VARS, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(*)                                  :: spc_param
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(INOUT) :: a

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip
  TINTEGER sign_pass, off_pass
  TREAL fi,fk,f,fcut_log,damp

! #######################################################################
  IF ( spc_param(1) .GT. 0 ) THEN;
     sign_pass = 1.   ! HIGHPASS
     off_pass  = 0.
  ELSE                ! spc_param(1) <= 0
     sign_pass =-1.   ! LOWPASS
     off_pass  = 1.
  ENDIF

  fcut_log = LOG(spc_param(1))
  DO k = 1,nz
#ifdef USE_MPI
     kglobal = k + ims_offset_k
#else
     kglobal = k
#endif
     IF ( kglobal .LE. g(3)%size/2+1 ) THEN; fk = M_REAL(kglobal-1)/g(3)%scale
     ELSE;                                   fk =-M_REAL(g(3)%size+1-kglobal)/g(3)%scale; ENDIF

     DO i = 1,nx/2+1
#ifdef USE_MPI
        iglobal = i + ims_offset_i/2
#else
        iglobal = i
#endif
        fi = M_REAL(iglobal-1)/g(1)%scale

        f = SQRT(fi**2 + fk**2)
        IF ( f .GT. 0 )  THEN;  damp = (ERF((LOG(f) - fcut_log)/spc_param(2)) + 1.)/2.
        ELSE ;                  damp = C_0_R;
        ENDIF

        ! Set to high- or low-pass
        ! high-pass: damp = 0.0 + damp
        ! low-pass:  damp = 1.0 - damp
        damp = off_pass + sign_pass*damp

        ! apply filter
        DO j = 1,ny
           ip = (j-1)*(nx/2+1) + i
           a(ip,k) = damp*a(ip,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE OPR_FILTER_ERF_2D

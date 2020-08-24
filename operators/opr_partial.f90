#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_PARTIAL1(nlines, bcs, g, u,result, wrk2d)

  USE DNS_TYPES, ONLY : grid_dt

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2),          INTENT(IN)    :: bcs    ! BCs at xmin (1) and xmax (2):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip

! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4P_RHS(g%size,nlines, u, result)

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_RHS(g%size,nlines, u, result)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_RHS(g%size,nlines, u, result)

     END SELECT

     CALL TRIDPSS(g%size,nlines, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5), result,wrk2d)

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C1N6_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_DIRECT   ) ! Not yet implemented
        CALL FDM_C1N6_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     END SELECT

     ip = (bcs(1) + bcs(2)*2)*3
     CALL TRIDSS(g%size,nlines, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3), result)

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL1

! ###################################################################
! ###################################################################
SUBROUTINE OPR_PARTIAL2(nlines, bcs, g, u,result, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2,*),        INTENT(IN)    :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines,g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines,g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(nlines,g%size), INTENT(INOUT) :: wrk3d  ! First derivative, in case needed

! -------------------------------------------------------------------
  TINTEGER ip

! ###################################################################
! Check whether to calculate 1. order derivative
  IF ( .NOT. g%uniform ) THEN
     IF ( g%mode_fdm .eq. FDM_COM4_JACOBIAN .OR. &
          g%mode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
          g%mode_fdm .eq. FDM_COM8_JACOBIAN      ) THEN
        CALL OPR_PARTIAL1(nlines, bcs, g, u,wrk3d, wrk2d)
     ENDIF
  ENDIF

! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_RHS(g%size,nlines, u, result)

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
       ! CALL FDM_C2N6P_RHS(g%size,nlines, u, result)
       CALL FDM_C2N6HP_RHS(g%size,nlines, u, result)

     CASE( FDM_COM8_JACOBIAN )                  ! Not yet implemented
        CALL FDM_C2N6P_RHS(g%size,nlines, u, result)

     END SELECT

     CALL TRIDPSS(g%size,nlines, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3),g%lu2(1,4),g%lu2(1,5), result,wrk2d)

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        IF ( g%uniform ) THEN
           CALL FDM_C2N4_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE ! Not yet implemented
        ENDIF
     CASE( FDM_COM6_JACOBIAN )
        IF ( g%uniform ) THEN
          ! CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
          CALL FDM_C2N6H_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
          ! CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
          CALL FDM_C2N6HNJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented; defaulting to 6. order
        IF ( g%uniform ) THEN
           CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
           CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF

     CASE( FDM_COM6_DIRECT   )
        CALL FDM_C2N6ND_RHS(g%size,nlines, g%lu2(1,4), u, result)

     END SELECT

     ip = (bcs(1,2) + bcs(2,2)*2)*3
     CALL TRIDSS(g%size,nlines, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3), result)

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL2

! ###################################################################
! ###################################################################
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Routines for different specific directions
!########################################################################
SUBROUTINE OPR_PARTIAL_X(type, nx,ny,nz, bcs, g, u, result, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: type      ! OPR_P1     1.order derivative
                                                         ! OPR_P2     2.order derivative
                                                         ! OPR_P2_P1  2. and 1.order derivatives (1. in tmp1)
  TINTEGER,                   INTENT(IN)    :: nx,ny,nz  ! array sizes
  TINTEGER, DIMENSION(2,*),   INTENT(IN)    :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(ny*nz),    INTENT(INOUT) :: wrk2d

  TARGET u, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nyz

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c, p_d

#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => result
     p_b => wrk3d
     p_c => result
     p_d => tmp1
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => u
     p_b => result
     IF ( type .EQ. OPR_P2_P1 ) THEN
        p_c => tmp1
        p_d => wrk3d
     ELSE
        p_c => wrk3d
        p_d => tmp1
     ENDIF
     nyz = ny*nz
#ifdef USE_MPI
  ENDIF
#endif

! -------------------------------------------------------------------
! Local transposition: make x-direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, g%size, g%size, nyz,    p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, g%size, nyz,    g%size, p_b, nyz)
#endif

! ###################################################################
  SELECT CASE( type )

  CASE( OPR_P2 )
     CALL OPR_PARTIAL2(nyz, bcs, g, p_b,p_c, wrk2d,p_d)

  CASE( OPR_P1 )
     CALL OPR_PARTIAL1(nyz, bcs, g, p_b,p_c, wrk2d    )

  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(nyz, bcs, g, p_b,p_c, wrk2d,p_d)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nyz, bcs, g, p_b, p_d, wrk2d)
     ENDIF

  END SELECT

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO       (p_c, nyz, nyz,    g%size, p_b, g%size)
#else
  CALL DNS_TRANSPOSE(p_c, nyz, g%size, nyz,    p_b, g%size)
#endif

  IF ( type .EQ. OPR_P2_P1 ) THEN
#ifdef USE_ESSL
  CALL DGETMO       (p_d, nyz, nyz,    g%size, p_c, g%size)
#else
  CALL DNS_TRANSPOSE(p_d, nyz, g%size, nyz,    p_c, g%size)
#endif
  ENDIF

#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     IF ( type .EQ. OPR_P2_P1 ) THEN ! only if you really want first derivative back
        CALL DNS_MPI_TRPB_I(p_c, tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     ENDIF
     CALL DNS_MPI_TRPB_I(p_b, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c,p_d)

  RETURN
END SUBROUTINE OPR_PARTIAL_X

!########################################################################
!########################################################################
SUBROUTINE OPR_PARTIAL_Z(type, nx,ny,nz, bcs, g, u, result, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: type      ! OPR_P1     1.order derivative
                                                         ! OPR_P2     2.order derivative
                                                         ! OPR_P2_P1  2. and 1.order derivatives (1. in tmp1)
  TINTEGER,                   INTENT(IN)    :: nx,ny,nz  ! array sizes
  TINTEGER, DIMENSION(2,*),   INTENT(IN)    :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*ny),    INTENT(INOUT) :: wrk2d

  TARGET u, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c

#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_K_PARTIAL
#endif

! ###################################################################
  IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
     result = C_0_R
     IF ( type .EQ. OPR_P2_P1 ) tmp1 = C_0_R

  ELSE
! ###################################################################
! -------------------------------------------------------------------
! MPI Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(u, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => result
     IF ( type .EQ. OPR_P2_P1 ) THEN
        p_b => tmp1
        p_c => wrk3d
     ELSE
        p_b => wrk3d
        p_c => tmp1
     ENDIF
     nxy = ims_size_k(id)
 ELSE
#endif
    p_a => u
    p_b => result
    p_c => tmp1
    nxy = nx*ny
#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
  SELECT CASE( type )

  CASE( OPR_P2 )
     CALL OPR_PARTIAL2(nxy, bcs, g, p_a,p_b, wrk2d,p_c)

  CASE( OPR_P1 )
     CALL OPR_PARTIAL1(nxy, bcs, g, p_a,p_b, wrk2d    )

  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(nxy, bcs, g, p_a,p_b, wrk2d,p_c)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxy, bcs, g, p_a,p_c, wrk2d)
     ENDIF

  END SELECT

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_b, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     IF ( type .EQ. OPR_P2_P1 ) THEN
        CALL DNS_MPI_TRPB_K(p_c, tmp1, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c)

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL_Z

!########################################################################
!########################################################################
SUBROUTINE OPR_PARTIAL_Y(type, nx,ny,nz, bcs, g, u, result, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: type      ! OPR_P1     1.order derivative
                                                         ! OPR_P2     2.order derivative
                                                         ! OPR_P2_P1  2. and 1.order derivatives (1. in tmp1)
  TINTEGER,                   INTENT(IN)    :: nx,ny,nz  ! array sizes
  TINTEGER, DIMENSION(2,*),   INTENT(IN)    :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*nz),    INTENT(INOUT) :: wrk2d

  TARGET u, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy, nxz
  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c

! ###################################################################
  IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
     result = C_0_R
     IF ( type .EQ. OPR_P2_P1 ) tmp1 = C_0_R

  ELSE
! ###################################################################
  nxy = nx*ny
  nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_a => u
     p_b => result
     p_c => tmp1
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (u, nxy, nxy, nz, result, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, result, nz)
#endif
     p_a => result
     IF ( type .EQ. OPR_P2_P1 ) THEN
        p_b => tmp1
        p_c => wrk3d
     ELSE
        p_b => wrk3d
        p_c => tmp1
     ENDIF
  ENDIF

! ###################################################################
  SELECT CASE( type )

  CASE( OPR_P2 )
     CALL OPR_PARTIAL2(nxz, bcs, g, p_a,p_b, wrk2d,p_c)

  CASE( OPR_P1 )
     CALL OPR_PARTIAL1(nxz, bcs, g, p_a,p_b, wrk2d    )

  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(nxz, bcs, g, p_a,p_b, wrk2d,p_c)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxz, bcs, g, p_a,p_c, wrk2d)
     ENDIF

  END SELECT

! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_b, nz, nz, nxy, result, nxy)
#else
     CALL DNS_TRANSPOSE(p_b, nz, nxy, nz, result, nxy)
#endif
     IF ( type .EQ. OPR_P2_P1 ) THEN
#ifdef USE_ESSL
        CALL DGETMO       (p_c, nz, nz, nxy, tmp1, nxy)
#else
        CALL DNS_TRANSPOSE(p_c, nz, nxy, nz, tmp1, nxy)
#endif
     ENDIF
  ENDIF

  NULLIFY(p_a,p_b,p_c)

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL_Y

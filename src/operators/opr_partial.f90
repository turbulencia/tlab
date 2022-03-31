#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_PARTIAL1(nlines, bcs, g, u,result, wrk2d)

  USE TLAB_TYPES, ONLY : grid_dt

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
        CALL FDM_C1N4P_RHS( g%size,nlines, u, result)

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_RHS( g%size,nlines, u, result)

     CASE( FDM_COM6_JACPENTA ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6MP_RHS(g%size,nlines, u, result)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_RHS( g%size,nlines, u, result)

     END SELECT

     IF (.NOT. (g%mode_fdm .EQ. FDM_COM6_JACPENTA)) THEN
        CALL TRIDPSS(  g%size,nlines, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5),                       result,wrk2d)
     ELSE
        CALL PENTADPSS(g%size,nlines, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5),g%lu1(1,6),g%lu1(1,7), result)
     ENDIF

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_RHS( g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C1N6_RHS( g%size,nlines, bcs(1),bcs(2), u, result)
     
     CASE( FDM_COM6_JACPENTA )
        CALL FDM_C1N6M_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_RHS( g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_DIRECT   ) ! Not yet implemented
        CALL FDM_C1N6_RHS( g%size,nlines, bcs(1),bcs(2), u, result)

     END SELECT

     ip = (bcs(1) + bcs(2)*2)*5
     IF (.NOT. (g%mode_fdm .EQ. FDM_COM6_JACPENTA)) THEN
        CALL TRIDSS(   g%size,nlines, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3),                             result)
     ELSE
        CALL PENTADSS2(g%size,nlines, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3),g%lu1(1,ip+4),g%lu1(1,ip+5), result)
   ENDIF

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL1
! ###################################################################
! ###################################################################
SUBROUTINE OPR_PARTIAL1_IBM(nlines, bcs, g, u,result, wrk2d, wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt
  USE DNS_IBM,    ONLY : fld_ibm
  USE DNS_IBM,    ONLY : nobi,    nobj,   nobk
  USE DNS_IBM,    ONLY : nobi_b,  nobj_b, nobk_b 
  USE DNS_IBM,    ONLY : nobi_e,  nobj_e, nobk_e 
  USE DNS_IBM,    ONLY : isize_nobi,    isize_nobj,    isize_nobk
  USE DNS_IBM,    ONLY : isize_nobi_be, isize_nobj_be, isize_nobk_be 
  USE DNS_IBM,    ONLY : ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

! ############################################# ! 
! DEBUG
#ifdef IBM_DEBUG
#ifdef USE_MPI
  use TLAB_MPI_VARS,   only : ims_pro
#endif
#endif
! ############################################# ! 
   
  IMPLICIT NONE

#include "integers.h"
   
  TINTEGER,                        INTENT(IN)    :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2,*),        INTENT(IN)    :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(nlines*g%size), INTENT(INOUT) :: wrk3d

  TINTEGER, PARAMETER                            :: is = i0 ! scalar index; if 0, then velocity

  ! -------------------------------------------------------------------

! ############################################# ! 
! debugging
#ifdef IBM_DEBUG
#ifdef USE_MPI
#else
  TINTEGER, parameter  ::  ims_pro=0  
#endif
#endif
! ############################################ ! 

  ! IBM not for scalar fields! (will be implemented later)
  ! modify incoming u fields (fill solids with spline functions, depending on direction)

  SELECT CASE (g%name)
   
  CASE('x')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_partial_', g%name ! debug
#endif
    IF (ims_pro_ibm_x) THEN
      CALL IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
      CALL OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL1(nlines, bcs, g, u,       result, wrk2d)  ! no splines needed  
    ENDIF

  CASE('y')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_partial_', g%name ! debug
#endif
    IF (ims_pro_ibm_y) THEN
      CALL IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
      CALL OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL1(nlines, bcs, g, u,       result, wrk2d)  ! no splines needed
    ENDIF

  CASE('z')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_partial_', g%name ! debug
#endif
    IF (ims_pro_ibm_z) THEN
      CALL IBM_SPLINE_XYZ(is, u, fld_ibm, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
      CALL OPR_PARTIAL1(nlines, bcs, g, fld_ibm, result, wrk2d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL1(nlines, bcs, g, u,       result, wrk2d)  ! no splines needed
    ENDIF
   
  END SELECT

  RETURN
END SUBROUTINE OPR_PARTIAL1_IBM
! ###################################################################
! ###################################################################
SUBROUTINE OPR_IBM(nlines, g, u,result, wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt
  USE DNS_IBM,    ONLY : nobi,    nobj,   nobk
  USE DNS_IBM,    ONLY : nobi_b,  nobj_b, nobk_b 
  USE DNS_IBM,    ONLY : nobi_e,  nobj_e, nobk_e 
  USE DNS_IBM,    ONLY : isize_nobi,    isize_nobj,    isize_nobk
  USE DNS_IBM,    ONLY : isize_nobi_be, isize_nobj_be, isize_nobk_be 
   
  IMPLICIT NONE

#include "integers.h"
   
  TINTEGER,                        INTENT(IN)    :: nlines
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines*g%size), INTENT(INOUT) :: wrk3d

  TINTEGER, PARAMETER                            :: is = i0 ! scalar index; if 0, then velocity
  
  ! -------------------------------------------------------------------
  ! IBM not for scalar fields! (will be implemented later)
  ! modify incoming u fields (fill solids with spline functions, depending on direction)
  
  SELECT CASE (g%name)
  CASE('x')
    CALL IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
  CASE('y')
    CALL IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
  CASE('z')
    CALL IBM_SPLINE_XYZ(is, u, result, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
  END SELECT

  RETURN
END SUBROUTINE OPR_IBM
! ###################################################################
! ###################################################################
SUBROUTINE OPR_PARTIAL2(nlines, bcs, g, u,result, wrk2d,wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt

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

  ! --- local declarations -------------------------------------------------------------------

  TINTEGER is !,ip
  PARAMETER(is=-1)

  ! ### end of declarations ##################################################################

  CALL OPR_PARTIAL2D(is,nlines, bcs, g, u,result, wrk2d,wrk3d)
! ! ###################################################################
! ! Check whether to calculate 1. order derivative
!   IF ( .NOT. g%uniform ) THEN
!      IF ( g%mode_fdm .eq. FDM_COM4_JACOBIAN .OR. &
!           g%mode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
!           g%mode_fdm .eq. FDM_COM8_JACOBIAN      ) THEN
!         CALL OPR_PARTIAL1(nlines, bcs, g, u,wrk3d, wrk2d)
!      ENDIF
!   ENDIF
!
! ! ###################################################################
!   IF ( g%periodic ) THEN
!      SELECT CASE( g%mode_fdm )
!
!      CASE( FDM_COM4_JACOBIAN )
!         CALL FDM_C2N4P_RHS(g%size,nlines, u, result)
!
!      CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA ) ! Direct = Jacobian because uniform grid
!        ! CALL FDM_C2N6P_RHS(g%size,nlines, u, result)
!        CALL FDM_C2N6HP_RHS(g%size,nlines, u, result)
!
!      CASE( FDM_COM8_JACOBIAN )                  ! Not yet implemented
!         CALL FDM_C2N6P_RHS(g%size,nlines, u, result)
!
!      END SELECT
!
!      CALL TRIDPSS(g%size,nlines, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3),g%lu2(1,4),g%lu2(1,5), result,wrk2d)
!
! ! -------------------------------------------------------------------
!   ELSE
!      SELECT CASE( g%mode_fdm )
!
!      CASE( FDM_COM4_JACOBIAN )
!         IF ( g%uniform ) THEN
!            CALL FDM_C2N4_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
!         ELSE ! Not yet implemented
!         ENDIF
!      CASE( FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA )
!         IF ( g%uniform ) THEN
!           CALL FDM_C2N6H_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
!         ELSE
!           CALL FDM_C2N6HNJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
!         ENDIF
!
!      CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented; defaulting to 6. order
!         IF ( g%uniform ) THEN
!            CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
!         ELSE
!            CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
!         ENDIF
!
!      CASE( FDM_COM6_DIRECT   )
!         CALL FDM_C2N6ND_RHS(g%size,nlines, g%lu2(1,4), u, result)
!
!      END SELECT

!      ip = (bcs(1,2) + bcs(2,2)*2)*3
!      CALL TRIDSS(g%size,nlines, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3), result)

!   ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL2
! ###################################################################################
! ### First Derivative
! ### Second Derivative includes
! ###       factor 1    (is=-1)
! ###       viscosity   (is= 0)
! ###       diffusivity (is= 1,inb_scal)
SUBROUTINE OPR_PARTIAL2D(is,nlines, bcs, g, u,result, wrk2d,wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: is     ! scalar index; if 0, then velocity
  TINTEGER,                        INTENT(IN)    :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2,*),        INTENT(IN)    :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines,g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines,g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(nlines,g%size), INTENT(INOUT) :: wrk3d  ! First derivative

  TREAL, DIMENSION(:,:), POINTER :: lu2_p


! -------------------------------------------------------------------
  TINTEGER ip

  ! ###################################################################
  ! always calculate first derivative as this routine is normally called from opr_burgers
  ! ###################################################################
  CALL OPR_PARTIAL1(nlines, bcs, g, u,wrk3d, wrk2d)

  IF ( is .GE. 0 ) THEN
     IF ( g%periodic ) THEN
        lu2_p => g%lu2d(:,is*5+1:)  ! periodic;     including diffusivity/viscosity
     ELSE
        lu2_p => g%lu2d(:,is*3+1:)  ! non-periodic; including diffusivity/viscosity
     ENDIF
  ELSE
     IF ( g%periodic ) THEN
        lu2_p => g%lu2(:,1:)        ! periodic;     plain derivative
     ELSE
        ip = (bcs(1,2)+bcs(2,2)*2)*3! non-periodic; plain derivative
        lu2_p => g%lu2(:,ip+1:)
     ENDIF
  ENDIF


  ! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_RHS(g%size,nlines, u, result)

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACPENTA ) ! Direct = Jacobian because uniform grid
        CALL FDM_C2N6HP_RHS(g%size,nlines, u, result)

     CASE( FDM_COM8_JACOBIAN )                  ! Not yet implemented
        CALL FDM_C2N6P_RHS(g%size,nlines, u, result)

     END SELECT

     CALL TRIDPSS(g%size,nlines, lu2_p(1,1),lu2_p(1,2),lu2_p(1,3),lu2_p(1,4),lu2_p(1,5), result,wrk2d)

  ! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )

     CASE( FDM_COM4_JACOBIAN )
        IF ( g%uniform ) THEN
           CALL FDM_C2N4_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE ! Not yet implemented
        ENDIF

     CASE( FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA )
        IF ( g%uniform ) THEN
          CALL FDM_C2N6H_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
           ! need first derivative from above
          CALL FDM_C2N6HNJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented; defaulting to 6. order
        IF ( g%uniform ) THEN
           CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
           ! Need first derivative from above
           CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF

     CASE( FDM_COM6_DIRECT   )
        CALL FDM_C2N6ND_RHS(g%size,nlines, g%lu2(1,4), u, result)

     END SELECT

     CALL TRIDSS(g%size,nlines, lu2_p(1,1),lu2_p(1,2),lu2_p(1,3), result)

  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL2D
! ###################################################################
! ###################################################################
SUBROUTINE OPR_PARTIAL2D_IBM(is, nlines, bcs, g, u, result, wrk2d, wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt
  USE DNS_IBM,    ONLY : fld_ibm
  USE DNS_IBM,    ONLY : nobi,    nobj,   nobk
  USE DNS_IBM,    ONLY : nobi_b,  nobj_b, nobk_b 
  USE DNS_IBM,    ONLY : nobi_e,  nobj_e, nobk_e 
  USE DNS_IBM,    ONLY : isize_nobi,    isize_nobj,    isize_nobk
  USE DNS_IBM,    ONLY : isize_nobi_be, isize_nobj_be, isize_nobk_be 
  USE DNS_IBM,    ONLY : ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z

! ############################################# ! 
! DEBUG
#ifdef IBM_DEBUG
#ifdef USE_MPI
  use TLAB_MPI_VARS,   only : ims_pro
#endif
#endif
! ############################################# ! 
   
  IMPLICIT NONE
   
  TINTEGER,                        INTENT(IN)             :: is     ! scalar index; if 0, then velocity
  TINTEGER,                        INTENT(IN)             :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2,*),        INTENT(IN)             :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
                                                                    !     0 biased, non-zero
                                                                    !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)             :: g
  TREAL, DIMENSION(nlines,g%size), INTENT(IN),    TARGET  :: u
  TREAL, DIMENSION(nlines,g%size), INTENT(OUT)            :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT)          :: wrk2d
  TREAL, DIMENSION(nlines,g%size), INTENT(INOUT)          :: wrk3d  ! First derivative 

  TREAL, DIMENSION(:,:),                          POINTER :: p_fld
  TREAL, DIMENSION(:),                            POINTER :: p_fld_ibm
  
  ! -------------------------------------------------------------------

! ############################################# ! 
! debugging
#ifdef IBM_DEBUG
#ifdef USE_MPI
#else
  TINTEGER, parameter  ::  ims_pro=0  
#endif
#endif
! ############################################ ! 

  ! pointer to field
  p_fld => u

  ! -------------------------------------------------------------------

  ! IBM not for scalar fields! (will be implemented later)
  ! modify incoming u fields (fill solids with spline functions, depending on direction)

  SELECT CASE (g%name)
   
  CASE('x')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_burgers_', g%name ! debug
#endif
    IF (ims_pro_ibm_x) THEN
      CALL IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, wrk3d)
      p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld,     result, wrk2d, wrk3d)  ! no splines needed  
    ENDIF

  CASE('y')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_burgers_', g%name ! debug
#endif
    IF (ims_pro_ibm_y) THEN
      CALL IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, wrk3d)
      p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld,     result, wrk2d, wrk3d)  ! no splines needed
    ENDIF

  CASE('z')
#ifdef IBM_DEBUG
    IF (ims_pro == 0) write(*,*) 'ibm_burgers_', g%name ! debug
#endif
    IF (ims_pro_ibm_z) THEN
      CALL IBM_SPLINE_XYZ(is, p_fld, fld_ibm, g, nlines, isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, wrk3d)
      p_fld_ibm => fld_ibm                                                     ! pointer to modified velocity
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld_ibm, result, wrk2d, wrk3d)  ! now with modified u fields
    ELSE
      CALL OPR_PARTIAL2D(is, nlines, bcs, g, p_fld,     result, wrk2d, wrk3d)  ! no splines needed
    ENDIF
   
  END SELECT

  ! -------------------------------------------------------------------

  NULLIFY(p_fld, p_fld_ibm)
   
  RETURN
END SUBROUTINE OPR_PARTIAL2D_IBM
! ###################################################################
#include "dns_error.h"
! ###################################################################

SUBROUTINE OPR_PARTIAL0_INT(dir, nlines, g, u,result, wrk2d)

  USE TLAB_TYPES,     ONLY : grid_dt
  USE TLAB_PROCS,     ONLY : TLAB_STOP, TLAB_WRITE_ASCII
  USE TLAB_CONSTANTS, ONLY : efile

  IMPLICIT NONE
 
  TINTEGER,                            INTENT(IN)    :: dir    ! scalar direction flag
                                                               !     0 'vp' --> vel. to pre. grid
                                                               !     1 'pv' --> pre. to vel. grid
  TINTEGER,                            INTENT(IN)    :: nlines ! number of lines to be solved
  TYPE(grid_dt),                       INTENT(IN)    :: g
  TREAL, DIMENSION(nlines,g%size),     INTENT(IN)    :: u
  TREAL, DIMENSION(nlines,g%size),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),            INTENT(INOUT) :: wrk2d
 
! ###################################################################
! Interpolation, direction 'vp': vel. --> pre. grid
  IF ( dir .EQ. 0 ) THEN 
    IF ( g%periodic ) THEN
      SELECT CASE( g%mode_fdm )        
      CASE DEFAULT
        CALL FDM_C0INTVP6P_RHS(g%size,nlines, u, result)
      END SELECT
      CALL TRIDPSS(g%size,nlines, g%lu0i(1,1),g%lu0i(1,2),g%lu0i(1,3),g%lu0i(1,4),g%lu0i(1,5), result,wrk2d)
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
      CALL TLAB_STOP(DNS_ERROR_NOTIMPL)
    ENDIF
! Interpolation, direction 'pv': pre. --> vel. grid
  ELSE IF ( dir .EQ. 1 ) THEN 
    IF ( g%periodic ) THEN
      SELECT CASE( g%mode_fdm )        
      CASE DEFAULT
        CALL FDM_C0INTPV6P_RHS(g%size,nlines, u, result)
      END SELECT
      CALL TRIDPSS(g%size,nlines, g%lu0i(1,1),g%lu0i(1,2),g%lu0i(1,3),g%lu0i(1,4),g%lu0i(1,5), result,wrk2d)
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL0_INT. Non-periodic case not implemented.')
      CALL TLAB_STOP(DNS_ERROR_NOTIMPL)
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL0_INT

! ###################################################################
! ###################################################################

SUBROUTINE OPR_PARTIAL1_INT(dir, nlines, g, u,result, wrk2d)

  USE TLAB_TYPES,     ONLY : grid_dt
  USE TLAB_PROCS,     ONLY : TLAB_STOP, TLAB_WRITE_ASCII
  USE TLAB_CONSTANTS, ONLY : efile

  IMPLICIT NONE
 
  TINTEGER,                            INTENT(IN)    :: dir    ! scalar direction flag
                                                               !     0 'vp' --> vel. to pre. grid
                                                               !     1 'pv' --> pre. to vel. grid
  TINTEGER,                            INTENT(IN)    :: nlines ! number of lines to be solved
  TYPE(grid_dt),                       INTENT(IN)    :: g
  TREAL, DIMENSION(nlines,g%size),     INTENT(IN)    :: u
  TREAL, DIMENSION(nlines,g%size),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),            INTENT(INOUT) :: wrk2d
  
! ###################################################################
! 1st interpolatory derivative, direction 'vp': vel. --> pre. grid
  IF ( dir .EQ. 0 ) THEN
    IF ( g%periodic ) THEN
      SELECT CASE( g%mode_fdm )
      CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C1INTVP6P_RHS(g%size,nlines, u, result)
      END SELECT
      CALL TRIDPSS(g%size,nlines, g%lu1i(1,1),g%lu1i(1,2),g%lu1i(1,3),g%lu1i(1,4),g%lu1i(1,5), result,wrk2d)
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
      CALL TLAB_STOP(DNS_ERROR_NOTIMPL)
    ENDIF
! 1st interpolatory derivative, direction 'pv': pre. --> vel. grid
  ELSE IF ( dir .EQ. 1 ) THEN
    IF ( g%periodic ) THEN
      SELECT CASE( g%mode_fdm )
      CASE( FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM8_JACOBIAN )
        CALL FDM_C1INTPV6P_RHS(g%size,nlines, u, result)
      END SELECT
      CALL TRIDPSS(g%size,nlines, g%lu1i(1,1),g%lu1i(1,2),g%lu1i(1,3),g%lu1i(1,4),g%lu1i(1,5), result,wrk2d)
    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'OPR_PARTIAL1_INT. Non-periodic case not implemented.')
      CALL TLAB_STOP(DNS_ERROR_NOTIMPL)
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE OPR_PARTIAL1_INT
 
! ###################################################################
! ###################################################################

#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Routines for different specific directions
!########################################################################
SUBROUTINE OPR_PARTIAL_X(type, nx,ny,nz, bcs, g, u, result, tmp1, wrk2d,wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt
  USE DNS_IBM,    ONLY : ibm_partial
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                   INTENT(IN)    :: type      ! OPR_P1           1.order derivative
                                                         ! OPR_P2           2.order derivative
                                                         ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
                                                         ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
                                                         ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
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
  TINTEGER, PARAMETER :: id = TLAB_MPI_I_PARTIAL
#endif

! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL TLAB_MPI_TRPF_I(u, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
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
     CALL OPR_PARTIAL2(     nyz, bcs, g, p_b,p_c, wrk2d,p_d)

  CASE( OPR_P1 )
    IF (ibm_partial) THEN
      CALL OPR_PARTIAL1_IBM(nyz, bcs, g, p_b,p_c, wrk2d,p_d)
    ELSE
      CALL OPR_PARTIAL1(    nyz, bcs, g, p_b,p_c, wrk2d    )
    ENDIF

  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(     nyz, bcs, g, p_b,p_c, wrk2d,p_d)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nyz, bcs, g, p_b, p_d, wrk2d)
     ENDIF
  
  CASE( OPR_P0_INT_VP )
     CALL OPR_PARTIAL0_INT(i0, nyz, g, p_b,p_c, wrk2d)

  CASE( OPR_P0_INT_PV )
     CALL OPR_PARTIAL0_INT(i1, nyz, g, p_b,p_c, wrk2d)

  CASE( OPR_P1_INT_VP )
     CALL OPR_PARTIAL1_INT(i0, nyz, g, p_b,p_c, wrk2d)

  CASE( OPR_P1_INT_PV )
     CALL OPR_PARTIAL1_INT(i1, nyz, g, p_b,p_c, wrk2d)

  CASE( OPR_P0_IBM )
     CALL OPR_IBM(             nyz, g, p_b,p_c, p_d)

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
        CALL TLAB_MPI_TRPB_I(p_c, tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     ENDIF
     CALL TLAB_MPI_TRPB_I(p_b, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c,p_d)

  RETURN
END SUBROUTINE OPR_PARTIAL_X

!########################################################################
!########################################################################
SUBROUTINE OPR_PARTIAL_Z(type, nx,ny,nz, bcs, g, u, result, tmp1, wrk2d,wrk3d)

  USE TLAB_TYPES,    ONLY : grid_dt
  USE DNS_IBM,       ONLY : ibm_partial
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                   INTENT(IN)    :: type      ! OPR_P1           1.order derivative
                                                         ! OPR_P2           2.order derivative
                                                         ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
                                                         ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
                                                         ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
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
  TINTEGER, PARAMETER :: id = TLAB_MPI_K_PARTIAL
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
     CALL TLAB_MPI_TRPF_K(u, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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
    IF (ibm_partial) THEN
      CALL OPR_PARTIAL1_IBM(nxy, bcs, g, p_a,p_b, wrk2d,p_c)
    ELSE
      CALL OPR_PARTIAL1(    nxy, bcs, g, p_a,p_b, wrk2d    )
    ENDIF
  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(     nxy, bcs, g, p_a,p_b, wrk2d,p_c)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxy, bcs, g, p_a,p_c, wrk2d)
     ENDIF
  
  CASE( OPR_P0_INT_VP )
     CALL OPR_PARTIAL0_INT(i0, nxy, g, p_a,p_b, wrk2d)
 
  CASE( OPR_P0_INT_PV )
     CALL OPR_PARTIAL0_INT(i1, nxy, g, p_a,p_b, wrk2d)

  CASE( OPR_P1_INT_VP )
     CALL OPR_PARTIAL1_INT(i0, nxy, g, p_a,p_b, wrk2d)

  CASE( OPR_P1_INT_PV )
     CALL OPR_PARTIAL1_INT(i1, nxy, g, p_a,p_b, wrk2d)

  CASE( OPR_P0_IBM )
    CALL OPR_IBM(              nxy, g, p_a,p_b, p_c)

  END SELECT

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL TLAB_MPI_TRPB_K(p_b, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     IF ( type .EQ. OPR_P2_P1 ) THEN
        CALL TLAB_MPI_TRPB_K(p_c, tmp1, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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

  USE TLAB_TYPES, ONLY : grid_dt
  USE DNS_IBM,    ONLY : ibm_partial
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                       INTENT(IN)    :: type      ! OPR_P1           1.order derivative
                                                             ! OPR_P2           2.order derivative
                                                             ! OPR_P2_P1        2. and 1.order derivatives (1. in tmp1)
                                                             ! OPR_P0_INT_VP/PV interpolation              (vel.<->pre.)
                                                             ! OPR_P1_INT_VP/PV 1.order int. derivative    (vel.<->pre.)
  TINTEGER,                       INTENT(IN)    :: nx,ny,nz  ! array sizes
  TINTEGER, DIMENSION(2,*),       INTENT(IN)    :: bcs       ! BCs at xmin (1,*) and xmax (2,*)
  TYPE(grid_dt),                  INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: u
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d 
  TREAL, DIMENSION(nx*nz),        INTENT(INOUT) :: wrk2d

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
    IF (ibm_partial) THEN
      CALL OPR_PARTIAL1_IBM(nxz, bcs, g, p_a,p_b, wrk2d,p_c)
    ELSE
      CALL OPR_PARTIAL1(    nxz, bcs, g, p_a,p_b, wrk2d    )
    ENDIF

  CASE( OPR_P2_P1 )
     CALL OPR_PARTIAL2(     nxz, bcs, g, p_a,p_b, wrk2d,p_c)

! Check whether we need to calculate the 1. order derivative
     IF ( g%uniform .OR. g%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxz, bcs, g, p_a,p_c, wrk2d)
     ENDIF
  
  CASE( OPR_P0_INT_VP )
     CALL OPR_PARTIAL0_INT(i0, nxz, g, p_a,p_b, wrk2d)
 
  CASE( OPR_P0_INT_PV )
     CALL OPR_PARTIAL0_INT(i1, nxz, g, p_a,p_b, wrk2d)
 
  CASE( OPR_P1_INT_VP )
     CALL OPR_PARTIAL1_INT(i0, nxz, g, p_a,p_b, wrk2d)
 
  CASE( OPR_P1_INT_PV )
     CALL OPR_PARTIAL1_INT(i1, nxz, g, p_a,p_b, wrk2d)

  CASE( OPR_P0_IBM )
     CALL OPR_IBM(             nxz, g, p_a,p_b, p_c)

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
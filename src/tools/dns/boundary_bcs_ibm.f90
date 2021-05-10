#include "types.h"
#include "dns_const.h"

! #include "dns_error.h"        ! needed if DNS_STOP is used, in this case: define specific error 
                                ! (e.g. number of objects do not fit with grid)
! #ifdef USE_MPI
! #include "dns_const_mpi.h"    ! most likely not needed, if no specific MPI calls are present
! #endif

!########################################################################
!# HISTORY
!#
!# 2021/05/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the boundary values for immersed object(s) in the flow  
!# domain. The epsilon field epsi can be seen as an indicator field:
!#    epsi(i,j,k) = 0   for 
!#    epsi(i,j,k) = 1   for  
!#
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################
module DNS_IBM

  ! use ...

  implicit none

  ! TINTEGER
  ! TREAL 

  ! all functions/subroutines are private by default, make neeeded ones public 
  private 
  public :: BOUNDARY_BCS_IBM

contains
 !########################################################################
  SUBROUTINE BOUNDARY_BCS_IBM(nx, ny, nz, g)  ! (ibc, nx,ny,nz, g, u, bcs_hb,bcs_ht, wrk1d,tmp1,tmp2)

    USE DNS_TYPES, ONLY     : grid_dt
    USE DNS_CONSTANTS, ONLY : lfile

    USE DNS_GLOBAL,ONLY     : imax,jmax,kmax
  ! MPI
#ifdef USE_MPI 
    USE DNS_MPI, ONLY       : ims_offset_i, ims_offset_j, ims_offset_k, ims_pro
    USE DNS_MPI, ONLY       : ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k
    USE DNS_MPI, ONLY       : ims_npro, ims_npro_i, ims_npro_j, ims_npro_k
#endif 

    IMPLICIT NONE

#include "integers.h"
  ! MPI
#ifdef USE_MPI
#include "mpif.h"
#else
    TINTEGER, parameter  :: ims_offset_i=0, ims_offset_j=0, ims_offset_k=0
    TINTEGER, parameter  :: ims_pro_i=0,    ims_pro_j=0,    ims_pro_k=0,    ims_pro=0
    TINTEGER, parameter  :: ims_npro_i=1,   ims_npro_j=1,   ims_npro_k=1,   ims_npro=0 
#endif

! initialize with intent in and out and inout 
    TINTEGER nx,ny,nz!, ibc
    TINTEGER :: istart, iend, jstart, jend, kstart, kend
    TYPE(grid_dt),                      INTENT(IN)  :: g
    ! TREAL, DIMENSION(nx*nz,ny), TARGET, INTENT(IN)  :: u         ! they are transposed below
    ! TREAL, DIMENSION(nx*nz,ny), TARGET              :: tmp1,tmp2 ! they are transposed below
    ! TREAL, DIMENSION(g%size,3), TARGET              :: wrk1d
    ! TREAL, DIMENSION(nx*nz),    TARGET, INTENT(OUT) :: bcs_hb,bcs_ht

  ! -------------------------------------------------------------------
    ! TINTEGER nxz, nxy

    ! TREAL, DIMENSION(:,:), POINTER :: p_org,p_dst
    ! TREAL, DIMENSION(:),   POINTER :: p_bcs_hb,p_bcs_ht
    ! TREAL, DIMENSION(:),   POINTER :: a,b,c
    CHARACTER*250 line1


  ! ###################################################################


    istart=ims_offset_i; iend=ims_offset_i + imax -1   
    jstart=ims_offset_j; jend=ims_offset_j + jmax -1
    kstart=ims_offset_k; kend=ims_offset_k + kmax -1

    ! for debugging
    IF ( ims_pro .EQ. 0 ) THEN 
      WRITE(*,*) '=== Initialization of Grid and Decomposition ==='
      WRITE(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax*ims_npro_j ,' x ', kmax*ims_npro_k 
      WRITE(*,*) 'DECOMP: ranks', ims_npro_i,' x ',ims_npro_j,' x ',ims_npro_k 
      WRITE(*,*) '        grid ', imax,      ' x ',jmax,      ' x ',kmax 
      WRITE(*,*) '=== Global Grid Information of each Task     ==='
    ENDIF
    write(*,*) 'Task_i number ', ims_pro_i, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend
    write(*,*) 'Task_j number ', ims_pro_j, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend
    write(*,*) 'Task_k number ', ims_pro_k, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend

      ! WRITE(line1,*) itime; line1 = 'Starting time integration at It'//TRIM(ADJUSTL(line1))//'.'
      ! CALL IO_WRITE_ASCII(lfile,line1)
  


     
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !   IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
  !   bcs_hb = C_0_R; bcs_ht = C_0_R 

  !   ELSE
  ! ! ###################################################################
  !   nxy = nx*ny 
  !   nxz = nx*nz

  !   a => wrk1d(:,1)
  !   b => wrk1d(:,2)
  !   c => wrk1d(:,3)

  ! ! -------------------------------------------------------------------
  ! ! Make y  direction the last one
  ! ! -------------------------------------------------------------------
  !   IF ( nz .EQ. 1 ) THEN
  !     p_org => u
  !     p_dst => tmp1
  !     p_bcs_hb => bcs_hb
  !     p_bcs_ht => bcs_ht
  !   ELSE
  ! #ifdef USE_ESSL
  !     CALL DGETMO(u, nxy, nxy, nz, tmp1, nz)
  ! #else
  !     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, tmp1, nz)
  ! #endif
  !     p_org => tmp1
  !     p_dst => tmp2
  !     p_bcs_hb => tmp1(:,1)
  !     p_bcs_ht => tmp1(:,2)
  !   ENDIF

  ! ! ###################################################################
  !   SELECT CASE( g%mode_fdm )
      
  !   CASE( FDM_COM4_JACOBIAN ) !not yet implemented
      
  !   CASE( FDM_COM6_JACOBIAN )
  !     CALL FDM_C1N6_BCS_LHS(ny,     ibc, g%jac, a,b,c)
  !     CALL FDM_C1N6_BCS_RHS(ny,nxz, ibc,        p_org,p_dst)

  !   CASE( FDM_COM6_DIRECT   ) !not yet implemented
  !     CALL FDM_C1N6_BCS_LHS(ny,     ibc, g%jac, a,b,c)
  !     CALL FDM_C1N6_BCS_RHS(ny,nxz, ibc,        p_org,p_dst)
      
  !   CASE( FDM_COM8_JACOBIAN ) !not yet implemented
      
  !   END SELECT
    
  ! ! -------------------------------------------------------------------
  !   IF      ( ibc .EQ. 1 ) THEN
  !     CALL TRIDFS(ny-1,     a(2),b(2),c(2))
  !     CALL TRIDSS(ny-1,nxz, a(2),b(2),c(2), p_dst(1,2))
  !     p_bcs_hb(:) = p_dst(:,1 ) + c(1) *p_dst(:,2)   

  !   ELSE IF ( ibc .EQ. 2 ) THEN
  !     CALL TRIDFS(ny-1,     a,b,c)
  !     CALL TRIDSS(ny-1,nxz, a,b,c, p_dst)
  !     p_bcs_ht(:) = p_dst(:,ny) + a(ny)*p_dst(:,ny-1)

  !   ELSE IF ( ibc .EQ. 3 ) THEN
  !     CALL TRIDFS(ny-2,     a(2),b(2),c(2))
  !     CALL TRIDSS(ny-2,nxz, a(2),b(2),c(2), p_dst(1,2))
  !     p_bcs_hb(:) = p_dst(:,1 ) + c(1) *p_dst(:,2)   
  !     p_bcs_ht(:) = p_dst(:,ny) + a(ny)*p_dst(:,ny-1)

  !   ENDIF

  ! ! ###################################################################
  ! ! -------------------------------------------------------------------
  ! ! Put bcs arrays in correct order
  ! ! -------------------------------------------------------------------
  !   IF ( nz .GT. 1 ) THEN
  ! #ifdef USE_ESSL
  !     IF ( ibc .EQ. 1 .OR. ibc .EQ. 3 ) CALL DGETMO(p_bcs_hb, nz, nz, nx, bcs_hb, nx)
  !     IF ( ibc .EQ. 2 .OR. ibc .EQ. 3 ) CALL DGETMO(p_bcs_ht, nz, nz, nx, bcs_ht, nx)
  ! #else
  !     IF ( ibc .EQ. 1 .OR. ibc .EQ. 3 ) CALL DNS_TRANSPOSE(p_bcs_hb, nz, nx, nz, bcs_hb, nx)
  !     IF ( ibc .EQ. 2 .OR. ibc .EQ. 3 ) CALL DNS_TRANSPOSE(p_bcs_ht, nz, nx, nz, bcs_ht, nx)
  ! #endif
  !   ENDIF
  !   NULLIFY(p_org,p_dst,p_bcs_hb,p_bcs_ht,a,b,c)

  !   ENDIF

    RETURN
  END SUBROUTINE BOUNDARY_BCS_IBM
  !########################################################################
end module DNS_IBM
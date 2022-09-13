#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_CHEM(chemistry, nx,ny,nz, is, s, source)

  USE TLAB_TYPES,  ONLY : term_dt, profiles_dt
  USE TLAB_VARS, ONLY : sbg, damkohler
  USE TLAB_VARS, ONLY : g

  IMPLICIT NONE

  TYPE(term_dt),                INTENT(IN)  :: chemistry
  TINTEGER,                     INTENT(IN)  :: nx,ny,nz, is
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: source

! -----------------------------------------------------------------------
  TREAL dummy, dummy2, PROFILES
  type(profiles_dt) prof_loc
  TINTEGER i,j,k
  external PROFILES

!########################################################################
  SELECT CASE( chemistry%type )

  CASE( EQNS_CHEM_LAYEREDRELAXATION )
     prof_loc%type = PROFILE_TANH
     prof_loc%ymean = sbg(is)%ymean
     prof_loc%thick =-chemistry%parameters(3) *C_05_R
     prof_loc%mean = C_05_R
     prof_loc%delta = C_1_R
     prof_loc%lslope=C_0_R
     prof_loc%uslope=C_0_R
     DO i=1,nx
        DO k=1,nz
           DO j=1,ny
              source(i+(j-1)*nx+(k-1)*nx*ny) = PROFILES(prof_loc, g(2)%nodes(j)-chemistry%parameters(2)) ! strength constant
           ENDDO
        ENDDO
     ENDDO
     
     dummy  =-damkohler(is) /chemistry%parameters(1)
     source = dummy *source *s(:,is)

  CASE( EQNS_CHEM_QUADRATIC  )
     dummy  = damkohler(is) *chemistry%parameters(is)
     source = dummy *s(:,2) *s(:,3)

  CASE( EQNS_CHEM_QUADRATIC3 )
     dummy  = damkohler(is) *chemistry%parameters(is)

     IF      ( is .GE. 1 .AND. is .LE. 3 ) THEN
        source = dummy *s(:,2) *s(:,3)
     ELSE IF ( is .GE. 4 .AND. is .LE. 6 ) THEN
        source = dummy *s(:,4) *s(:,5)
     ELSE IF ( is .GE. 7 .AND. is .LE. 9 ) THEN
        source = dummy *s(:,7) *s(:,8)
     ENDIF
     
  CASE( EQNS_CHEM_OZONE )
     dummy  = damkohler(is)
     IF ( is .EQ. 4 ) dummy =-dummy
     
     source =-chemistry%parameters(1) /( C_1_R +chemistry%parameters(2)*s(:,1) )
     source = EXP(source)

     IF ( is .EQ. 4 ) THEN
        dummy2 = C_1_R + chemistry%parameters(3)
        source = dummy *( dummy2 *s(:,4) - source *s(:,2) *s(:,3) )
     ELSE
        source = dummy *(         s(:,4) - source *s(:,2) *s(:,3) )
     ENDIF
     
  END SELECT
  
  RETURN
END SUBROUTINE FI_CHEM

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2011/04/18 - Alberto de Lozar
!#              Implementing inertial effects
!# 2011/10/01 - J.P. Mellado
!#              Moving the buffer to be part of the pressure calculation
!# 2011/11/28 - C. Ansorge
!#              Combining BLAS and OMP
!#
!########################################################################
!# DESCRIPTION
!#
!# Branching among different formulations of the RHS.
!# 
!# Be careful to define here the pointers and to enter the RHS routines 
!# with individual fields. Definition of pointer inside of RHS routines
!# decreased performance considerably (at least in JUGENE)
!#
!########################################################################
!# ARGUMENTS 
!#
!# txc     Aux   3D array size 6
!#
!########################################################################
SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT(dte,etime, &
     x,y,z,dx,dy,dz, q,hq,s,hs, txc, vaux, wrk1d,wrk2d,wrk3d, &
     l_q, l_hq, l_txc, l_tags, l_comm)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field
  USE DNS_GLOBAL, ONLY : inb_flow,inb_scal, inb_scal_array
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : iadvection, idiffusion, ibodyforce, iviscous
  USE DNS_GLOBAL, ONLY : p_init, damkohler
  USE DNS_GLOBAL, ONLY : body_param
  USE DNS_GLOBAL, ONLY : icalc_particle, isize_particle, inb_particle
  USE THERMO_GLOBAL, ONLY : imixture
  USE LAGRANGE_GLOBAL, ONLY: particle_number, particle_bumper

  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : VA_BCS_HT, VA_BCS_HB, VA_BCS_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type
  USE DNS_LOCAL,  ONLY : imode_rhs

  IMPLICIT NONE

  TREAL dte,etime
  TREAL, DIMENSION(*)                 :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field, *)    :: q,hq, s,hs
  TREAL, DIMENSION(isize_txc_field,*) :: txc
  TREAL, DIMENSION(*)                 :: wrk1d,wrk2d,wrk3d, vaux

  TREAL, DIMENSION(isize_particle,*)  :: l_q, l_hq
  TREAL, DIMENSION(*)                 :: l_comm, l_txc
  INTEGER(8), DIMENSION(*)            :: l_tags

  TARGET :: q

! -----------------------------------------------------------------------
  TINTEGER ij, is
  TREAL dummy
  TINTEGER srt,end,siz    !  Variables for OpenMP Partitioning 

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w

#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT')
#endif

#ifdef USE_BLAS
  ilen = isize_field
#endif

! Define pointers
  u => q(:,1)
  v => q(:,2)
  w => q(:,3)

! #######################################################################
! Evaluate standard RHS of incompressible equations
! #######################################################################

! -----------------------------------------------------------------------
  IF      ( iadvection .EQ. EQNS_DIVERGENCE .AND. &
            iviscous   .EQ. EQNS_EXPLICIT   .AND. & 
            idiffusion .EQ. EQNS_EXPLICIT         ) THEN
     CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3(dte, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3),s, &
          q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
          vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), vaux, &
          wrk1d,wrk2d,wrk3d)
     
     DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3(is, dx,dy,dz, u,v,w,s(1,is),hs(1,is), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
     ENDDO
        
! -----------------------------------------------------------------------
  ELSE IF ( iadvection .EQ. EQNS_SKEWSYMMETRIC .AND. &
            iviscous   .EQ. EQNS_EXPLICIT      .AND. & 
            idiffusion .EQ. EQNS_EXPLICIT            ) THEN
     CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2(dte, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3),s, &
          q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
          vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), vaux, &
          wrk1d,wrk2d,wrk3d)
     
     DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2(is, dx,dy,dz, u,v,w,s(1,is), hs(1,is), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
     ENDDO
     
! -----------------------------------------------------------------------
  ELSE IF ( iadvection .EQ. EQNS_CONVECTIVE .AND. &
            iviscous   .EQ. EQNS_EXPLICIT   .AND. & 
            idiffusion .EQ. EQNS_EXPLICIT         ) THEN
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_SUPSAT ) THEN
        CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_ALWATER(dte,etime, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3),s, &
           txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
           vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), &
           wrk1d,wrk2d,wrk3d)

        IF ( imixture .EQ. MIXT_TYPE_SUPSAT ) THEN
           CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_SUPSAT(dte, dx,dy,dz, u,v,w, hq(1,1),hq(1,2),hq(1,3), s, hs, &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk1d,wrk2d,wrk3d)
        ELSE ! imixture .EQ. MIXT_TYPE_AIRWATER
           DO is = 1,inb_scal
              CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_ALWATER(is, dte, dx,dy,dz, u,v,w,s(1,1), hs(1,is), &
                 txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d,hq(1,1),hq(1,2),hq(1,3))
           ENDDO
        ENDIF

     ELSE
        IF      ( imode_rhs .EQ. EQNS_RHS_SPLIT       ) THEN 
           CALL RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1(dte,etime, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3),s, &
                q,hq, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
                vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), vaux, &
                wrk1d,wrk2d,wrk3d)
           CALL FI_SOURCES_SCAL(y,dy, s, hs, txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d)
           DO is = 1,inb_scal
              CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1(is, dte, dx,dy,dz, u,v,w,s(1,is),hs(1,is), s,&
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
           ENDDO

        ELSE IF ( imode_rhs .EQ. EQNS_RHS_COMBINED    ) THEN 
           CALL RHS_GLOBAL_INCOMPRESSIBLE_1(dte,etime, x,y,z,dx,dy,dz, u,v,w,hq(1,1),hq(1,2),hq(1,3), &
                q,hq, s,hs, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
                vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), vaux, &
                wrk1d,wrk2d,wrk3d)
        ELSE IF ( imode_rhs .EQ. EQNS_RHS_NONBLOCKING ) THEN 
#ifdef USE_PSFFT 
           CALL RHS_GLOBAL_INCOMPRESSIBLE_NBC(dte,etime, x,y,z,dx,dy,dz,&
                u(1),v(1),w(1),s(1,1),&
                txc(1,1), txc(1,2), &
                txc(1,3), txc(1,4), txc(1,5), txc(1,6), txc(1,7), txc(1,8),txc(1,9),txc(1,10), &
                txc(1,11),txc(1,12),txc(1,13),txc(1,14),&
                hq(1,1),hq(1,2),hq(1,3), hs(1,1), &
                vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux(vindex(VA_BCS_VI)), vaux, &
                wrk1d,wrk2d,wrk3d)
#else
           CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT. Need compiling flag -DUSE_PSFFT.')
           CALL DNS_STOP(DNS_ERROR_PSFFT)
#endif
        ENDIF
     ENDIF
     
! -----------------------------------------------------------------------
  ELSE
     CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT. Undeveloped formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     
  ENDIF


! #######################################################################
! Call RHS particle algorithm
! #######################################################################
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL RHS_PARTICLE_GLOBAL(x,y,z,dx,dy,dz,q,s,wrk1d,wrk2d,wrk3d,txc,l_q,l_hq, &
          l_tags, l_comm)
     
!    CALL FIELD_TO_PARTICLE &
!    (q(:,1), wrk1d, wrk2d, wrk3d,x ,y, z, l_txc, l_tags, l_hq, l_q)
     
  END IF

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     ! Flow part needs to be taken into account in the pressure
     DO is = 1,inb_scal
     CALL BOUNDARY_BUFFER_RELAXATION_SCAL(is,&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), x,y, q,s(1,is),hs(1,is))
     ENDDO
  ENDIF

! #######################################################################
! Perform the time stepping for incompressible equations
! #######################################################################
#ifdef USE_OPENMP
!$omp parallel default(shared) &
#ifdef USE_BLAS
!$omp private (ilen,is,srt,end,siz)
#else
!$omp private (ij,  is,srt,end,siz)
#endif 
#endif 

     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
#ifdef USE_BLAS 
     ilen = siz
#endif 
     DO is = 1,inb_flow

#ifdef USE_BLAS
        CALL DAXPY(ilen, dte, hq(srt,is), 1, q(srt,is), 1)
#else
        DO ij = srt,end
           q(ij,is) = q(ij,is) + dte*hq(ij,is)
        ENDDO
#endif
     ENDDO

     DO is = 1,inb_scal
#ifdef BLAS
        CALL DAXPY(ilen, dte, hs(srt,is), 1, s(srt,is), 1)
#else
        DO ij = srt,end 
           s(ij,is) = s(ij,is) + dte*hs(ij,is)
        ENDDO
#endif
     ENDDO
#ifdef USE_OPENMP
!$omp end parallel 
#endif 
        
! ######################################################################
! Particle POSTION UPDATED and  SEND/RECV TO THE NEW PROCESSOR
! ######################################################################
  IF ( icalc_particle .EQ. 1 ) THEN 
    CALL PARTICLE_TIME_SUBSTEP(dte, x, z, l_q, l_hq,l_tags, l_comm )    
  END IF 

! ###################################################################
! Calculate other intensive thermodynamic variables
! ###################################################################
  IF      ( (imixture .EQ. MIXT_TYPE_AIRWATER) .OR. &
            ( (imixture .EQ. MIXT_TYPE_SUPSAT) .AND. ( damkohler(1) .LE. C_0_R) ) ) THEN
     dummy = p_init/C_1_R ! MRATIO = 1
     CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(:,2), dummy, s(:,1))

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
     CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array), wrk3d)

  ENDIF

! ###################################################################
! Simulation control
! ###################################################################
  CALL DNS_CONTROL_SCAL(s)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_SUBSTEP_INCOMPRESSIBLE')
#endif

  RETURN
END SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_EXPLICIT

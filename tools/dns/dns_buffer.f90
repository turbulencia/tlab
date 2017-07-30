#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Note that with the energy formulation the forcing terms have been
!# made proportional to differences on the conserved variables,
!# i.e. rho, rho*V, rho*(e+v^2/2). 
!#
!########################################################################
SUBROUTINE BOUNDARY_BUFFER_INITIALIZE(q,s, txc, wrk3d)

  USE DNS_CONSTANTS, ONLY : tag_flow,tag_scal, lfile
  USE DNS_GLOBAL,    ONLY : imode_eqns
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, isize_field
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : itime
  USE DNS_GLOBAL,    ONLY : mach
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,     ONLY : BuffFlowJmin,BuffFlowJmax,BuffFlowImin,BuffFlowImax, buff_load
  USE DNS_LOCAL,     ONLY : BuffScalJmin,BuffScalJmax,BuffScalImin,BuffScalImax
#ifdef USE_MPI
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)    :: q, s
  TREAL, DIMENSION(imax*jmax*kmax,2), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk3d
  
! -------------------------------------------------------------------
  TINTEGER j,jloc, i,iloc, iq,is, idummy,io_sizes(5)
  TREAL dummy, var_max,var_min, prefactor

  CHARACTER*32 name, str, varname(inb_flow+inb_scal)
  CHARACTER*128 line 

! ###################################################################
! Prepare data for subarrays
  DO iq = 1,MAX(inb_flow,inb_scal)
     WRITE(varname(iq),*) iq; varname(iq) = TRIM(ADJUSTL(varname(iq)))
  ENDDO

#ifdef USE_MPI
! I/O routines not yet developed for this particular case
  IF ( ims_npro_i .GT. 1 .AND. ( BuffFlowImin%size .GT. 0 .OR. BuffFLowImax%size .GT. 0 ) ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_INIT. I/O routines undeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)     
  ENDIF
#endif

! ###################################################################
! Define profiles of relaxation times
  DO iq = 1,inb_flow
     
     IF ( BuffFlowJmin%size .GT. 0 ) THEN
        dummy = C_1_R /( g(2)%nodes(BuffFlowJmin%size) - g(2)%nodes(1) ) ! Inverse of segment length
        DO jloc = 1,BuffFlowJmin%size
           j = jloc
           BuffFlowJmin%tau(jloc,iq) = BuffFlowJmin%strength(iq) *( (g(2)%nodes(BuffFlowJmin%size)-g(2)%nodes(j)) *dummy ) **BuffFlowJmin%sigma(iq)
        ENDDO
     ENDIF
     
     IF ( BuffFlowJmax%size .GT. 0 ) THEN
        dummy = C_1_R /( g(2)%nodes(g(2)%size) - g(2)%nodes(g(2)%size-BuffFlowJmax%size+1) ) ! Inverse of segment length
        DO jloc = 1,BuffFlowJmax%size
           j = g(2)%size -BuffFlowJmax%size +jloc
           BuffFlowJmax%tau(jloc,iq) = BuffFlowJmax%strength(iq) *( (g(2)%nodes(j)-g(2)%nodes(g(2)%size-BuffFlowJmax%size+1)) *dummy ) **BuffFlowJmax%sigma(iq)
        ENDDO
     ENDIF

     IF ( BuffFlowImin%size .GT. 0 ) THEN
        dummy = C_1_R /( g(1)%nodes(BuffFlowImin%size) - g(1)%nodes(1) ) ! Inverse of segment length
        DO iloc = 1,BuffFlowImin%size
           i = iloc
           BuffFlowImin%tau(iloc,iq) = BuffFlowImin%strength(iq) *( (g(1)%nodes(BuffFlowImin%size)-g(1)%nodes(i)) *dummy ) **BuffFlowImin%sigma(iq)
        ENDDO
     ENDIF
     
     IF ( BuffFlowImax%size .GT. 0 ) THEN
        dummy = C_1_R /( g(1)%nodes(g(1)%size) - g(1)%nodes(g(1)%size-BuffFlowImax%size+1) ) ! Inverse of segment length
        DO iloc = 1,BuffFlowImax%size
           i = g(1)%size -BuffFlowImax%size +iloc
           BuffFlowImax%tau(iloc,iq) = BuffFlowImax%strength(iq) *( (g(1)%nodes(i)-g(1)%nodes(g(1)%size-BuffFlowImax%size+1)) *dummy ) **BuffFlowImax%sigma(iq)
        ENDDO
     ENDIF

  ENDDO
  
  DO is = 1,inb_scal
     
     IF ( BuffScalJmin%size .GT. 0 ) THEN
        dummy = C_1_R /( g(2)%nodes(BuffScalJmin%size) - g(2)%nodes(1) ) ! Inverse of segment length
        DO jloc = 1,BuffScalJmin%size
           j = jloc
           BuffScalJmin%tau(jloc,is) = BuffScalJmin%strength(is) *( (g(2)%nodes(BuffScalJmin%size)-g(2)%nodes(j)) *dummy ) **BuffScalJmin%sigma(is)
        ENDDO
     ENDIF
     
     IF ( BuffScalJmax%size .GT. 0 ) THEN
        dummy = C_1_R /( g(2)%nodes(g(2)%size) - g(2)%nodes(g(2)%size-BuffScalJmax%size+1) ) ! Inverse of segment length
        DO jloc = 1,BuffScalJmax%size
           j = g(2)%size -BuffScalJmax%size +jloc
           BuffScalJmax%tau(jloc,is) = BuffScalJmax%strength(is) *( (g(2)%nodes(j)-g(2)%nodes(g(2)%size-BuffScalJmax%size+1)) *dummy ) **BuffScalJmax%sigma(is)
        ENDDO
     ENDIF

     IF ( BuffScalImin%size .GT. 0 ) THEN
        dummy = C_1_R /( g(1)%nodes(BuffScalImin%size) - g(1)%nodes(1) ) ! Inverse of segment length
        DO iloc = 1,BuffScalImin%size
           i = iloc
           BuffScalImin%tau(iloc,is) = BuffScalImin%strength(is) *( (g(1)%nodes(BuffScalImin%size)-g(1)%nodes(i)) *dummy ) **BuffScalImin%sigma(is)
        ENDDO
     ENDIF
     
     IF ( BuffScalImax%size .GT. 0 ) THEN
        dummy = C_1_R /( g(1)%nodes(g(1)%size) - g(1)%nodes(g(1)%size-BuffScalImax%size+1) ) ! Inverse of segment length
        DO iloc = 1,BuffScalImax%size
           i = g(1)%size -BuffScalImax%size +iloc
           BuffScalImax%tau(iloc,is) = BuffScalImax%strength(is) *( (g(1)%nodes(i)-g(1)%nodes(g(1)%size-BuffScalImax%size+1)) *dummy ) **BuffScalImax%sigma(is)
        ENDDO
     ENDIF

  ENDDO
  
! ###################################################################
! Read buffer zones
  IF ( buff_load .EQ. 1 ) THEN

     IF ( BuffFlowJmin%size .GT. 0 ) THEN
        name = TRIM(ADJUSTL(tag_flow))//'bcs.jmin'
        idummy = imax*BuffFlowJmin%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_READ_SUBARRAY8(i4, name, varname, BuffFlowJmin%ref, io_sizes, wrk3d)
     ENDIF
     IF ( BuffScalJmin%size .GT. 0 ) THEN
        name = TRIM(ADJUSTL(tag_scal))//'bcs.jmin'
        idummy = imax*BuffScalJmin%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_scal/)
        CALL IO_READ_SUBARRAY8(i4, name, varname, BuffScalJmin%ref, io_sizes, wrk3d)
     ENDIF

     IF ( BuffFlowJmax%size .GT. 0 ) THEN 
        name = TRIM(ADJUSTL(tag_flow))//'bcs.jmax'
        idummy = imax*BuffFlowJmax%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_READ_SUBARRAY8(i5, name, varname, BuffFlowJmax%ref, io_sizes, wrk3d)
     ENDIF
     IF ( BuffFlowJmax%size .GT. 0 ) THEN
        name = TRIM(ADJUSTL(tag_scal))//'bcs.jmax'
        idummy = imax*BuffScalJmax%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_scal/)
        CALL IO_READ_SUBARRAY8(i5, name, varname, BuffScalJmax%ref, io_sizes, wrk3d)
     ENDIF

     IF ( BuffFlowImin%size .GT. 0 ) THEN 
        name = TRIM(ADJUSTL(tag_flow))//'bcs.imin'
        CALL DNS_READ_FIELDS(name, i0, BuffFlowImin%size,jmax,kmax, inb_flow, i0, isize_field, BuffFlowImin%ref, wrk3d)
     ENDIF
     IF ( BuffFlowImin%size .GT. 0 ) THEN
        name = TRIM(ADJUSTL(tag_scal))//'bcs.imin'
        CALL DNS_READ_FIELDS(name, i0, BuffScalImin%size,jmax,kmax, inb_scal, i0, isize_field, BuffScalImin%ref, wrk3d)
     ENDIF

     IF ( BuffFlowImax%size .GT. 0 ) THEN 
        name = TRIM(ADJUSTL(tag_flow))//'bcs.imax'
        CALL DNS_READ_FIELDS(name, i0, BuffFlowImax%size,jmax,kmax, inb_flow, i0, isize_field, BuffFlowImax%ref, wrk3d)
     ENDIF
     IF ( BuffFlowImax%size .GT. 0 ) THEN
        name = TRIM(ADJUSTL(tag_scal))//'bcs.imax'
        CALL DNS_READ_FIELDS(name, i0, BuffScalImax%size,jmax,kmax, inb_scal, i0, isize_field, BuffScalImax%ref, wrk3d)
     ENDIF

  ELSE
! ###################################################################
! Set and save buffer zones
     IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
        prefactor = C_05_R *(gama0-C_1_R)*mach*mach
        txc(:,2) = q(:,4) + prefactor *( q(:,1)*q(:,1) +q(:,2)*q(:,2) +q(:,3)*q(:,3) )
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     ELSE
        txc(:,1) = C_1_R
     ENDIF

     IF ( BuffFlowJmin%size .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_HB(q,s, txc, BuffFlowJmin%ref, BuffScalJmin%ref)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.jmin.'//TRIM(ADJUSTL(name))
        idummy = imax*BuffFlowJmin%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_WRITE_SUBARRAY8(i4, name, varname, BuffFlowJmin%ref, io_sizes, wrk3d)
     ENDIF
     IF ( BuffScalJmin%size .GT. 0 ) THEN 
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.jmin.'//TRIM(ADJUSTL(name))
        idummy = imax*BuffScalJmin%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_scal/)
        CALL IO_WRITE_SUBARRAY8(i4, name, varname, BuffScalJmin%ref, io_sizes, wrk3d)
     ENDIF

     IF ( BuffFlowJmax%size .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_HT(q,s, txc, BuffFlowJmax%ref, BuffScalJmax%ref)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.jmax.'//TRIM(ADJUSTL(name))
        idummy = imax*BuffFlowJmax%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_flow/)
        CALL IO_WRITE_SUBARRAY8(i5, name, varname, BuffFlowJmax%ref,                   io_sizes, wrk3d)
     ENDIF
     IF ( BuffScalJmax%size .GT. 0 ) THEN 
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.jmax.'//TRIM(ADJUSTL(name))
        idummy = imax*BuffScalJmax%size*kmax; io_sizes = (/idummy,1,idummy,1,inb_scal/)
        CALL IO_WRITE_SUBARRAY8(i5, name, varname, BuffScalJmax%ref, io_sizes, wrk3d)
     ENDIF

     IF ( BuffFlowImin%size .GT. 0 ) THEN 
        CALL BOUNDARY_INIT_VI(q,s, txc, BuffFlowImin%ref, BuffScalImin%ref)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.imin.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, BuffFlowImin%size,jmax,kmax, inb_flow, isize_field, BuffFlowImin%ref, wrk3d)
     ENDIF
     IF ( BuffScalImin%size .GT. 0 ) THEN 
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.imin.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, BuffScalImin%size,jmax,kmax, inb_scal, isize_field, BuffScalImin%ref, wrk3d)
     ENDIF

     IF ( BuffFlowImax%size .GT. 0 ) THEN
        CALL BOUNDARY_INIT_VO(q,s, txc, BuffFlowImax%ref, BuffScalImax%ref)

        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_flow))//'bcs.imax.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, BuffFlowImax%size,jmax,kmax, inb_flow, isize_field, BuffFlowImax%ref, wrk3d)
     ENDIF
     IF ( BuffScalImax%size .GT. 0 ) THEN
        WRITE(name, *) itime; name = TRIM(ADJUSTL(tag_scal))//'bcs.imax.'//TRIM(ADJUSTL(name))
        CALL DNS_WRITE_FIELDS(name, i0, BuffScalImax%size,jmax,kmax, inb_scal, isize_field, BuffScalImax%ref, wrk3d)
     ENDIF

  ENDIF

! ###################################################################
! Control
  IF ( BuffFlowJmax%size .GT. 0 ) THEN
     DO iq = 1, inb_flow
        var_max = MAXVAL(BuffFlowJmax%ref(:,:,:,iq)); var_min = MINVAL(BuffFlowJmax%ref(:,:,:,iq))
#ifdef USE_MPI
        CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        var_max = dummy
        CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        var_min = dummy
#endif
        WRITE(line,10) var_max
        WRITE(str, 10) var_min
        line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
        WRITE(str,*) iq
        line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at upper buffer are '//TRIM(ADJUSTL(line))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDDO
  ENDIF

  IF ( BuffScalJmax%size .GT. 0 ) THEN
     DO iq = 1, inb_scal
        var_max = MAXVAL(BuffScalJmax%ref(:,:,:,iq)); var_min = MINVAL(BuffScalJmax%ref(:,:,:,iq))
#ifdef USE_MPI
        CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        var_max = dummy
        CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        var_min = dummy
#endif
        WRITE(line,10) var_max
        WRITE(str, 10) var_min
        line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
        WRITE(str,*) iq
        line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at upper buffer are '//TRIM(ADJUSTL(line))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDDO
  ENDIF

  IF ( BuffFlowJmin%size .GT. 0 ) THEN 
     DO iq = 1, inb_flow
        var_max = MAXVAL(BuffFlowJmin%ref(:,:,:,iq)); var_min = MINVAL(BuffFlowJmin%ref(:,:,:,iq))
#ifdef USE_MPI
        CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        var_max = dummy
        CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        var_min = dummy
#endif
        WRITE(line,10) var_max
        WRITE(str, 10) var_min
        line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
        WRITE(str,*) iq
        line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at lower buffer are '//TRIM(ADJUSTL(line))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDDO
  ENDIF

  IF ( BuffScalJmin%size .GT. 0 ) THEN 
     DO iq = 1, inb_scal
        var_max = MAXVAL(BuffScalJmin%ref(:,:,:,iq)); var_min = MINVAL(BuffScalJmin%ref(:,:,:,iq))
#ifdef USE_MPI
        CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        var_max = dummy
        CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        var_min = dummy
#endif
        WRITE(line,10) var_max
        WRITE(str, 10) var_min
        line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
        WRITE(str,*) iq
        line = 'Bounds of field '//TRIM(ADJUSTL(str))//' at lower buffer are '//TRIM(ADJUSTL(line))
        CALL IO_WRITE_ASCII(lfile,line)
     ENDDO
  ENDIF

  RETURN
10 FORMAT(G_FORMAT_R)
END SUBROUTINE BOUNDARY_BUFFER_INITIALIZE

! ###################################################################
! ###################################################################
SUBROUTINE BOUNDARY_BUFFER_RELAXATION_FLOW(q, hq)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : imode_eqns, inb_flow
  USE DNS_GLOBAL, ONLY : mach
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,  ONLY : BuffFlowJmin, BuffFlowJmax, BuffFlowImin, BuffFlowImax

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax,inb_flow), INTENT(IN)  :: q
  TREAL, DIMENSION(imax,jmax,kmax,inb_flow), INTENT(OUT) :: hq

  TARGET q
  
! -------------------------------------------------------------------
  TINTEGER i,j, iloc,jloc, iq, iq_max
  TREAL prefactor

  TREAL, DIMENSION(:,:,:), POINTER :: rho
  
! ###################################################################
  SELECT CASE( imode_eqns )

  CASE ( DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL )
     prefactor = C_05_R *(gama0-C_1_R) *mach*mach
     rho => q(:,:,:,5)

     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN; iq_max = 3;       ! The energy needs special treatment
     ELSE;                                       iq_max = 4; ENDIF
        
     DO iq = 1,iq_max
           
        DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
           j = jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( rho(:,j,:) *q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
        ENDDO
        
        DO jloc = 1,BuffFlowJmax%size ! Top boundary
           j = jmax -BuffFlowJmax%size +jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( rho(:,j,:) *q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImin%size ! Inflow boundary
           i = iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( rho(i,:,:) *q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImax%size ! Outflow boundary
           i = imax -BuffFlowImax%size +iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( rho(i,:,:) *q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
        ENDDO

     ENDDO

     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        iq = 4

        DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
           j = jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( rho(:,j,:) *( q(:,j,:,iq) &
                +prefactor *( q(:,j,:,1)*q(:,j,:,1) +q(:,j,:,2)*q(:,j,:,2) +q(:,j,:,3)*q(:,j,:,3) ) ) -BuffFlowJmin%Ref(:,jloc,:,iq))
        ENDDO
        
        DO jloc = 1,BuffFlowJmax%size ! Top boundary
           j = jmax -BuffFlowJmax%size +jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( rho(:,j,:) *( q(:,j,:,iq) &
                +prefactor *( q(:,j,:,1)*q(:,j,:,1) +q(:,j,:,2)*q(:,j,:,2) +q(:,j,:,3)*q(:,j,:,3) ) ) -BuffFlowJmax%Ref(:,jloc,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImin%size ! Inflow boundary
           i = iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( rho(i,:,:) *( q(i,:,:,iq)&
                +prefactor *( q(i,:,:,1)*q(i,:,:,1) +q(i,:,:,2)*q(i,:,:,2) +q(i,:,:,3)*q(i,:,:,3) ) ) -BuffFlowImin%Ref(iloc,:,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImax%size ! Outflow boundary
           i = imax -BuffFlowImax%size +iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( rho(i,:,:) *( q(i,:,:,iq)&
                +prefactor *( q(i,:,:,1)*q(i,:,:,1) +q(i,:,:,2)*q(i,:,:,2) +q(i,:,:,3)*q(i,:,:,3) ) ) -BuffFlowImax%Ref(iloc,:,:,iq))
        ENDDO
     ENDIF
     
     iq = 5 ! Density
     DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
        j = jloc
        hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
     ENDDO
     
     DO jloc = 1,BuffFlowJmax%size ! Top boundary
        j = jmax -BuffFlowJmax%size +jloc
        hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
     ENDDO
     
     DO iloc = 1,BuffFlowImin%size ! Inflow boundary
        i = iloc
        hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
     ENDDO
     
     DO iloc = 1,BuffFlowImax%size ! Outflow boundary
        i = imax -BuffFlowImax%size +iloc
        hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
     ENDDO
     
  CASE ( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
     DO iq = 1,inb_flow
           
        DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
           j = jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
        ENDDO
        
        DO jloc = 1,BuffFlowJmax%size ! Top boundary
           j = jmax -BuffFlowJmax%size +jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImin%size ! Inflow boundary
           i = iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
        ENDDO
        
        DO iloc = 1,BuffFlowImax%size ! Outflow boundary
           i = imax -BuffFlowImax%size +iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
        ENDDO
        
     ENDDO
     
  END SELECT
     
  RETURN
END SUBROUTINE BOUNDARY_BUFFER_RELAXATION_FLOW

! ###################################################################
! ###################################################################
SUBROUTINE BOUNDARY_BUFFER_RELAXATION_SCAL(is, rho,s, hs)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : imode_eqns
  USE DNS_LOCAL,  ONLY : BuffScalJmin, BuffScalJmax, BuffScalImin, BuffScalImax

  IMPLICIT NONE

  TINTEGER,                         INTENT(IN)  :: is
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(IN)  :: rho, s
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT) :: hs

! -------------------------------------------------------------------
  TINTEGER i,j, iloc,jloc

! ###################################################################
  SELECT CASE( imode_eqns )

  CASE ( DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL )     
     DO jloc = 1,BuffScalJmin%size ! Bottom boundary
        j = jloc
        hs(:,j,:) = hs(:,j,:) -BuffScalJmin%Tau(jloc,is) *( rho(:,j,:) *s(:,j,:) -BuffScalJmin%Ref(:,jloc,:,is))
     ENDDO
     
     DO jloc = 1,BuffScalJmax%size ! Top boundary
        j = jmax -BuffScalJmax%size +jloc
        hs(:,j,:) = hs(:,j,:) -BuffScalJmax%Tau(jloc,is) *( rho(:,j,:) *s(:,j,:) -BuffScalJmax%Ref(:,jloc,:,is))
     ENDDO
     
     DO iloc = 1,BuffScalImin%size ! Inflow boundary
        i = iloc
        hs(i,:,:) = hs(i,:,:) -BuffScalImin%Tau(iloc,is) *( rho(i,:,:) *s(i,:,:) -BuffScalImin%Ref(iloc,:,:,is))
     ENDDO
     
     DO iloc = 1,BuffScalImax%size ! Outflow boundary
        i = imax -BuffScalImax%size +iloc
        hs(i,:,:) = hs(i,:,:) -BuffScalImax%Tau(iloc,is) *( rho(i,:,:) *s(i,:,:) -BuffScalImax%Ref(iloc,:,:,is))
     ENDDO
     
  CASE ( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
     DO jloc = 1,BuffScalJmin%size ! Bottom boundary
        j = jloc
        hs(:,j,:) = hs(:,j,:) -BuffScalJmin%Tau(jloc,is) *( s(:,j,:) -BuffScalJmin%Ref(:,jloc,:,is))
     ENDDO
     
     DO jloc = 1,BuffScalJmax%size ! Top boundary
        j = jmax -BuffScalJmax%size +jloc
        hs(:,j,:) = hs(:,j,:) -BuffScalJmax%Tau(jloc,is) *( s(:,j,:) -BuffScalJmax%Ref(:,jloc,:,is))
     ENDDO
     
     DO iloc = 1,BuffScalImin%size ! Inflow boundary
        i = iloc
        hs(i,:,:) = hs(i,:,:) -BuffScalImin%Tau(iloc,is) *( s(i,:,:) -BuffScalImin%Ref(iloc,:,:,is))
     ENDDO
     
     DO iloc = 1,BuffScalImax%size ! Outflow boundary
        i = imax -BuffScalImax%size +iloc
        hs(i,:,:) = hs(i,:,:) -BuffScalImax%Tau(iloc,is) *( s(i,:,:) -BuffScalImax%Ref(iloc,:,:,is))
     ENDDO
     
  END SELECT
     
  RETURN
END SUBROUTINE BOUNDARY_BUFFER_RELAXATION_SCAL


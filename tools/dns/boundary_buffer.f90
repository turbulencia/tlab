#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

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
  USE DNS_LOCAL,     ONLY : BuffFlowJmin,BuffFlowJmax,BuffFlowImin,BuffFlowImax, BuffLoad
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
  IF ( BuffLoad ) THEN

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

        IF ( BuffFlowJmin%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
              j = jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( rho(:,j,:) *q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowJmax%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmax%size ! Top boundary
              j = jmax -BuffFlowJmax%size +jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( rho(:,j,:) *q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImin%active(iq) ) THEN
           DO iloc = 1,BuffFlowImin%size ! Inflow boundary
              i = iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( rho(i,:,:) *q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImax%active(iq) ) THEN
           DO iloc = 1,BuffFlowImax%size ! Outflow boundary
              i = imax -BuffFlowImax%size +iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( rho(i,:,:) *q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
        
     ENDDO

     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        iq = 4

        IF ( BuffFlowJmin%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
              j = jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( rho(:,j,:) *( q(:,j,:,iq) &
                   +prefactor *( q(:,j,:,1)*q(:,j,:,1) +q(:,j,:,2)*q(:,j,:,2) +q(:,j,:,3)*q(:,j,:,3) ) ) -BuffFlowJmin%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowJmax%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmax%size ! Top boundary
              j = jmax -BuffFlowJmax%size +jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( rho(:,j,:) *( q(:,j,:,iq) &
                   +prefactor *( q(:,j,:,1)*q(:,j,:,1) +q(:,j,:,2)*q(:,j,:,2) +q(:,j,:,3)*q(:,j,:,3) ) ) -BuffFlowJmax%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImin%active(iq) ) THEN
           DO iloc = 1,BuffFlowImin%size ! Inflow boundary
              i = iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( rho(i,:,:) *( q(i,:,:,iq)&
                   +prefactor *( q(i,:,:,1)*q(i,:,:,1) +q(i,:,:,2)*q(i,:,:,2) +q(i,:,:,3)*q(i,:,:,3) ) ) -BuffFlowImin%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImax%active(iq) ) THEN
           DO iloc = 1,BuffFlowImax%size ! Outflow boundary
              i = imax -BuffFlowImax%size +iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( rho(i,:,:) *( q(i,:,:,iq)&
                   +prefactor *( q(i,:,:,1)*q(i,:,:,1) +q(i,:,:,2)*q(i,:,:,2) +q(i,:,:,3)*q(i,:,:,3) ) ) -BuffFlowImax%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
     ENDIF
     
     iq = 5 ! Density
     IF ( BuffFlowJmin%active(iq) ) THEN
        DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
           j = jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
        ENDDO
     ENDIF
     
     IF ( BuffFlowJmax%active(iq) ) THEN
        DO jloc = 1,BuffFlowJmax%size ! Top boundary
           j = jmax -BuffFlowJmax%size +jloc
           hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
        ENDDO
     ENDIF
     
     IF ( BuffFlowImin%active(iq) ) THEN
        DO iloc = 1,BuffFlowImin%size ! Inflow boundary
           i = iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
        ENDDO
     ENDIF
     
     IF ( BuffFlowImax%active(iq) ) THEN
        DO iloc = 1,BuffFlowImax%size ! Outflow boundary
           i = imax -BuffFlowImax%size +iloc
           hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
        ENDDO
     ENDIF
     
  CASE ( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
     DO iq = 1,inb_flow
        
        IF ( BuffFlowJmin%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmin%size ! Bottom boundary
              j = jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmin%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmin%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowJmax%active(iq) ) THEN
           DO jloc = 1,BuffFlowJmax%size ! Top boundary
              j = jmax -BuffFlowJmax%size +jloc
              hq(:,j,:,iq) = hq(:,j,:,iq) -BuffFlowJmax%Tau(jloc,iq) *( q(:,j,:,iq) -BuffFlowJmax%Ref(:,jloc,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImin%active(iq) ) THEN
           DO iloc = 1,BuffFlowImin%size ! Inflow boundary
              i = iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImin%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImin%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
        
        IF ( BuffFlowImax%active(iq) ) THEN
           DO iloc = 1,BuffFlowImax%size ! Outflow boundary
              i = imax -BuffFlowImax%size +iloc
              hq(i,:,:,iq) = hq(i,:,:,iq) -BuffFlowImax%Tau(iloc,iq) *( q(i,:,:,iq) -BuffFlowImax%Ref(iloc,:,:,iq))
           ENDDO
        ENDIF
        
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
     IF ( BuffScalJmin%active(is) ) THEN
        DO jloc = 1,BuffScalJmin%size ! Bottom boundary
           j = jloc
           hs(:,j,:) = hs(:,j,:) -BuffScalJmin%Tau(jloc,is) *( rho(:,j,:) *s(:,j,:) -BuffScalJmin%Ref(:,jloc,:,is))
        ENDDO
     ENDIF
     
     IF ( BuffScalJmax%active(is) ) THEN
        DO jloc = 1,BuffScalJmax%size ! Top boundary
           j = jmax -BuffScalJmax%size +jloc
           hs(:,j,:) = hs(:,j,:) -BuffScalJmax%Tau(jloc,is) *( rho(:,j,:) *s(:,j,:) -BuffScalJmax%Ref(:,jloc,:,is))
        ENDDO
     ENDIF
     
     IF ( BuffScalImin%active(is) ) THEN
        DO iloc = 1,BuffScalImin%size ! Inflow boundary
           i = iloc
           hs(i,:,:) = hs(i,:,:) -BuffScalImin%Tau(iloc,is) *( rho(i,:,:) *s(i,:,:) -BuffScalImin%Ref(iloc,:,:,is))
        ENDDO
     ENDIF
     
     IF ( BuffScalImax%active(is) ) THEN
        DO iloc = 1,BuffScalImax%size ! Outflow boundary
           i = imax -BuffScalImax%size +iloc
           hs(i,:,:) = hs(i,:,:) -BuffScalImax%Tau(iloc,is) *( rho(i,:,:) *s(i,:,:) -BuffScalImax%Ref(iloc,:,:,is))
        ENDDO
     ENDIF
     
  CASE ( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC )
     IF ( BuffScalJmin%active(is) ) THEN
        DO jloc = 1,BuffScalJmin%size ! Bottom boundary
           j = jloc
           hs(:,j,:) = hs(:,j,:) -BuffScalJmin%Tau(jloc,is) *( s(:,j,:) -BuffScalJmin%Ref(:,jloc,:,is))
        ENDDO
     ENDIF
     
     IF ( BuffScalJmax%active(is) ) THEN
        DO jloc = 1,BuffScalJmax%size ! Top boundary
           j = jmax -BuffScalJmax%size +jloc
           hs(:,j,:) = hs(:,j,:) -BuffScalJmax%Tau(jloc,is) *( s(:,j,:) -BuffScalJmax%Ref(:,jloc,:,is))
        ENDDO
     ENDIF
     
     IF ( BuffScalImin%active(is) ) THEN
        DO iloc = 1,BuffScalImin%size ! Inflow boundary
           i = iloc
           hs(i,:,:) = hs(i,:,:) -BuffScalImin%Tau(iloc,is) *( s(i,:,:) -BuffScalImin%Ref(iloc,:,:,is))
        ENDDO
     ENDIF
  
     IF ( BuffScalImax%active(is) ) THEN
        DO iloc = 1,BuffScalImax%size ! Outflow boundary
           i = imax -BuffScalImax%size +iloc
           hs(i,:,:) = hs(i,:,:) -BuffScalImax%Tau(iloc,is) *( s(i,:,:) -BuffScalImax%Ref(iloc,:,:,is))
        ENDDO
     ENDIF
     
  END SELECT
     
  RETURN
END SUBROUTINE BOUNDARY_BUFFER_RELAXATION_SCAL

!########################################################################
!########################################################################
SUBROUTINE BOUNDARY_INIT_HT(q,s, txc, buffer_q, buffer_s)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : imode_sim, imode_eqns
  USE DNS_GLOBAL,    ONLY : inb_scal, area
  USE DNS_LOCAL,     ONLY : BuffFlowJmax, BuffScalJmax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)  :: q, s, txc
  TREAL, DIMENSION(imax,BuffFlowJmax%size,kmax,*) :: buffer_q
  TREAL, DIMENSION(imax,BuffScalJmax%size,kmax,*) :: buffer_s

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG_IK, AVG1V1D, COV2V2D, COV2V1D
  TINTEGER i, j, is, iq, jloc

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Temporally evolving shear layer
! ###################################################################
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

! Flow variables
     DO iq = 1,3
        DO j = 1, BuffFlowJmax%size
           jloc = jmax -BuffFlowJmax%size +j
           IF ( .NOT. BuffFlowJmax%hard ) &
                BuffFlowJmax%hardvalues(iq) = COV2V2D(imax,jmax,kmax, jloc, r_loc,q(1,iq))
           buffer_q(:,j,:,iq) = BuffFlowJmax%hardvalues(iq)
        ENDDO
     ENDDO

! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1, BuffFlowJmax%size
           jloc = jmax -BuffFlowJmax%size +j
           IF ( .NOT. BuffFlowJmax%hard ) &
                BuffFlowJmax%hardvalues(4) = COV2V2D(imax,jmax,kmax, jloc, r_loc,e_loc)
           buffer_q(:,j,:,4) = BuffFlowJmax%hardvalues(4)
           IF ( .NOT. BuffFlowJmax%hard ) &
                BuffFlowJmax%hardvalues(5) = AVG_IK(imax,jmax,kmax, jloc, r_loc, g(1)%jac,g(3)%jac, area)
           buffer_q(:,j,:,5) = BuffFlowJmax%hardvalues(5)
        ENDDO
     ENDIF
     
! Scalars
     DO is = 1,inb_scal
        DO j = 1,BuffScalJmax%size
           jloc = jmax -BuffScalJmax%size +j
           IF ( .NOT. BuffScalJmax%hard ) &
                BuffScalJmax%hardvalues(is) = COV2V2D(imax,jmax,kmax, jloc, r_loc,s(1,is))
           buffer_s(:,j,:,is) = BuffScalJmax%hardvalues(is)
        ENDDO
     ENDDO
     
! ###################################################################
! Spatially evolving jet
! ###################################################################
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN  

! Flow variables
     DO iq = 1,3
        DO j = 1, BuffFlowJmax%size
           jloc = jmax -BuffFlowJmax%size +j
           DO i = 1,imax
              IF ( .NOT. BuffFlowJmax%hard ) &
                   BuffFlowJmax%hardvalues(iq) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,q(1,iq))
              buffer_q(i,j,:,iq) = BuffFlowJmax%hardvalues(iq)
           ENDDO
        ENDDO
     ENDDO
        
! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1, BuffFlowJmax%size
           jloc = jmax -BuffFlowJmax%size +j
           DO i = 1,imax
              IF ( .NOT. BuffFlowJmax%hard ) &
                   BuffFlowJmax%hardvalues(4) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,e_loc)
              buffer_q(i,j,:,4) = BuffFlowJmax%hardvalues(4)
              IF ( .NOT. BuffFlowJmax%hard ) &
                   BuffFlowJmax%hardvalues(5) = AVG1V1D(imax,jmax,kmax, i,jloc, i1, r_loc)
              buffer_q(i,j,:,5) = BuffFlowJmax%hardvalues(5)
           ENDDO
        ENDDO
     ENDIF
        
! Scalars
     DO is = 1,inb_scal
        DO j = 1,BuffScalJmax%size
           jloc = jmax -BuffScalJmax%size +j
           DO i = 1,imax
              IF ( .NOT. BuffScalJmax%hard ) &
                   BuffScalJmax%hardvalues(is) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,s(1,is))
              buffer_s(i,j,:,is) = BuffScalJmax%hardvalues(is)
           ENDDO
        ENDDO
     ENDDO
  
  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_INIT_HT

!########################################################################
!########################################################################
SUBROUTINE BOUNDARY_INIT_HB(q,s, txc, buffer_q, buffer_s)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : imode_sim, imode_eqns
  USE DNS_GLOBAL,    ONLY : inb_scal, area
  USE DNS_LOCAL,     ONLY : BuffFlowJmin, BuffScalJmin

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)  :: q, s, txc
  TREAL, DIMENSION(imax,BuffFlowJmin%size,kmax,*) :: buffer_q
  TREAL, DIMENSION(imax,BuffScalJmin%size,kmax,*) :: buffer_s

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG_IK, AVG1V1D, COV2V2D, COV2V1D
  TINTEGER i, j, is, iq

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Temporally evolving shear layer
! ###################################################################
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

! Flow variables
     DO iq = 1,3
        DO j = 1, BuffFlowJmin%size
           IF ( .NOT. BuffFlowJmin%hard ) &
                BuffFlowJmin%hardvalues(iq) = COV2V2D(imax,jmax,kmax, j, r_loc,q(1,iq))
           buffer_q(:,j,:,iq) = BuffFlowJmin%hardvalues(iq)
        ENDDO
     ENDDO
     
! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1, BuffFlowJmin%size
           IF ( .NOT. BuffFlowJmin%hard ) &
                BuffFlowJmin%hardvalues(4) = COV2V2D(imax,jmax,kmax, j, r_loc,e_loc)
           buffer_q(:,j,:,4) = BuffFlowJmin%hardvalues(4)
           IF ( .NOT. BuffFlowJmin%hard ) &
                BuffFlowJmin%hardvalues(5) = AVG_IK(imax,jmax,kmax, j, r_loc, g(1)%jac,g(3)%jac, area)
           buffer_q(:,j,:,5) = BuffFlowJmin%hardvalues(5)
        ENDDO
     ENDIF

! Scalars
     DO is = 1,inb_scal
        DO j = 1,BuffScalJmin%size
           IF ( .NOT. BuffScalJmin%hard ) &
                BuffScalJmin%hardvalues(is) = COV2V2D(imax,jmax,kmax, j, r_loc,s(1,is))
           buffer_s(:,j,:,is) = BuffScalJmin%hardvalues(is)
        ENDDO
     ENDDO

! ###################################################################
! Spatially evolving jet
! ###################################################################
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN  

! Flow variables
     DO iq = 1,3
        DO j = 1, BuffFlowJmin%size
           DO i = 1,imax
              IF ( .NOT. BuffFlowJmin%hard ) &
                BuffFlowJmin%hardvalues(iq) = COV2V1D(imax,jmax,kmax, i,j, r_loc,q(1,iq))
              buffer_q(i,j,:,iq) = BuffFlowJmin%hardvalues(iq)
           ENDDO
        ENDDO
     ENDDO
     
! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        DO j = 1, BuffFlowJmin%size
           DO i = 1,imax
              IF ( .NOT. BuffFlowJmin%hard ) &
                   BuffFlowJmin%hardvalues(4) = COV2V1D(imax,jmax,kmax, i,j, r_loc,e_loc)
              buffer_q(i,j,:,4) = BuffFlowJmin%hardvalues(4)
              IF ( .NOT. BuffFlowJmin%hard ) &
                   BuffFlowJmin%hardvalues(5) = AVG1V1D(imax,jmax,kmax, i,j, i1, r_loc)
              buffer_q(i,j,:,5) = BuffFlowJmin%hardvalues(5)
           ENDDO
        ENDDO
     ENDIF

! Scalars
     DO is = 1,inb_scal
        DO j = 1,BuffScalJmin%size
           DO i = 1,imax
              IF ( .NOT. BuffScalJmin%hard ) &
                   BuffScalJmin%hardvalues(is) = COV2V1D(imax,jmax,kmax, i,j, r_loc,s(1,is))
              buffer_s(i,j,:,is) = BuffScalJmin%hardvalues(is)
           ENDDO
        ENDDO
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_INIT_HB

!########################################################################
!########################################################################
!It only makes sense in spatially evolving cases.
SUBROUTINE BOUNDARY_INIT_VI(q,s, txc, buffer_q, buffer_s)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : imode_eqns, inb_scal, inb_vars
  USE DNS_LOCAL,  ONLY : BuffFlowImin, BuffScalImin

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)  :: q, s, txc
  TREAL, DIMENSION(BuffFlowImin%size,jmax,kmax,*) :: buffer_q
  TREAL, DIMENSION(BuffScalImin%size,jmax,kmax,*) :: buffer_s

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG1V1D, COV2V1D
  TREAL dbuff(inb_vars)
  TINTEGER i, j, iq, is

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Shear layer and jet profile
! ###################################################################
  DO iq = 1,3
     DO j = 1,jmax
        DO i = 1,BuffFlowImin%size
           dbuff(iq) = COV2V1D(imax,jmax,kmax, i,j, r_loc,q(1,iq))
           buffer_q(i,j,:,iq) = dbuff(iq)
        ENDDO
     ENDDO
  ENDDO

! if compressible
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1,jmax
        DO i = 1,BuffFlowImin%size
           dbuff(4) = COV2V1D(imax,jmax,kmax, i,j, r_loc,e_loc)
           buffer_q(i,j,:,4) = dbuff(4)
           dbuff(5) = AVG1V1D(imax,jmax,kmax, i,j, i1, r_loc)
           buffer_q(i,j,:,5) = dbuff(5)
        ENDDO
     ENDDO
  ENDIF
     
  DO is = 1,inb_scal
     DO j = 1,jmax
        DO i = 1,BuffScalImin%size
           dbuff(is) = COV2V1D(imax,jmax,kmax, i,j, r_loc,s(1,is))
           buffer_s(i,j,:,is) = dbuff(is)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE BOUNDARY_INIT_VI

!########################################################################
!########################################################################
!It only makes sense in spatially evolving cases.
SUBROUTINE BOUNDARY_INIT_VO(q,s, txc, buffer_q, buffer_s)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : imode_eqns, inb_scal, inb_vars
  USE DNS_LOCAL,  ONLY : BuffFlowImax, BuffScalImax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)  :: q, s, txc
  TREAL, DIMENSION(BuffFlowImax%size,jmax,kmax,*) :: buffer_q
  TREAL, DIMENSION(BuffScalImax%size,jmax,kmax,*) :: buffer_s

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG1V1D, COV2V1D
  TREAL dbuff(inb_vars)
  TINTEGER i, j, iq, is, iloc

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Shear layer and jet profile
! ###################################################################
  DO iq = 1,3
     DO j = 1,jmax
        DO i = 1,BuffFlowImax%size
           iloc = imax -BuffFlowImax%size +i
           dbuff(iq) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,q(1,iq))
           buffer_q(i,j,:,iq) = dbuff(iq)
        ENDDO
     ENDDO
  ENDDO
  
! if compressible
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1,jmax
        DO i = 1,BuffFlowImax%size
           iloc = imax -BuffFlowImax%size +i
           dbuff(4) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,e_loc)
           buffer_q(i,j,:,4) = dbuff(4)
           dbuff(5) = AVG1V1D(imax,jmax,kmax, iloc,j, i1, r_loc)
           buffer_q(i,j,:,5) = dbuff(5)
        ENDDO
     ENDDO
  ENDIF
  
  DO is = 1,inb_scal
     DO j = 1,jmax
        DO i = 1,BuffScalImax%size
           iloc = imax -BuffScalImax%size +i
           dbuff(is) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,s(1,is))
           buffer_s(i,j,:,is) = dbuff(is)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE BOUNDARY_INIT_VO

!########################################################################
!########################################################################
SUBROUTINE BOUNDARY_BUFFER_FILTER(rho,u,v,w,e,z1, txc1,txc2,txc3,txc4,txc5, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : icalc_scal, inb_scal
  USE DNS_LOCAL,  ONLY : BuffFlowImax, BuffFlowImax, BuffFlowJmin, BuffFlowJmax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax  )        :: rho, u, v, w, e
  TREAL, DIMENSION(imax,jmax,kmax,*)        :: z1
  TREAL, DIMENSION(BuffFlowImax%size,jmax,kmax) :: txc1, txc2, txc3, txc4, txc5 
  TREAL, DIMENSION(*)                       :: wrk1d,wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER id, iloc, i, j, k, is, ibc_x(4), ibc_y(4), ibc_z(4), buff_imax
  TREAL eta, delta, amp, ampr, rho_ratio

! ###################################################################
!!! Routines OPR_FILTER have been changed. This routine needs to be updates !!!
  CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
  CALL DNS_STOP(DNS_ERROR_UNDEVELOP)


! BCs for the filters (see routine FILTER)
  ! ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
  ! ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ! ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

! ###################################################################
! Bottom boundary
! ###################################################################
  IF ( BuffFlowJmin%size .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! ###################################################################
! Top boundary
! ###################################################################
  IF ( BuffFlowJmax%size .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! ###################################################################
! Outflow boundary
! ###################################################################
  IF ( BuffFlowImax%size .GT. 1 ) THEN
     id = DNS_MPI_K_OUTBCS
     buff_imax = imax -BuffFlowImax%size +iloc
! -------------------------------------------------------------------
! Flow
! -------------------------------------------------------------------
     DO k = 1,kmax
        DO j = 1,jmax
           DO iloc = 1,BuffFlowImax%size ! Outflow boundary
              i = imax -BuffFlowImax%size +iloc
              txc1(iloc,j,k) = rho(i,j,k)
              txc2(iloc,j,k) = rho(i,j,k)*u(i,j,k)
              txc3(iloc,j,k) = rho(i,j,k)*v(i,j,k)
              txc4(iloc,j,k) = rho(i,j,k)*w(i,j,k)
              txc5(iloc,j,k) = rho(i,j,k)*e(i,j,k)
!              txc2(iloc,j,k) = u(i,j,k)
!              txc3(iloc,j,k) = v(i,j,k)
!              txc4(iloc,j,k) = w(i,j,k)
!              txc5(iloc,j,k) = e(i,j,k)
           ENDDO
        ENDDO
     ENDDO

     CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc1, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc3, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc4, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc5, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
     delta = (g(1)%nodes(imax)-g(1)%nodes(buff_imax))/C_16_R
     DO i = buff_imax,imax
        iloc = i-buff_imax+1

        eta = g(1)%nodes(i) - C_05_R*( g(1)%nodes(imax)+g(1)%nodes(buff_imax) )
        amp = C_05_R*(C_1_R+TANH(C_05_R*eta/delta))
        ampr= C_1_R-amp

        DO k = 1,kmax
           DO j = 1,jmax
              rho_ratio = rho(i,j,k)
              rho(i,j,k) = ampr*rho(i,j,k) + amp*txc1(iloc,j,k)
              rho_ratio = rho_ratio/rho(i,j,k)

              u(i,j,k) = ampr*rho_ratio*u(i,j,k) + amp*txc2(iloc,j,k)/rho(i,j,k)
              v(i,j,k) = ampr*rho_ratio*v(i,j,k) + amp*txc3(iloc,j,k)/rho(i,j,k)
              w(i,j,k) = ampr*rho_ratio*w(i,j,k) + amp*txc4(iloc,j,k)/rho(i,j,k)
              e(i,j,k) = ampr*rho_ratio*e(i,j,k) + amp*txc5(iloc,j,k)/rho(i,j,k)

! store rho_ratio to be used in scalar section
              txc1(iloc,j,k) = rho_ratio

!              rho(i,j,k) = ampr*rho(i,j,k) + amp*txc1(iloc,j,k)
!              u(i,j,k)   = ampr*u(i,j,k)   + amp*txc2(iloc,j,k)
!              v(i,j,k)   = ampr*v(i,j,k)   + amp*txc3(iloc,j,k)
!              w(i,j,k)   = ampr*w(i,j,k)   + amp*txc4(iloc,j,k)
!              e(i,j,k)   = ampr*e(i,j,k)   + amp*txc5(iloc,j,k)

           ENDDO
        ENDDO
     ENDDO

! -------------------------------------------------------------------
! Scalar
! -------------------------------------------------------------------
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal

           DO k = 1,kmax
              DO j = 1,jmax
                 DO iloc = 1,BuffFlowImax%size ! Outflow boundary
                    i = imax -BuffFlowImax%size +iloc
                    txc2(iloc,j,k) = rho(i,j,k)*z1(i,j,k,is)
!                    txc2(iloc,j,k) = z1(i,j,k,is)
                 ENDDO
              ENDDO
           ENDDO

           CALL OPR_FILTER(i2, BuffFlowImax%size,jmax,kmax, ibc_x,ibc_y,ibc_z, id, &
                txc2, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)

! thickness \delta_\theta s.t. 2\delta_w = L_buffer/2
           delta = (g(1)%nodes(imax)-g(1)%nodes(buff_imax))/C_16_R
           DO i = buff_imax,imax
              iloc = i-buff_imax+1

              eta = g(1)%nodes(i) - C_05_R*( g(1)%nodes(imax)+g(1)%nodes(buff_imax) )
              amp = C_05_R*(C_1_R+TANH(C_05_R*eta/delta))
              ampr= C_1_R-amp

              DO k = 1,kmax
                 DO j = 1,jmax
                    z1(i,j,k,is) = ampr*txc1(iloc,j,k)*z1(i,j,k,is) + &
                         amp*txc2(iloc,j,k)/rho(i,j,k)
!                    z1(i,j,k,is) = ampr*z1(i,j,k,is) + amp*txc2(iloc,j,k)
                 ENDDO
              ENDDO

           ENDDO

        ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BUFFER_FILTER

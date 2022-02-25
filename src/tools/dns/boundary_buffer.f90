#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# Implementation of relaxation terms in buffer regions of the computational domain
!#
!# Data is read in dns_read_local
!#
!# Note that with the energy formulation the forcing terms have been
!# made proportional to differences on the conserved variables, i.e. rho, rho*V, rho*(e+v^2/2).
!#
!# The buffer files need to written without header information, so that
!# the header global variables (like itime) are not updated within this routine
!#
!########################################################################
MODULE BOUNDARY_BUFFER

  USE TLAB_TYPES,     ONLY : filter_dt

  USE TLAB_CONSTANTS, ONLY : tag_flow,tag_scal, wfile,lfile,efile, MAX_VARS
  USE TLAB_VARS,    ONLY : imode_eqns, imode_sim
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, isize_field
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : itime
  USE TLAB_VARS,    ONLY : mach
  USE TLAB_VARS,    ONLY : io_aux
  USE TLAB_PROCS
  USE THERMO_VARS, ONLY : gama0
  USE IO_FIELDS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_VARS, ONLY : ims_comm_z
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k, ims_pro_i, ims_pro_k
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: FORM_POWER_MIN = 1
  TINTEGER, PARAMETER :: FORM_POWER_MAX = 2

  TYPE buffer_dt
    SEQUENCE
    TINTEGER TYPE                                  ! relaxation, filter...
    TINTEGER size                                  ! # points in buffer layer
    TINTEGER offset                                ! position in absolute grid
    TINTEGER nfields                               ! number of fields to apply tosition in absolute grid
    LOGICAL active(MAX_VARS), hard
    TREAL strength(MAX_VARS), sigma(MAX_VARS)      ! parameters in relaxation term
    TINTEGER form                                  ! form of function of relaxation term
    TREAL hardvalues(MAX_VARS)                     ! Fixed reference values
    TREAL, ALLOCATABLE, DIMENSION(:,:)     :: tau  ! relaxation timescale for each field
    TREAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ref  ! reference fields
    TINTEGER id_subarray
  END TYPE buffer_dt

  TINTEGER,        PUBLIC :: BuffType
  LOGICAL,         PUBLIC :: BuffLoad
  TYPE(buffer_dt), PUBLIC :: BuffFlowImin,BuffFlowImax,BuffFlowJmin,BuffFlowJmax
  TYPE(buffer_dt), PUBLIC :: BuffScalImin,BuffScalImax,BuffScalJmin,BuffScalJmax
  !  TYPE(filter_dt), DIMENSION(3) :: FilterBuffer
  ! BufferFilter should then be a block in dns.ini as [Filter], which is read in io_read_global.

  PUBLIC :: BOUNDARY_BUFFER_READBLOCK
  PUBLIC :: BOUNDARY_BUFFER_INITIALIZE
  PUBLIC :: BOUNDARY_BUFFER_RELAX_FLOW
  PUBLIC :: BOUNDARY_BUFFER_RELAX_SCAL
  PUBLIC :: BOUNDARY_BUFFER_RELAX_SCAL_I
  PUBLIC :: BOUNDARY_BUFFER_FILTER

  PRIVATE

  TINTEGER j,jloc, i,iloc, iq,is, idummy
  TREAL dummy

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE BOUNDARY_BUFFER_READBLOCK(bakfile,inifile,tag,variable,nfields)
    IMPLICIT NONE

    CHARACTER(LEN=*) bakfile,inifile,tag
    TYPE(buffer_dt) variable
    TINTEGER nfields

    TREAL dummies(nfields+1)
    CHARACTER(LEN=512) sRes
    CHARACTER(LEN=4) str

    ! ###################################################################
    variable%active(:) = .FALSE.; variable%hard = .FALSE.
    IF ( variable%size .GT. 0 ) THEN
      CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'Parameters'//TRIM(ADJUSTL(tag)), 'void', sRes)
      IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN
        str = TRIM(ADJUSTL(tag))
        CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'Parameters'//str(1:1), '1.0,2.0', sRes)
      ENDIF
      is = nfields+1; CALL LIST_REAL(sRes, is, dummies)
      IF      ( is .EQ. 1 ) THEN
        variable%strength(:) = dummies(1)
        variable%sigma(:) = C_2_R
      ELSE IF ( is .EQ. 2 ) THEN
        variable%strength(:) = dummies(1)
        variable%sigma(:) = dummies(2)
      ELSE IF ( is .EQ. nfields+1 ) THEN
        variable%strength(1:nfields) = dummies(1:nfields)
        variable%sigma(:) = dummies(nfields+1)
      ELSE
        CALL TLAB_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Wrong number of values in BufferZone.ParametersUImin.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
      ENDIF
      DO is = 1,nfields
        IF ( variable%strength(is) .NE. C_0_R ) variable%active(is) = .TRUE.
      ENDDO

      CALL SCANINICHAR(bakfile, inifile, 'BufferZone', 'HardValues'//TRIM(ADJUSTL(tag)), 'void', sRes)
      IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
        is = nfields; CALL LIST_REAL(sRes, is, variable%hardvalues)
        IF ( is .EQ. nfields ) THEN
          variable%hard = .TRUE.
        ELSE
          CALL TLAB_WRITE_ASCII(wfile, 'DNS_READ_LOCAL. Wrong number of values in BufferZone.HardValues.'//TRIM(ADJUSTL(tag))//'.')
          CALL TLAB_STOP(DNS_ERROR_OPTION)
        ENDIF
      ENDIF

    ENDIF

  END SUBROUTINE BOUNDARY_BUFFER_READBLOCK

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE BOUNDARY_BUFFER_INITIALIZE(q,s, txc, wrk3d)
    IMPLICIT NONE

    TREAL, DIMENSION(isize_field,*), INTENT(IN   ) :: q, s
    TREAL, DIMENSION(isize_field,6), INTENT(INOUT) :: txc
    TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk3d

    ! ###################################################################
    BuffFlowImin%offset = 0; BuffFlowImax%offset = g(1)%size -BuffFlowImax%size
    BuffFlowJmin%offset = 0; BuffFlowJmax%offset = g(2)%size -BuffFlowJmax%size
    BuffScalImin%offset = 0; BuffScalImax%offset = g(1)%size -BuffScalImax%size
    BuffScalJmin%offset = 0; BuffScalJmax%offset = g(2)%size -BuffScalJmax%size

    BuffFlowImin%nfields = inb_flow; BuffFlowImax%nfields = inb_flow
    BuffFlowJmin%nfields = inb_flow; BuffFlowJmax%nfields = inb_flow
    BuffScalImin%nfields = inb_scal; BuffScalImax%nfields = inb_scal
    BuffScalJmin%nfields = inb_scal; BuffScalJmax%nfields = inb_scal

    BuffFlowImin%form = FORM_POWER_MIN; BuffFlowImax%form = FORM_POWER_MAX
    BuffFlowJmin%form = FORM_POWER_MIN; BuffFlowJmax%form = FORM_POWER_MAX
    BuffScalImin%form = FORM_POWER_MIN; BuffScalImax%form = FORM_POWER_MAX
    BuffScalJmin%form = FORM_POWER_MIN; BuffScalJmax%form = FORM_POWER_MAX

    ! ###################################################################
    txc(:,2:inb_flow+1) = q(:,1:inb_flow)
    IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
      txc(:,1) = q(:,5) ! Density
      dummy = C_05_R *(gama0-C_1_R)*mach*mach
      txc(:,5) = q(:,4) + dummy *( q(:,1)*q(:,1) +q(:,2)*q(:,2) +q(:,3)*q(:,3) )
      txc(:,6) = C_1_R
    ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
      txc(:,1) = q(:,5) ! Density
      txc(:,6) = C_1_R
    ELSE
      txc(:,1) = C_1_R ! Density
    ENDIF

    IF ( BuffFlowImin%size > 0 ) CALL INI_BLOCK(1, TRIM(ADJUSTL(tag_flow))//'bcs.imin', BuffFlowImin, txc(1,1), txc(1,2), wrk3d )
    IF ( BuffFlowImax%size > 0 ) CALL INI_BLOCK(1, TRIM(ADJUSTL(tag_flow))//'bcs.imax', BuffFlowImax, txc(1,1), txc(1,2), wrk3d )
    IF ( BuffFlowJmin%size > 0 ) CALL INI_BLOCK(2, TRIM(ADJUSTL(tag_flow))//'bcs.jmin', BuffFlowJmin, txc(1,1), txc(1,2), wrk3d )
    IF ( BuffFlowJmax%size > 0 ) CALL INI_BLOCK(2, TRIM(ADJUSTL(tag_flow))//'bcs.jmax', BuffFlowJmax, txc(1,1), txc(1,2), wrk3d )

    IF ( BuffScalImin%size > 0 ) CALL INI_BLOCK(1, TRIM(ADJUSTL(tag_scal))//'bcs.imin', BuffScalImin, txc(1,1), s, wrk3d )
    IF ( BuffScalImax%size > 0 ) CALL INI_BLOCK(1, TRIM(ADJUSTL(tag_scal))//'bcs.imax', BuffScalImax, txc(1,1), s, wrk3d )
    IF ( BuffScalJmin%size > 0 ) CALL INI_BLOCK(2, TRIM(ADJUSTL(tag_scal))//'bcs.jmin', BuffScalJmin, txc(1,1), s, wrk3d )
    IF ( BuffScalJmax%size > 0 ) CALL INI_BLOCK(2, TRIM(ADJUSTL(tag_scal))//'bcs.jmax', BuffScalJmax, txc(1,1), s, wrk3d )

    RETURN
  END SUBROUTINE BOUNDARY_BUFFER_INITIALIZE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE INI_BLOCK(idir, tag,item, rho, a, wrk3d)
#ifdef USE_MPI
    USE MPI
#endif

    IMPLICIT NONE

    TINTEGER,         INTENT(IN   ) :: idir         ! Wall-normal direction of buffer zone
    CHARACTER(LEN=*), INTENT(IN   ) :: tag          ! File name information
    TYPE(buffer_dt),  INTENT(INOUT) :: item
    TREAL,            INTENT(IN   ) :: rho(isize_field), a(isize_field,item%nfields)
    TREAL,            INTENT(INOUT) :: wrk3d(isize_field)

    ! -------------------------------------------------------------------
    TINTEGER io_sizes(5), id
    TREAL var_max,var_min

    CHARACTER*32 str, varname(item%nfields)
    CHARACTER*128 line

    TREAL COV2V1D, COV2V2D

! #ifdef USE_MPI
!     TINTEGER                :: ndims
!     TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
! #endif
!
    ! ###################################################################
    ! Reference fields
    IF ( idir == 1 ) ALLOCATE( item%ref(item%size,jmax,kmax,item%nfields) )
    IF ( idir == 2 ) ALLOCATE( item%ref(imax,item%size,kmax,item%nfields) )

    DO iq = 1,item%nfields; WRITE(varname(iq),*) iq; ENDDO

    SELECT CASE ( idir )
    CASE( 1 )
      id = IO_SUBARRAY_BUFFER_ZOY
      io_aux(id)%offset = 0
#ifdef USE_MPI
      IF ( ims_pro_i ==  ( item%offset /imax) ) io_aux(id)%active = .TRUE. ! info must belong to only 1 PE
      io_aux(id)%communicator = ims_comm_z
      io_aux(id)%subarray = IO_CREATE_SUBARRAY_ZOY( jmax*item%size,kmax, MPI_REAL8 )

      ! ndims = 2
      ! sizes(1)   = jmax *item%size; sizes(2)   = kmax *ims_npro_k
      ! locsize(1) = jmax *item%size; locsize(2) = kmax
      ! offset(1)  = 0;               offset(2)  = ims_offset_k
      !
      ! CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
      !     MPI_ORDER_FORTRAN, MPI_REAL8, io_aux(id)%subarray, ims_err)
      ! CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)
#endif
      idummy = item%size*jmax*kmax; io_sizes = (/idummy,1,idummy,1,item%nfields/)

    CASE( 2 )
      id = IO_SUBARRAY_BUFFER_XOZ
      io_aux(id)%offset = 0
#ifdef USE_MPI
      io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = MPI_COMM_WORLD
      io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ( imax,item%size,kmax, MPI_REAL8 )
      !
      ! ndims = 3
      ! sizes(1)  =imax *ims_npro_i; sizes(2)   = item%size; sizes(3)   = kmax *ims_npro_k
      ! locsize(1)=imax;             locsize(2) = item%size; locsize(3) = kmax
      ! offset(1) =ims_offset_i;     offset(2)  = 0;         offset(3)  = ims_offset_k
#endif
      idummy = imax*item%size*kmax; io_sizes = (/idummy,1,idummy,1,item%nfields/)

    END SELECT

! #ifdef USE_MPI
!     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
!         MPI_ORDER_FORTRAN, MPI_REAL8, io_aux(id)%subarray, ims_err)
!     CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)
! #endif

    IF ( BuffLoad ) THEN
      CALL IO_READ_SUBARRAY8(id, tag, varname, item%ref, io_sizes, wrk3d)

    ELSE

      SELECT CASE ( idir )
      CASE( 1 )
        DO iq = 1,item%nfields
          DO j = 1,jmax
            DO iloc = 1, item%size
              i = item%offset +iloc
              IF ( .NOT. item%hard ) item%hardvalues(iq) = COV2V1D(imax,jmax,kmax, i,j, rho, a(1,iq))
              item%ref(iloc,j,:,iq) = item%hardvalues(iq)
            ENDDO
          ENDDO
        ENDDO

      CASE( 2 )
        IF ( imode_sim == DNS_MODE_TEMPORAL ) THEN
          DO iq = 1,item%nfields
            DO jloc = 1, item%size
              j = item%offset +jloc
              IF ( .NOT. item%hard ) item%hardvalues(iq) = COV2V2D(imax,jmax,kmax, j, rho, a(1,iq))
              item%ref(:,jloc,:,iq) = item%hardvalues(iq)
            ENDDO
          ENDDO
        ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
          DO iq = 1,item%nfields
            DO jloc = 1, item%size
              j = item%offset +jloc
              DO i = 1,imax
                IF ( .NOT. item%hard ) item%hardvalues(iq) = COV2V1D(imax,jmax,kmax, i,j, rho, a(1,iq))
                item%ref(i,jloc,:,iq) = item%hardvalues(iq)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      END SELECT

      WRITE(str, *) itime; str = TRIM(ADJUSTL(tag))//'.'//TRIM(ADJUSTL(str))
      CALL IO_WRITE_SUBARRAY8(id, str, varname, item%ref, io_sizes, wrk3d)

    ENDIF

    DO iq = 1, item%nfields ! Control
      var_max = MAXVAL(item%ref(:,:,:,iq)); var_min = MINVAL(item%ref(:,:,:,iq))
#ifdef USE_MPI
      CALL MPI_ALLREDUCE(var_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
      var_max = dummy
      CALL MPI_ALLREDUCE(var_min, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
      var_min = dummy
#endif
      WRITE(line,10) var_max
      WRITE(str, 10) var_min
      line = TRIM(ADJUSTL(str))//' and '//TRIM(ADJUSTL(line))
      line = 'Checking bounds of field '//TRIM(ADJUSTL(tag))//'.'//TRIM(ADJUSTL(varname(iq)))//': '//TRIM(ADJUSTL(line))
      CALL TLAB_WRITE_ASCII(lfile,line)
    ENDDO

    ! -----------------------------------------------------------------------
    ! Strength of the telaxation terms
    ALLOCATE( item%tau(item%size,item%nfields) )
    DO iq = 1,item%nfields
      dummy = C_1_R /( g(idir)%nodes( item%offset+item%size ) - g(idir)%nodes( item%offset+1 ) ) ! Inverse of segment length
      DO jloc = 1,item%size
        j = item%offset +jloc
        IF ( item%form == FORM_POWER_MAX ) item%tau(jloc,iq) = item%strength(iq) *( ( g(idir)%nodes(j) -g(idir)%nodes(item%offset+1) ) *dummy ) **item%sigma(iq)
        IF ( item%form == FORM_POWER_MIN ) item%tau(jloc,iq) = item%strength(iq) *( ( g(idir)%nodes(item%offset+item%size) -g(idir)%nodes(j)) *dummy ) **item%sigma(iq)
      ENDDO
    ENDDO

    ! -----------------------------------------------------------------------
    ! Filters at boundaries; nseeds to be checked
#ifdef USE_MPI
    IF ( item%type == DNS_BUFFER_FILTER ) THEN
      SELECT CASE ( idir )
      CASE( 1 )
        CALL TLAB_WRITE_ASCII(lfile,'Initialize MPI types for Ox BCs explicit filter.')
        id     = TLAB_MPI_K_OUTBCS
        idummy = item%size*jmax
        CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, idummy, 1,1,1,1, &
            ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

      CASE( 2 )
        CALL TLAB_WRITE_ASCII(lfile,'Initialize MPI types for Oy BCs explicit filter.')
        id     = TLAB_MPI_K_TOPBCS
        idummy = imax*item%size
        CALL TLAB_MPI_TYPE_K(ims_npro_k, kmax, idummy, 1,1,1,1, &
            ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

      END SELECT
    ENDIF
#endif

    RETURN
10  FORMAT(G_FORMAT_R)
  END SUBROUTINE INI_BLOCK

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE BOUNDARY_BUFFER_RELAX_SCAL(s,hs, q)
    IMPLICIT NONE

    TREAL, INTENT(IN   ) :: s(isize_field,inb_scal)
    TREAL, INTENT(  OUT) :: hs(isize_field,inb_scal)
    TREAL, INTENT(IN   ) :: q(isize_field,inb_flow)   ! In case I need the density

    ! ###################################################################
    SELECT CASE( imode_eqns )
    CASE ( DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL )
      IF ( BuffScalImin%size > 0 ) CALL RELAX_BLOCK_RHO( 1, BuffScalImin, s,hs, q(:,5) )
      IF ( BuffScalImax%size > 0 ) CALL RELAX_BLOCK_RHO( 1, BuffScalImax, s,hs, q(:,5) )
      IF ( BuffScalJmin%size > 0 ) CALL RELAX_BLOCK_RHO( 2, BuffScalJmin, s,hs, q(:,5) )
      IF ( BuffScalJmax%size > 0 ) CALL RELAX_BLOCK_RHO( 2, BuffScalJmax, s,hs, q(:,5) )
    CASE DEFAULT
      IF ( BuffScalImin%size > 0 ) CALL RELAX_BLOCK( 1, BuffScalImin, s,hs )
      IF ( BuffScalImax%size > 0 ) CALL RELAX_BLOCK( 1, BuffScalImax, s,hs )
      IF ( BuffScalJmin%size > 0 ) CALL RELAX_BLOCK( 2, BuffScalJmin, s,hs )
      IF ( BuffScalJmax%size > 0 ) CALL RELAX_BLOCK( 2, BuffScalJmax, s,hs )
    END SELECT

    RETURN
  END SUBROUTINE BOUNDARY_BUFFER_RELAX_SCAL

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE BOUNDARY_BUFFER_RELAX_SCAL_I(is, s,hs)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: is ! field to which relaxation is applied
    TREAL,    INTENT(IN   ) :: s(isize_field)
    TREAL,    INTENT(  OUT) :: hs(isize_field)

    ! ###################################################################
    IF ( BuffScalImin%size > 0 ) CALL RELAX_BLOCK_I( 1, BuffScalImin, s,hs, is )
    IF ( BuffScalImax%size > 0 ) CALL RELAX_BLOCK_I( 1, BuffScalImax, s,hs, is )
    IF ( BuffScalJmin%size > 0 ) CALL RELAX_BLOCK_I( 2, BuffScalJmin, s,hs, is )
    IF ( BuffScalJmax%size > 0 ) CALL RELAX_BLOCK_I( 2, BuffScalJmax, s,hs, is )

    RETURN
  END SUBROUTINE BOUNDARY_BUFFER_RELAX_SCAL_I

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE BOUNDARY_BUFFER_RELAX_FLOW(q,hq)
    IMPLICIT NONE

    TREAL, INTENT(IN   ) :: q(isize_field,inb_flow)
    TREAL, INTENT(  OUT) :: hq(isize_field,inb_flow)

    ! ###################################################################
    SELECT CASE( imode_eqns )
    CASE ( DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL )
      IF ( BuffFlowImin%size > 0 ) CALL RELAX_BLOCK_CF( 1, BuffFlowImin, q,hq )
      IF ( BuffFlowImax%size > 0 ) CALL RELAX_BLOCK_CF( 1, BuffFlowImax, q,hq )
      IF ( BuffFlowJmin%size > 0 ) CALL RELAX_BLOCK_CF( 2, BuffFlowJmin, q,hq )
      IF ( BuffFlowJmax%size > 0 ) CALL RELAX_BLOCK_CF( 2, BuffFlowJmax, q,hq )
    CASE DEFAULT
      IF ( BuffFlowImin%size > 0 ) CALL RELAX_BLOCK( 1, BuffFlowImin, q,hq )
      IF ( BuffFlowImax%size > 0 ) CALL RELAX_BLOCK( 1, BuffFlowImax, q,hq )
      IF ( BuffFlowJmin%size > 0 ) CALL RELAX_BLOCK( 2, BuffFlowJmin, q,hq )
      IF ( BuffFlowJmax%size > 0 ) CALL RELAX_BLOCK( 2, BuffFlowJmax, q,hq )
    END SELECT

    RETURN
  END SUBROUTINE BOUNDARY_BUFFER_RELAX_FLOW

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE RELAX_BLOCK( idir, item, s, hs )
    IMPLICIT NONE

    TINTEGER,         INTENT(IN   ) :: idir
    TYPE(buffer_dt),  INTENT(IN   ) :: item
    TREAL,            INTENT(IN   ) :: s(imax,jmax,kmax,item%nfields)
    TREAL,            INTENT(  OUT) :: hs(imax,jmax,kmax,item%nfields)

    ! ###################################################################
    SELECT CASE ( idir )
    CASE( 1 )
      DO iq = 1,item%nfields
        DO jloc = 1, item%size
          j = item%offset +jloc
          hs(j,:,:,iq) = hs(j,:,:,iq) -item%Tau(jloc,iq) *( s(j,:,:,iq) -item%Ref(jloc,:,:,iq))
        ENDDO
      ENDDO

    CASE( 2 )
      DO iq = 1,item%nfields
        DO jloc = 1, item%size
          j = item%offset +jloc
          hs(:,j,:,iq) = hs(:,j,:,iq) -item%Tau(jloc,iq) *( s(:,j,:,iq) -item%Ref(:,jloc,:,iq))
        ENDDO
      ENDDO

    END SELECT

    RETURN
  END SUBROUTINE RELAX_BLOCK

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE RELAX_BLOCK_RHO( idir, item, s, hs, rho )
    IMPLICIT NONE

    TINTEGER,         INTENT(IN   ) :: idir
    TYPE(buffer_dt),  INTENT(IN   ) :: item
    TREAL,            INTENT(IN   ) :: s(imax,jmax,kmax,item%nfields)
    TREAL,            INTENT(  OUT) :: hs(imax,jmax,kmax,item%nfields)
    TREAL,            INTENT(IN   ) :: rho(imax,jmax,kmax)

    ! ###################################################################
    SELECT CASE ( idir )
    CASE( 1 )
      DO iq = 1,item%nfields
        DO jloc = 1, item%size
          j = item%offset +jloc
          hs(j,:,:,iq) = hs(j,:,:,iq) -item%Tau(jloc,iq) *( rho(j,:,:) *s(j,:,:,iq) -item%Ref(jloc,:,:,iq))
        ENDDO
      ENDDO

    CASE( 2 )
      DO iq = 1,item%nfields
        DO jloc = 1, item%size
          j = item%offset +jloc
          hs(:,j,:,iq) = hs(:,j,:,iq) -item%Tau(jloc,iq) *( rho(:,j,:) *s(:,j,:,iq) -item%Ref(:,jloc,:,iq))
        ENDDO
      ENDDO

    END SELECT

    RETURN
  END SUBROUTINE RELAX_BLOCK_RHO

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE RELAX_BLOCK_I( idir, item, s, hs, iq )
    IMPLICIT NONE

    TINTEGER,         INTENT(IN   ) :: idir,iq
    TYPE(buffer_dt),  INTENT(IN   ) :: item
    TREAL,            INTENT(IN   ) :: s(imax,jmax,kmax)
    TREAL,            INTENT(  OUT) :: hs(imax,jmax,kmax)

    ! ###################################################################
    SELECT CASE ( idir )
    CASE( 1 )
      DO jloc = 1, item%size
        j = item%offset +jloc
        hs(j,:,:) = hs(j,:,:) -item%Tau(jloc,iq) *( s(j,:,:) -item%Ref(jloc,:,:,iq) )
      ENDDO

    CASE( 2 )
      DO jloc = 1, item%size
        j = item%offset +jloc
        hs(:,j,:) = hs(:,j,:) -item%Tau(jloc,iq) *( s(:,j,:) -item%Ref(:,jloc,:,iq) )
      ENDDO

    END SELECT

    RETURN
  END SUBROUTINE RELAX_BLOCK_I

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE RELAX_BLOCK_CF( idir, item, q,hq )  ! Compressible flow case
    IMPLICIT NONE

    TINTEGER,         INTENT(IN   ) :: idir
    TYPE(buffer_dt),  INTENT(IN   ) :: item
    TREAL,            INTENT(IN   ) :: q(imax,jmax,kmax,item%nfields)
    TREAL,            INTENT(  OUT) :: hq(imax,jmax,kmax,item%nfields)

    ! ###################################################################
    IF ( imode_eqns == DNS_EQNS_TOTAL ) THEN
      dummy = C_05_R *(gama0-C_1_R) *mach*mach
    END IF

    SELECT CASE ( idir )
    CASE( 1 )
      DO iq = 1,item%nfields-2
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(j,:,:,iq) = hq(j,:,:,iq) -item%Tau(jloc,iq) *( q(j,:,:,5) *q(j,:,:,iq) -item%Ref(jloc,:,:,iq))
        ENDDO
      ENDDO
      iq = item%nfields-1       ! Energy variable
        IF ( imode_eqns == DNS_EQNS_TOTAL ) THEN
          DO jloc = 1, item%size
            j = item%offset +jloc
            hq(j,:,:,iq) = hq(j,:,:,iq) -item%Tau(jloc,iq) *( q(j,:,:,5) *( q(j,:,:,iq) &
                          +dummy *( q(j,:,:,1)*q(j,:,:,1) +q(j,:,:,2)*q(j,:,:,2) +q(j,:,:,3)*q(j,:,:,3) ) ) &
                          -item%Ref(jloc,:,:,iq) )
          ENDDO
        ELSE
          DO jloc = 1, item%size
            j = item%offset +jloc
            hq(j,:,:,iq) = hq(j,:,:,iq) -item%Tau(jloc,iq) *( q(j,:,:,5) *q(j,:,:,iq) -item%Ref(jloc,:,:,iq))
          ENDDO
        ENDIF
      iq = item%nfields         ! Density, should be iq = 5
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(j,:,:,iq) = hq(j,:,:,iq) -item%Tau(jloc,iq) *( q(j,:,:,iq) -item%Ref(jloc,:,:,iq))
        ENDDO

    CASE( 2 )
      DO iq = 1,item%nfields-2
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(:,j,:,iq) = hq(:,j,:,iq) -item%Tau(jloc,iq) *( q(:,j,:,5) *q(:,j,:,iq) -item%Ref(:,jloc,:,iq))
        ENDDO
      ENDDO
      iq = item%nfields-1       ! Energy variable
      IF ( imode_eqns == DNS_EQNS_TOTAL ) THEN
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(:,j,:,iq) = hq(:,j,:,iq) -item%Tau(jloc,iq) *( q(:,j,:,5) *( q(:,j,:,iq) &
                        +dummy *( q(:,j,:,1)*q(:,j,:,1) +q(:,j,:,2)*q(:,j,:,2) +q(:,j,:,3)*q(:,j,:,3) ) ) &
                        -item%Ref(:,jloc,:,iq) )
        END DO
      ELSE
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(:,j,:,iq) = hq(:,j,:,iq) -item%Tau(jloc,iq) *( q(:,j,:,5) *q(:,j,:,iq) -item%Ref(:,jloc,:,iq))
        ENDDO
      ENDIF
      iq = item%nfields         ! Density, should be iq = 5
        DO jloc = 1, item%size
          j = item%offset +jloc
          hq(:,j,:,iq) = hq(:,j,:,iq) -item%Tau(jloc,iq) *( q(:,j,:,iq) -item%Ref(:,jloc,:,iq))
        ENDDO

    END SELECT

    RETURN
  END SUBROUTINE RELAX_BLOCK_CF

  !########################################################################
  !########################################################################
  SUBROUTINE BOUNDARY_BUFFER_FILTER(rho,u,v,w,e,z1, txc1,txc2,txc3,txc4,txc5, wrk1d,wrk2d,wrk3d)

    IMPLICIT NONE

#include "integers.h"

    TREAL, DIMENSION(imax,jmax,kmax  )        :: rho, u, v, w, e
    TREAL, DIMENSION(imax,jmax,kmax,*)        :: z1
    TREAL, DIMENSION(BuffFlowImax%size,jmax,kmax) :: txc1, txc2, txc3, txc4, txc5
    TREAL, DIMENSION(*)                       :: wrk1d,wrk2d,wrk3d

    ! -------------------------------------------------------------------
    TINTEGER id, k, ibc_x(4), ibc_y(4), ibc_z(4), buff_imax
    TREAL eta, delta, amp, ampr, rho_ratio

    ! ###################################################################
!!! Routines OPR_FILTER have been changed. This routine needs to be updates !!!
    CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
    CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)


    ! BCs for the filters (see routine FILTER)
    ! ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
    ! ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0
    ! ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0

    ! ###################################################################
    ! Bottom boundary
    ! ###################################################################
    IF ( BuffFlowJmin%size .GT. 1 ) THEN
      CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
      CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
    ENDIF

    ! ###################################################################
    ! Top boundary
    ! ###################################################################
    IF ( BuffFlowJmax%size .GT. 1 ) THEN
      CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Filter not yet implemented.')
      CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
    ENDIF

    ! ###################################################################
    ! Outflow boundary
    ! ###################################################################
    IF ( BuffFlowImax%size .GT. 1 ) THEN
      id = TLAB_MPI_K_OUTBCS
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
      !     IF ( icalc_scal .EQ. 1 ) THEN
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
      !     ENDIF

    ENDIF

    RETURN
  END SUBROUTINE BOUNDARY_BUFFER_FILTER
END MODULE BOUNDARY_BUFFER

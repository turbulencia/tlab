#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

MODULE PLANES

  USE DNS_CONSTANTS, ONLY : lfile, efile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field, inb_scal_array, inb_flow_array
  USE DNS_GLOBAL, ONLY : rbackground, g
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : io_aux
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE
  SAVE
  PRIVATE

  TINTEGER, PARAMETER, PUBLIC :: MAX_SAVEPLANES = 20

  TYPE planes_dt
    SEQUENCE
    TINTEGER n
    TINTEGER nodes(MAX_SAVEPLANES)
    TINTEGER size
    TINTEGER io(5)
  END TYPE planes_dt

  TYPE(planes_dt), PUBLIC :: iplanes,jplanes,kplanes
  CHARACTER*32 varname(1)
  TINTEGER idummy

  PUBLIC :: PLANES_INITIALIZE, PLANES_SAVE

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE PLANES_INITIALIZE()

#ifdef USE_MPI
    USE DNS_MPI
#endif

    IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"

    ! -----------------------------------------------------------------------
    TINTEGER                :: ndims, id
    TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
#endif

    ! ###################################################################
    iplanes%size = (inb_flow_array +inb_scal_array +1) *iplanes%n           ! Flow and scal variables, pressure
    jplanes%size = (inb_flow_array +inb_scal_array +1) *jplanes%n
    kplanes%size = (inb_flow_array +inb_scal_array +1) *kplanes%n
    IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) jplanes%size = jplanes%size +2  ! Add LWP and intgral of TWP

    IF ( iplanes%size > imax ) THEN
      CALL IO_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size imax is is insufficient.')
      CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
    END IF
    IF ( jplanes%size > jmax ) THEN
      CALL IO_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size jmax is is insufficient.')
      CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
    END IF
    IF ( kplanes%size > kmax ) THEN
      CALL IO_WRITE_ASCII(efile, 'PLANES_INITIALIZE. Array size kmax is is insufficient.')
      CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
    END IF

    ! Info for IO routines: total size, lower bound, upper bound, stride, # variables
    idummy = jmax* kmax *iplanes%size; iplanes%io = [ idummy, 1, idummy, 1, 1 ]
    idummy = kmax *imax *jplanes%size; jplanes%io = [ idummy, 1, idummy, 1, 1 ]
    idummy = imax *jmax *kplanes%size; kplanes%io = [ idummy, 1, idummy, 1, 1 ]
    varname = [ '' ]

    io_aux(:)%offset = 0        ! defaults

#ifdef USE_MPI
    io_aux(:)%active = .FALSE.  ! defaults

    id = IO_SUBARRAY_PLANES_XOY
    IF ( kplanes%n > 0 ) THEN ! Saving full vertical xOy planes; writing only info of PE containing the first plane
      IF ( ims_pro_k == ( kplanes%nodes(1) /kmax) ) io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = ims_comm_x

      ndims = 2
      sizes(1)   = imax *ims_npro_i; sizes(2)   = jmax *kplanes%size
      locsize(1) = imax;             locsize(2) = jmax *kplanes%size
      offset(1)  = ims_offset_i;     offset(2)  = 0

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF

    id = IO_SUBARRAY_PLANES_ZOY
    IF ( iplanes%n > 0 ) THEN ! Saving full vertical zOy planes; writing only info of PE containing the first plane
      IF ( ims_pro_i ==  ( iplanes%nodes(1) /imax) ) io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = ims_comm_z

      ndims = 2
      sizes(1)   = jmax *iplanes%size;  sizes(2)   = kmax *ims_npro_k
      locsize(1) = jmax *iplanes%size;  locsize(2) = kmax
      offset(1)  = 0;                   offset(2)  = ims_offset_k

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF

    id = IO_SUBARRAY_PLANES_XOZ
    IF ( jplanes%n > 0 ) THEN ! Saving full blocks xOz planes for prognostic variables
      io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = MPI_COMM_WORLD

      ndims = 3 ! Subarray for the output of the 2D data
      sizes(1)  =imax *ims_npro_i;  sizes(2)   = jplanes%size;  sizes(3)   = kmax *ims_npro_k
      locsize(1)=imax;              locsize(2) = jplanes%size;  locsize(3) = kmax
      offset(1) =ims_offset_i;      offset(2)  = 0;             offset(3)  = ims_offset_k

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF
#endif

    RETURN
  END SUBROUTINE PLANES_INITIALIZE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE PLANES_SAVE(q,s, p, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)
    IMPLICIT NONE

    TREAL, INTENT(IN   ) :: q(imax,jmax,kmax,inb_flow_array)
    TREAL, INTENT(IN   ) :: s(imax,jmax,kmax,inb_scal_array)
    TREAL, INTENT(INOUT) :: p(imax,jmax,kmax)             ! larger arrays for the Poisson solver,
    TREAL, INTENT(INOUT) :: tmp1(imax,jmax,kplanes%size)  ! but local shapes are used
    TREAL, INTENT(INOUT) :: tmp2(imax,jplanes%size,kmax)  ! to arrange plane data.
    TREAL, INTENT(INOUT) :: wrk3d(jmax,iplanes%size,kmax) ! In this array we transpose data
    TREAL, INTENT(INOUT) :: tmp3(imax,jmax,kmax,3)
    TREAL, INTENT(INOUT) :: wrk1d(*)
    TREAL, INTENT(INOUT) :: wrk2d(imax,kmax,*)

    ! -------------------------------------------------------------------
    TINTEGER offset, j, k
    CHARACTER*32 fname
    CHARACTER*250 line1

    ! ###################################################################
    WRITE(fname,100) rtime
    WRITE(line1,*) itime; line1 = 'Writing planes at It'//TRIM(ADJUSTL(line1))//' and time '//TRIM(ADJUSTL(fname))//'.'
    CALL IO_WRITE_ASCII(lfile,line1)
    100 FORMAT(G_FORMAT_R)

    CALL FI_PRESSURE_BOUSSINESQ(q,s, p, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

    IF ( kplanes%n > 0 ) THEN
      offset = 0
      DO idummy = 1,inb_flow_array
        tmp1(:,:,1+offset:kplanes%n+offset) = q(:,:,kplanes%nodes(1:kplanes%n),idummy)
        offset = offset +kplanes%n
      END DO
      DO idummy = 1,inb_scal_array
        tmp1(:,:,1+offset:kplanes%n+offset) = s(:,:,kplanes%nodes(1:kplanes%n),idummy)
        offset = offset +kplanes%n
      END DO
      tmp1(:,:,1+offset:kplanes%n+offset) = p(:,:,kplanes%nodes(1:kplanes%n))
      offset = offset +kplanes%n
      WRITE(fname,*) itime; fname = 'planesK.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOY, fname, varname, tmp1, kplanes%io, wrk3d)
    END IF

    IF ( jplanes%n > 0 ) THEN
      offset = 0
      DO idummy = 1,inb_flow_array
        tmp2(:,1+offset:jplanes%n+offset,:) = q(:,jplanes%nodes(1:jplanes%n),:,idummy)
        offset = offset +jplanes%n
      END DO
      DO idummy = 1,inb_scal_array
        tmp2(:,1+offset:jplanes%n+offset,:) = s(:,jplanes%nodes(1:jplanes%n),:,idummy)
        offset = offset +jplanes%n
      END DO
      tmp2(:,1+offset:jplanes%n+offset,:) = p(:,jplanes%nodes(1:jplanes%n),:)
      offset = offset +jplanes%n
      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN    ! Add LWP and intgral of TWP
        CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, s(1,1,1,inb_scal_array  ), wrk2d, wrk1d,wrk3d)
        tmp2(:,1 +offset,:) = wrk2d(:,:,1)
        offset = offset +1
        CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, s(1,1,1,inb_scal_array-1), wrk2d, wrk1d,wrk3d)
        tmp2(:,1 +offset,:) = wrk2d(:,:,1)
        offset = offset +1
    END IF
      WRITE(fname,*) itime; fname = 'planesJ.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOZ, fname, varname, tmp2, jplanes%io, wrk3d)
    END IF

    IF ( iplanes%n > 0 ) THEN       ! We transpose to make j-lines together in memory
      offset = 0
      DO idummy = 1,inb_flow_array
        DO k = 1,kmax; DO j = 1,jmax
          wrk3d(j,1+offset:iplanes%n+offset,k) = q(iplanes%nodes(1:iplanes%n),j,k,idummy)
        END DO; END DO
        offset = offset +iplanes%n
      END DO
      DO idummy = 1,inb_scal_array
        DO k = 1,kmax; DO j = 1,jmax
          wrk3d(j,1+offset:iplanes%n+offset,k) = s(iplanes%nodes(1:iplanes%n),j,k,idummy)
        END DO; END DO
        offset = offset +iplanes%n
      END DO
      DO k = 1,kmax; DO j = 1,jmax
        wrk3d(j,1+offset:iplanes%n+offset,k) = p(iplanes%nodes(1:iplanes%n),j,k)
      END DO; END DO
      offset = offset +iplanes%n
      WRITE(fname,*) itime; fname = 'planesI.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_ZOY, fname, varname, wrk3d, iplanes%io, tmp1)
    END IF

    RETURN
  END SUBROUTINE PLANES_SAVE

END MODULE PLANES

#include "types.h"
#include "dns_const.h"

MODULE PLANES

  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field, inb_scal_array, inb_flow_array
  USE DNS_GLOBAL, ONLY : rbackground, g
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : io_aux
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE
  SAVE
  PRIVATE

  TINTEGER, PARAMETER,                 PUBLIC :: MAX_SAVEPLANES = 20
  TINTEGER,                            PUBLIC :: nplanes_i, nplanes_j, nplanes_k
  TINTEGER,                            PUBLIC :: pplanes_j, nplanes_j_aux ! to be removed
  TINTEGER, DIMENSION(MAX_SAVEPLANES), PUBLIC :: planes_i,  planes_j,  planes_k

  TINTEGER idummy, splanes_i(5), splanes_j(5), splanes_k(5), splanes_jp(5)
  CHARACTER*32 varname(1)

  PUBLIC :: PLANES_INITIALIZE, PLANES_SAVE

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE PLANES_INITIALIZE()
    IMPLICIT NONE

#ifdef USE_MPI
    CALL DNS_MPIO_AUX
#else
    io_aux(:)%offset = 0
#endif

    ! Sizes information for saving planes
    idummy     = kmax*jmax*nplanes_i*(inb_flow_array+inb_scal_array)
    splanes_i  = (/idummy,1,idummy,1,1/)
    idummy     = imax*jmax*nplanes_k*(inb_flow_array+inb_scal_array)
    splanes_k  = (/idummy,1,idummy,1,1/)
    idummy     = imax*kmax*nplanes_j*(inb_flow_array+inb_scal_array)
    splanes_j  = (/idummy,1,idummy,1,1/)
    idummy     = imax*kmax*nplanes_j
    splanes_jp = (/idummy,1,idummy,1,1/)
    varname    = (/''/)

    RETURN
  END SUBROUTINE PLANES_INITIALIZE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE PLANES_SAVE(q,s, hq,txc, wrk1d,wrk2d,wrk3d)
    IMPLICIT NONE

    TREAL, DIMENSION(isize_field,*), INTENT(IN   ) :: q,s
    TREAL, DIMENSION(*),             INTENT(INOUT) :: hq,txc, wrk1d,wrk2d,wrk3d

    ! -------------------------------------------------------------------
    TINTEGER is, iq, ip
    CHARACTER*250 line1
    CHARACTER*32 fname

    ! ###################################################################
    WRITE(fname,100) rtime
    WRITE(line1,*) itime; line1 = 'Writing planes at It'//TRIM(ADJUSTL(line1))//' and time '//TRIM(ADJUSTL(fname))//'.'
    CALL IO_WRITE_ASCII(lfile,line1)
100 FORMAT(G_FORMAT_R)

    IF ( nplanes_k .GT. 0 ) THEN
      CALL REDUCE_Z_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, nplanes_k,planes_k, txc)
      WRITE(fname,*) itime; fname = 'planesK.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOY, fname, varname, txc, splanes_k, hq)
    ENDIF

    IF ( nplanes_j .GT. 0 ) THEN
      IF ( nplanes_j_aux .GT. 0 ) THEN ! Calculate integrals
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
          ip = 1
          DO iq = 1,inb_flow_array
            CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, q(1,iq), wrk3d(ip), wrk1d,hq) !hq is aux variable
            ip = ip + imax*kmax
          ENDDO
          DO is = 1,inb_scal_array
            CALL THERMO_ANELASTIC_LWP(imax,jmax,kmax, g(2), rbackground, s(1,is), wrk3d(ip), wrk1d,hq)
            ip = ip + imax*kmax
          ENDDO
        ELSE
          wrk3d(1:(inb_flow_array+inb_scal_array)*imax*kmax) = C_0_R
        ENDIF
      ENDIF
      CALL REDUCE_Y_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, wrk3d, nplanes_j,nplanes_j_aux,planes_j, txc)
      WRITE(fname,*) itime; fname = 'planesJ.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOZ, fname, varname, txc, splanes_j, hq)

      IF ( pplanes_j .EQ. 1 ) THEN   !calculate and write pressure for xOz planes
        CALL FI_PRESSURE_BOUSSINESQ(q,s,txc(1), &
            txc(1+isize_txc_field),txc(1+2*isize_txc_field),txc(1+3*isize_txc_field), &
            wrk1d,wrk2d,wrk3d)
        CALL REDUCE_Y_ALL(imax,jmax,kmax, 1 ,txc(1), 0,s, wrk3d, nplanes_j,nplanes_j_aux,planes_j, txc)
        WRITE(fname,*) itime; fname = 'pressrJ.'//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_XOZ_P, fname, varname, txc, splanes_jp, hq)
      ENDIF
    ENDIF

    IF ( nplanes_i .GT. 0 ) THEN
      CALL REDUCE_X_ALL(imax,jmax,kmax, inb_flow_array,q, inb_scal_array,s, nplanes_i,planes_i, txc)
      WRITE(fname,*) itime; fname = 'planesI.'//TRIM(ADJUSTL(fname))
      CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_PLANES_ZOY, fname, varname, txc, splanes_i, hq)
    ENDIF

    RETURN
  END SUBROUTINE PLANES_SAVE

  ! ###################################################################
  ! ###################################################################
#ifdef USE_MPI

  SUBROUTINE DNS_MPIO_AUX()
    USE DNS_MPI

    IMPLICIT NONE

#include "mpif.h"

    ! -----------------------------------------------------------------------
    TINTEGER                :: ndims, idummy, id
    TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

    ! #######################################################################
    io_aux(:)%active = .FALSE. ! defaults
    io_aux(:)%offset = 0

    ! ###################################################################
    ! Subarray information to write planes
    ! ###################################################################
    idummy = inb_flow_array +inb_scal_array

    id = IO_SUBARRAY_PLANES_XOY
    IF ( nplanes_k .GT. 0 ) THEN ! Saving full vertical xOy planes; writing only info of PE containing the first plane
      IF ( ims_pro_k .EQ. ( planes_k(1) /kmax) ) io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = ims_comm_x

      ndims = 2
      sizes(1)   = imax *ims_npro_i; sizes(2)   = jmax *nplanes_k *idummy
      locsize(1) = imax;             locsize(2) = jmax *nplanes_k *idummy
      offset(1)  = ims_offset_i;     offset(2)  = 0

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF

    id = IO_SUBARRAY_PLANES_ZOY
    IF ( nplanes_i .GT. 0 ) THEN ! Saving full vertical zOy planes; writing only info of PE containing the first plane
      IF ( ims_pro_i .EQ.  ( planes_i(1) /imax) ) io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = ims_comm_z

      ndims = 2
      sizes(1)   = jmax *nplanes_i *idummy; sizes(2)   = kmax *ims_npro_k
      locsize(1) = jmax *nplanes_i *idummy; locsize(2) = kmax
      offset(1)  = 0;                       offset(2)  = ims_offset_k

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF

    id = IO_SUBARRAY_PLANES_XOZ
    IF ( nplanes_j .GT. 0 ) THEN ! Saving full blocks xOz planes for prognostic variables
      io_aux(id)%active = .TRUE.
      io_aux(id)%communicator = MPI_COMM_WORLD

      ndims = 3 ! Subarray for the output of the 2D data
      sizes(1)  =imax *ims_npro_i; sizes(2)   = nplanes_j*idummy; sizes(3)   = kmax *ims_npro_k
      locsize(1)=imax;             locsize(2) = nplanes_j*idummy; locsize(3) = kmax
      offset(1) =ims_offset_i;     offset(2)  = 0;                offset(3)  = ims_offset_k

      CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
      CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ENDIF

    IF ( pplanes_j .EQ. 1 ) THEN

      id = IO_SUBARRAY_PLANES_XOZ_P
      IF ( nplanes_j .GT. 0 ) THEN ! Saving full blocks xOz planes for pressure
        io_aux(id)%active = .TRUE.
        io_aux(id)%communicator = MPI_COMM_WORLD

        ndims = 3 ! Subarray for the output of the 2D data
        sizes(1)  =imax*ims_npro_i; sizes(2)   = nplanes_j; sizes(3)   = kmax *ims_npro_k
        locsize(1)=imax;            locsize(2) = nplanes_j; locsize(3) = kmax
        offset(1) =ims_offset_i;    offset(2)  = 0;         offset(3)  = ims_offset_k

        CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
            MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
        CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE DNS_MPIO_AUX

#endif

END MODULE PLANES

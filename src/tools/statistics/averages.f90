#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

#define C_FILE_LOC "AVERAGES"

PROGRAM AVERAGES

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : ifile,efile,lfile,gfile, tag_flow,tag_scal,tag_part
  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_err
  USE DNS_MPI, ONLY : ims_npro_i, ims_npro_k
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
  USE TLAB_MPI_PROCS
#endif
  USE THERMO_GLOBAL, ONLY : imixture
  USE LAGRANGE_GLOBAL
  USE LAGRANGE_ARRAYS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  ! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  =  20
  TINTEGER, PARAMETER :: igate_size_max =   8
  TINTEGER, PARAMETER :: params_size_max=   2

  ! -------------------------------------------------------------------
  ! Additional local arrays
  TREAL,      ALLOCATABLE, SAVE         :: mean(:), y_aux(:)
  INTEGER(1), ALLOCATABLE, SAVE         :: gate(:)
  TYPE(pointers_dt)                     :: vars(16)
  TREAL,      ALLOCATABLE, SAVE         :: surface(:,:,:)       ! Gate envelopes

  ! -------------------------------------------------------------------
  ! Local variables
  ! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, bakfile
  CHARACTER*32 varname(16)

  TINTEGER opt_main, opt_block, opt_order
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER nfield, ifield, is, ij, k, bcs(2,2)
  TREAL eloc1, eloc2, eloc3, cos1, cos2, cos3, dummy
  TINTEGER jmax_aux, iread_flow, iread_scal, idummy

  ! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  TREAL gate_threshold(igate_size_max)
  INTEGER(1) gate_level

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TREAL opt_vec(iopt_size_max)
  TREAL opt_vec2(iopt_size_max)

  TINTEGER params_size
  TREAL params(params_size_max)

  TINTEGER io_sizes(5)
#ifdef USE_MPI
  TINTEGER                :: ndims, id
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
#endif

  ! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u, v, w

  !########################################################################
  !########################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL TLAB_START()

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part == 1 ) THEN
    CALL PARTICLE_READ_GLOBAL(ifile)
  END IF

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  ! -------------------------------------------------------------------
  ! File names
  ! -------------------------------------------------------------------
#include "dns_read_times.h"

  ! -------------------------------------------------------------------
  ! Additional options
  ! -------------------------------------------------------------------
  opt_main  =-1 ! default values
  opt_block = 1
  gate_level= 0
  opt_order = 1

  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamAverages', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes == '-1' ) THEN
#ifdef USE_MPI
#else
    WRITE(*,*) 'Option ?'
    WRITE(*,*) ' 1. Conventional averages'
    WRITE(*,*) ' 2. Intermittency or gate function'
    WRITE(*,*) ' 3. Momentum equation'
    WRITE(*,*) ' 4. Main variables'
    WRITE(*,*) ' 5. Enstrophy W_iW_i/2 equation'
    WRITE(*,*) ' 6. Strain 2S_ijS_ij/2 equation'
    WRITE(*,*) ' 7. Scalar gradient G_iG_i/2 equation'
    WRITE(*,*) ' 8. Velocity gradient invariants'
    WRITE(*,*) ' 9. Scalar gradient components'
    WRITE(*,*) '10. Eigenvalues of rate-of-strain tensor'
    WRITE(*,*) '11. Eigenframe of rate-of-strain tensor'
    WRITE(*,*) '12. Longitudinal velocity derivatives'
    WRITE(*,*) '13. Vertical fluxes'
    WRITE(*,*) '14. Pressure partition'
    WRITE(*,*) '15. Dissipation'
    WRITE(*,*) '16. Third-order scalar covariances'
    WRITE(*,*) '17. Potential vorticity'
    READ(*,*) opt_main

    WRITE(*,*) 'Planes block size ?'
    READ(*,*) opt_block

    IF ( opt_main > 2 ) THEN
      WRITE(*,*) 'Gate level to be used ?'
      READ(*,*) gate_level
      WRITE(*,*) 'Number of moments ?'
      READ(*,*) opt_order
    END IF

#endif
  ELSE
    opt_main = INT(opt_vec(1))
    IF ( iopt_size >= 2 ) opt_block = INT(opt_vec(2))
    IF ( iopt_size >= 3 ) gate_level= INT(opt_vec(3),KIND=1)
    IF ( iopt_size >= 3 ) opt_order = INT(opt_vec(4))

  END IF

  IF ( opt_main < 0 ) THEN ! Check
    CALL TLAB_WRITE_ASCII(efile, 'AVERAGES. Missing input [ParamAverages] in dns.ini.')
    CALL TLAB_STOP(DNS_ERROR_INVALOPT)
  END IF

  IF ( opt_block < 1 ) THEN
    CALL TLAB_WRITE_ASCII(efile, 'AVERAGES. Invalid value of opt_block.')
    CALL TLAB_STOP(DNS_ERROR_INVALOPT)
  END IF

  ! -------------------------------------------------------------------
  iread_flow = 0
  iread_scal = 0
  inb_txc    = 0
  nfield     = 2

  SELECT CASE ( opt_main )
  CASE ( 1 )
    inb_txc = MAX(inb_txc,9)
    iread_flow = icalc_flow; iread_scal = icalc_scal
  CASE ( 2 )
    ifourier = 0
  CASE ( 3 )
    nfield = 14
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,12)
  CASE ( 4 )
    nfield = 6 +inb_scal
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,3)
    IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) inb_txc = MAX(inb_txc,6)
  CASE ( 5 ) ! enstrophy
    nfield = 7
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,8)
  CASE ( 6 )
    nfield = 5
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,8)
  CASE ( 7 ) ! scalar gradient
    nfield = 5
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,6)
  CASE ( 8 )
    nfield = 3

    iread_flow = 1;                 inb_txc = MAX(inb_txc,6)
  CASE ( 9 )
    nfield = 5
    iread_scal = 1; inb_txc = MAX(inb_txc,4)
  CASE (10 ) ! eigenvalues
    nfield = 3
    iread_flow = 1;                 inb_txc = MAX(inb_txc,9)
  CASE (11 ) ! eigenframe
    nfield = 6
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,10)
  CASE (12 ) ! longitudinal velocity derivatives
    nfield = 3
    iread_flow = 1; iread_scal = 0; inb_txc = MAX(inb_txc,3)
  CASE (13 ) ! Vertical flux
    nfield = 2*(3+inb_scal_array)
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(MAX(inb_txc,3+inb_scal_array),4)
  CASE (14 ) ! pressure partition
    nfield = 3
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,7)
  CASE (15 ) ! dissipation partition
    nfield = 1
    iread_flow = 1; iread_scal = 0; inb_txc = MAX(inb_txc,6)
  CASE (16 ) ! third-order scalar covariances
    nfield = 3
    iread_flow = 0; iread_scal = 1; inb_txc = MAX(inb_txc,3)
  CASE (17 ) ! potential vorticity
    nfield = 2
    iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,6)
  END SELECT

  ! -------------------------------------------------------------------
  ! Defining gate levels for conditioning
  ! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
  opt_cond_relative = 0
  igate_size    = 0

  IF ( opt_main == 2 .OR. gate_level /= 0 ) THEN
#include "dns_read_partition.h"
    IF ( opt_cond > 1 ) inb_txc = MAX(inb_txc,5)
  END IF

  ! -------------------------------------------------------------------
  ! Allocating memory space
  ! -------------------------------------------------------------------
  ALLOCATE(gate(isize_field))

  ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
  jmax_aux = g(2)%size/opt_block
  ALLOCATE(y_aux(jmax_aux)) ! Reduced vertical grid

  ! Subarray information to read and write envelopes
  IF ( opt_main == 2 .AND. opt_cond > 1 .AND. igate_size /= 0 ) THEN
    io_sizes = (/imax*2*igate_size*kmax,1,imax*2*igate_size*kmax,1,1/)

#ifdef USE_MPI
    id = IO_SUBARRAY_ENVELOPES

    io_aux(id)%active = .TRUE.
    io_aux(id)%communicator = MPI_COMM_WORLD

    ndims = 3
    sizes(1)  =imax *ims_npro_i; sizes(2)   = igate_size *2; sizes(3)   = kmax *ims_npro_k
    locsize(1)=imax;             locsize(2) = igate_size *2; locsize(3) = kmax
    offset(1) =ims_offset_i;     offset(2)  = 0;             offset(3)  = ims_offset_k

    CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
        MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
    CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

#else
    io_aux(:)%offset = 0
#endif

    ALLOCATE(surface(imax,2*igate_size,kmax))

  END IF

  IF ( opt_main == 1 ) THEN
    ALLOCATE(mean(jmax*MAX_AVG_TEMPORAL))
  ELSE IF ( opt_main == 2 ) THEN
    ALLOCATE(mean(igate_size*(jmax_aux+1)))
  ELSE
    ALLOCATE(mean(opt_order*nfield*(jmax_aux+1)))
  END IF

  isize_wrk3d = MAX(isize_field,opt_order*nfield*jmax)
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)
  IF ( icalc_part == 1) THEN
    isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  END IF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL PARTICLE_ALLOCATE(C_FILE_LOC)

  ! -------------------------------------------------------------------
  ! Initialize
  ! -------------------------------------------------------------------
  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  IF ( ifourier == 1 ) THEN         ! For Poisson solver
    CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  END IF

  IF ( iread_flow == 1 ) THEN       ! We need array space
    CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  END IF

  CALL FI_PROFILES_INITIALIZE(wrk1d)  ! Initialize thermodynamic quantities

  y_aux(:) = 0                        ! Reduced vertical grid
  DO ij = 1,jmax_aux*opt_block
    is = (ij-1)/opt_block + 1
    y_aux(is) = y_aux(is) + y(ij,1)/M_REAL(opt_block)
  END DO

  IF ( iread_flow == 1 ) THEN
    u => q(:,1)
    v => q(:,2)
    w => q(:,3)
  END IF

  ! ###################################################################
  ! Postprocess given list of files
  ! ###################################################################
  DO it = 1,itime_size
    itime = itime_vec(it)

    WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
    CALL TLAB_WRITE_ASCII(lfile, sRes)

    IF ( iread_scal == 1 ) THEN
      WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
      CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)
    END IF

    IF ( iread_flow == 1 ) THEN
      WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
      CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, inb_flow,i0, isize_wrk3d, q, wrk3d)
    END IF

    CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

    ! -------------------------------------------------------------------
    ! Calculate intermittency
    ! -------------------------------------------------------------------
    IF      ( opt_cond == 1 ) THEN ! External file
      WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname)); params_size = 2
      CALL IO_READ_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
      igate_size = INT(params(2))

      IF ( opt_main == 2 ) rtime = params(1)

    ELSE IF ( opt_cond > 1 ) THEN
      opt_cond_scal = 1
      IF ( imixture == MIXT_TYPE_AIRWATER .OR. imixture == MIXT_TYPE_AIRWATER_LINEAR ) THEN
        opt_cond_scal = inb_scal_array
      END IF

      CALL TLAB_WRITE_ASCII(lfile,'Calculating partition...')
      CALL FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
          imax,jmax,kmax, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)

      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE_INT1(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, gate, wrk1d)
      END IF

    END IF

    ! -------------------------------------------------------------------
    ! Type of averages
    ! -------------------------------------------------------------------
    SELECT CASE ( opt_main )

      ! ###################################################################
      ! Conventional statistics
      ! ###################################################################
    CASE ( 1 )
      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,9), txc(1,1),txc(1,2), txc(1,4), wrk1d,wrk2d,wrk3d)
      END IF

      IF ( icalc_scal == 1 ) THEN
        DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
          txc(1:isize_field,6) = txc(1:isize_field,9) ! Pass the pressure in tmp6
          CALL AVG_SCAL_XZ(is, q,s, s(1,is), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
        END DO

        ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
        IF ( buoyancy%TYPE == EQNS_BOD_QUADRATIC   .OR. buoyancy%TYPE == EQNS_BOD_BILINEAR    .OR. &
            imixture == MIXT_TYPE_AIRWATER        .OR. imixture == MIXT_TYPE_AIRWATER_LINEAR ) THEN
          IF ( buoyancy%TYPE == EQNS_EXPLICIT ) THEN
            CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,7))
          ELSE
            wrk1d(1:jmax,1) = C_0_R
            CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,7), wrk1d)
          END IF
          dummy = C_1_R/froude
          txc(1:isize_field,7) = txc(1:isize_field,7) *dummy

          txc(1:isize_field,6) = txc(1:isize_field,9) ! Pass the pressure in tmp6
          CALL AVG_SCAL_XZ(is, q,s, txc(1,7), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)

        END IF

        IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
          IF ( imixture == MIXT_TYPE_AIRWATER ) THEN
            is = is + 1
            CALL THERMO_ANELASTIC_THETA_L(imax,jmax,kmax, s, epbackground,pbackground, txc(1,7))
            !                 CALL THERMO_ANELASTIC_STATIC_CONSTANTCP(imax,jmax,kmax, s, epbackground, txc(1,7))
            txc(1:isize_field,6) = txc(1:isize_field,9) ! Pass the pressure in tmp6
            CALL AVG_SCAL_XZ(is, q,s, txc(1,7), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
          END IF
        END IF

      END IF

      ! Lagrange Liquid and Liquid without diffusion
      IF ( icalc_part == 1 ) THEN
        IF ( ilagrange == LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange == LAG_TYPE_BIL_CLOUD_4 ) THEN
          WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
          CALL IO_READ_PARTICLE(fname, l_g, l_q)

          l_txc = C_1_R; ! We want density
          CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,7), wrk2d,wrk3d)

          txc(:,7) = txc(:,7) + C_SMALL_R
          idummy = inb_part - 3 ! # scalar properties solved in the lagrangian
          DO is = inb_scal_array +1 +1, inb_scal_array+1 +idummy
            sbg(is)%mean = C_1_R; sbg(is)%delta = C_0_R; sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)
            l_txc(:,1)=l_q(:,3+is-inb_scal_array-1) !!! DO WE WANT l_txc(:,is) ???
            CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,8), wrk2d,wrk3d)
            txc(:,8) = txc(:,8) /txc(:,7)
            sbg(is)%mean  = sbg(inb_scal_array)%mean
            sbg(is)%delta = sbg(inb_scal_array)%delta
            sbg(is)%ymean = sbg(inb_scal_array)%ymean
            txc(1:isize_field,6) = txc(1:isize_field,9) ! Pass the pressure in tmp6
            CALL AVG_SCAL_XZ(is, q,s, txc(1,8), &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
          END DO
        END IF
      END IF

      IF ( icalc_flow == 1 ) THEN
        txc(1:isize_field,3) = txc(1:isize_field,9) ! Pass the pressure in tmp3
        CALL AVG_FLOW_XZ(q,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
            txc(1,7),txc(1,8),txc(1,9), mean, wrk1d,wrk2d,wrk3d)
      END IF

      ! ###################################################################
      ! Partition of field
      ! ###################################################################
    CASE ( 2 )
      DO is = 1,igate_size
        WRITE(varname(is),*) is; varname(is) = 'Partition'//TRIM(ADJUSTL(varname(is)))
      END DO
      WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
      CALL INTER_N_XZ(fname, itime,rtime, imax,jmax,kmax, igate_size, varname, gate, y, mean)

      IF ( opt_cond > 1 ) THEN ! write only if the gate information has not been read
        WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname))
        params(1) = rtime; params(2) = M_REAL(igate_size); params_size = 2
        CALL IO_WRITE_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)

        DO is = 1,igate_size
          gate_level = INT(is,KIND=1)
          CALL BOUNDARY_LOWER_INT1(imax,jmax,kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1,2))
          DO k = 1,kmax ! rearranging
            ij = (k-1) *imax +1
            surface(1:imax,is           ,k) = wrk2d(ij:ij+imax-1,1)
          END DO
          CALL BOUNDARY_UPPER_INT1(imax,jmax,kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1,2))
          DO k = 1,kmax ! rearranging
            ij = (k-1) *imax +1
            surface(1:imax,is+igate_size,k) = wrk2d(ij:ij+imax-1,1)
          END DO
        END DO
        varname = ''
        WRITE(fname,*) itime; fname = 'envelopesJ.'//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_ENVELOPES, fname, varname, surface, io_sizes, wrk3d)

      END IF

      ! ###################################################################
      ! Momentum equation
      ! ###################################################################
    CASE ( 3 )
      WRITE(fname,*) itime; fname='avgMom'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      ifield = ifield+1; vars(ifield)%field => q(:,1);   vars(ifield)%tag = 'U'
      ifield = ifield+1; vars(ifield)%field => q(:,3);   vars(ifield)%tag = 'W'

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'Uy'
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'Uyy'
      CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,1), txc(1,2), txc(1,1), wrk2d,wrk3d)
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'Wy'
      ifield = ifield+1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'Wyy'
      CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,3), txc(1,4), txc(1,3), wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'VU)y'
      txc(1:isize_field,6) = q(1:isize_field,2) *q(1:isize_field,1)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,6), txc(1,5), wrk3d, wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,6); vars(ifield)%tag = 'VUy'
      txc(1:isize_field,6) = q(1:isize_field,2) *txc(1:isize_field,1)
      ifield = ifield+1; vars(ifield)%field => txc(:,7); vars(ifield)%tag = 'UUx'
      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,7), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,7) = q(1:isize_field,1) *txc(1:isize_field,7)
      ifield = ifield+1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'WUz'
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,1), txc(1,8), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,8) = q(1:isize_field,3) *txc(1:isize_field,8)

      ifield = ifield+1; vars(ifield)%field => txc(:,9); vars(ifield)%tag = 'WV)y'
      txc(1:isize_field,10) = q(1:isize_field,2) *q(1:isize_field,3)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,10), txc(1,9), wrk3d, wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,10); vars(ifield)%tag = 'VWy'
      txc(1:isize_field,10) = q(1:isize_field,2) *txc(1:isize_field,3)
      ifield = ifield+1; vars(ifield)%field => txc(:,11); vars(ifield)%tag = 'UWx'
      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,3), txc(1,11), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,11) = q(1:isize_field,1) *txc(1:isize_field,11)
      ifield = ifield+1; vars(ifield)%field => txc(:,12); vars(ifield)%tag = 'WWz'
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,12), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,12) = q(1:isize_field,3) *txc(1:isize_field,12)

      ! ###################################################################
      ! Main variables
      ! ###################################################################
    CASE ( 4 )
      WRITE(fname,*) itime; fname='avgMain'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      ifield = ifield+1; vars(ifield)%field => u(:); vars(ifield)%tag = 'U'
      ifield = ifield+1; vars(ifield)%field => v(:); vars(ifield)%tag = 'V'
      ifield = ifield+1; vars(ifield)%field => w(:); vars(ifield)%tag = 'W'

      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
        ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'P'
      ELSE
        ifield = ifield+1; vars(ifield)%field => q(:,5); vars(ifield)%tag = 'R'
        ifield = ifield+1; vars(ifield)%field => q(:,6); vars(ifield)%tag = 'P'
        ifield = ifield+1; vars(ifield)%field => q(:,7); vars(ifield)%tag = 'T'
      END IF

      IF ( icalc_scal == 1 ) THEN
        DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
          ifield = ifield+1; vars(ifield)%field => s(:,is); WRITE(vars(ifield)%tag,*) is; vars(ifield)%tag = 'Scalar'//TRIM(ADJUSTL(vars(ifield)%tag))
        END DO
      END IF

      ! ###################################################################
      ! Enstrophy equation
      ! ###################################################################
    CASE ( 5 )
      WRITE(fname,*) itime; fname='avgW2'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      ! result vector in txc4, txc5, txc6
      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        IF ( buoyancy%TYPE == EQNS_NONE ) THEN
          txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
        ELSE
          IF ( buoyancy%TYPE == EQNS_EXPLICIT ) THEN
            CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, wrk3d)
          ELSE
            wrk1d(1:jmax,1) = C_0_R
            CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
          END IF
          DO ij = 1,isize_field
            s(ij,1) = wrk3d(ij)*buoyancy%vector(2)
          END DO

          CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(:,4) =-txc(:,4)
          txc(:,5) = C_0_R
          CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,6), wrk3d, wrk2d,wrk3d)
        END IF

      ELSE
        CALL FI_VORTICITY_BAROCLINIC(imax,jmax,kmax, q(1,5),q(1,6), txc(1,4), txc(1,3),txc(1,7), wrk2d,wrk3d)
      END IF

      CALL FI_CURL(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3), txc(1,7), wrk2d,wrk3d)
      txc(1:isize_field,8) = txc(1:isize_field,1)*txc(1:isize_field,4) &
          + txc(1:isize_field,2)*txc(1:isize_field,5) + txc(1:isize_field,3)*txc(1:isize_field,6)

      CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, u,v,w, txc(1,1),&
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

      CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, u,v,w, txc(1,2),&
          txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
      txc(1:isize_field,2) = visc *txc(1:isize_field,2)

      CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4),txc(1,5), wrk2d,wrk3d)  ! Enstrophy
      CALL FI_INVARIANT_P(imax,jmax,kmax, u,v,w, txc(1,4),txc(1,5), wrk2d,wrk3d)  ! Dilatation

      txc(1:isize_field,6) = txc(1:isize_field,4) *txc(1:isize_field,3) ! -w^2 div(u)
      txc(1:isize_field,5) = txc(1:isize_field,1) /txc(1:isize_field,3) ! production rate
      txc(1:isize_field,4) = LOG(txc(1:isize_field,3))                  ! ln(w^2)

      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'EnstrophyW_iW_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'LnEnstrophyW_iW_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'ProductionW_iW_jS_ij'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'DiffusionNuW_iLapW_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,6); vars(ifield)%tag = 'DilatationMsW_iW_iDivU'
      ifield = ifield +1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'Baroclinic'
      ifield = ifield +1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'RateAN_iN_jS_ij'

      ! ###################################################################
      ! Strain equation
      ! ###################################################################
    CASE ( 6 )
      WRITE(fname,*) itime; fname='avgS2'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
        CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, u,v,w, txc(1,1), &
            txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      ELSE
        CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, u,v,w, q(1,6), &
            txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      END IF
      txc(1:isize_field,1) = C_2_R *txc(1:isize_field,2)

      CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, u,v,w, &
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
      txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

      CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, u,v,w, &
          txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk2d,wrk3d)
      txc(1:isize_field,3) = C_2_R *visc *txc(1:isize_field,3)

      CALL FI_STRAIN(imax,jmax,kmax, u,v,w, txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      txc(1:isize_field,4) = C_2_R *txc(1:isize_field,4)
      txc(1:isize_field,5) = LOG( txc(1:isize_field,4) )

      ifield = ifield +1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'Strain2S_ijS_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'LnStrain2S_ijS_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'ProductionMs2S_ijS_jkS_ki'
      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'DiffusionNuS_ijLapS_ij'
      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'Pressure2S_ijP_ij'

      ! ###################################################################
      ! Scalar gradient equation
      ! ###################################################################
    CASE ( 7 )
      WRITE(fname,*) itime; fname='avgG2'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s, u,v,w, &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

      CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s, &   ! array u used as auxiliar
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),u, wrk2d,wrk3d)
      txc(1:isize_field,2) = txc(1:isize_field,2) *visc /schmidt(inb_scal)

      CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,3), txc(1,4), wrk2d,wrk3d)
      txc(1:isize_field,5) = txc(1:isize_field,1) /txc(1:isize_field,3)
      txc(1:isize_field,4) = LOG(txc(1:isize_field,3))

      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'GradientG_iG_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'LnGradientG_iG_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'ProductionMsG_iG_jS_ij'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'DiffusionNuG_iLapG_i'
      ifield = ifield +1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'StrainAMsN_iN_jS_ij'

      ! ###################################################################
      ! Velocity gradient invariants
      ! ###################################################################
    CASE ( 8 )
      WRITE(fname,*) itime; fname='avgInv'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_INVARIANT_R(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_INVARIANT_Q(imax,jmax,kmax, u,v,w, txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
      CALL FI_INVARIANT_P(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4), wrk2d,wrk3d)

      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'InvariantP'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'InvariantQ'
      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'InvariantR'

      ! ###################################################################
      ! Scalar gradient components
      ! ###################################################################
    CASE ( 9 )
      WRITE(fname,*) itime; fname='avgGi'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,1), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,2), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,3), wrk3d, wrk2d,wrk3d)
      DO ij = 1,isize_field                       ! Angles; s array is overwritten to save space
        dummy = txc(ij,2) /SQRT(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
        txc(ij,4) = ASIN(dummy)                  ! with Oy
        s(ij,1)  = ATAN2(txc(ij,3),txc(ij,1))    ! with Ox in plane xOz
      END DO

      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'GradientX'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'GradientY'
      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'GradientZ'
      ifield = ifield +1; vars(ifield)%field => s(:,1);   vars(ifield)%tag = 'Theta'
      ifield = ifield +1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'Phi'

      ! ###################################################################
      ! eigenvalues of rate-of-strain tensor
      ! ###################################################################
    CASE ( 10 )
      WRITE(fname,*) itime; fname='avgEig'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_STRAIN_TENSOR(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))

      ifield = ifield +1; vars(ifield)%field => txc(:,7); vars(ifield)%tag = 'Lambda1'
      ifield = ifield +1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'Lambda2'
      ifield = ifield +1; vars(ifield)%field => txc(:,9); vars(ifield)%tag = 'Lambda3'

      ! ###################################################################
      ! eigenframe of rate-of-strain tensor
      ! ###################################################################
    CASE ( 11 )
      WRITE(fname,*) itime; fname='avgCos'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_STRAIN_TENSOR(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))  ! txc7-txc9
      CALL FI_TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1,1), txc(1,7))   ! txc1-txc6

      CALL FI_CURL(imax,jmax,kmax, u,v,w, txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk2d,wrk3d)
      DO ij = 1,isize_field                                             ! local direction cosines of vorticity vector
        dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
        u(ij) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
        v(ij) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
        eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
        eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
        eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
        w(ij) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
      END DO

      ifield = ifield +1; vars(ifield)%field => u; vars(ifield)%tag = 'cos(w,lambda1)'
      ifield = ifield +1; vars(ifield)%field => v; vars(ifield)%tag = 'cos(w,lambda2)'
      ifield = ifield +1; vars(ifield)%field => w; vars(ifield)%tag = 'cos(w,lambda3)'

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,7), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,8), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,9), wrk3d, wrk2d,wrk3d)
      DO ij = 1,isize_field                                             ! local direction cosines of scalar gradient vector
        dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
        cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
        cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
        eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
        eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
        eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
        cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3
      END DO

      ifield = ifield +1; vars(ifield)%field => txc(:,7); vars(ifield)%tag = 'cos(G,lambda1)'
      ifield = ifield +1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'cos(G,lambda2)'
      ifield = ifield +1; vars(ifield)%field => txc(:,9); vars(ifield)%tag = 'cos(G,lambda3)'

      ! ###################################################################
      ! longitudinal velocity derivatives
      ! ###################################################################
    CASE ( 12 )
      WRITE(fname,*) itime; fname='avgDer'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, txc(1,1), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, txc(1,2), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, txc(1,3), wrk3d, wrk2d,wrk3d)

      ifield = ifield +1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'dudx'
      ifield = ifield +1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'dvdy'
      ifield = ifield +1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'dwdz'

      ! ###################################################################
      ! Vertical fluxes
      ! ###################################################################
    CASE ( 13 )
      ifield = 0
      WRITE(fname,*) itime; fname='avgFluxY'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, txc(:,1), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, txc(:,2), wrk3d, wrk2d,wrk3d)
      txc(:,1) = ( txc(:,1) + txc(:,2) ) *visc
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'tauyx'

      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, txc(:,2), wrk3d, wrk2d,wrk3d)
      txc(:,2) =   txc(:,2) *C_2_R       *visc
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'tauyy'

      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, txc(:,3), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, txc(:,4), wrk3d, wrk2d,wrk3d)
      txc(:,3) = ( txc(:,3) + txc(:,4) ) *visc
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'tauyz'

      DO is = 1,inb_scal_array
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(:,is), txc(:,3+is), wrk3d, wrk2d,wrk3d)
        txc(:,3+is) =   txc(:,3+is) *visc /schmidt(inb_scal)
        ifield = ifield+1; vars(ifield)%field => txc(:,3+is); WRITE(vars(ifield)%tag,*) is; vars(ifield)%tag = 'tauy'//TRIM(ADJUSTL(vars(ifield)%tag))
      END DO

      u = u*v
      ifield = ifield+1; vars(ifield)%field => u; vars(ifield)%tag = 'vu'
      ! I need v below
      ifield = ifield+1; vars(ifield)%field => v; vars(ifield)%tag = 'vv'
      w = w*v
      ifield = ifield+1; vars(ifield)%field => w; vars(ifield)%tag = 'vw'
      DO is = 1,inb_scal_array
        s(:,is) = s(:,is)*v
        ifield = ifield+1; vars(ifield)%field => s(:,is); WRITE(vars(ifield)%tag,*) is; vars(ifield)%tag = 'v'//TRIM(ADJUSTL(vars(ifield)%tag))
      END DO
      v = v*v ! I need v above for the scalar fluxes

      ! ###################################################################
      ! Hydrostatic and dynamic pressure
      ! ###################################################################
    CASE ( 14 )
      WRITE(fname,*) itime; fname='avgP'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'P'

      q = C_0_R
      CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,2), txc(1,3),txc(1,4), txc(1,5), wrk1d,wrk2d,wrk3d)
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'Psta'

      txc(:,3) = txc(:,1) - txc(:,2)
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'Pdyn'

      ! ###################################################################
      ! Dissipation
      ! ###################################################################
    CASE ( 15 )
      WRITE(fname,*) itime; fname='avgEps'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_DISSIPATION(i1,imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'Eps'

      ! ###################################################################
      ! Covariances among scalars
      ! ###################################################################
    CASE ( 16 )
      WRITE(fname,*) itime; fname='avgSiCov'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, s(1,1))
      CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, s(1,2))

      txc(1:isize_field,1) = s(1:isize_field,1)   *s(1:isize_field,2)
      txc(1:isize_field,2) = txc(1:isize_field,1) *s(1:isize_field,1)
      txc(1:isize_field,3) = txc(1:isize_field,1) *s(1:isize_field,2)

      ifield = 0
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 's1s2'
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 's1s2s1'
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 's1s2s2'

      ! ###################################################################
      ! Potential vorticity
      ! ###################################################################
    CASE ( 17 )
      WRITE(fname,*) itime; fname='avgPV'//TRIM(ADJUSTL(fname))
      CALL TLAB_WRITE_ASCII(lfile,'Computing '//TRIM(ADJUSTL(fname))//'...')
      ifield = 0

      CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
      txc(1:isize_field,6) = txc(1:isize_field,1)*txc(1:isize_field,1) &
          + txc(1:isize_field,2)*txc(1:isize_field,2) &
          + txc(1:isize_field,3)*txc(1:isize_field,3) ! Enstrophy
      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,1) =                       txc(1:isize_field,1)*txc(1:isize_field,4)
      txc(1:isize_field,5) =                       txc(1:isize_field,4)*txc(1:isize_field,4) ! norm grad b
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,2)*txc(1:isize_field,4)
      txc(1:isize_field,5) = txc(1:isize_field,5) +txc(1:isize_field,4)*txc(1:isize_field,4) ! norm grad b
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
      txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,3)*txc(1:isize_field,4)
      txc(1:isize_field,5) = txc(1:isize_field,5) +txc(1:isize_field,4)*txc(1:isize_field,4) ! norm grad b

      txc(1:isize_field,5) = SQRT( txc(1:isize_field,5) +C_SMALL_R)
      txc(1:isize_field,6) = SQRT( txc(1:isize_field,6) +C_SMALL_R)
      txc(1:isize_field,2) = txc(1:isize_field,1) /( txc(1:isize_field,5) *txc(1:isize_field,6) ) ! Cosine of angle between 2 vectors

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'PV'
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'Cos'

    END SELECT

    IF ( opt_main > 2 ) THEN
      IF ( nfield < ifield ) THEN ! Check
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
        CALL TLAB_STOP(DNS_ERROR_WRKSIZE)
      END IF

      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
      END IF

      CALL AVG_N_XZ(fname, itime, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_order, vars, gate_level,gate, y_aux, mean)

    END IF

  END DO

  CALL TLAB_STOP(0)
END PROGRAM AVERAGES

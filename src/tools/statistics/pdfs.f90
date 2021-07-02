#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "PDFS"

PROGRAM PDFS

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : ifile,efile,lfile ,gfile, tag_flow,tag_scal
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  =  20
  TINTEGER, PARAMETER :: igate_size_max =   8
  TINTEGER, PARAMETER :: params_size_max =  2

  TREAL,      ALLOCATABLE, SAVE   :: x(:,:), y(:,:), z(:,:)
  TREAL,      ALLOCATABLE, TARGET :: s(:,:), q(:,:), txc(:,:)
  TREAL,      ALLOCATABLE         :: wrk2d(:,:)
  TREAL,      ALLOCATABLE         :: pdf(:), y_aux(:), wrk1d(:), wrk3d(:)
  INTEGER(1), ALLOCATABLE         :: gate(:)
  TYPE(pointers_dt)               :: vars(16)

  ! -------------------------------------------------------------------
  ! Local variables
  ! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, bakfile
  CHARACTER*64 str, line

  TINTEGER opt_main, opt_block, opt_bins(2)
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER nfield, ifield, isize_wrk3d, ij, is, bcs(2,2), isize_pdf
  TREAL dummy, eloc1, eloc2, eloc3, cos1, cos2, cos3
  TINTEGER jmax_aux, iread_flow, iread_scal, ierr, idummy
  TINTEGER ibc(16)
  TREAL vmin(16), vmax(16)
  LOGICAL reduce_data

  ! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  TREAL gate_threshold(igate_size_max)
  INTEGER(1) gate_level

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TINTEGER opt_vec(iopt_size_max)
  TREAL opt_vec2(iopt_size_max)

  TINTEGER params_size
  TREAL params(params_size_max)

  !########################################################################
  !########################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(ifile)

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
  opt_bins  =16

  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_INTEGER(sRes, iopt_size, opt_vec)

  IF ( sRes == '-1' ) THEN
#ifdef USE_MPI
#else
    WRITE(*,*) 'Option ?'
    WRITE(*,*) ' 1. Main variables'
    WRITE(*,*) ' 2. Scalar gradient G_iG_i/2 equation'
    WRITE(*,*) ' 3. Enstrophy W_iW_i/2 equation'
    WRITE(*,*) ' 4. Strain 2S_ijS_ij/2 equation'
    WRITE(*,*) ' 5. Velocity gradient invariants'
    WRITE(*,*) ' 6. Flamelet equation'
    WRITE(*,*) ' 7. Joint enstrophy and strain'
    WRITE(*,*) ' 9. Joint scalar and scalar gradient'
    WRITE(*,*) '10. Scalar gradient components'
    WRITE(*,*) '11. Eigenvalues of rate-of-strain tensor'
    WRITE(*,*) '12. Eigenframe of rate-of-strain tensor'
    WRITE(*,*) '13. Longitudinal velocity derivatives'
    WRITE(*,*) '14. Potential vorticity'
    WRITE(*,*) '15. Analysis of B and V'
    READ(*,*) opt_main

    WRITE(*,*) 'Planes block size ?'
    READ(*,*) opt_block

    WRITE(*,*) 'Gate level to be used ?'
    READ(*,*) gate_level

    WRITE(*,*) 'Number of PDF bins ?'
    READ(*,'(A)') sRes
    idummy = 2
    CALL LIST_INTEGER(sRes, idummy, opt_bins)

#endif
  ELSE
    opt_main = opt_vec(1)
    IF ( iopt_size >= 2 ) opt_block = opt_vec(2)
    IF ( iopt_size >= 3 ) gate_level= INT(opt_vec(3),KIND=1)
    IF ( iopt_size >= 4 ) opt_bins  = opt_vec(4:5)

  END IF

  IF ( opt_main < 0 ) THEN ! Check
    CALL IO_WRITE_ASCII(efile, 'PDFS. Missing input [ParamPdfs] in dns.ini.')
    CALL DNS_STOP(DNS_ERROR_INVALOPT)
  END IF

  IF ( opt_block < 1 ) THEN
    CALL IO_WRITE_ASCII(efile, 'PDFS. Invalid value of opt_block.')
    CALL DNS_STOP(DNS_ERROR_INVALOPT)
  END IF

  ! -------------------------------------------------------------------
  iread_flow = 0
  iread_scal = 0
  IF      ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN; inb_txc = 6;
  ELSE;                                                                                         inb_txc = 1
  END IF
  IF ( ifourier == 1 ) inb_txc = MAX(inb_txc,1)
  nfield     = 2

  SELECT CASE ( opt_main )
  CASE( 1 )
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,6)
    nfield = 4 +inb_scal_array;     isize_pdf = opt_bins(1) +2
    IF ( imode_eqns == DNS_EQNS_INTERNAL .OR. imode_eqns == DNS_EQNS_TOTAL ) nfield = nfield +2
  CASE( 2 ) ! Scalar gradient equation
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,6)
    nfield = 5;                     isize_pdf = opt_bins(1) +2
  CASE( 3 ) ! Enstrophy equation
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,8)
    nfield = 7;                     isize_pdf = opt_bins(1) +2
  CASE( 4 ) ! Strain equation
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,8)
    nfield = 5;                     isize_pdf = opt_bins(1) +2
  CASE( 5 ) ! Invariants
    iread_flow = 1; inb_txc = MAX(inb_txc,6)
    nfield = 3;                     isize_pdf = opt_bins(1)*opt_bins(2) +2 +2*opt_bins(1)
  CASE( 6 ) ! Chi-flamelet
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,6)
    nfield = 2;                     isize_pdf = opt_bins(1) +2
  CASE( 7 )
    iread_flow = 1; inb_txc = MAX(inb_txc,4)
    nfield = 2;                     isize_pdf = opt_bins(1)*opt_bins(2) +2 +2*opt_bins(1)
  CASE( 9 )
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,3)
    nfield = 2;                     isize_pdf = opt_bins(1)*opt_bins(2) +2 +2*opt_bins(1)
  CASE( 10 )
    iread_scal = 1;                 inb_txc = MAX(inb_txc,4)
    nfield = 5;                     isize_pdf = opt_bins(1)*opt_bins(2) +2 +2*opt_bins(1)
  CASE( 11 ) ! eigenvalues
    iread_flow = 1; inb_txc = MAX(inb_txc,9)
    nfield = 3;                     isize_pdf = opt_bins(1) +2
  CASE( 12 ) ! eigenframe
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,9)
    nfield = 6;                     isize_pdf = opt_bins(1) +2
  CASE( 13 ) ! longitudinal velocity derivatives
    iread_flow = 1; inb_txc = MAX(inb_txc,3)
    nfield = 3;                     isize_pdf = opt_bins(1) +2
  CASE( 14 ) ! potential vorticity
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,6)
    nfield = 2;                     isize_pdf = opt_bins(1) +2
  CASE( 15 ) ! joint s and v
    iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,8)
    nfield = 2;                     isize_pdf = opt_bins(1)*opt_bins(2) +2 +2*opt_bins(1)
  END SELECT

  ! -------------------------------------------------------------------
  ! Defining gate levels for conditioning
  ! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
  opt_cond_relative = 0
  igate_size    = 0

  IF ( gate_level /= 0 ) THEN
#include "dns_read_partition.h"
    IF ( opt_cond > 1 ) inb_txc = MAX(inb_txc,5)
  END IF

  ! -------------------------------------------------------------------
  ! Allocating memory space
  ! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  inb_wrk2d = MAX(inb_wrk2d,4)
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(gate(isize_field))

  ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
  jmax_aux = g(2)%size/opt_block
  ALLOCATE(y_aux(jmax_aux)) ! Reduced vertical grid

  ! Space for the min and max of sampling variable at opt_bins+1,opt_bins+2
  ! Space for the 3D pdf at jmax_aux+1
  ALLOCATE( pdf( isize_pdf *(jmax_aux+1) *nfield ) )

  isize_wrk3d = MAX(isize_field,isize_txc_field)
#include "dns_alloc_arrays.h"

  ! -------------------------------------------------------------------
  ! Initialize
  ! -------------------------------------------------------------------
#include "dns_read_grid.h"

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

  ibc(1:nfield) = 1

  ! ###################################################################
  ! Postprocess given list of files
  ! ###################################################################
  DO it=1, itime_size
    itime = itime_vec(it)

    WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
    CALL IO_WRITE_ASCII(lfile, sRes)

    IF ( iread_scal == 1 ) THEN
      WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
      CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)

      IF      ( imixture == MIXT_TYPE_AIRWATER .AND. damkohler(3) <= C_0_R ) THEN ! Calculate q_l
        CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2),s(1,1), epbackground,pbackground)
      ELSE IF ( imixture == MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
      END IF

    END IF

    IF ( iread_flow == 1 ) THEN
      WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
      CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)
    END IF

    ! -------------------------------------------------------------------
    ! Calculate intermittency
    ! -------------------------------------------------------------------
    IF      ( opt_cond == 1 ) THEN ! External file
      WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname)); params_size = 2
      CALL IO_READ_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
      igate_size = INT(params(2))

    ELSE IF ( opt_cond > 1 ) THEN
      opt_cond_scal = 1 ! Scalar field to use for the conditioning
      IF ( imixture == MIXT_TYPE_AIRWATER .OR. imixture == MIXT_TYPE_AIRWATER_LINEAR ) THEN
        opt_cond_scal = inb_scal_array
      END IF

      CALL IO_WRITE_ASCII(lfile,'Calculating partition...')
      CALL FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
          imax,jmax,kmax, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)

      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE_INT1(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, gate, wrk1d)
      END IF

    END IF

    ! -------------------------------------------------------------------
    ! Type of PDFs
    ! -------------------------------------------------------------------
    ifield = 0
    reduce_data = .TRUE.

    SELECT CASE ( opt_main )

      ! ###################################################################
      ! Main variable 2D-PDF
      ! ###################################################################
    CASE( 1 )
      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)

      ELSE
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density

        CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1),txc(1,2),txc(1,3), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1))

      END IF

      ifield = ifield+1; vars(ifield)%field => q(:,1);   vars(ifield)%tag = 'u'
      ifield = ifield+1; vars(ifield)%field => q(:,2);   vars(ifield)%tag = 'v'
      ifield = ifield+1; vars(ifield)%field => q(:,3);   vars(ifield)%tag = 'w'
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'p'
      IF ( imode_eqns == DNS_EQNS_INTERNAL .OR. imode_eqns == DNS_EQNS_TOTAL ) THEN
        ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'r'
        ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 't'
      END IF

      DO is = 1,inb_scal_array
        ifield = ifield+1; vars(ifield)%field => s(:,is); vars(ifield)%tag = 's'
        WRITE(str,*) is; vars(ifield)%tag=TRIM(ADJUSTL(vars(ifield)%tag))//TRIM(ADJUSTL(str))
      END DO

      DO is = 1,ifield ! In case we want same interval for all heights
        IF ( ibc(is) == 0 ) CALL MINMAX(imax,jmax,kmax, vars(is)%field, vmin(is),vmax(is))
      END DO

      ! ###################################################################
      ! Scalar gradient equation
      ! ###################################################################
    CASE ( 2 )
      CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient equation...')

      CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s, & ! array q used as auxiliar
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),q(1,1), wrk2d,wrk3d)
      txc(1:isize_field,2) = txc(1:isize_field,2) *visc /schmidt(inb_scal)
      CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,3), txc(1,4), wrk2d,wrk3d)
      txc(1:isize_field,5) = txc(1:isize_field,1) /txc(1:isize_field,3)
      txc(1:isize_field,4) = LOG(txc(1:isize_field,3))

      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'GiGi';                 ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'LnGiGi';               ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'ProductionMsGiGjSij';  ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'DiffusionNuGiLapGi';   ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'StrainAMsNiNjSij';     ibc(ifield) = 2

      ! ###################################################################
      ! Enstrophy equation
      ! ###################################################################
    CASE ( 3 )
      CALL IO_WRITE_ASCII(lfile,'Computing enstrophy equation...')

      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        IF ( buoyancy%TYPE == EQNS_NONE ) THEN
          txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
        ELSE
          IF ( buoyancy%TYPE == EQNS_EXPLICIT ) THEN
            CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, wrk3d)
          ELSE
            wrk1d(1:jmax) = C_0_R
            CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
          END IF
          s(1:isize_field,1) = wrk3d(1:isize_field) *buoyancy%vector(2)

          CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(:,4) =-txc(:,4)
          txc(:,5) = C_0_R
          CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,6), wrk3d, wrk2d,wrk3d)
        END IF

      ELSE
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density

        CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1)) ! pressure in txc1
        ! result vector in txc4, txc5, txc6
        CALL FI_VORTICITY_BAROCLINIC(imax,jmax,kmax, txc(1,2),txc(1,1), txc(1,4), txc(1,3),txc(1,7), wrk2d,wrk3d)
      END IF
      ! result vector in txc1, txc2, txc3
      CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), txc(1,7), wrk2d,wrk3d)
      ! scalar product, store in txc8
      txc(1:isize_field,8) = txc(1:isize_field,1)*txc(1:isize_field,4) + txc(1:isize_field,2)*txc(1:isize_field,5) + txc(1:isize_field,3)*txc(1:isize_field,6)

      CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

      CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2),&
          txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
      txc(1:isize_field,2) = visc *txc(1:isize_field,2)

      CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,3), txc(1,4),txc(1,5), wrk2d,wrk3d)

      CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,4), txc(1,5), wrk2d,wrk3d)

      txc(1:isize_field,5) = txc(1:isize_field,4) *txc(1:isize_field,3) ! -w^2 div(u)
      txc(1:isize_field,4) = txc(1:isize_field,1) /txc(1:isize_field,3) ! production rate
      txc(1:isize_field,6) = LOG(txc(1:isize_field,3))                  ! ln(w^2)

      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'WiWi';                 ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,6); vars(ifield)%tag = 'LnWiWi';               ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'ProductionWiWjSij';    ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'DiffusionNuWiLapWi';   ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,5); vars(ifield)%tag = 'DilatationMsWiWiDivU'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'Baroclinic';           ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'RateANiNjSij';         ibc(ifield) = 2

      ! ###################################################################
      ! Strain equation
      ! ###################################################################
    CASE ( 4 )
      CALL IO_WRITE_ASCII(lfile,'Computing strain equation...')

      IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)

      ELSE
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, q(1,4), wrk3d)! energy
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, q(1,5), wrk3d)! density
        CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), q(1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), q(1,7), txc(1,1))         ! pressure in txc1

      END IF
      CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      txc(1:isize_field,1) = C_2_R *txc(1:isize_field,2)

      CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
      txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

      CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
          txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk2d,wrk3d)
      txc(1:isize_field,3) = C_2_R *visc *txc(1:isize_field,3)

      CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      txc(1:isize_field,4) = C_2_R *txc(1:isize_field,4)
      txc(1:isize_field,5) = LOG( txc(1:isize_field,4) )

      ifield = ifield+1; vars(1)%field => txc(:,4); vars(ifield)%tag = '2SijSij';                  ibc(ifield) = 2
      ifield = ifield+1; vars(2)%field => txc(:,5); vars(ifield)%tag = 'Ln2SijSij';                ibc(ifield) = 2
      ifield = ifield+1; vars(3)%field => txc(:,2); vars(ifield)%tag = 'ProductionMs2SijSjkS_ki';  ibc(ifield) = 2
      ifield = ifield+1; vars(4)%field => txc(:,3); vars(ifield)%tag = 'DiffusionNuSijLapSij';     ibc(ifield) = 2
      ifield = ifield+1; vars(5)%field => txc(:,1); vars(ifield)%tag = 'Pressure2SijPij';          ibc(ifield) = 2

      ! ###################################################################
      ! Velocity gradient invariants
      ! ###################################################################
    CASE ( 5 )
      CALL IO_WRITE_ASCII(lfile,'Computing velocity gradient invariants...')

      CALL FI_INVARIANT_R(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_INVARIANT_Q(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
      CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,3), txc(1,4), wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'InvP'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'InvQ'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'InvR'; ibc(ifield) = 2

      IF (  jmax_aux*opt_block /= g(2)%size .AND. reduce_data ) THEN ! I already need it here
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
        reduce_data = .FALSE.
      END IF

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.RQ'
      CALL PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, y_aux, txc(1,1),txc(1,2), pdf, wrk2d )

      ! ###################################################################
      ! Chi flamelet equation PDF
      ! ###################################################################
    CASE ( 6 )
      CALL IO_WRITE_ASCII(lfile,'Computing flamelet equation...')

      CALL FI_STRAIN_A(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'StrainAGiGi'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'StrainA';     ibc(ifield) = 2

      ! ###################################################################
      ! Joint PDF W^2 and 2S^2
      ! ###################################################################
    CASE ( 7 )
      CALL IO_WRITE_ASCII(lfile,'Computing enstrophy-strain pdf...')

      CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
      CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3), txc(1,4), wrk2d,wrk3d)
      txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

      IF (  jmax_aux*opt_block /= g(2)%size .AND. reduce_data ) THEN ! I already need it here
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
        reduce_data = .FALSE.
      END IF

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.WS'
      CALL PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, y_aux, txc(1,1),txc(1,2), pdf, wrk2d )

      ! ###################################################################
      ! Joint PDF Scalar and Scalar Gradient
      ! ###################################################################
    CASE ( 9 )
      CALL IO_WRITE_ASCII(lfile,'Computing scalar-scalar--gradient pdf...')

      CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,1), txc(1,2), wrk2d,wrk3d)
      txc(1:isize_field,2) = LOG(txc(1:isize_field,1))

      ifield = ifield+1; vars(1)%field => s(:,1);   vars(ifield)%tag = 's';      ibc(ifield) = 1
      ifield = ifield+1; vars(2)%field => txc(:,1); vars(ifield)%tag = 'GiGi';   ibc(ifield) = 2
      ifield = ifield+1; vars(3)%field => txc(:,2); vars(ifield)%tag = 'LnGiGi'; ibc(ifield) = 3

      IF (  jmax_aux*opt_block /= g(2)%size .AND. reduce_data ) THEN ! I already need it here
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
        reduce_data = .FALSE.
      END IF

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.SLnG'
      CALL PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, s(1,1),txc(1,2), y_aux, pdf, wrk2d )

      WRITE(fname,*) itime; fname='cavgGiGi'//TRIM(ADJUSTL(fname))
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          1, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,1), y_aux, pdf, wrk1d)

      WRITE(fname,*) itime; fname='cavgLnGiGi'//TRIM(ADJUSTL(fname))
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          1, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,2), y_aux, pdf, wrk1d)

      ! ###################################################################
      ! Scalar gradient components
      ! ###################################################################
    CASE ( 10 )
      CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient components...')

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,1), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,2), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,3), wrk3d, wrk2d,wrk3d)
      ! Angles; s array is overwritten to save space
      DO ij = 1,isize_field
        dummy = txc(ij,2)/SQRT(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
        txc(ij,4) = ASIN(dummy)                 ! with Oy
        s(ij,1)  = ATAN2(txc(ij,3),txc(ij,1))  ! with Ox in plane xOz
      END DO

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'Gx';     ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'Gy';     ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'Gz';     ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => s(:,1);   vars(ifield)%tag = 'Gtheta'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,4); vars(ifield)%tag = 'Gphi';   ibc(ifield) = 2

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.GphiS'
      CALL PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, s(1,1),txc(1,4), y_aux, pdf, wrk2d )

      ! ###################################################################
      ! eigenvalues of rate-of-strain tensor
      ! ###################################################################
    CASE ( 11 )
      CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues of Sij...')

      CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7)) ! txc7-txc9

      ifield = ifield+1; vars(ifield)%field => txc(:,7); vars(ifield)%tag = 'Lambda1'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'Lambda2'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,9); vars(ifield)%tag = 'Lambda3'; ibc(ifield) = 2

      ! ###################################################################
      ! eigenframe of rate-of-strain tensor
      ! ###################################################################
    CASE ( 12 )
      CALL IO_WRITE_ASCII(lfile,'Computing eigenframe of Sij...')

      CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
      CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7)) ! txc7-txc9
      CALL FI_TENSOR_EIGENFRAME(imax,jmax,kmax,  txc(1,1), txc(1,7)) ! txc1-txc6
      CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk2d,wrk3d)

      DO ij = 1,isize_field ! local direction cosines of vorticity vector
        dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
        q(ij,1) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
        q(ij,2) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
        eloc1   = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
        eloc2   = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
        eloc3   = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
        q(ij,3) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
      END DO

      ifield = ifield+1; vars(ifield)%field => q(:,1); vars(ifield)%tag = 'cos(w,lambda1)'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => q(:,2); vars(ifield)%tag = 'cos(w,lambda2)'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => q(:,3); vars(ifield)%tag = 'cos(w,lambda3)'; ibc(ifield) = 2

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,7), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,8), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,9), wrk3d, wrk2d,wrk3d)

      DO ij = 1,isize_field ! local direction cosines of scalar gradient vector
        dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
        cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
        cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
        eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
        eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
        eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
        cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3
      END DO

      ifield = ifield+1; vars(ifield)%field => txc(:,7); vars(ifield)%tag = 'cos(G,lambda1)'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,8); vars(ifield)%tag = 'cos(G,lambda2)'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,9); vars(ifield)%tag = 'cos(G,lambda3)'; ibc(ifield) = 2

      ! ###################################################################
      ! Longitudinal velocity derivatives
      ! ###################################################################
    CASE ( 13 )
      CALL IO_WRITE_ASCII(lfile,'Computing longitudinal velocity derivatives...')

      CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,1), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), q(1,2), txc(1,2), wrk3d, wrk2d,wrk3d)
      CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,3), wrk3d, wrk2d,wrk3d)

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'Sxx'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'Syy'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,3); vars(ifield)%tag = 'Szz'; ibc(ifield) = 2

      ! ###################################################################
      ! Potential vorticity
      ! ###################################################################
    CASE ( 14 )
      CALL IO_WRITE_ASCII(lfile,'Computing potntial vorticity...')

      CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
      txc(1:isize_field,6) = txc(1:isize_field,1)**2 +txc(1:isize_field,2)**2 +txc(1:isize_field,3)**2 ! Enstrophy
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

      txc(1:isize_field,1) = txc(1:isize_field,1)*txc(1:isize_field,1) ! Squared of the potential voticity
      txc(1:isize_field,1) = LOG(txc(1:isize_field,1)+C_SMALL_R)

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'LnPotentialEnstrophy'; ibc(ifield) = 2
      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'CosPotentialEnstrophy';ibc(ifield) = 2

      ! ###################################################################
      ! Analysis of B and V
      ! ###################################################################
    CASE ( 15 )
      CALL IO_WRITE_ASCII(lfile,'Computing analysis of B and V...')

      ifield = ifield+1; vars(ifield)%field => txc(:,1); vars(ifield)%tag = 'b'; ibc(ifield) = 1
      IF ( buoyancy%TYPE == EQNS_EXPLICIT ) THEN
        CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,1))
      ELSE
        wrk1d(1:jmax) = C_0_R
        CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
      END IF
      dummy =  C_1_R /froude
      txc(1:isize_field,1) = txc(1:isize_field,1) *dummy

      ifield = ifield+1; vars(ifield)%field => txc(:,2); vars(ifield)%tag = 'v'; ibc(ifield) = 1
      txc(1:isize_field,2) = q(1:isize_field,2)

      ! I need tmp1 w/o reduction to calculate derivatives
      CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), txc(1,1),txc(1,5), txc(1,6), wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,4), txc(1,6), wrk2d,wrk3d)
      CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), txc(1,1),txc(1,3), txc(1,6), wrk2d,wrk3d)
      txc(1:isize_field,3) = txc(1:isize_field,3) +txc(1:isize_field,4) +txc(1:isize_field,5)

      ! -------------------------------------------------------------------
      IF (  jmax_aux*opt_block /= g(2)%size .AND. reduce_data ) THEN ! I already need it here
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
        reduce_data = .FALSE.
      END IF

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.bv'
      CALL PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), y_aux, pdf, wrk2d )

      ! -------------------------------------------------------------------
      WRITE(fname,*) itime; fname='cavgB'//TRIM(ADJUSTL(fname))
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,1), y_aux, pdf, wrk1d)

      WRITE(fname,*) itime; fname='cavgBii'//TRIM(ADJUSTL(fname))
      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
      END IF
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

      ! -------------------------------------------------------------------
      WRITE(fname,*) itime; fname='cavgU'//TRIM(ADJUSTL(fname))
      txc(1:isize_field,3) = q(1:isize_field,1)
      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
      END IF
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

      WRITE(fname,*) itime; fname='cavgW'//TRIM(ADJUSTL(fname))
      txc(1:isize_field,3) = q(1:isize_field,3)
      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
      END IF
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

      WRITE(fname,*) itime; fname='cavgVii'//TRIM(ADJUSTL(fname))
      CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), q(1,2),txc(1,5), txc(1,6), wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), q(1,2),txc(1,4), txc(1,6), wrk2d,wrk3d)
      CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), q(1,2),txc(1,3), txc(1,6), wrk2d,wrk3d)
      txc(1:isize_field,3) = txc(1:isize_field,3) +txc(1:isize_field,4) +txc(1:isize_field,5)
      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
      END IF
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

      ! -------------------------------------------------------------------
      bbackground = C_0_R
      CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,3), txc(1,4),txc(1,5), txc(1,6), wrk1d,wrk2d,wrk3d)
      CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,3),txc(1,4), wrk3d, wrk2d,wrk3d)
      IF (  jmax_aux*opt_block /= g(2)%size ) THEN
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
        CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,4), wrk1d)
      END IF

      WRITE(fname,*) itime; fname='cavgP'//TRIM(ADJUSTL(fname))
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

      WRITE(fname,*) itime; fname='cavgPy'//TRIM(ADJUSTL(fname))
      CALL CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, txc(1,4), y_aux, pdf, wrk1d)
      fname = TRIM(ADJUSTL(fname))//'.bv'
      CALL CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,4), y_aux, pdf, wrk2d )

    END SELECT

    ! ###################################################################
    IF ( ifield > 0 ) THEN
      IF ( nfield < ifield ) THEN
        CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
        CALL DNS_STOP(DNS_ERROR_WRKSIZE)
      END IF

      IF (  jmax_aux*opt_block /= g(2)%size .AND. reduce_data ) THEN
        DO is = 1,ifield
          CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
        END DO
      END IF

      WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
      CALL PDF1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
          ifield, opt_bins(1), ibc, vmin,vmax,vars, gate_level,gate, y_aux, pdf, wrk1d)

    END IF

  END DO

  CALL DNS_STOP(0)
END PROGRAM PDFS

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "PDFS"

PROGRAM PDFS

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile,lfile ,gfile, tag_flow,tag_scal
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  ! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  =  20
  TINTEGER, PARAMETER :: igate_size_max =   8
  TINTEGER, PARAMETER :: params_size_max =  2

  ! Arrays declarations
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE   :: x,y,z
  TREAL,      DIMENSION(:,:), ALLOCATABLE, TARGET :: s, q, txc
  TREAL,      DIMENSION(:,:), ALLOCATABLE         :: wrk2d
  TREAL,      DIMENSION(:),   ALLOCATABLE         :: pdf, y_aux, wrk1d, wrk3d
  INTEGER(1), DIMENSION(:),   ALLOCATABLE         :: gate

  TYPE(pointers_dt), DIMENSION(16) :: DATA

  ! -------------------------------------------------------------------
  ! Local variables
  ! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, inifile, bakfile, varname(16)
  CHARACTER*64 str, line

  TINTEGER opt_main, opt_block, opt_bins(2)
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER nfield, isize_wrk3d, ij, is, bcs(2,2), isize_pdf
  TREAL dummy, eloc1, eloc2, eloc3, cos1, cos2, cos3
  TINTEGER jmax_aux, iread_flow, iread_scal, ierr, idummy
  TINTEGER ibc(16)
  TREAL amin(16), amax(16)

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

  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)

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

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_INTEGER(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
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
     IF ( iopt_size .GE. 2 ) opt_block = opt_vec(2)
     IF ( iopt_size .GE. 3 ) gate_level= INT(opt_vec(3),KIND=1)
     IF ( iopt_size .GE. 4 ) opt_bins  = opt_vec(4:5)

  ENDIF

  IF ( opt_main .LT. 0 ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'PDFS. Missing input [ParamPdfs] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  IF ( opt_block .LT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile, 'PDFS. Invalid value of opt_block.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  ! -------------------------------------------------------------------
  iread_flow = 0
  iread_scal = 0
  IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN; inb_txc = 6;
  ELSE;                                                                                             inb_txc = 1
  ENDIF
  IF ( ifourier .EQ. 1 ) inb_txc = MAX(inb_txc,1)
  nfield     = 2

  SELECT CASE ( opt_main )
  CASE( 1 )
     iread_scal = 1; iread_flow = 1; inb_txc = MAX(inb_txc,6)
     nfield = 4 +inb_scal_array;     isize_pdf = opt_bins(1) +2
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) nfield = nfield +2
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

  IF ( gate_level .NE. 0 ) THEN
#include "dns_read_partition.h"
     IF ( opt_cond .GT. 1 ) inb_txc = MAX(inb_txc,5)
  ENDIF

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

  IF ( ifourier .EQ. 1 ) THEN         ! For Poisson solver
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( iread_flow .EQ. 1 ) THEN       ! We need array space
     CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

  CALL FI_PROFILES_INITIALIZE(wrk1d)  ! Initialize thermodynamic quantities

  y_aux(:) = 0                        ! Reduced vertical grid
  DO ij = 1,jmax_aux*opt_block
     is = (ij-1)/opt_block + 1
     y_aux(is) = y_aux(is) + y(ij,1)/M_REAL(opt_block)
  ENDDO

  ibc(1:nfield) = 1

  ! ###################################################################
  ! Postprocess given list of files
  ! ###################################################################
  DO it=1, itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)

     IF ( iread_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)

        IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Calculate q_l
           CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2),s(1,1), epbackground,pbackground)
        ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
           CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
        ENDIF

     ENDIF

     IF ( iread_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)
     ENDIF

     ! -------------------------------------------------------------------
     ! Calculate intermittency
     ! -------------------------------------------------------------------
     IF      ( opt_cond .EQ. 1 ) THEN ! External file
        WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname)); params_size = 2
        CALL IO_READ_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
        igate_size = INT(params(2))

     ELSE IF ( opt_cond .GT. 1 ) THEN
        opt_cond_scal = 1 ! Scalar field to use for the conditioning
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
           opt_cond_scal = inb_scal_array
        ENDIF

        CALL IO_WRITE_ASCII(lfile,'Calculating partition...')
        CALL FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
             imax,jmax,kmax, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE_INT1(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, gate, wrk1d)
        ENDIF

     ENDIF

     ! -------------------------------------------------------------------
     ! Type of PDFs
     ! -------------------------------------------------------------------
     SELECT CASE ( opt_main )

        ! ###################################################################
        ! Main variable 2D-PDF
        ! ###################################################################
     CASE( 1 )
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)

        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density

           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1),txc(1,2),txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1))

        ENDIF

        nfield = 0
        nfield = nfield+1; DATA(nfield)%field => q(:,1);   varname(nfield) = 'u'
        nfield = nfield+1; DATA(nfield)%field => q(:,2);   varname(nfield) = 'v'
        nfield = nfield+1; DATA(nfield)%field => q(:,3);   varname(nfield) = 'w'
        nfield = nfield+1; DATA(nfield)%field => txc(:,1); varname(nfield) = 'p'
        IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           nfield = nfield+1; DATA(nfield)%field => txc(:,2); varname(nfield) = 'r'
           nfield = nfield+1; DATA(nfield)%field => txc(:,3); varname(nfield) = 't'
        ENDIF

        DO is = 1,inb_scal_array
           nfield = nfield+1; DATA(nfield)%field => s(:,is); varname(nfield) = 's'
           WRITE(str,*) is; varname(nfield)=TRIM(ADJUSTL(varname(nfield)))//TRIM(ADJUSTL(str))
        ENDDO

        DO is = 1,nfield ! In case we want same interval for all heights
           IF ( ibc(is) .EQ. 0 ) CALL MINMAX(imax,jmax,kmax, DATA(is)%field, amin(is),amax(is))
        ENDDO

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Scalar gradient equation
        ! ###################################################################
     CASE ( 2 )
        CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s, & ! array q used as auxiliar
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),q(1,1), wrk2d,wrk3d)
        dummy = visc /schmidt(inb_scal)
        txc(1:isize_field,2) = dummy *txc(1:isize_field,2)
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,3), txc(1,4), wrk2d,wrk3d)
        txc(1:isize_field,5) = txc(1:isize_field,1) /txc(1:isize_field,3)
        txc(1:isize_field,4) = LOG(txc(1:isize_field,3))

        DATA(1)%field => txc(:,3); varname(1) = 'GiGi';                 ibc(1) = 2
        DATA(2)%field => txc(:,4); varname(2) = 'LnGiGi';               ibc(2) = 2
        DATA(3)%field => txc(:,1); varname(3) = 'ProductionMsGiGjSij';  ibc(3) = 2
        DATA(4)%field => txc(:,2); varname(4) = 'DiffusionNuGiLapGi';   ibc(4) = 2
        DATA(5)%field => txc(:,5); varname(5) = 'StrainAMsNiNjSij';     ibc(5) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Enstrophy equation
        ! ###################################################################
     CASE ( 3 )
        CALL IO_WRITE_ASCII(lfile,'Computing baroclinic term...')
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           IF ( buoyancy%TYPE .EQ. EQNS_NONE ) THEN
              txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
           ELSE
              IF ( buoyancy%TYPE .EQ. EQNS_EXPLICIT ) THEN
                 CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, wrk3d)
              ELSE
                 wrk1d(1:jmax) = C_0_R
                 CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
              ENDIF
              DO ij = 1,isize_field
                 s(ij,1) = wrk3d(ij)*buoyancy%vector(2)
              ENDDO

              CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,4), wrk3d, wrk2d,wrk3d)
              txc(:,4) =-txc(:,4)
              txc(:,5) = C_0_R
              CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,6), wrk3d, wrk2d,wrk3d)
           ENDIF

        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density

           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1)) ! pressure in txc1
           ! result vector in txc4, txc5, txc6
           CALL FI_VORTICITY_BAROCLINIC(imax,jmax,kmax, txc(1,2),txc(1,1), txc(1,4), txc(1,3),txc(1,7), wrk2d,wrk3d)
        ENDIF
        ! result vector in txc1, txc2, txc3
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), txc(1,7), wrk2d,wrk3d)
        ! scalar product, store in txc8
        DO ij = 1,isize_field
           txc(ij,8) = txc(ij,1)*txc(ij,4) + txc(ij,2)*txc(ij,5) + txc(ij,3)*txc(ij,6)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
        CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
        CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2),&
             txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        txc(1:isize_field,2) = visc *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy...')
        CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,3), txc(1,4),txc(1,5), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing dilatation term...')
        CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,4), txc(1,5), wrk2d,wrk3d)

        txc(1:isize_field,5) = txc(1:isize_field,4) *txc(1:isize_field,3) ! -w^2 div(u)
        txc(1:isize_field,4) = txc(1:isize_field,1) /txc(1:isize_field,3) ! production rate
        txc(1:isize_field,6) = LOG(txc(1:isize_field,3))                  ! ln(w^2)

        DATA(1)%field => txc(:,3); varname(1) = 'WiWi';                 ibc(1) = 2
        DATA(2)%field => txc(:,6); varname(2) = 'LnWiWi';               ibc(2) = 2
        DATA(3)%field => txc(:,1); varname(3) = 'ProductionWiWjSij';    ibc(3) = 2
        DATA(4)%field => txc(:,2); varname(4) = 'DiffusionNuWiLapWi';   ibc(4) = 2
        DATA(5)%field => txc(:,5); varname(5) = 'DilatationMsWiWiDivU'; ibc(5) = 2
        DATA(6)%field => txc(:,8); varname(6) = 'Baroclinic';           ibc(6) = 2
        DATA(7)%field => txc(:,4); varname(7) = 'RateANiNjSij';         ibc(7) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Strain equation
        ! ###################################################################
     CASE ( 4 )
        CALL IO_WRITE_ASCII(lfile,'Computing strain pressure...')
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)

        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, q(1,4), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, q(1,5), wrk3d)! density
           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), q(1,7), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), q(1,7), txc(1,1))         ! pressure in txc1

        ENDIF
        CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,1) = C_2_R *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
        CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
        CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
             txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk2d,wrk3d)
        txc(1:isize_field,3) = C_2_R *visc *txc(1:isize_field,3)

        CALL IO_WRITE_ASCII(lfile,'Computing strain...')
        CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,4) = C_2_R *txc(1:isize_field,4)
        txc(1:isize_field,5) = LOG( txc(1:isize_field,4) )

        DATA(1)%field => txc(:,4); varname(1) = '2SijSij';                  ibc(1) = 2
        DATA(2)%field => txc(:,5); varname(2) = 'Ln2SijSij';                ibc(2) = 2
        DATA(3)%field => txc(:,2); varname(3) = 'ProductionMs2SijSjkS_ki';  ibc(3) = 2
        DATA(4)%field => txc(:,3); varname(4) = 'DiffusionNuSijLapSij';     ibc(4) = 2
        DATA(5)%field => txc(:,1); varname(5) = 'Pressure2SijPij';          ibc(5) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Velocity gradient invariants
        ! ###################################################################
     CASE ( 5 )
        CALL FI_INVARIANT_R(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL FI_INVARIANT_Q(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,3), txc(1,4), wrk2d,wrk3d)

        DATA(1)%field => txc(:,3); varname(1) = 'InvP'; ibc(1) = 2
        DATA(2)%field => txc(:,2); varname(2) = 'InvQ'; ibc(2) = 2
        DATA(3)%field => txc(:,1); varname(3) = 'InvR'; ibc(3) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        WRITE(fname,*) itime; fname='jpdf'//TRIM(ADJUSTL(fname))//'.RQ'
        CALL PDF2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, y_aux, DATA(1)%field,DATA(2)%field, pdf, wrk2d )

        ! ###################################################################
        ! Chi flamelet equation PDF
        ! ###################################################################
     CASE ( 6 )
        CALL FI_STRAIN_A(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        DATA(1)%field => txc(:,1); varname(1) = 'StrainAGiGi'; ibc(1) = 2
        DATA(2)%field => txc(:,2); varname(2) = 'StrainA';     ibc(2) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Joint PDF W^2 and 2S^2
        ! ###################################################################
     CASE ( 7 )
        CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
        CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3), txc(1,4), wrk2d,wrk3d)
        txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

        WRITE(fname,*) itime; fname='jpdf'//TRIM(ADJUSTL(fname))//'.WS'
        ! We pass ny=1 and it only calculates 3D pdfs (twice, but it allows us to reuse existing routines)
        CALL PDF2V(fname, imax*jmax, 1, kmax, opt_bins, y_aux, txc(1,1),txc(1,2), pdf, wrk2d )

        ! ###################################################################
        ! Joint PDF Scalar and Scalar Gradient
        ! ###################################################################
     CASE ( 9 )
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,1), txc(1,2), wrk2d,wrk3d)
        txc(1:isize_field,1) = LOG(txc(1:isize_field,1))

        WRITE(fname,*) itime; fname='jpdf'//TRIM(ADJUSTL(fname))//'.SG'
        ! We pass ny=1 and it only calculates 3D pdfs (twice, but it allows us to reuse existing routines)
        CALL PDF2V(fname, imax*jmax, 1, kmax, opt_bins, y_aux, s(1,1),txc(1,1), pdf, wrk2d )

        ! ###################################################################
        ! Scalar gradient components
        ! ###################################################################
     CASE ( 10 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,3), wrk3d, wrk2d,wrk3d)
        ! Angles; s array is overwritten to save space
        DO ij = 1,imax*jmax*kmax
           dummy = txc(ij,2)/SQRT(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
           txc(ij,4) = ASIN(dummy)                 ! with Oy
           s(ij,1)  = ATAN2(txc(ij,3),txc(ij,1))  ! with Ox in plane xOz
        ENDDO

        DATA(1)%field => txc(:,1); varname(1) = 'Gx';     ibc(1) = 2
        DATA(2)%field => txc(:,2); varname(2) = 'Gy';     ibc(2) = 2
        DATA(3)%field => txc(:,3); varname(3) = 'Gz';     ibc(3) = 2
        DATA(4)%field => s(:,1);   varname(4) = 'Gtheta'; ibc(4) = 2
        DATA(5)%field => txc(:,4); varname(5) = 'Gphi';   ibc(5) = 2

        WRITE(fname,*) itime; fname='jpdfGi'//TRIM(ADJUSTL(fname))
        ! We pass ny=1 and it only calculates 3D pdfs (twice, but it allows us to reuse existing routines)
        CALL PDF2V(fname, imax*jmax, 1, kmax, opt_bins, y_aux, s(1,1),txc(1,4), pdf, wrk2d )

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! eigenvalues of rate-of-strain tensor
        ! ###################################################################
     CASE ( 11 )
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7)) ! txc7-txc9

        DATA(1)%field => txc(:,7); varname(1) = 'Lambda1'; ibc(1) = 2
        DATA(2)%field => txc(:,8); varname(2) = 'Lambda2'; ibc(2) = 2
        DATA(3)%field => txc(:,9); varname(3) = 'Lambda3'; ibc(3) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! eigenframe of rate-of-strain tensor
        ! ###################################################################
     CASE ( 12 )
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7)) ! txc7-txc9
        CALL FI_TENSOR_EIGENFRAME(imax,jmax,kmax,  txc(1,1), txc(1,7)) ! txc1-txc6
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk2d,wrk3d)

        ! local direction cosines of vorticity vector
        DO ij = 1,imax*jmax*kmax
           dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           q(ij,1) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           q(ij,2) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1   = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2   = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3   = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           q(ij,3) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        ENDDO

        DATA(1)%field => q(:,1); varname(1) = 'cos(w,lambda1)'; ibc(1) = 2
        DATA(2)%field => q(:,2); varname(2) = 'cos(w,lambda2)'; ibc(2) = 2
        DATA(3)%field => q(:,3); varname(3) = 'cos(w,lambda3)'; ibc(3) = 2

        ! local direction cosines of scalar gradient vector
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...') ! txc7-txc9
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,7), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,8), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,9), wrk3d, wrk2d,wrk3d)

        DO ij = 1,imax*jmax*kmax
           dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
           txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3
        ENDDO

        DATA(4)%field => txc(:,7); varname(4) = 'cos(G,lambda1)'; ibc(4) = 2
        DATA(5)%field => txc(:,8); varname(5) = 'cos(G,lambda2)'; ibc(5) = 2
        DATA(6)%field => txc(:,9); varname(6) = 'cos(G,lambda3)'; ibc(6) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Longitudinal velocity derivatives
        ! ###################################################################
     CASE ( 13 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), q(1,2), txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,3), wrk3d, wrk2d,wrk3d)

        DATA(1)%field => txc(:,1); varname(1) = 'Sxx'; ibc(1) = 2
        DATA(2)%field => txc(:,2); varname(2) = 'Syy'; ibc(2) = 2
        DATA(3)%field => txc(:,3); varname(3) = 'Szz'; ibc(3) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! Potential vorticity
        ! ###################################################################
     CASE ( 14 )
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

        DATA(1)%field => txc(:,1); varname(1) = 'LnPotentialEnstrophy'; ibc(1) = 2
        DATA(2)%field => txc(:,2); varname(2) = 'CosPotentialEnstrophy';ibc(2) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)

        ! ###################################################################
        ! joint scalar and vertical velocity
        ! ###################################################################
     CASE ( 15 )
        nfield = 0

        nfield = nfield+1; DATA(nfield)%field => txc(:,1); varname(nfield) = 'b'; ibc(nfield) = 1
        IF ( buoyancy%TYPE .EQ. EQNS_EXPLICIT ) THEN
           CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,1))
        ELSE
           wrk1d(1:jmax) = C_0_R
           CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
        ENDIF
        dummy =  C_1_R /froude
        txc(1:isize_field,1) = txc(1:isize_field,1) *dummy

        nfield = nfield+1; DATA(nfield)%field => txc(:,2); varname(nfield) = 'v'; ibc(nfield) = 1
        txc(1:isize_field,2) = q(1:isize_field,2)

        ! I need tmp1 w/o reduction to calculate derivatives
        CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), txc(1,1),txc(1,5), txc(1,6), wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,4), txc(1,6), wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), txc(1,1),txc(1,3), txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,3) = txc(1:isize_field,3) +txc(1:isize_field,4) +txc(1:isize_field,5)

        ! -------------------------------------------------------------------
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, DATA(is)%field, wrk1d)
           ENDDO
        ENDIF
        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, y_aux, pdf, wrk1d)
        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))//'.bv'
        CALL PDF2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), y_aux, pdf, wrk2d )

        ! -------------------------------------------------------------------
        WRITE(fname,*) itime; fname='cavgB'//TRIM(ADJUSTL(fname))
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,1), y_aux, pdf, wrk1d)

        WRITE(fname,*) itime; fname='cavgBii'//TRIM(ADJUSTL(fname))
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
        ENDIF
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

        ! -------------------------------------------------------------------
        WRITE(fname,*) itime; fname='cavgU'//TRIM(ADJUSTL(fname))
        txc(1:isize_field,3) = q(1:isize_field,1)
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
        ENDIF
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

        WRITE(fname,*) itime; fname='cavgW'//TRIM(ADJUSTL(fname))
        txc(1:isize_field,3) = q(1:isize_field,3)
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
        ENDIF
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

        WRITE(fname,*) itime; fname='cavgVii'//TRIM(ADJUSTL(fname))
        CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), q(1,2),txc(1,5), txc(1,6), wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), q(1,2),txc(1,4), txc(1,6), wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), q(1,2),txc(1,3), txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,3) = txc(1:isize_field,3) +txc(1:isize_field,4) +txc(1:isize_field,5)
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
        ENDIF
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

        ! -------------------------------------------------------------------
        bbackground = C_0_R
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,3), txc(1,4),txc(1,5), txc(1,6), wrk1d,wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,3),txc(1,4), wrk3d, wrk2d,wrk3d)
        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,3), wrk1d)
           CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, txc(1,4), wrk1d)
        ENDIF

        WRITE(fname,*) itime; fname='cavgP'//TRIM(ADJUSTL(fname))
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,3), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,3), y_aux, pdf, wrk2d )

        WRITE(fname,*) itime; fname='cavgPy'//TRIM(ADJUSTL(fname))
        CALL CAVG1V_N(fname, varname, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins(1), ibc, amin,amax,DATA, gate_level,gate, txc(1,4), y_aux, pdf, wrk1d)
        fname = TRIM(ADJUSTL(fname))//'.bv'
        CALL CAVG2V(fname, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1,1),txc(1,2), txc(1,4), y_aux, pdf, wrk2d )

     END SELECT
  ENDDO

  CALL DNS_END(0)

  STOP

END PROGRAM PDFS

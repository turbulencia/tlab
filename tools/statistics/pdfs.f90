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
#ifdef USE_MPI
#include "mpif.h"
#endif

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

  TYPE(pointers_dt), DIMENSION(16) :: data

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 varname(16)
  CHARACTER*64 str, line

  TINTEGER opt_main, opt_block, opt_bins, opt_bcs, opt_bins_2
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER nfield, isize_wrk3d, ij, is, bcs(2,2)
  TREAL diff, dummy
  TREAL eloc1, eloc2, eloc3, cos1, cos2, cos3
  TINTEGER jmax_aux, iread_flow, iread_scal, ierr, idummy
  TINTEGER npdf_size, ibc(16)
  TREAL amin(16), amax(16)

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

  inb_wrk2d = MAX(inb_wrk2d,4)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(gate(isize_field))

  ALLOCATE(y_aux(g(2)%size)) ! Reduced vertical grid

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Additional options
! -------------------------------------------------------------------
  opt_main  =-1 ! default values
  opt_block = 1
  gate_level  = 0
  opt_bins  = 16
  opt_bcs   = 1

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'Option ?'
     WRITE(*,*) ' 1. Main variables'
     WRITE(*,*) ' 2. Scalar'
     WRITE(*,*) ' 3. Scalar gradient G_iG_i/2 equation'
     WRITE(*,*) ' 4. Enstrophy W_iW_i/2 equation'
     WRITE(*,*) ' 5. Strain 2S_ijS_ij/2 equation'
     WRITE(*,*) ' 6. Velocity gradient invariants'
     WRITE(*,*) ' 7. Flamelet equation'
     WRITE(*,*) ' 9. Joint WS 3D-PDF'
     WRITE(*,*) '10. Joint scalar scalar-gradient 3D-PDF'
     WRITE(*,*) '11. Conditional scalar gradient G_iG_i 3D-PDF'
     WRITE(*,*) '12. Scalar gradient components'
     WRITE(*,*) '13. Gradient trajectories data'
     WRITE(*,*) '14. Eigenvalues of rate-of-strain tensor'
     WRITE(*,*) '15. Eigenframe of rate-of-strain tensor'
     WRITE(*,*) '16. Longitudinal velocity derivatives'
     WRITE(*,*) '17. Potential vorticity'
     READ(*,*) opt_main

     WRITE(*,*) 'Planes block size ?'
     READ(*,*) opt_block  

     WRITE(*,*) 'Gate level to be used ?'
     READ(*,*) gate_level

     WRITE(*,*) 'Number of PDF bins ?'
     READ(*,*) opt_bins
     WRITE(*,*) 'Local interval (1) or homogeneous (0) ?'
     READ(*,*) opt_bcs

     IF ( opt_main .EQ. 10 ) THEN ! conditional 3D-PDFs on scalars
        WRITE(*,*) 'Number of scalar bins ?'
        READ(*,*) opt_bins_2
     ENDIF

#endif
  ELSE
     opt_main = INT(opt_vec(1))
     IF ( iopt_size .GE. 2 ) opt_block = INT(opt_vec(2))
     IF ( iopt_size .GE. 3 ) gate_level  = INT(opt_vec(3),KIND=1)
     IF ( iopt_size .GE. 4 ) opt_bins  = INT(opt_vec(4))
     IF ( iopt_size .GE. 5 ) opt_bcs   = INT(opt_vec(5))
     IF ( opt_main .EQ. 10 ) THEN
        IF ( iopt_size .GE. 6 ) opt_bins_2  = INT(opt_vec(6))
     ENDIF

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
  ELSE;                                                                                             inb_txc = 1; ENDIF
  IF ( ifourier .EQ. 1 ) inb_txc = MAX(inb_txc,1)
  nfield     = 2 

  SELECT CASE ( opt_main )

  CASE( 1 )
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,6)
     nfield = 4 +inb_scal
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) nfield = nfield +2
  CASE( 2 )
     iread_scal = 1
     nfield = inb_scal
  CASE( 3 ) ! Scalar gradient equation
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,6)
     nfield = 5
  CASE( 4 ) ! Enstrophy equation
     iread_flow = 1
     iread_scal = 1
     inb_txc = MAX(inb_txc,8)
     nfield = 7
  CASE( 5 ) ! Strain equation
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,8)
     nfield = 5
  CASE( 6 ) ! Invariants
     iread_flow = 1
     inb_txc = MAX(inb_txc,6)
     nfield = 3
  CASE( 7 ) ! Chi-flamelet 
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,6)
     nfield = 2
  CASE( 9 )
     iread_flow = 1
     inb_txc = MAX(inb_txc,4)
  CASE( 10 )
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,3)
  CASE( 11 )
     iread_scal = 1
     iread_flow = 1
     inb_txc = MAX(inb_txc,3)
  CASE( 12 )
     iread_scal = 1
     inb_txc = MAX(inb_txc,4)
     nfield = 5
  CASE( 13 )
!     inb_txc = MAX(inb_txc,2)
!     nfield = 2
     inb_txc = MAX(inb_txc,1)
     nfield = 1
  CASE( 14 ) ! eigenvalues
     iread_flow = 1
     inb_txc = MAX(inb_txc,9)
     nfield = 3
  CASE( 15 ) ! eigenframe
     iread_flow = 1
     iread_scal = 1
     inb_txc = MAX(inb_txc,9)
     nfield = 6
  CASE( 16 ) ! longitudinal velocity derivatives
     iread_flow = 1
     iread_scal = 0
     inb_txc = MAX(inb_txc,3)
     nfield = 3
  CASE (17 ) ! potential vorticity
     nfield = 1
     inb_txc = MAX(inb_txc,4)
     iread_flow = 1
     iread_scal = 1
  END SELECT

! -------------------------------------------------------------------
! Defining gate levels for conditioning
! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
  opt_cond_scal = 1
  opt_cond_relative = 0
  igate_size    = 0

  IF ( gate_level .NE.0 ) THEN

#include "dns_read_partition.h"

  ENDIF

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
! in case g(2)%size is not divisible by opt_block, drop the upper most planes
  jmax_aux = g(2)%size/opt_block
  
! Space for the min and max of sampling variable at opt_bins+1,opt_bins+2
! Space for the 3D pdf at jmax_aux+1
  npdf_size = (jmax_aux+1) *(opt_bins+2) *nfield 

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------      
  ALLOCATE(pdf(npdf_size))

! -------------------------------------------------------------------
  isize_txc   = isize_txc_field*inb_txc
  isize_wrk3d = MAX(isize_field,isize_txc_field)

#include "dns_alloc_arrays.h"

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! ------------------------------------------------------------------------
! Define size of blocks
! ------------------------------------------------------------------------
  y_aux(:) = 0
  do ij = 1,jmax
     is = (ij-1)/opt_block + 1 
     y_aux(is) = y_aux(is) + y(ij,1)/M_REAL(opt_block)
  enddo

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( iread_flow .EQ. 1 ) THEN ! We need array space
     CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

! -------------------------------------------------------------------
! Initialize thermodynamic quantities
! -------------------------------------------------------------------
  CALL FI_PROFILES_INITIALIZE(wrk1d)
  
! ###################################################################
! Calculating statistics
! ###################################################################
  DO it=1, itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)

     IF ( iread_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)
     ENDIF

     IF ( iread_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)
     ENDIF
     IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
     ELSE;                                  diff = visc/schmidt(inb_scal); ENDIF

! -------------------------------------------------------------------
! Calculate intermittency
! -------------------------------------------------------------------
     IF ( opt_cond .GT. 0 ) THEN
        IF      ( opt_cond .EQ. 1 ) THEN ! External file
           WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname)); params_size = 2
           CALL IO_READ_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
           igate_size = INT(params(2))
           
        ELSE
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
              opt_cond_scal = inb_scal_array
           ENDIF
           
           CALL IO_WRITE_ASCII(lfile,'Calculating partition...')
           CALL FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                imax,jmax,kmax, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)
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
        nfield = nfield+1; data(nfield)%field => q(:,1);   varname(nfield) = 'Pu'
        nfield = nfield+1; data(nfield)%field => q(:,2);   varname(nfield) = 'Pv'
        nfield = nfield+1; data(nfield)%field => q(:,3);   varname(nfield) = 'Pw'
        nfield = nfield+1; data(nfield)%field => txc(:,1); varname(nfield) = 'Pp'
        IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           nfield = nfield+1; data(nfield)%field => txc(:,2); varname(nfield) = 'Pr'
           nfield = nfield+1; data(nfield)%field => txc(:,3); varname(nfield) = 'Pt'
        ENDIF

        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal
              nfield = nfield+1; data(nfield)%field => s(:,is); varname(nfield) = 'P'
              WRITE(str,*) is; varname(nfield)=TRIM(ADJUSTL(varname(nfield)))//TRIM(ADJUSTL(str))
           ENDDO
        ENDIF

        ibc(1:nfield) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Scalar 2D-PDF
! ###################################################################
     CASE ( 2 )
        nfield = 0
        DO is = 1,inb_scal
           nfield = nfield+1; data(is)%field => s(:,is); varname(is) = 'Scalar '
           WRITE(str,*) is; varname(is)=TRIM(ADJUSTL(varname(is)))//TRIM(ADJUSTL(str))
           IF ( opt_bcs .EQ. 0 ) THEN
              ibc(is) = 0
              CALL MINMAX(imax,jmax,kmax, s(1,is), amin(is),amax(is))
           ELSE
              ibc(is) = 2
           ENDIF
        ENDDO

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdfSc'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Scalar gradient equation
! ###################################################################
     CASE ( 3 )
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
        CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

! array u used as auxiliar
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
        CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s, &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),q(1,1), wrk2d,wrk3d)
        txc(1:isize_field,2) = diff *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,3), txc(1,4), wrk2d,wrk3d)
        txc(1:isize_field,5) = txc(1:isize_field,1) /txc(1:isize_field,3)
        txc(1:isize_field,4) = log(txc(1:isize_field,3))

        data(1)%field => txc(:,3); varname(1) = 'GradientG_iG_i'        ; ibc(1) = 2
        data(2)%field => txc(:,4); varname(2) = 'LnGradientG_iG_i'      ; ibc(2) = 2
        data(3)%field => txc(:,1); varname(3) = 'ProductionMsG_iG_jS_ij'; ibc(3) = 2
        data(4)%field => txc(:,2); varname(4) = 'DiffusionNuG_iLapG_i'  ; ibc(4) = 2
        data(5)%field => txc(:,5); varname(5) = 'StrainAMsN_iN_jS_ij'   ; ibc(5) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Enstrophy equation
! ###################################################################
     CASE ( 4 )
        CALL IO_WRITE_ASCII(lfile,'Computing baroclinic term...')
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           IF ( buoyancy%type .EQ. EQNS_NONE ) THEN
              txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
           ELSE
! calculate buoyancy vector along Oy
              IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
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
           
           CALL THERMO_CALORIC_TEMPERATURE&
                (imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE&
                (imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1)) ! pressure in txc1
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
        txc(1:isize_field,6) = log(txc(1:isize_field,3))                  ! ln(w^2)

        data(1)%field => txc(:,3); varname(1) = 'EnstrophyW_iW_i'     ;   ibc(1) = 2
        data(2)%field => txc(:,6); varname(2) = 'LnEnstrophyW_iW_i'   ;   ibc(2) = 2
        data(3)%field => txc(:,1); varname(3) = 'ProductionW_iW_jS_ij';   ibc(3) = 2
        data(4)%field => txc(:,2); varname(4) = 'DiffusionNuW_iLapW_i';   ibc(4) = 2
        data(5)%field => txc(:,5); varname(5) = 'DilatationMsW_iW_iDivU'; ibc(5) = 2
        data(6)%field => txc(:,8); varname(6) = 'Baroclinic';             ibc(6) = 2
        data(7)%field => txc(:,4); varname(7) = 'RateAN_iN_jS_ij'     ;   ibc(7) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Strain equation
! ###################################################################
     CASE ( 5 )
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
        txc(1:isize_field,5) = log( txc(1:isize_field,4) )
              
        data(1)%field => txc(:,4); varname(1) = 'Strain2S_ijS_i'           ; ibc(1) = 2
        data(2)%field => txc(:,5); varname(2) = 'LnStrain2S_ijS_i'         ; ibc(2) = 2
        data(3)%field => txc(:,2); varname(3) = 'ProductionMs2S_ijS_jkS_ki'; ibc(3) = 2
        data(4)%field => txc(:,3); varname(4) = 'DiffusionNuS_ijLapS_ij'   ; ibc(4) = 2
        data(5)%field => txc(:,1); varname(5) = 'Pressure2S_ijP_ij'        ; ibc(5) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Velocity gradient invariants
! ###################################################################
     CASE ( 6 )
        CALL IO_WRITE_ASCII(lfile,'Computing third invariant R...') ! txc1 contains R
        CALL FI_INVARIANT_R(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing second invariant Q...') ! txc2 contains Q
        CALL FI_INVARIANT_Q(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing first invariant P...')
        CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,3), txc(1,4), wrk2d,wrk3d)

        data(1)%field => txc(:,3); varname(1) = 'InvariantP'; ibc(1) = 2
        data(2)%field => txc(:,2); varname(2) = 'InvariantQ'; ibc(2) = 2
        data(3)%field => txc(:,1); varname(3) = 'InvariantR'; ibc(3) = 2

        WRITE(fname,*) itime; fname='jpdfRQ'//TRIM(ADJUSTL(fname))
        CALL JPDF3D(fname, i0, gate_level, i0, imax, jmax, kmax, i0, i0,&
             gate, txc(1,2), txc(1,1), opt_bins, opt_bins, wrk2d(1,1), wrk2d(1,2), wrk2d(1,3), wrk1d)

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Chi flamelet equation PDF
! ###################################################################
     CASE ( 7 )
        CALL FI_STRAIN_A(imax,jmax,kmax, s, q(1,1),q(1,2),q(1,3), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        data(1)%field => txc(:,1); varname(1) = 'StrainAG_iG_i'; ibc(1) = 2
        data(2)%field => txc(:,2); varname(2) = 'StrainA';       ibc(2) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Joint PDF W^2 and 2S^2
! ###################################################################
     CASE ( 9 )
! txc1 contains w^2
        CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
! in case we want to condition on w
        DO ij = 1,imax*jmax*kmax
           s(ij,1) = SQRT(txc(ij,1))
        ENDDO
! txc2 contains 2s^2
        CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,2), txc(1,3), txc(1,4), wrk2d,wrk3d)
        DO ij = 1,imax*jmax*kmax
           txc(ij,2) = C_2_R*txc(ij,2)
        ENDDO

        WRITE(fname,*) itime; fname='jpdfWS'//TRIM(ADJUSTL(fname))
        CALL JPDF3D(fname, i0, gate_level, i0, imax, jmax, kmax, i0, i0,&
             gate, txc(1,2), txc(1,1), opt_bins, opt_bins, wrk2d(1,1), wrk2d(1,2), wrk2d(1,3), wrk1d)

! ###################################################################
! Conditional scalar gradient 3D-PDFs 
! ###################################################################
     CASE ( 10 )
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,1), txc(1,2), wrk2d,wrk3d)
        varname(1) = 'ScalarGradient'
        txc(1:isize_field,1) = log(txc(1:isize_field,1))

        WRITE(fname,*) itime; fname='cpdf'//TRIM(ADJUSTL(fname))
        CALL CPDF3D_N(fname, varname, gate_level, i1, i1, rtime, imax, jmax, kmax,&
             i1, opt_bins_2, opt_bins, gate, s, txc, pdf, wrk1d)

! ###################################################################
! Joint PDF Scalar and Scalar Gradient 
! ###################################################################
     CASE ( 11 )
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,1), txc(1,2), wrk2d,wrk3d)
        txc(1:isize_field,1) = log(txc(1:isize_field,1))

        WRITE(fname,*) itime; fname='jpdfXiZ'//TRIM(ADJUSTL(fname))
        CALL JPDF3D(fname, i1, gate_level, i1, imax, jmax, kmax, i0, i0,&
             gate, s, txc(1,1), opt_bins, opt_bins, wrk2d(1,1), wrk2d(1,2), wrk2d(1,3), wrk1d)

! ###################################################################
! Scalar gradient components
! ###################################################################
     CASE ( 12 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,3), wrk3d, wrk2d,wrk3d)
! Angles; s array is overwritten to save space
        DO ij = 1,imax*jmax*kmax
           dummy = txc(ij,2)/sqrt(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
           txc(ij,4) = asin(dummy)                 ! with Oy
           s(ij,1)  = atan2(txc(ij,3),txc(ij,1))  ! with Ox in plane xOz
        ENDDO
        
        data(1)%field => txc(:,1); varname(1) = 'GradientX'; ibc(1) = 2
        data(2)%field => txc(:,2); varname(2) = 'GradientY'; ibc(2) = 2
        data(3)%field => txc(:,3); varname(3) = 'GradientZ'; ibc(3) = 2
        data(4)%field => s(:,1);   varname(4) = 'Theta';     ibc(4) = 2
        data(5)%field => txc(:,4); varname(5) = 'Phi';       ibc(5) = 2

        WRITE(fname,*) itime; fname='jpdfGi'//TRIM(ADJUSTL(fname))
        CALL JPDF3D(fname, i1, gate_level, i1, imax, jmax, kmax, i0, i0,&
             gate, s, txc(1,4), opt_bins, opt_bins, wrk2d(1,1), wrk2d(1,2), wrk2d(1,3), wrk1d)

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Gradient trajectory data from external file
! ###################################################################
     CASE ( 13 )
        WRITE(fname,*) itime; fname = 'de'//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i0, imax,jmax,kmax, i1,i0, isize_wrk3d, txc, wrk3d)

        data(1)%field => txc(:,1); varname(1) = 'ScalarDifference'
!        data(2)%field => txc(:,2); varname(2) = 'ScalarMean'

        IF ( opt_bcs .EQ. 0 ) THEN
           ibc(1) = 0; amin(1) = 0.; amax(1) = 1.
!           ibc(2) = 0; amin(2) = 0.; amax(2) = 1.
        ELSE
           ibc(1) = 2
!           ibc(2) = 2
        ENDIF

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdfDe'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

!        WRITE(fname,*) itime; fname='jpdfDe'//TRIM(ADJUSTL(fname))
!        CALL JPDF3D(fname, i1, gate_level, i1, imax, jmax, kmax, i0, i0,&
!             gate, txc(1,2), txc(1,1), opt_bins, opt_bins, wrk2d(1,1), wrk2d(1,2), wrk2d(1,3), wrk1d)

! ###################################################################
! eigenvalues of rate-of-strain tensor
! ###################################################################
     CASE ( 14 )
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...') ! txc6-txc9
        CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7))

        data(1)%field => txc(:,7); varname(1) = 'Lambda1'; ibc(1) = 2
        data(2)%field => txc(:,8); varname(2) = 'Lambda2'; ibc(2) = 2
        data(3)%field => txc(:,9); varname(3) = 'Lambda3'; ibc(3) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! eigenframe of rate-of-strain tensor
! ###################################################################
     CASE ( 15 )
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...') ! txc7-txc9
        CALL FI_TENSOR_EIGENVALUES(imax,jmax,kmax, txc(1,1), txc(1,7))

        CALL IO_WRITE_ASCII(lfile,'Computing eigenframe...') ! txc1-txc6
        CALL FI_TENSOR_EIGENFRAME(imax,jmax,kmax,  txc(1,1), txc(1,7))

! local direction cosines of vorticity vector
        CALL IO_WRITE_ASCII(lfile,'Computing vorticity vector...') ! txc7-txc9
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk2d,wrk3d)

        DO ij = 1,imax*jmax*kmax
           dummy = sqrt(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           q(ij,1) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           q(ij,2) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1   = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2   = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3   = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           q(ij,3) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        ENDDO

        data(1)%field => q(:,1); varname(1) = 'cos(w,lambda1)'; ibc(1) = 2
        data(2)%field => q(:,2); varname(2) = 'cos(w,lambda2)'; ibc(2) = 2
        data(3)%field => q(:,3); varname(3) = 'cos(w,lambda3)'; ibc(3) = 2

! local direction cosines of scalar gradient vector
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...') ! txc7-txc9
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,7), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,8), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,9), wrk3d, wrk2d,wrk3d)

        DO ij = 1,imax*jmax*kmax
           dummy = sqrt(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
           txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3 
        ENDDO

        data(4)%field => txc(:,7); varname(4) = 'cos(G,lambda1)'; ibc(4) = 2
        data(5)%field => txc(:,8); varname(5) = 'cos(G,lambda2)'; ibc(5) = 2
        data(6)%field => txc(:,9); varname(6) = 'cos(G,lambda3)'; ibc(6) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Scalar gradient components
! ###################################################################
     CASE ( 16 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), q(1,2), txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,3), wrk3d, wrk2d,wrk3d)

        data(1)%field => txc(:,1); varname(1) = 'dudx'
        data(2)%field => txc(:,2); varname(2) = 'dvdy'
        data(3)%field => txc(:,3); varname(3) = 'dwdz'

        ibc(1:nfield) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

! ###################################################################
! Potential vorticity
! ###################################################################
     CASE ( 17 )
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) =                       txc(1:isize_field,1)*txc(1:isize_field,4) 
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,2)*txc(1:isize_field,4) 
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s(1,1), txc(1,4), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1) +txc(1:isize_field,3)*txc(1:isize_field,4)

        txc(1:isize_field,1) = txc(1:isize_field,1)*txc(1:isize_field,1)
        txc(1:isize_field,1) = LOG(txc(1:isize_field,1)+C_SMALL_R)

        data(1)%field => txc(:,1); varname(1) = 'LnPotentialEnstrophy'; ibc(1) = 2

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='pdf'//TRIM(ADJUSTL(fname))
        CALL PDF2D_N(fname, varname, gate_level, rtime, &
             imax*opt_block, jmax_aux, kmax, nfield, ibc, amin, amax, gate, &
             data, opt_bins, npdf_size, pdf, wrk1d)

     END SELECT
  ENDDO

  CALL DNS_END(0)

  STOP

100 FORMAT(G_FORMAT_R)

END PROGRAM PDFS

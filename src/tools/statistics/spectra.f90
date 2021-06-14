#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define SPEC_SINGLE  0
#define SPEC_AVERAGE 1

#define C_FILE_LOC "SPECTRA"

!########################################################################
!# Tool/Library SPECTRA
!#
!########################################################################
!# HISTORY
!#
!# 2012/04/26 - C. Ansorge
!#              Created
!# 2014/01/01 - J.P. Mellado
!#              Adding correlation, and cross terms
!# 2015/01/01 - J.P. Mellado
!#              Parallelizing the 1D spectra; radial spectre not yet
!#
!########################################################################
!# DESCRIPTION
!#
!# Postprocessing tool to compute spectra in temporal mode
!#
!########################################################################
PROGRAM SPECTRA

  USE DNS_TYPES, ONLY : pointers_dt, subarray_dt
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif

#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

#include "integers.h"

#ifdef USE_MPI
#include "mpif.h"
#endif

! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  =  20

! Arrays declarations
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q, s
  TREAL, DIMENSION(:,:), ALLOCATABLE :: wrk1d,wrk2d, txc
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk3d, p_aux, y_aux, samplesize

  TREAL, DIMENSION(:,:), ALLOCATABLE :: out2d, outx,outz,outr

  TYPE(pointers_dt), DIMENSION(16) :: vars

  TARGET q, s, p_aux

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 varname(16)
  CHARACTER*64 str, line
  CHARACTER*8  tag_file, tag_name, tag_var(16)
  TINTEGER p_pairs(16,2)

  TINTEGER opt_main, opt_ffmt, opt_time, opt_block, flag_buoyancy
  TINTEGER flag_mode, iread_flow, iread_scal, ierr
  TINTEGER isize_wrk3d, isize_out2d, isize_aux, sizes(5)
  TINTEGER nfield, nfield_ref
  TINTEGER is, iv, iv_offset, iv1, iv2, ip, j
  TINTEGER jmax_aux, kxmax,kymax,kzmax
  TINTEGER icalc_radial
  TREAL norm, dummy

  TINTEGER kx_total,ky_total,kz_total, kr_total, isize_spec2dr

  TREAL AVG1V2D, COV2V2D

  TINTEGER inb_scal_min, inb_scal_max ! Iterval of scalars to calculate, to be able reduce memory constraints (hard coded)

! Reading variables
  CHARACTER*512 sRes

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TREAL opt_vec(iopt_size_max)

#ifdef USE_MPI
  TINTEGER id
#endif

!########################################################################
!########################################################################
  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d))

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
  opt_ffmt  = 0
  opt_time  = 0

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamSpectra', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'Option ?'
     WRITE(*,*) '1. Main variables 2D spectra'
     WRITE(*,*) '2. Main variables 2D cross-spectra'
     WRITE(*,*) '3. Main variables 2D correlation'
     WRITE(*,*) '4. Main variables 2D cross-correlation'
     WRITE(*,*) '5. Main variables 3D spectra'
     READ(*,*) opt_main

     IF ( opt_main .LT. 5 ) THEN
        WRITE(*,*) 'Planes block size ?'
        READ(*,*) opt_block
     ENDIF
     WRITE(*,*) 'Save full spectra fields to disk (1-yes/0-no) ?'
     READ(*,*) opt_ffmt
     WRITE(*,*) 'Average over time (1-yes/0-no) ?'
     READ(*,*) opt_time
#endif
  ELSE
     opt_main = INT(opt_vec(1))
     IF ( iopt_size .GT. 1 ) opt_block= INT(opt_vec(2))
     IF ( iopt_size .GT. 2 ) opt_ffmt = INT(opt_vec(3))
     IF ( iopt_size .GT. 3 ) opt_time = INT(opt_vec(4))
  ENDIF

  IF ( opt_main .LT. 0 ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'SPECTRA. Missing input [ParamSpectra] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  IF ( opt_block .LT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile, 'SPECTRA. Invalid value of opt_block.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  IF ( opt_time .NE. SPEC_SINGLE .AND. opt_time .NE. SPEC_AVERAGE ) THEN
     CALL IO_WRITE_ASCII(efile, 'SPECTRA. Invalid value of opt_time.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
! in case g(2)%size is not divisible by opt_block, drop the upper most planes
  jmax_aux = g(2)%size/opt_block

  flag_buoyancy = 0 ! default

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
! in case we need the buoyancy statistics
     IF ( buoyancy%type .EQ. EQNS_BOD_QUADRATIC   .OR. &
          buoyancy%type .EQ. EQNS_BOD_BILINEAR    .OR. &
          imixture .EQ. MIXT_TYPE_AIRWATER        .OR. &
          imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        flag_buoyancy = 1
        inb_scal_array= inb_scal_array+1             ! space for the buoyancy field
     ENDIF
  ENDIF

  IF      ( opt_main .EQ. 1 ) THEN; flag_mode = 1 ! spectra
  ELSE IF ( opt_main .EQ. 2 ) THEN; flag_mode = 1
  ELSE IF ( opt_main .EQ. 3 ) THEN; flag_mode = 2 ! correlations
  ELSE IF ( opt_main .EQ. 4 ) THEN; flag_mode = 2
  ELSE IF ( opt_main .EQ. 5 ) THEN; flag_mode = 1 ! spectra
  ENDIF

  IF      ( flag_mode .EQ. 1 ) THEN                 ! spectra
     kxmax = imax/2; kymax = jmax/2; kzmax = kmax/2
  ELSE                                              ! correlation
     kxmax = imax;   kymax = jmax;   kzmax = kmax
  ENDIF
  isize_out2d = imax*jmax_aux*kmax                  ! accumulation of 2D data

! -------------------------------------------------------------------
!  maximum wavenumber & length lag; radial data is not really parallelized yet
  kx_total = MAX(g(1)%size/2,1); ky_total = MAX(g(2)%size/2,1); kz_total = MAX(g(3)%size/2,1)

  IF ( opt_main .EQ. 4 ) THEN    ! Cross-correlations need the full length
     kx_total = MAX(g(1)%size,1); ky_total = MAX(g(2)%size,1); kz_total = MAX(g(3)%size,1)
  ENDIF

  IF ( opt_main .GE. 5 ) THEN ! 3D spectrum
!     kr_total =  INT(SQRT(M_REAL( (kx_total-1)**2 + (kz_total-1)**2 + (ky_total-1)**2))) + 1 ! Use if need to check Parseval's in output data
     kr_total = MIN(kx_total,MIN(ky_total,kz_total))
  ELSE
!     kr_total =  INT(SQRT(M_REAL( (kx_total-1)**2 + (kz_total-1)**2))) + 1 ! Use if need to check Parseval's in output data
     kr_total = MIN(kx_total,kz_total)
  ENDIF

  IF ( opt_main .GE. 5 ) THEN; isize_spec2dr = kr_total            ! 3D spectrum
  ELSE;                        isize_spec2dr = kr_total *jmax_aux; ENDIF

! -------------------------------------------------------------------
! Define MPI-type for writing spectra
! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL SPECTRA_MPIO_AUX(opt_main, opt_block)
#else
  io_aux(:)%offset = 0
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
  iread_flow = icalc_flow
  iread_scal = icalc_scal

  inb_scal_min = 1              ! Change this values if you want to reduce the number of scalars to process
  inb_scal_max = inb_scal_array ! and thereby reduced memory requirements
  ! inb_scal_min = 4
  ! inb_scal_max = 4

  nfield_ref  = 0     ! defining the number of accesible fields
  IF ( icalc_flow .EQ. 1 ) nfield_ref = nfield_ref + inb_flow + 1 ! pressure
  IF ( icalc_scal .EQ. 1 ) nfield_ref = nfield_ref + inb_scal_array

  nfield = 0          ! defining the number of accessed fields
  IF      ( opt_main .EQ. 1 .OR. opt_main .EQ. 3 ) THEN ! Auto-spectra & correlations
     IF ( icalc_flow .EQ. 1 ) nfield = nfield + inb_flow + 1 ! pressure
     IF ( icalc_scal .EQ. 1 ) nfield = nfield +(inb_scal_max - inb_scal_min + 1)

  ELSE IF ( opt_main .EQ. 2 .OR. opt_main .EQ. 4 ) THEN ! cross-spectra and cross-correlations
     IF ( icalc_flow .EQ. 1 ) nfield = nfield + 3
     IF ( icalc_scal .EQ. 1 ) nfield = nfield + 3 *(inb_scal_max - inb_scal_min + 1)

  ELSE
     nfield = nfield_ref

  ENDIF

  inb_txc = 5 ! default

  isize_aux = jmax_aux
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     IF ( MOD(jmax_aux,ims_npro_k) .NE. 0 ) THEN
        isize_aux = ims_npro_k *(jmax_aux/ims_npro_k+1)
     ENDIF

     CALL IO_WRITE_ASCII(lfile,'Initialize MPI type 2 for Oz spectra integration.')
     id = DNS_MPI_K_AUX2
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_aux, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

  ENDIF
#endif

  isize_wrk2d = MAX(isize_wrk2d,isize_aux*kmax); inb_wrk2d = MAX(inb_wrk2d,6)
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))

  ALLOCATE(outx(kxmax*jmax_aux, nfield))
  ALLOCATE(outz(kzmax*jmax_aux, nfield))

  WRITE(str,*) nfield; line = 'Allocating array outr  of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_spec2dr; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(outr(isize_spec2dr, nfield),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'SPECTRA. Not enough memory for spectral data.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  IF ( flag_mode .EQ. 2 ) THEN
     ALLOCATE( samplesize(kr_total) )
  ENDIF

  IF ( opt_ffmt .EQ. 1 ) THEN ! need additional space for 2d spectra
     WRITE(str,*) nfield; line = 'Allocating array out2d of size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) isize_out2d; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(out2d(isize_out2d, nfield),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'SPECTRA. Not enough memory for spectral data.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     WRITE(str,*) isize_txc_field; line = 'Allocating array p_aux of size '//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(p_aux(isize_txc_field),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'SPECTRA. Not enough memory for p_aux.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

! extend array by complex nyquist frequency in x (+1 TCOMPLEX = +2 TREAL)
!              by boundary conditions in y       (+1 TCOMPLEX = +2 TREAL)

  isize_wrk3d = isize_txc_field                ! default
  isize_wrk3d = MAX(isize_wrk3d,isize_spec2dr) ! space needed in INTEGRATE_SPECTRUM

#include "dns_alloc_arrays.h"

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
#include "dns_read_grid.h"

  icalc_radial = 0
  IF ( flag_mode .EQ. 1 .AND. g(1)%size     .EQ. g(3)%size     ) icalc_radial = 1 ! Calculate radial spectra
  IF ( flag_mode .EQ. 2 .AND. g(1)%jac(1,1) .EQ. g(3)%jac(1,1) ) icalc_radial = 1 ! Calculate radial correlations

! ------------------------------------------------------------------------
! Define size of blocks
! ------------------------------------------------------------------------
  y_aux(:) = 0
  DO j = 1,jmax
     is = (j-1)/opt_block + 1
     y_aux(is) = y_aux(is) + y(j,1)/M_REAL(opt_block)
  ENDDO

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

! -------------------------------------------------------------------
! Initialize thermodynamic quantities
! -------------------------------------------------------------------
  CALL FI_PROFILES_INITIALIZE(wrk1d)

! -------------------------------------------------------------------
! Initialize
! -------------------------------------------------------------------
  outx = C_0_R; outz = C_0_R; outr = C_0_R
  IF ( opt_ffmt .EQ. 1 ) out2d = C_0_R

! Normalization
  IF ( opt_main .GE. 5 ) THEN ! 3D spectra
     norm = C_1_R / M_REAL(g(1)%size*g(3)%size*g(2)%size)
  ELSE
     norm = C_1_R / M_REAL(g(1)%size*g(3)%size)
  ENDIF

! Define tags
  IF      ( flag_mode .EQ. 1 ) THEN; tag_file = 'sp';  tag_name = 'E' ! spectra
  ELSE IF ( flag_mode .EQ. 2 ) THEN; tag_file = 'cr';  tag_name = 'C' ! correlations
  ENDIF

! Define reference pointers and tags
  iv = 0
  IF ( icalc_flow .EQ. 1 ) THEN
     iv = iv+1; vars(iv)%field => q(:,1); tag_var(iv) = 'u'
     iv = iv+1; vars(iv)%field => q(:,2); tag_var(iv) = 'v'
     iv = iv+1; vars(iv)%field => q(:,3); tag_var(iv) = 'w'
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        iv = iv+1; vars(iv)%field => q(:,6); tag_var(iv) = 'p'
        iv = iv+1; vars(iv)%field => q(:,5); tag_var(iv) = 'r'
        iv = iv+1; vars(iv)%field => q(:,7); tag_var(iv) = 't'
     ELSE
        iv = iv+1; vars(iv)%field => p_aux;  tag_var(iv) = 'p'
     ENDIF
  ENDIF
  iv_offset = iv

  IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal_array
        WRITE(sRes,*) is
        iv = iv+1; vars(iv)%field => s(:,is); tag_var(iv) = TRIM(ADJUSTL(sRes))
     ENDDO
  ENDIF

  IF ( nfield_ref .NE. iv ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'SPECTRA. Array space nfield_ref incorrect.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

! Define pairs
  iv = 0
  IF      ( opt_main .EQ. 1 .OR. opt_main .EQ. 3 ) THEN ! Auto-spectra & correlations
     DO ip = 1,iv_offset
        iv = iv+1; p_pairs(iv,1) = iv; p_pairs(iv,2) = iv
     ENDDO
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = inb_scal_min,inb_scal_max
           ip = is+iv_offset
           iv = iv+1; p_pairs(iv,1) = ip; p_pairs(iv,2) = ip
        ENDDO
     ENDIF

  ELSE IF ( opt_main .EQ. 2 .OR. opt_main .EQ. 4 ) THEN ! Cross-spectra & correlations
     IF ( icalc_flow .EQ. 1 ) THEN
        iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 2
        iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 3
        iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = 3
        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = inb_scal_min,inb_scal_max
              ip = is+iv_offset
              iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = ip
              iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = ip
              iv = iv+1; p_pairs(iv,1) = 3; p_pairs(iv,2) = ip
           ENDDO
        ENDIF
! block to calculate the pressure-velocity and triple-velocity correlation terms in turbulent transport
        ! iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 4
        ! iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = 4
        ! iv = iv+1; p_pairs(iv,1) = 3; p_pairs(iv,2) = 4
        ! IF ( icalc_scal .EQ. 1 ) THEN ! aux array for u_iu_i/2
        !    s(:,1) = C_05_R*( q(:,1)*q(:,1) + q(:,2)*q(:,2) + q(:,3)*q(:,3) ); tag_var(5) = 'q'
        !    iv = iv+1; p_pairs(iv,1) = 1; p_pairs(iv,2) = 5
        !    iv = iv+1; p_pairs(iv,1) = 2; p_pairs(iv,2) = 5
        !    iv = iv+1; p_pairs(iv,1) = 3; p_pairs(iv,2) = 5
        ! ENDIF
     ELSE
        CALL IO_WRITE_ASCII(efile, 'SPECTRA. Cross-spectra needs flow fields.')
        CALL DNS_STOP(DNS_ERROR_INVALOPT)
     ENDIF

  ENDIF

  IF ( nfield .NE. iv ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'SPECTRA. Array space nfield incorrect.')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

  DO iv = 1,nfield ! define variable names
     varname(iv) = tag_name(1:1)//TRIM(ADJUSTL(tag_var(p_pairs(iv,1))))//TRIM(ADJUSTL(tag_var(p_pairs(iv,2))))
  ENDDO

! ###################################################################
! Calculating statistics
! ###################################################################
  DO it = 1,itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)

! Read data
     IF ( iread_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_flow,i0, isize_wrk3d, q, wrk3d)
     ENDIF

     IF ( iread_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)
     ENDIF

! Calculate diagnostic quantities to be processed
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL FI_PRESSURE_BOUSSINESQ(q,s, p_aux, txc(1,1),txc(1,2), txc(1,3), wrk1d,wrk2d,wrk3d)
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
           CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal+1))
        ENDIF
        IF ( flag_buoyancy .EQ. 1 ) THEN
           IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
              CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, s(1,inb_scal_array))
           ELSE
              wrk1d(1:jmax,1) = C_0_R
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, s(1,inb_scal_array), wrk1d)
           ENDIF
           dummy = C_1_R /froude
           s(:,inb_scal_array) = s(:,inb_scal_array)*dummy
        ENDIF
     ELSE
        CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), q(1,7), wrk3d)
        CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), q(1,7), q(1,6))
     ENDIF

! Remove fluctuation
     IF ( opt_main .GE. 5 ) THEN ! 3D spectra
        DO iv = 1,nfield_ref
          dummy = AVG1V2D(imax*jmax,i1,kmax, i1, i1, vars(iv)%field)  ! 3D average
           vars(iv)%field =  vars(iv)%field - dummy
        ENDDO
     ELSE
        DO iv = 1,nfield_ref
           CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, vars(iv)%field)
        ENDDO
     ENDIF

! reset if needed
     IF( opt_time .EQ. SPEC_SINGLE ) THEN
        outx = C_0_R; outz = C_0_R; outr = C_0_R
        IF ( opt_ffmt .EQ. 1 ) out2d = C_0_R
     ENDIF

! ###################################################################
! 2D Spectra & Correlations
! ###################################################################
     IF ( opt_main .LE. 4 ) THEN

! -------------------------------------------------------------------
! Calculate 2d spectra into array out2d and 1d spectra into arrays outX
! -------------------------------------------------------------------
        DO iv = 1,nfield
           iv1 = p_pairs(iv,1); iv2 = p_pairs(iv,2)

           wrk1d(:,1:3) = C_0_R ! variance to normalize and check Parseval's relation
           DO j=1,jmax
              wrk1d(j,1) = COV2V2D(imax,jmax,kmax,j,vars(iv1)%field,vars(iv2)%field)
              wrk1d(j,2) = COV2V2D(imax,jmax,kmax,j,vars(iv1)%field,vars(iv1)%field)
              wrk1d(j,3) = COV2V2D(imax,jmax,kmax,j,vars(iv2)%field,vars(iv2)%field)
           ENDDO

           txc(1:isize_field,1) =  vars(iv1)%field(1:isize_field)
           IF ( iv2 .EQ. iv1 ) THEN
              CALL OPR_FOURIER_CONVOLUTION_FXZ(i1, flag_mode, imax,jmax,kmax, &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
           ELSE
              txc(1:isize_field,2) = vars(iv2)%field(1:isize_field)
              CALL OPR_FOURIER_CONVOLUTION_FXZ(i2, flag_mode, imax,jmax,kmax, &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
           ENDIF

           IF      ( flag_mode .EQ. 1 ) THEN ! Spectra
              txc(:,1) = txc(:,1) *norm*norm

! Reduce 2D spectra into array wrk3d
              wrk3d = C_0_R
              CALL REDUCE_SPECTRUM(imax,jmax,kmax, opt_block, &
                   txc(1,1),wrk3d, txc(1,3),wrk1d(1,4))

! Calculate and accumulate 1D spectra; only the half of wrk3d with the power data is necessary
              CALL INTEGRATE_SPECTRUM(imax/2, jmax_aux, kmax, kr_total, isize_aux, &
                   wrk3d, outx(1,iv),outz(1,iv),outr(1,iv), wrk2d(1,1), wrk2d(1,3), wrk2d(1,5))

           ELSE IF ( flag_mode .EQ. 2 ) THEN  ! Correlations
              txc(:,2) = txc(:,2) *norm*norm

! Reduce 2D correlation into array wrk3d and accumulate 1D correlation
              wrk3d = C_0_R
              CALL REDUCE_CORRELATION(imax,jmax,kmax, opt_block, kr_total, &
                   txc(1,2), wrk3d, outx(1,iv),outz(1,iv),outr(1,iv), wrk1d(1,2),wrk1d(1,4),icalc_radial)
           ENDIF

! Check Parseval's relation
           ip = g(2)%size - MOD(g(2)%size,opt_block)  ! Drop the uppermost ny%nblock
           WRITE(line,100) MAXVAL( ABS(wrk1d(1:ip,4) - wrk1d(1:ip,1)) )
           WRITE(str, *  ) MAXLOC( ABS(wrk1d(1:ip,4) - wrk1d(1:ip,1)) )
           line = 'Checking Parseval: Maximum residual '//TRIM(ADJUSTL(line))//' at level '//TRIM(ADJUSTL(str))//'.'
           CALL IO_WRITE_ASCII(lfile, line)

! Accumulate 2D information, if needed
           IF ( opt_ffmt .EQ. 1 ) out2d(1:isize_out2d,iv) = out2d(1:isize_out2d,iv) + wrk3d(1:isize_out2d)

        ENDDO

        IF ( flag_mode .EQ. 2 .AND. icalc_radial .EQ. 1 ) THEN  ! Calculate sampling size for radial correlation
           samplesize = C_0_R
           CALL RADIAL_SAMPLESIZE(imax,kmax, kr_total, samplesize)
        ENDIF

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
        IF ( opt_time .EQ. SPEC_SINGLE .OR. it .EQ. itime_size ) THEN

! Normalizing accumulated spectra
           ip = opt_block
           IF( opt_time .EQ. SPEC_AVERAGE ) ip = ip*itime_size
           dummy = C_1_R/M_REAL(ip)
           IF ( ip .GT. 1 ) THEN
              outx = outx *dummy; outz = outz *dummy; outr = outr *dummy
              IF ( opt_ffmt .EQ. 1 ) out2d= out2d*dummy
           ENDIF

! Reducing radial data
#ifdef USE_MPI
           DO iv = 1,nfield
              CALL MPI_Reduce(outr(1,iv), wrk3d, isize_spec2dr, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
              IF ( ims_pro .EQ. 0 ) outr(1:isize_spec2dr,iv) = wrk3d(1:isize_spec2dr)
           ENDDO

           IF ( flag_mode .EQ. 2 .AND. icalc_radial .EQ. 1 ) THEN ! Calculate sampling size for radial correlation
              CALL MPI_Reduce(samplesize, wrk3d, kr_total, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
              IF ( ims_pro .EQ. 0 ) samplesize(1:kr_total) = wrk3d(1:kr_total)
           ENDIF
#endif

! Normalize radial correlation
#ifdef USE_MPI
           IF ( ims_pro .EQ. 0 ) THEN
#endif
           IF ( flag_mode .EQ. 2 .AND. icalc_radial .EQ. 1 ) THEN
              DO iv1 = 1, kr_total
                 IF ( samplesize(iv1) .GT. C_0_R ) samplesize(iv1) = C_1_R /samplesize(iv1)
              ENDDO

              DO iv = 1,nfield
                 DO j = 1,jmax_aux
                    DO iv1 = 1, kr_total
                       ip = iv1 + (j-1)*kr_total
                       outr(ip,iv) = outr(ip,iv) *samplesize(iv1)
                    ENDDO
                 ENDDO
              ENDDO

           ENDIF
#ifdef USE_MPI
           ENDIF
#endif

! Saving 1D fields
           IF ( opt_time .EQ. SPEC_AVERAGE ) THEN
             WRITE(str,*) itime; WRITE(fname,*) itime_vec(1); str = TRIM(ADJUSTL(fname))//'-'//TRIM(ADJUSTL(str))
           ELSE
             WRITE(str,*) itime;
           ENDIF
           fname = 'x'//TRIM(ADJUSTL(tag_file))//TRIM(ADJUSTL(str))
           sizes(1) = kxmax*jmax_aux; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
           CALL IO_WRITE_SUBARRAY4(i1, fname, varname, outx, sizes, wrk3d)

           IF ( g(3)%size .GT. 1 ) THEN
              fname = 'z'//TRIM(ADJUSTL(tag_file))//TRIM(ADJUSTL(str))
              sizes(1) = kzmax*jmax_aux; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
              CALL IO_WRITE_SUBARRAY4(i2, fname, varname, outz, sizes, wrk3d)
           ENDIF

           IF ( icalc_radial .EQ. 1 ) THEN
              fname = 'r'//TRIM(ADJUSTL(tag_file))//TRIM(ADJUSTL(str))
              CALL WRITE_SPECTRUM1D(fname, varname, kr_total*jmax_aux, nfield, outr)
           ENDIF

! Saving 2D fields
           IF ( opt_ffmt .EQ. 1 ) THEN
              IF ( flag_mode .EQ. 2 ) THEN ! correlations
                 fname = 'cor'//TRIM(ADJUSTL(str))
                 sizes(1) = isize_out2d; sizes(2) = 1; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                 CALL IO_WRITE_SUBARRAY4(i3, fname, varname, out2d, sizes, wrk3d)

              ELSE                         ! spectra
                 fname = 'pow'//TRIM(ADJUSTL(str))
                 sizes(1) = isize_out2d; sizes(2) = 1; sizes(3) = sizes(1) /2; sizes(4) = 1; sizes(5) = nfield
                 CALL IO_WRITE_SUBARRAY4(i3, fname, varname, out2d, sizes, wrk3d)

                 fname = 'pha'//TRIM(ADJUSTL(str))
                 sizes(1) = isize_out2d; sizes(2) = 1+sizes(1) /2; sizes(3) = sizes(1); sizes(4) = 1; sizes(5) = nfield
                 CALL IO_WRITE_SUBARRAY4(i3, fname, varname, out2d, sizes, wrk3d)

              ENDIF

           ENDIF

        ENDIF

! ###################################################################
! 3D Spectra
! ###################################################################
     ELSE IF ( opt_main .EQ. 5 ) THEN

        DO iv = 1,nfield
           txc(1:isize_field,1) = vars(iv)%field(1:isize_field)
           CALL OPR_FOURIER_F(i3, imax,jmax,kmax, txc(1,1),txc(1,2), txc(1,3),wrk2d,wrk3d)

           CALL OPR_FOURIER_SPECTRA_3D(imax,jmax,kmax, isize_spec2dr, txc(1,2), outr(1,iv), wrk3d)
        ENDDO

        outr = outr *norm*norm

        WRITE(fname,*) itime; fname = 'rsp'//TRIM(ADJUSTL(fname))
        CALL WRITE_SPECTRUM1D(fname, varname, kr_total, nfield, outr)

     ENDIF

  ENDDO ! Loop in itime

  100 FORMAT(G_FORMAT_R)
  CALL DNS_STOP(0)
END PROGRAM SPECTRA

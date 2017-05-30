#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "VISUALS"

!########################################################################
!# DESCRIPTION
!#
!# Creating data blocks for visualizations. Derived from ensight.f
!# Partition to be incorporated via a MASK routine before VISUALS_WRITE
!#
!########################################################################
PROGRAM VISUALS_MAIN

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE THERMO_GLOBAL, ONLY : NSP, THERMO_SPNAME
  USE LAGRANGE_GLOBAL
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
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: q, s, txc
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_txc
  INTEGER(8), DIMENSION(:),   ALLOCATABLE, SAVE :: l_tags
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: wrk2d,wrk3d
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: wrk1d
  
  INTEGER(1), DIMENSION(:),   ALLOCATABLE, SAVE :: gate
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: surface

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 flow_file, scal_file, part_file, aux_file, plot_file, time_str
  CHARACTER*64 str, line

  TINTEGER opt_format, flag_buoyancy
  TINTEGER opt_cond, opt_cond_scal, opt_threshold
  TINTEGER isize_wrk3d, ij, is, n, bcs(2,2)
  TINTEGER iscal_offset, iread_flow, iread_scal, iread_part, idummy, ierr, MaskSize
  TREAL diff, umin, umax, dummy
  TINTEGER subdomain(6)

! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  INTEGER(1) opt_gate, igate_vec(igate_size_max)
  TREAL gate_threshold(igate_size_max)

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size, iv
  TINTEGER opt_vec(iopt_size_max)
  TREAL opt_vec2(iopt_size_max)

  TINTEGER params_size
  TREAL params(params_size_max)

#ifdef USE_MPI
  INTEGER icount
#endif

!########################################################################
!########################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) &
     CALL PARTICLE_READ_GLOBAL(inifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

  ALLOCATE(gate(isize_field))

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
  opt_format = 2 ! default values
  opt_gate   = 0

  IF    ( imixture .EQ. MIXT_TYPE_NONE ) THEN; iscal_offset = 9    ! define iscal_offset
  ELSE;                                        iscal_offset = 9+NSP; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamVisuals', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_INTEGER(sRes, iopt_size, opt_vec)
  
  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,'(A)') 'Option?'
     WRITE(*,'(A)') ' 0. Grid'
     WRITE(*,'(A)') ' 1. VelocityX'
     WRITE(*,'(A)') ' 2. VelocityY'
     WRITE(*,'(A)') ' 3. VelocityZ'
     WRITE(*,'(A)') ' 4. VelocityVector'
     WRITE(*,'(A)') ' 5. Velocity V_iV_i'
     WRITE(*,'(A)') ' 6. Density'
     WRITE(*,'(A)') ' 7. Temperature'
     WRITE(*,'(A)') ' 8. Pressure'
     WRITE(*,'(A)') ' 9. Scalars'
     IF ( imixture .NE. MIXT_TYPE_NONE ) THEN
        DO is = 1,NSP
           WRITE(*,'(I2,A)') 9+is,'. '//THERMO_SPNAME(is)
        ENDDO
     ENDIF

     WRITE(*,'(I2,A)') iscal_offset+ 1,'. ScalarGradientVector'
     WRITE(*,'(I2,A)') iscal_offset+ 2,'. ScalarGradient G_iG_i (Ln)'
     WRITE(*,'(I2,A)') iscal_offset+ 3,'. ScalarGradientEquation'
     WRITE(*,'(I2,A)') iscal_offset+ 4,'. VorticityVector'
     WRITE(*,'(I2,A)') iscal_offset+ 5,'. Enstrophy W_iW_i (Ln)'
     WRITE(*,'(I2,A)') iscal_offset+ 6,'. EnstrophyEquation'
     WRITE(*,'(I2,A)') iscal_offset+ 7,'. StrainTensor (not yet)'
     WRITE(*,'(I2,A)') iscal_offset+ 8,'. Strain 2S_ijS_ij (Ln)'
     WRITE(*,'(I2,A)') iscal_offset+ 9,'. StrainEquation'
     WRITE(*,'(I2,A)') iscal_offset+10,'. VelocityGradientInvariants'
     WRITE(*,'(I2,A)') iscal_offset+11,'. Space partition'
     WRITE(*,'(I2,A)') iscal_offset+12,'. Buoyancy'
     WRITE(*,'(I2,A)') iscal_offset+13,'. Envelopes'
     WRITE(*,'(I2,A)') iscal_offset+14,'. HorizontalDivergence'
     WRITE(*,'(I2,A)') iscal_offset+15,'. Turbulent quantities'
     WRITE(*,'(I2,A)') iscal_offset+16,'. Radiative forcing'
     WRITE(*,'(I2,A)') iscal_offset+17,'. Relative humidity'
     WRITE(*,'(I2,A)') iscal_offset+18,'. Particle Density'
     WRITE(*,'(I2,A)') iscal_offset+19,'. Thermodynamic quantities'
     READ(*,'(A512)') sRes        
#endif
  ENDIF
  iopt_size = iopt_size_max
  CALL LIST_INTEGER(sRes, iopt_size, opt_vec)

  IF ( opt_vec(1) .LT. 0 ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'VISUALS. Missing input [PostProcessing.ParamVisuals] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT) 
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'Subdomain', '-1', sRes)
  
  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'Subdomain limits ?'
     READ(*,'(A64)') sRes
#endif
  ENDIF
  idummy = 6
  CALL LIST_INTEGER(sRes, idummy, subdomain)
  
  IF ( idummy .LT. 6 ) THEN ! default
     subdomain(1) = 1; subdomain(2) = g(1)%size
     subdomain(3) = 1; subdomain(4) = g(2)%size
     subdomain(5) = 1; subdomain(6) = g(3)%size
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'Format', '-1', sRes)
  
  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'File Format ?'
     WRITE(*,*) ' 0. General restart format'
     WRITE(*,*) ' 1. Ensight'
     WRITE(*,*) ' 2. Raw, single precision, no header'
     READ(*,'(A64)') sRes
#endif
  ENDIF
  IF ( LEN_TRIM(ADJUSTL(sRes)) .GT. 0 ) THEN
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'general' ) THEN; opt_format = 0
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'ensight' ) THEN; opt_format = 1
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'single'  ) THEN; opt_format = 2
     ELSE
        READ(sRes,*) opt_format
     ENDIF
  ENDIF
  
  IF ( opt_format .LT. 0 ) opt_format = 2 ! default is single precission, no header
  
! -------------------------------------------------------------------
! Defining gate levels for conditioning
! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
  opt_cond_scal = 1
  igate_size    = 0
  opt_threshold = 0

  DO iv = 1,iopt_size
     IF ( opt_vec(iv) .EQ. iscal_offset+11 ) opt_gate = 1 ! partition
  ENDDO
  IF ( opt_gate .EQ. 1 ) THEN

#include "dns_read_partition.h"

  ENDIF

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
  MaskSize    = 6

  flag_buoyancy = 0 ! default
  
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN 
! in case we need the buoyancy statistics
     IF ( buoyancy%type .EQ. EQNS_BOD_QUADRATIC   .OR. &
          buoyancy%type .EQ. EQNS_BOD_BILINEAR    .OR. &       
          imixture .EQ. MIXT_TYPE_AIRWATER        .OR. &
          imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        flag_buoyancy = 1
     ENDIF
  ENDIF
  
! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
  iread_flow = 0
  iread_scal = 0
  iread_part = 0
  inb_txc    = 0

  IF      ( opt_cond .EQ. 2 ) THEN
     iread_scal = 1
     inb_txc    = 1
  ELSE IF ( opt_cond .EQ. 3 ) THEN
     inb_txc = MAX(inb_txc,3)
     iread_scal = 1
  ELSE IF ( opt_cond .EQ. 4 ) THEN
     inb_txc = MAX(inb_txc,3)
     iread_flow = 1
  ELSE IF ( opt_cond .EQ. 5 ) THEN
     iread_flow = 1
     inb_txc    = 1
  ENDIF

  DO iv = 1,iopt_size
     IF ( opt_vec(iv) .EQ. 1              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
     IF ( opt_vec(iv) .EQ. 2              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
     IF ( opt_vec(iv) .EQ. 3              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
     IF ( opt_vec(iv) .EQ. 4              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. 5              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
     IF ( opt_vec(iv) .EQ. 6              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. 7              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. 8              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,6); ENDIF
     IF ( opt_vec(iv) .EQ. 9              ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,1); ENDIF
     IF ( opt_vec(iv) .GT. 9 .AND. &
          opt_vec(iv) .LE. iscal_offset   ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,4); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+1 ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+2 ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+3 ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,6); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+4 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,4); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+5 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+6 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+7 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,4); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+8 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+9 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+10) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+12) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,4); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+14) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+15) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,5); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+16) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+17) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+18) THEN; iread_part = 1;                 inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. iscal_offset+19) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF

  ENDDO

  isize_txc   = isize_txc_field *inb_txc
  isize_wrk3d = isize_txc_field
#ifdef USE_MPI
  isize_wrk3d = isize_wrk3d + isize_field ! more space in wrk3d array needed in IO_WRITE_VISUALS
#endif
  IF ( icalc_particle .eq. 1) THEN
     isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  END IF

#include "dns_alloc_arrays.h"

  IF ( iread_part .EQ. 1 ) THEN ! Particle variables
     inb_particle_txc = 1 ! so far, not needed
#include "dns_alloc_larrays.h"
  ENDIF

  DO iv = 1,iopt_size
     IF ( opt_vec(iv) .EQ. iscal_offset+13 ) idummy = 1 ! envelopes/superlayers
  ENDDO
  IF ( idummy .EQ. 1 ) ALLOCATE(surface(imax*kmax,inb_scal))

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1 .AND. inb_txc .GE. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF
 
  IF ( iread_flow .EQ. 1 .AND. inb_txc .GE. 3 ) THEN ! We need array space
     CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

! -------------------------------------------------------------------
! Initialize MPI/IO subarray information
! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL VISUALS_MPIO_AUX(opt_format, subdomain)
#else
  io_aux(:)%offset = 0
#endif

! -------------------------------------------------------------------
! Initialize thermodynamic quantities
! -------------------------------------------------------------------
  CALL FI_PROFILES(wrk1d)
  
! ###################################################################
! Grid
! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     DO iv = 1,iopt_size
        IF ( opt_vec(iv) .EQ. 0 ) THEN
           CALL ENSIGHT_GRID('grid.ensight', imax,jmax,kmax_total, subdomain, x,y,z)
        ENDIF
     ENDDO
#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Postprocess given list of files
! ###################################################################
  DO it=1, itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))//'.'
     CALL IO_WRITE_ASCII(lfile, sRes)
     
     IF ( iread_scal .EQ. 1 ) THEN ! Scalar variables
        WRITE(scal_file,*) itime; scal_file = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(scal_file))
        CALL DNS_READ_FIELDS(scal_file, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)

        IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Calculate q_l
           CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2),s(1,1), epbackground,pbackground)
           
        ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
           CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
           
        ENDIF

     ENDIF

     IF ( iread_flow .EQ. 1 ) THEN ! Flow variables
        WRITE(flow_file,*) itime; flow_file = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(flow_file))
        CALL DNS_READ_FIELDS(flow_file, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)
     ENDIF

     IF ( iread_part .EQ. 1 ) THEN ! Particle variables
        WRITE(part_file,*) itime; part_file="particle."//TRIM(ADJUSTL(part_file))
     ENDIF

     WRITE(sRes,100) rtime; sRes = 'Physical time '//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)
     
! -------------------------------------------------------------------
! Calculate intermittency
! -------------------------------------------------------------------
     IF ( opt_gate .EQ. 1 ) THEN
#include "dns_calc_partition.h"
     ENDIF

! -------------------------------------------------------------------
! define time string
! -------------------------------------------------------------------
     DO ij=MaskSize,1,-1
        time_str(ij:ij)='0'
     ENDDO
     WRITE(plot_file,'(I10)') itime
     time_str(MaskSize-LEN_TRIM(ADJUSTL(plot_file))+1:Masksize)=TRIM(ADJUSTL(plot_file))

! -------------------------------------------------------------------
! Loop over options
! -------------------------------------------------------------------
     DO iv = 1,iopt_size

! ###################################################################
! Velocities
! ###################################################################
        IF      ( opt_vec(iv) .EQ. 1 ) THEN
           txc(1:isize_field,1) = q(1:isize_field,1)

           plot_file = 'VelocityX'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 2 ) THEN
           txc(1:isize_field,1) = q(1:isize_field,2)

           plot_file = 'VelocityY'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 3 ) THEN
           txc(1:isize_field,1) = q(1:isize_field,3)

           plot_file = 'VelocityZ'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 4 ) THEN
           txc(1:isize_field,1:3) = q(1:isize_field,1:3)

           plot_file = 'VelocityVector'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i3, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 5 ) THEN
           txc(1:isize_field,1) = q(1:isize_field,1)*q(1:isize_field,1) &
                                + q(1:isize_field,2)*q(1:isize_field,2) &
                                + q(1:isize_field,3)*q(1:isize_field,3)

           plot_file = 'VelocityMagnitude'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ENDIF

! ###################################################################
! Thermodynamic state
! ###################################################################
        IF ( opt_vec(iv) .EQ. 6 .OR. opt_vec(iv) .EQ. 7 .OR. opt_vec(iv) .EQ. 8 ) THEN
! -------------------------------------------------------------------
! Incompressible
! -------------------------------------------------------------------
           IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
              IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density 
                 IF    ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
                    CALL THERMO_ANELASTIC_DENSITY(imax,jmax,kmax, s, epbackground,pbackground, txc(1,1))
                    
                 ELSE
                    wrk1d(1:jmax,1) = C_0_R
                    CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
                    dummy = C_1_R/froude
                    txc(1:isize_field,1) = txc(1:isize_field,1)*dummy + C_1_R
                 ENDIF
                 
                 plot_file = 'Density'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

              ELSE IF ( opt_vec(iv) .EQ. 7 ) THEN ! temperature
                 IF      ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
                    CALL THERMO_ANELASTIC_TEMPERATURE(imax,jmax,kmax, s, epbackground, txc(1,1))
                    
                    plot_file = 'Temperature'//time_str(1:MaskSize)
                    CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                    
                    IF ( damkohler(1) .GT. C_0_R ) THEN ! Supersaturated liquid
                       txc(1:isize_field,1:2) = s(1:isize_field,1:2)
                       CALL THERMO_AIRWATER_PH(imax,jmax,kmax, txc(1,2),txc(1,1), epbackground,pbackground)
                       txc(1:isize_field,3) = (s(1:isize_field,3)-txc(1:isize_field,3))/s(1,3)
                       
                       plot_file = 'Supsat'//time_str(1:MaskSize)
                       CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,3), wrk3d)

                    ENDIF

                 ENDIF

              ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
                 CALL IO_WRITE_ASCII(lfile,'Computing pressure...')
                 CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
                 
                 CALL IO_WRITE_ASCII(lfile,'Computing pressure gradient vector...')
                 CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc(1,1),txc(1,2), wrk3d, wrk2d,wrk3d)
                 CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,3), wrk3d, wrk2d,wrk3d)
                 CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc(1,1),txc(1,4), wrk3d, wrk2d,wrk3d)
                 

                 txc(1:isize_field,2) = txc(1:isize_field,2)*q(1:isize_field,1) & 
                                      + txc(1:isize_field,3)*q(1:isize_field,2) &
                                      + txc(1:isize_field,4)*q(1:isize_field,3)
                 txc(1:isize_field,2) = -txc(1:isize_field,2)

                 plot_file = 'Pressure'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

                 plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

              ENDIF

! -------------------------------------------------------------------
! compressible
! -------------------------------------------------------------------
           ELSE
              IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density
                 txc(1:isize_field,1) = q(1:isize_field,5)

                 plot_file = 'Density'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                 
              ELSE IF ( opt_vec(iv) .EQ. 7 ) THEN ! temperature
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), txc(1,1), wrk3d)

                 plot_file = 'Temperature'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                 
              ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), txc(1,1), wrk3d)
                 CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), txc(1,1), txc(1,2))

                 plot_file = 'Pressure'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
              ENDIF
              
           ENDIF

        ENDIF

! ###################################################################
! N Scalars
! ###################################################################
! -------------------------------------------------------------------
! All prognostic scalars
! -------------------------------------------------------------------
        IF      ( opt_vec(iv) .EQ. 9 ) THEN
           DO is = 1,inb_scal_array
              WRITE(str,*) is; str = 'Scalar'//TRIM(ADJUSTL(str))
              txc(1:isize_field,1) = s(1:isize_field,is)

              plot_file = TRIM(ADJUSTL(str))//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           ENDDO
          
! -------------------------------------------------------------------
! Individual and diagnostic scalars 
! -------------------------------------------------------------------
        ELSE IF ( opt_vec(iv) .GT. 9 .AND. opt_vec(iv) .LE. iscal_offset ) THEN
           
! -------------------------------------------------------------------
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! s(1,inb_scal+1) contains liquid mass fraction

              IF      ( opt_vec(iv) .EQ. 10 ) THEN ! vapor water mass fraction
                 s(:,1) = s(:,inb_scal) - s(:,inb_scal+1)
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(1)))//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s, wrk3d)
                 
              ELSE IF ( opt_vec(iv) .EQ. 11 ) THEN ! air mass fraction
                 s(:,1) = C_1_R - s(:,inb_scal) 
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(2)))//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s, wrk3d)
                 
              ELSE IF ( opt_vec(iv) .EQ. 12 ) THEN ! liquid mass fraction
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(3)))//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,inb_scal+1), wrk3d)
                 
              ENDIF
              
! -------------------------------------------------------------------
           ELSE ! Plot the chosen species
              is = opt_vec(iv) - 9
              plot_file = TRIM(ADJUSTL(THERMO_SPNAME(is)))//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,is), wrk3d)
              
           ENDIF

        ENDIF
        
! ###################################################################
! Scalar Derivatives
! ###################################################################
        IF ( opt_vec(iv) .GE. iscal_offset+1 .AND. opt_vec(iv) .LE. iscal_offset+3 ) THEN
           DO is = 1,inb_scal_array
              WRITE(str,*) is; str = 'Scalar'//TRIM(ADJUSTL(str))

! Scalar gradient vector
              IF ( opt_vec(iv) .EQ. iscal_offset+1 ) THEN
                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...')
                 CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s(1,is), txc(1,1), wrk3d, wrk2d,wrk3d)
                 CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,is), txc(1,2), wrk3d, wrk2d,wrk3d)
                 CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s(1,is), txc(1,3), wrk3d, wrk2d,wrk3d)
                 
                 plot_file = TRIM(ADJUSTL(str))//'GradientVector'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i3, subdomain, txc(1,1), wrk3d)
              ENDIF
              
! Scalar gradient magnitude
              IF ( opt_vec(iv) .EQ. iscal_offset+2 .OR. opt_vec(iv) .EQ. iscal_offset+3 ) THEN
                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
                 CALL FI_GRADIENT(imax,jmax,kmax, s(1,is), txc(1,1),txc(1,2), wrk2d,wrk3d)
                 
                 plot_file = TRIM(ADJUSTL(str))//'Gradient'//time_str(1:MaskSize)
                 
                 IF ( opt_vec(iv) .EQ. iscal_offset+2 ) THEN 
                    ! IF ( igate .GT. 0 ) THEN                              ! Conditioning
                    !    CALL MINMAX(imax,jmax,kmax, txc(1,1), umin,umax)
                    !    DO ij = 1,isize_field
                    !       IF ( MaskFile(ij) .EQ. 1 ) THEN; txc(ij,1)=log(txc(ij,1))
                    !       ELSE;                            txc(ij,1)=log(umin); ENDIF
                    !    ENDDO
                    ! ELSE                                                  ! Natural log
                       txc(1:isize_field,1) = LOG(txc(1:isize_field,1)+C_SMALL_R)
                       plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
                    ! ENDIF
                 ENDIF
                 
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
              ENDIF

! Scalar gradient equation
              IF ( opt_vec(iv) .EQ. iscal_offset+3 .AND. is .LE. inb_scal ) THEN
                 IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
                 ELSE;                                  diff = visc/schmidt(is); ENDIF

                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
                 CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s(1,is), q(1,1),q(1,2),q(1,3), &
                      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
                 plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                 
                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
                 CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s(1,is), &
                      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
                 txc(1:isize_field,1) = diff *txc(1:isize_field,1)
                 
                 plot_file = TRIM(ADJUSTL(str))//'GradientDiffusion'//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                 
              ENDIF
              
           ENDDO

        ENDIF
        
! ###################################################################
! Velocity Derivatives
! ###################################################################

! -------------------------------------------------------------------
! Vorticity 
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+4 ) THEN ! VorticityVector
           CALL IO_WRITE_ASCII(lfile,'Computing vorticity vector...')
           CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
           
           plot_file = 'VorticityVector'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i3, subdomain, txc(1,1), wrk3d)
        ENDIF
           
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+5 .OR. opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! Enstrophy
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy...')
           CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
           
           plot_file = 'Enstrophy'//time_str(1:MaskSize)

           IF ( opt_vec(iv) .EQ. iscal_offset+5 ) THEN ! Natural log
              txc(1:isize_field,1)=LOG(txc(1:isize_field,1)+C_SMALL_R)
              plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
           ENDIF
           
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
        ENDIF
        
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! EnstrophyEquation
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
           CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
                txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
           
           plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
           CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
                txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
           txc(1:isize_field,1)=visc*txc(1:isize_field,1)
           
           plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           
        ENDIF
        
! -------------------------------------------------------------------
! Strain
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+7 ) THEN ! Strain Tensor
           CALL IO_WRITE_ASCII(efile,'VISUALS. Strain tensor undevelop.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF

! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+8 .OR. opt_vec(iv) .EQ. iscal_offset+9 ) THEN ! Strain
           CALL IO_WRITE_ASCII(lfile,'Computing strain...')
           CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)
           txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)
           
           plot_file = 'Strain'//time_str(1:MaskSize)

           IF ( opt_vec(iv) .EQ. iscal_offset+8 ) THEN ! Natural log
              txc(1:isize_field,1)=LOG(txc(1:isize_field,1)+C_SMALL_R)
              plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
           ENDIF
           
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
        ENDIF

! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+9 ) THEN ! StrainEquation (I need the pressure)
              CALL IO_WRITE_ASCII(lfile,'Computing strain pressure...')
              IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
                 CALL IO_WRITE_ASCII(efile,'VISUALS. Strain eqn for incompressible undeveloped.')
                 CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
              ELSE
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), txc(1,1), wrk3d)
                 CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), txc(1,1), q(1,4)) ! pressure in q4
              ENDIF
              CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, q(1,1),q(1,2),q(1,3),q(1,4), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)

              plot_file = 'StrainPressure'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

              CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
              CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)

              plot_file = 'StrainProduction'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

              CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
              CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*visc*txc(1:isize_field,1)

              plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           ENDIF
           
! -------------------------------------------------------------------
! Invariants
! -------------------------------------------------------------------
           IF ( opt_vec(iv) .EQ. iscal_offset+10 ) THEN
              CALL IO_WRITE_ASCII(lfile,'Computing first invariant P...')
              CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), wrk2d,wrk3d)

              plot_file = 'InvariantP'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

              CALL IO_WRITE_ASCII(lfile,'Computing second invariant Q...')
              CALL FI_INVARIANT_Q(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)

              plot_file = 'InvariantQ'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

              CALL IO_WRITE_ASCII(lfile,'Computing third invariant R...')
              CALL FI_INVARIANT_R(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
                   txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

              plot_file = 'InvariantR'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

           ENDIF

! ###################################################################
! Partition
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+11 ) THEN           
           CALL IO_WRITE_ASCII(efile,'VISUALS. Partition undevelop.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF

! ###################################################################
! Buoyancy
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+12 ) THEN
           IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
              CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,1))
           ELSE
              wrk1d(1:jmax,1) = C_0_R
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
           ENDIF
           dummy =  C_1_R/froude
           txc(1:isize_field,1) = txc(1:isize_field,1) *dummy
           
           plot_file = 'Buoyancy'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

! buoyancy flux along Oy
           txc(1:isize_field,2) = txc(1:isize_field,1) *q(1:isize_field,2)

           plot_file = 'Cvb'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

! buoyancy fluctuation
           CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(1,1))

           plot_file = 'bPrime'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

! buoyancy source
           IF ( flag_buoyancy .EQ. 1 ) THEN
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
                 CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, txc(1,1),txc(1,2),txc(1,3))
                 CALL FI_GRADIENT(imax,jmax,kmax, txc(1,1),txc(1,2), txc(1,4), wrk2d,wrk3d)
                 
                 dummy = buoyancy%parameters(inb_scal_array)
                 txc(1:isize_field,2) = txc(1:isize_field,2) *txc(1:isize_field,3) *dummy
                 
              ELSE
                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
                 CALL FI_GRADIENT(imax,jmax,kmax, s, txc(1,1),txc(1,2), wrk2d,wrk3d)

                 CALL FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, txc(1,1), txc(1,2))
                 
              ENDIF

              dummy =  visc/schmidt(1) /froude

              txc(1:isize_field,1) = txc(1:isize_field,2) *dummy
              txc(1:isize_field,1) = LOG(ABS(txc(1:isize_field,1))+C_SMALL_R)
              
              plot_file = 'LnBuoyancySource'//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
              
           ENDIF

        ENDIF

! ###################################################################
! envelopes/superlayers
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+13 ) THEN
! we use the # scalars as a surrogate for the # envelopes to have control through dns.ini
           WRITE(aux_file,*) itime; aux_file='lower.'//TRIM(ADJUSTL(aux_file))
           CALL DNS_READ_FIELDS(aux_file, i1, imax,i1,kmax, inb_scal,i0, isize_wrk2d, surface, wrk2d)

           DO is = 1,inb_scal 
              subdomain(3) = 1; subdomain(4) = 1
              WRITE(str,*) is; plot_file = 'LowerEnvelope'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,i1,kmax, i1, subdomain, surface(1,is), wrk3d)
           ENDDO

           WRITE(aux_file,*) itime; aux_file='upper.'//TRIM(ADJUSTL(aux_file))
           CALL DNS_READ_FIELDS(aux_file, i1, imax,i1,kmax, inb_scal,i0, isize_wrk2d, surface, wrk2d)

           DO is = 1,inb_scal 
              subdomain(3) = 1; subdomain(4) = 1
              WRITE(str,*) is; plot_file = 'UpperEnvelope'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
              CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,i1,kmax, i1, subdomain, surface(1,is), wrk3d)

           ENDDO
        ENDIF

! ###################################################################
! Horizontal divergence
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+14 ) THEN
           CALL IO_WRITE_ASCII(lfile,'Computing horizontal divergence...')
           CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,2), wrk3d, wrk2d,wrk3d)
           CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,1), wrk3d, wrk2d,wrk3d)
           txc(1:isize_field,1) = txc(1:isize_field,1) + txc(1:isize_field,2)

           plot_file = 'HorizontalDivergence'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ENDIF

! ###################################################################
! Turbulent quantities k and \epsilon
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+15 ) THEN

! turbulent dissipation rate
           CALL IO_WRITE_ASCII(lfile,'Computing dissipation rate...')
           CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
                txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
           txc(1:isize_field,1)=LOG(txc(1:isize_field,1)+C_SMALL_R)

           plot_file = 'LnDissipation'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

! turbulent kinetic energy
           txc(1:isize_field,1) = q(1:isize_field,1); CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,1))
           txc(1:isize_field,2) = q(1:isize_field,2); CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,2))
           txc(1:isize_field,3) = q(1:isize_field,3); CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,3))

           txc(1:isize_field,1) = txc(1:isize_field,1)*txc(1:isize_field,1) &
                                + txc(1:isize_field,2)*txc(1:isize_field,2) &
                                + txc(1:isize_field,3)*txc(1:isize_field,3)

           plot_file = 'Tke'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           
        ENDIF

! ###################################################################
! Radiation
! ###################################################################
        IF (  opt_vec(iv) .EQ. iscal_offset+16) THEN    
           
           DO is = 1,inb_scal
              IF ( radiation%active(is) ) THEN
                 IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
                    CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, s(1,radiation%scalar(is)), txc(1,2))
                    CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), txc(1,2),                 txc(1,1), wrk1d,wrk3d)
                    CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, ribackground, txc(1,1))
                    
                 ELSE
                    CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(1)), txc(1,1), wrk1d,wrk3d)
                    
                 ENDIF
                 
                 WRITE(str,*) is; plot_file = 'Radiation'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
                 
              ENDIF
           ENDDO
           
        ENDIF
        
! ###################################################################
! Relative humidity
! ###################################################################
        IF (  opt_vec(iv) .EQ. iscal_offset+17) THEN    
           
           CALL THERMO_ANELASTIC_RELATIVEHUMIDITY(imax,jmax,kmax, s, epbackground,pbackground, wrk3d, txc(1,1))

           plot_file = 'RelativeHumidity'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           
        ENDIF
        
! ###################################################################
! Particle density
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+18 ) THEN
           CALL DNS_READ_PARTICLE(part_file,l_q)
           l_txc = C_1_R; ! We want density
           CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,1))
           str = 'ParticleDensity'
           plot_file = TRIM(ADJUSTL(str))//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           txc(:,1) = txc(:,1) + 0.00000001
           IF (inb_particle .GT. 3 ) THEN
              IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 )  THEN
                 DO is=1,2
                    l_txc(:,1)=l_q(:,3+is) !!! DO WE WANT l_txc(:,is) ???
                    CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,2))   
                    txc(:,2) = txc(:,2)/txc(:,1)
                    plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(is)))//time_str(1:MaskSize)
                    CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
                 END DO
              END IF
              IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
                 l_txc(:,1)=l_q(:,inb_particle) !inb_particle is the last component -> residence times in bil_cloud_4
                 CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,2))   
                 plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(3)))//time_str(1:MaskSize)
                 CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
              ENDIF
           END IF

        ENDIF
        
! ###################################################################
! Thermodynamic quantities
! ###################################################################
        IF (  opt_vec(iv) .EQ. iscal_offset+19) THEN    
           
           CALL VISUALS_FUNCTION1(imax,jmax,kmax, s,txc)

           plot_file = 'Enthalpy'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,1), wrk3d)
           plot_file = 'TotalWater'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,2), wrk3d)
           plot_file = 'LiquidWater'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,3), wrk3d)
           plot_file = 'Temperature'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
           plot_file = 'Density'//time_str(1:MaskSize)
           CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
           
        ENDIF

     ENDDO

  ENDDO
  
  CALL DNS_END(0)

  STOP

100 FORMAT(G_FORMAT_R)

END PROGRAM VISUALS_MAIN

!########################################################################
!# DESCRIPTION
!#
!# Calculate thermodynamic information for THERMO_AIRWATER_LINEAR
!#
!########################################################################
SUBROUTINE VISUALS_FUNCTION1(nx,ny,nz, s,txc)

  USE DNS_GLOBAL, ONLY : isize_txc_field
  USE DNS_GLOBAL, ONLY : epbackground, pbackground
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*)        :: s
  TREAL, DIMENSION(isize_txc_field,*) :: txc
  
! -----------------------------------------------------------------------
  TREAL qt_0,qt_1, h_0,h_1, p, T_0, C_0, PsiRef

! #######################################################################
  imixture = MIXT_TYPE_AIRWATER

  qt_0 = 9.0d-3;     qt_1 = 1.5d-3
  h_0  = 0.955376d0; h_1  = 0.981965d0
  p    = 0.940d0
  T_0  = 0.952181d0 ! 283.75 / TREF
  C_0  = 1.0089
  PsiRef= 6.57d-4
  
  s(:,3) = h_0  + s(:,1)*(h_1 -h_0 ) + s(:,2) *C_0 *T_0 *PsiRef ! enthalpy
  s(:,2) = qt_0 + s(:,1)*(qt_1-qt_0)                            ! total water, space for q_l
  s(:,1) = s(:,3)
  
  epbackground = C_0_R                                    ! potential energy
  pbackground  = p                                        ! pressure

  CALL THERMO_AIRWATER_PH(nx,ny,nz, s(1,2), s(1,1), epbackground,pbackground)
  CALL THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, s(1,1), epbackground, txc(1,1))
  CALL THERMO_ANELASTIC_DENSITY(nx,ny,nz, s(1,1), epbackground,pbackground, txc(1,2))

  RETURN
END SUBROUTINE VISUALS_FUNCTION1

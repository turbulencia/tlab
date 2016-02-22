#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "VISUALS"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2014/01/01 - J.P. Mellado
!#              Created
!#
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
#ifdef PARALLEL
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef PARALLEL
#include "mpif.h"
#endif

! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  =  20
  TINTEGER, PARAMETER :: igate_size_max =   8
  TINTEGER, PARAMETER :: params_size_max =  2

! Arrays declarations
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: x,y,z, dx,dy,dz
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: q, s, txc
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_txc
  INTEGER(8), DIMENSION(:),   ALLOCATABLE, SAVE :: l_tags
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: wrk2d,wrk3d
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: wrk1d
  
  INTEGER(1), DIMENSION(:),   ALLOCATABLE, SAVE :: gate
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: surface
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: tmp_mpi

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 flow_file, scal_file, part_file, aux_file, plot_file, time_str
  CHARACTER*64 str, line

  TINTEGER opt_format, flag_buoyancy
  TINTEGER opt_cond, opt_threshold
  TINTEGER isize_wrk3d, ij, is, n
  TINTEGER iscal_offset, iread_flow, iread_scal, iread_part, idummy, ierr, MaskSize
  TREAL diff, umin, umax, dummy, dummy2
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

#ifdef PARALLEL
  TINTEGER itmp_mpi
  INTEGER icount
#endif

!########################################################################
!########################################################################
  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef PARALLEL
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

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
  opt_format = 1 ! default values
  opt_gate   = 0

  IF    ( imixture .EQ. MIXT_TYPE_NONE ) THEN; iscal_offset = 9    ! define iscal_offset
  ELSE;                                        iscal_offset = 9+NSP; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamVisuals', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_INTEGER(sRes, iopt_size, opt_vec)
  
  IF ( sRes .EQ. '-1' ) THEN
#ifdef PARALLEL
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
     WRITE(*,'(I2,A)') iscal_offset+17,'. Particle Density'
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
#ifdef PARALLEL
#else
     WRITE(*,*) 'Subdomain limits ?'
     READ(*,'(A64)') sRes
#endif
  ENDIF
  idummy = 6
  CALL LIST_INTEGER(sRes, idummy, subdomain)
  
  IF ( idummy .LT. 6 ) THEN ! default
     subdomain(1) = 1; subdomain(2) = imax_total
     subdomain(3) = 1; subdomain(4) = jmax_total
     subdomain(5) = 1; subdomain(6) = kmax_total
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'Format', '-1', sRes)
  
  IF ( sRes .EQ. '-1' ) THEN
#ifdef PARALLEL
#else
     WRITE(*,*) 'File Format ?'
     WRITE(*,*) ' 0. General restart format'
     WRITE(*,*) ' 1. Ensight'
     WRITE(*,*) ' 2. Raw, single precision, no header'
     READ(*,'(A64)') sRes
#endif
  ENDIF
  READ(sRes,*) opt_format
  
  IF ( opt_format .LT. 0 ) opt_format = 1 ! default
  
! -------------------------------------------------------------------
! Defining gate levels for conditioning
! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
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

! in case we need the buoyancy statistics
  IF ( ibodyforce .EQ. EQNS_BOD_BILINEAR           .OR. &
       ibodyforce .EQ. EQNS_BOD_QUADRATIC          .OR. &
       ibodyforce .EQ. EQNS_BOD_PIECEWISE_LINEAR   .OR. &
       ibodyforce .EQ. EQNS_BOD_PIECEWISE_BILINEAR ) THEN
     flag_buoyancy = 1
  ELSE 
     flag_buoyancy = 0   
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
     IF ( opt_vec(iv) .EQ. 1              )       iread_flow = 1
     IF ( opt_vec(iv) .EQ. 2              )       iread_flow = 1
     IF ( opt_vec(iv) .EQ. 3              )       iread_flow = 1
     IF ( opt_vec(iv) .EQ. 4              )       iread_flow = 1
     IF ( opt_vec(iv) .EQ. 5              )       iread_flow = 1
     IF ( opt_vec(iv) .EQ. 6              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
     IF ( opt_vec(iv) .EQ. 7              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
     IF ( opt_vec(iv) .EQ. 8              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,4); ENDIF
     IF ( opt_vec(iv) .EQ. 9              )                       iread_scal = 1
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
     IF ( opt_vec(iv) .EQ. iscal_offset+16)                       iread_scal = 1
     IF ( opt_vec(iv) .EQ. iscal_offset+17) THEN; iread_part = 1;                 inb_txc=MAX(inb_txc,2); ENDIF ! Alberto check 2 or 1?

  ENDDO

  isize_txc   = isize_txc_field *inb_txc
  isize_wrk3d = isize_txc_field
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

! auxiliar array for in parallel mode
#ifdef PARALLEL
  itmp_mpi=isize_field*2
  WRITE(str,*) itmp_mpi; line = 'Allocating array tmp_mpi. Size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(tmp_mpi(itmp_mpi),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'VISUALS. Not enough memory for tmp_mpi.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
#endif

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
     CALL DNS_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Grid
! ###################################################################
#ifdef PARALLEL
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     DO iv = 1,iopt_size
        IF ( opt_vec(iv) .EQ. 0 ) THEN
           CALL ENSIGHT_GRID('grid.ensight', imax,jmax,kmax_total, subdomain, x,y,z)
        ENDIF
     ENDDO
#ifdef PARALLEL
  ENDIF
#endif

! ###################################################################
! Postprocess given list of files
! ###################################################################
  DO it=1, itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)
     
     IF ( iread_scal .EQ. 1 ) THEN ! Scalar variables
        WRITE(scal_file,*) itime; scal_file = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(scal_file))
        CALL DNS_READ_FIELDS(scal_file, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
     ENDIF

     IF ( iread_flow .EQ. 1 ) THEN ! Flow variables
        WRITE(flow_file,*) itime; flow_file = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(flow_file))
        CALL DNS_READ_FIELDS(flow_file, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)
     ENDIF

     IF ( iread_part .EQ. 1 ) THEN ! Particle variables
        WRITE(part_file,*) itime; part_file="particle."//TRIM(ADJUSTL(part_file))
     ENDIF

     IF      ( (imixture .EQ. MIXT_TYPE_AIRWATER) .OR. &                                     ! calculate liquid
               ( (imixture .EQ. MIXT_TYPE_SUPSAT) .AND. ( damkohler(1) .LE. C_0_R) ) ) THEN
        CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)
        CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(:,2), p_init, s(:,1))    
        
     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(:,inb_scal_array), wrk3d)
        
     ENDIF

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
           plot_file = 'VelocityX'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, q(:,1), tmp_mpi)

        ELSE IF ( opt_vec(iv) .EQ. 2 ) THEN
           plot_file = 'VelocityY'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, q(:,2), tmp_mpi)

        ELSE IF ( opt_vec(iv) .EQ. 3 ) THEN
           plot_file = 'VelocityZ'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, q(:,3), tmp_mpi)

        ELSE IF ( opt_vec(iv) .EQ. 4 ) THEN
           plot_file = 'VelocityVector'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i1, opt_format, imax,jmax,kmax, subdomain, q, tmp_mpi)

        ELSE IF ( opt_vec(iv) .EQ. 5 ) THEN
           wrk3d(1:isize_field) = q(1:isize_field,1)*q(1:isize_field,1) &
                                + q(1:isize_field,2)*q(1:isize_field,2) &
                                + q(1:isize_field,3)*q(1:isize_field,3)

           plot_file = 'VelocityMagnitude'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)

        ENDIF

! ###################################################################
! Thermodynamic state
! ###################################################################
        IF ( opt_vec(iv) .EQ. 6 .OR. opt_vec(iv) .EQ. 7 .OR. opt_vec(iv) .EQ. 8 ) THEN
! -------------------------------------------------------------------
! Incompressible
! -------------------------------------------------------------------
           IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
                imode_eqns .EQ. DNS_EQNS_ANELASTIC           ) THEN

              IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density 
                 wrk1d(1:jmax,1) = C_0_R
                 CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
                 dummy = C_1_R/froude
                 wrk3d(1:isize_field) = wrk3d(1:isize_field)*dummy + C_1_R

                 plot_file = 'Density'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)

              ELSE IF ( opt_vec(iv) .EQ. 7 ) THEN ! temperature

                 IF      ( imixture .EQ. MIXT_TYPE_AIRWATER )  THEN
                    CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho) ! calculate liquid
                    CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(:,2), p_init, s(:,1))    
                    CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1), wrk3d)
                    
                    plot_file = 'Temperature'//time_str(1:MaskSize)
                    CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)
                    
                 ELSE IF ( imixture .EQ. MIXT_TYPE_SUPSAT   ) THEN  ! liquid is already in the third scalar
                    CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1), wrk3d)
                    txc(1:isize_field,1) = (txc(1:isize_field,1) - txc(1,1))/(txc(isize_field,1) - txc(1,1)) 
                    
                    plot_file = 'Temperature'//time_str(1:MaskSize)
                    CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)
                    
! Buoyancy
                    CALL THERMO_THERMAL_DENSITY_HP_ALWATER(imax,jmax,kmax, s(1,2),s(1,1),p_init,txc(1,1))
                    txc(1:isize_field,1) = body_vector(2)*(txc(1:isize_field,1) - mean_rho)/mean_rho
                    
                    plot_file = 'Buoyancy'//time_str(1:MaskSize)
                    CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
                    
! Supersaturated liquid
                    txc(1:isize_field,1:2) = s(1:isize_field,1:2)
                    CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, txc(:,2), p_init, txc(:,1))
                    txc(1:isize_field,3) = (s(1:isize_field,3)-txc(1:isize_field,3))/s(1,3)
                     
                    plot_file = 'Supsat'//time_str(1:MaskSize)
                    CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,3), tmp_mpi)
                    
                 ENDIF

              ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
                 CALL IO_WRITE_ASCII(lfile,'Computing pressure...')
                 CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, q(:,1),q(:,2),q(:,3), s, &
                      txc(:,1),txc(:,2),txc(:,3),txc(:,4), wrk1d,wrk2d,wrk3d)
                 
                 CALL IO_WRITE_ASCII(lfile,'Computing pressure gradient vector...')
                 CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, txc(:,1),txc(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
                 CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, txc(:,1),txc(1,3), i0,i0, wrk1d,wrk2d,wrk3d)
                 CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, txc(:,1),txc(1,4), i0,i0, wrk1d,wrk2d,wrk3d)
                 

                 wrk3d(1:isize_field) = txc(1:isize_field,2)*q(1:isize_field,1) & 
                                      + txc(1:isize_field,3)*q(1:isize_field,2) &
                                      + txc(1:isize_field,4)*q(1:isize_field,3)
                 wrk3d(1:isize_field) = -wrk3d(1:isize_field)

                 plot_file = 'Pressure'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)

                 plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)

              ENDIF

! -------------------------------------------------------------------
! compressible
! -------------------------------------------------------------------
           ELSE
              IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density
                 plot_file = 'Density'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, q(:,5), tmp_mpi)
                 
              ELSE IF ( opt_vec(iv) .EQ. 7 ) THEN ! temperature
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(:,4), q(:,5), txc(:,1), wrk3d)

                 plot_file = 'Temperature'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)
                 
              ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(:,4), q(:,5), txc(:,1), wrk3d)
                 CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(:,5), txc(:,1), txc(:,2))

                 plot_file = 'Pressure'//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,2), tmp_mpi)
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
           DO is = 1,inb_scal
              WRITE(str,*) is; str = 'Scalar'//TRIM(ADJUSTL(str))
              
              plot_file = TRIM(ADJUSTL(str))//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s(1,is), tmp_mpi)
           ENDDO
          
! -------------------------------------------------------------------
! Individual and diagnostic scalars 
! -------------------------------------------------------------------
        ELSE IF ( opt_vec(iv) .GT. 9 .AND. opt_vec(iv) .LE. iscal_offset ) THEN
           
! -------------------------------------------------------------------
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! s(1,inb_scal+1) contains liquid mass fraction

              IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
                   imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN
                 CALL THERMO_AIRWATER_PHAL(i1, i1,i1, mean_i(2), p_init, mean_i(1))
                 CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)
                 CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(:,2), p_init, s(:,1))    

              ELSE ! Compressible
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(:,4), q(:,5), txc(:,1), wrk3d)
                 
              ENDIF

              IF      ( opt_vec(iv) .EQ. 10 ) THEN ! vapor water mass fraction
                 s(:,1) = s(:,inb_scal) - s(:,inb_scal+1)
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(1)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s, tmp_mpi)
                 
              ELSE IF ( opt_vec(iv) .EQ. 11 ) THEN ! air mass fraction
                 s(:,1) = C_1_R - s(:,inb_scal) 
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(2)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s, tmp_mpi)
                 
              ELSE IF ( opt_vec(iv) .EQ. 12 ) THEN ! liquid mass fraction
                 
                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(3)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s(1,inb_scal+1), tmp_mpi)
                 
              ENDIF
              
! -------------------------------------------------------------------
           ELSE IF  ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN

              IF      ( opt_vec(iv) .EQ. 10 ) THEN ! total water vapor water mass fraction

                 plot_file = TRIM(ADJUSTL('Totalwater'))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s, tmp_mpi)

              ELSE IF ( opt_vec(iv) .EQ. 11 ) THEN ! Enthalpy TO DO

                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(2)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s(1,2), tmp_mpi)

              ELSE IF ( opt_vec(iv) .EQ. 12 ) THEN ! liquid mass fraction
                 CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array))

                 plot_file = TRIM(ADJUSTL(THERMO_SPNAME(inb_scal_array)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s(1,inb_scal_array), tmp_mpi)

              ENDIF

! -------------------------------------------------------------------
           ELSE ! Plot the chosen species
              is = opt_vec(iv) - 9
              plot_file = TRIM(ADJUSTL(THERMO_SPNAME(is)))//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, s(1,is), tmp_mpi)
              
           ENDIF

        ENDIF
        
! ###################################################################
! Scalar Derivatives
! ###################################################################
        IF ( opt_vec(iv) .GE. iscal_offset+1 .AND. opt_vec(iv) .LE. iscal_offset+3 ) THEN
           DO is = 1,inb_scal
           WRITE(str,*) is; str = 'Scalar'//TRIM(ADJUSTL(str))

           IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
           ELSE;                                  diff = visc/schmidt(is); ENDIF

! Scalar gradient vector
           IF ( opt_vec(iv) .EQ. iscal_offset+1 ) THEN
              CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...')
              CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, s(1,is), txc(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
              CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s(1,is), txc(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
              CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, s(1,is), txc(1,3), i0,i0, wrk1d,wrk2d,wrk3d)
              
              plot_file = TRIM(ADJUSTL(str))//'GradientVector'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i1, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           ENDIF

! Scalar gradient
           IF ( opt_vec(iv) .EQ. iscal_offset+2 .OR. opt_vec(iv) .EQ. iscal_offset+3 ) THEN
              CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
              CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, s(1,is), &
                   txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)
              
              plot_file = TRIM(ADJUSTL(str))//'Gradient'//time_str(1:MaskSize)

              IF ( opt_vec(iv) .EQ. iscal_offset+2 ) THEN 
                 ! IF ( igate .GT. 0 ) THEN                              ! Conditioning
                 !    CALL MINMAX(imax,jmax,kmax, txc(1,1), umin,umax)
                 !    DO ij = 1,isize_field
                 !       IF ( MaskFile(ij) .EQ. 1 ) THEN; txc(ij,1)=log(txc(ij,1))
                 !       ELSE;                            txc(ij,1)=log(umin); ENDIF
                 !    ENDDO
                 ! ELSE                                                  ! Natural log
                    txc(1:isize_field,1) = log(txc(1:isize_field,1))
                    plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
                ! ENDIF
              ENDIF

              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           ENDIF

! Scalar gradient equation
           IF ( opt_vec(iv) .EQ. iscal_offset+3 ) THEN
              CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
              CALL FI_GRADIENT_PRODUCTION(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, s(1,is), &
                   q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
              plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

              CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
              CALL FI_GRADIENT_DIFFUSION&
                   (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, s(1,is), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
              txc(1:isize_field,1)=diff*txc(1:isize_field,1)

              plot_file = TRIM(ADJUSTL(str))//'GradientDiffusion'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

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
           CALL FI_CURL(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, &
                   q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d)
           
           plot_file = 'VorticityVector'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i1, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
        ENDIF
           
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+5 .OR. opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! Enstrophy
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy...')
           CALL FI_VORTICITY(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
                dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)
           
           plot_file = 'Enstrophy'//time_str(1:MaskSize)

           IF ( opt_vec(iv) .EQ. iscal_offset+5 ) THEN ! Natural log
              txc(1:isize_field,1)=log(txc(1:isize_field,1))
              plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
           ENDIF
           
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
        ENDIF
        
! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! EnstrophyEquation
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
           CALL FI_VORTICITY_PRODUCTION(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
           
           plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           
           CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
           CALL FI_VORTICITY_DIFFUSION&
                (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
           txc(1:isize_field,1)=visc*txc(1:isize_field,1)
           
           plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           
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
           CALL FI_STRAIN(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)
           txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)
           
           plot_file = 'Strain'//time_str(1:MaskSize)

           IF ( opt_vec(iv) .EQ. iscal_offset+8 ) THEN ! Natural log
              txc(1:isize_field,1)=log(txc(1:isize_field,1))
              plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
           ENDIF
           
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
        ENDIF

! -------------------------------------------------------------------
        IF ( opt_vec(iv) .EQ. iscal_offset+9 ) THEN ! StrainEquation (I need the pressure)
              CALL IO_WRITE_ASCII(lfile,'Computing strain pressure...')
              IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
                   imode_eqns .EQ. DNS_EQNS_ANELASTIC     )THEN
                 CALL IO_WRITE_ASCII(efile,'VISUALS. Strain eqn for incompressible undeveloped.')
                 CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
              ELSE
                 CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), txc(1,1), wrk3d)
                 CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), txc(1,1), q(1,4)) ! pressure in q4
              ENDIF
              CALL FI_STRAIN_PRESSURE&
                   (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3),q(1,4), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)

              plot_file = 'StrainPressure'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

              CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
              CALL FI_STRAIN_PRODUCTION(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)

              plot_file = 'StrainProduction'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

              CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
              CALL FI_STRAIN_DIFFUSION&
                   (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), &
                   txc(1,1), txc(1,2), txc(1,3), txc(1,4), txc(1,5), txc(1,6), wrk1d,wrk2d,wrk3d)
              txc(1:isize_field,1)=C_2_R*visc*txc(1:isize_field,1)

              plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           ENDIF
           
! -------------------------------------------------------------------
! Invariants
! -------------------------------------------------------------------
           IF ( opt_vec(iv) .EQ. iscal_offset+10 ) THEN
              CALL IO_WRITE_ASCII(lfile,'Computing first invariant P...')
              CALL FI_INVARIANT_P(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2), wrk1d,wrk2d,wrk3d)

              plot_file = 'InvariantP'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

              CALL IO_WRITE_ASCII(lfile,'Computing second invariant Q...')
              CALL FI_INVARIANT_Q(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d)

              plot_file = 'InvariantQ'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

              CALL IO_WRITE_ASCII(lfile,'Computing third invariant R...')
              CALL FI_INVARIANT_R(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                   dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1), &
                   txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)

              plot_file = 'InvariantR'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

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

           wrk1d(1:jmax,1) = C_0_R
           IF ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN ! update the liquid field in case of air water
              CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array))
           ENDIF
           CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
           dummy =  C_1_R/froude
           wrk3d(1:isize_field) = wrk3d(1:isize_field) *dummy
           
           plot_file = 'Buoyancy'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)

! buoyancy flux along Oy
           txc(1:isize_field,1) = wrk3d(1:isize_field) *q(1:isize_field,2)

           plot_file = 'Cvb'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)

! buoyancy fluctuation
           CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, wrk3d)

           plot_file = 'bPrime'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)

! buoyancy source
           IF ( flag_buoyancy .EQ. 1 ) THEN
              IF ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN
!                 CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array)) !Update the liquid function

!  Create xi intermediate array tmp1
                 dummy = ( body_param(1)*body_param(3)-body_param(2) ) /( C_1_R-body_param(3) )/body_param(3)
                 IF ( dummy .GT. C_SMALL_R ) THEN; dummy2 = body_param(5)*body_param(6)/dummy
                 ELSE;                             dummy2 = C_0_R; ENDIF ! Radiation only

                 txc(1:isize_field,4) = body_param(3) - s(1:isize_field,1) - dummy2*s(1:isize_field,2);

                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
                 CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                      dx,dy,dz, txc(1,4), txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)
               
                 CALL FI_BUOYANCY_SOURCE(ibodyforce, isize_field, body_param, txc(1,4), txc(1,1), wrk3d)
                 
              ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
                 CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array), txc(1,4)) ! calculate xi in tmp1
                 CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, txc(1,4), &
                      dx,dy,dz, txc(1,4), txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)
                 
                 CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, txc(1,2),txc(1,3), wrk3d)
                 
                 dummy = diff *body_param(3) /froude
                 txc(1:isize_field,1) = txc(1:isize_field,1) *txc(1:isize_field,3) *dummy ! evaporation source
                 
              ELSE
                 CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
                 CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                      dx,dy,dz, s, txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)

                 CALL FI_BUOYANCY_SOURCE(ibodyforce, isize_field, body_param, s, txc(1,1), wrk3d)
                 
              ENDIF

              dummy =  visc/schmidt(1) /froude

              txc(1:isize_field,1) = wrk3d(1:isize_field) *dummy
              txc(1:isize_field,1) = LOG(ABS(txc(1:isize_field,1)))
              
              plot_file = 'LnBuoyancySource'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(:,1), tmp_mpi)
              
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
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,i1,kmax, subdomain, surface(1,is), tmp_mpi)
           ENDDO

           WRITE(aux_file,*) itime; aux_file='upper.'//TRIM(ADJUSTL(aux_file))
           CALL DNS_READ_FIELDS(aux_file, i1, imax,i1,kmax, inb_scal,i0, isize_wrk2d, surface, wrk2d)

           DO is = 1,inb_scal 
              subdomain(3) = 1; subdomain(4) = 1
              WRITE(str,*) is; plot_file = 'UpperEnvelope'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,i1,kmax, subdomain, surface(1,is), tmp_mpi)

           ENDDO
        ENDIF

! ###################################################################
! Horizontal divergence
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+14 ) THEN
           CALL IO_WRITE_ASCII(lfile,'Computing horizontal divergence...')
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, q(1,1), txc(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, q(1,3), txc(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
           txc(1:isize_field,1) = txc(1:isize_field,1) + txc(1:isize_field,2)

           plot_file = 'HorizontalDivergence'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

        ENDIF

! ###################################################################
! Turbulent quantities k and \epsilon
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+15 ) THEN

! Initialize the density field in arrat tmp5 if needed
           IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
              txc(:,5) = C_1_R
           ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN ! not yet implemented
           ENDIF
                 
! turbulent dissipation rate into txc4 because I do not need the energy
           CALL IO_WRITE_ASCII(lfile,'Computing dissipation rate...')
           CALL FI_DISSIPATION(i1, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                area, visc, dx,dy,dz, txc(1,5),q(1,1),q(1,2),q(1,3), txc(1,1), &
                txc(1,2),txc(1,3),txc(1,4), wrk1d, wrk1d(1,6),wrk2d,wrk3d)
           txc(1:isize_field,1)=LOG(txc(1:isize_field,1))

           plot_file = 'LnDissipation'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)

! turbulent kinetic energy
           txc(1:isize_field,1) = q(1:isize_field,1); CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, txc(:,1))
           txc(1:isize_field,2) = q(1:isize_field,2); CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, txc(:,2))
           txc(1:isize_field,3) = q(1:isize_field,3); CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, txc(:,3))

           wrk3d(1:isize_field) = txc(1:isize_field,1)*txc(1:isize_field,1) &
                                + txc(1:isize_field,2)*txc(1:isize_field,2) &
                                + txc(1:isize_field,3)*txc(1:isize_field,3)

           plot_file = 'Tke'//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, wrk3d, tmp_mpi)
           
        ENDIF

! ###################################################################
! Radiation
! ###################################################################
        IF (  opt_vec(iv) .EQ. iscal_offset+16) THEN    
           
           IF ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN 
              CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array)) ! Update the liquid field
           ENDIF
           
           IF ( iradiation .NE. EQNS_NONE ) THEN
              CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), txc(:,1), wrk1d,wrk3d)
              
              plot_file = 'Radiation'//time_str(1:MaskSize)
              CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
              
           ENDIF

         ENDIF

! ###################################################################
! Particle density
! ###################################################################
        IF ( opt_vec(iv) .EQ. iscal_offset+17 ) THEN
           CALL DNS_READ_PARTICLE(part_file,l_q)
           l_txc = C_1_R; ! We want density
           CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,1))
           str = 'ParticleDensity'
           plot_file = TRIM(ADJUSTL(str))//time_str(1:MaskSize)
           CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,1), tmp_mpi)
           txc(:,1) = txc(:,1) + 0.00000001
           IF (inb_particle .GT. 3 ) THEN
              IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2 & 
                   .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
!                 DO is=1,inb_lag_aux_field
                 DO is=1,2
                    l_txc(:,1)=l_q(:,3+is) !!! DO WE WANT l_txc(:,is) ???
                    CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,2))   
                    txc(:,2) = txc(:,2)/txc(:,1)
                    plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(is)))//time_str(1:MaskSize)
                    CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,2), tmp_mpi)
                 END DO
              END IF
              IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
                 l_txc(:,1)=l_q(:,inb_particle) !inb_particle is the last component -> residence times in bil_cloud_4
                 CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,2))   
                 plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(3)))//time_str(1:MaskSize)
                 CALL VISUALS_WRITE(plot_file, i0, opt_format, imax,jmax,kmax, subdomain, txc(1,2), tmp_mpi)
              ENDIF
           END IF
           
        ENDIF

     ENDDO

! ###################################################################
! ###################################################################
! print out the time sequence for animation
     WRITE(sRes,100) rtime
     CALL IO_WRITE_ASCII('times.log','Time '//TRIM(ADJUSTL(sRes)))

  ENDDO
  
  CALL DNS_END(0)

  STOP

100 FORMAT(G_FORMAT_R)

END PROGRAM VISUALS_MAIN

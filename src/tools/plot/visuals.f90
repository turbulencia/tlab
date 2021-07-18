#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "VISUALS"

!########################################################################
!#
!# Creating data blocks for visualizations. Derived from ensight.f
!# Partition to be incorporated via a MASK routine before VISUALS_WRITE
!#
!########################################################################
PROGRAM VISUALS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
  USE THERMO_GLOBAL, ONLY : imixture
  USE THERMO_GLOBAL, ONLY : NSP, THERMO_SPNAME
  USE LAGRANGE_GLOBAL
  USE LAGRANGE_ARRAYS
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  ! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 3000
  TINTEGER, PARAMETER :: iopt_size_max  =  20
  TINTEGER, PARAMETER :: igate_size_max =   8
  TINTEGER, PARAMETER :: params_size_max =  2

  ! Additional local arrays
  INTEGER(1), DIMENSION(:),   ALLOCATABLE, SAVE :: gate

  ! -------------------------------------------------------------------
  ! Local variables
  ! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, bakfile
  CHARACTER*32 flow_file, scal_file, part_file, plot_file, time_str
  CHARACTER*64 str

  TINTEGER opt_format
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER ij, is, bcs(2,2)
  TINTEGER iscal_offset, iread_flow, iread_scal, iread_part, idummy, MaskSize
  TREAL diff, dummy
  TINTEGER subdomain(6)

  ! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  TREAL gate_threshold(igate_size_max)

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size, iv
  TINTEGER opt_vec(iopt_size_max)
  TREAL opt_vec2(iopt_size_max)

  TINTEGER params_size
  TREAL params(params_size_max)

  !########################################################################
  !########################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL DNS_START()

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part .EQ. 1 ) CALL PARTICLE_READ_GLOBAL(ifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  ! -------------------------------------------------------------------
  ! File names
  ! -------------------------------------------------------------------
#include "dns_read_times.h"

  ! -------------------------------------------------------------------
  ! Read local options
  ! -------------------------------------------------------------------
  opt_format = 2 ! default values

  IF    ( imixture .EQ. MIXT_TYPE_NONE ) THEN; iscal_offset = 9    ! define iscal_offset
  ELSE;                                        iscal_offset = 9+NSP
  ENDIF

  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamVisuals', '-1', sRes)
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
    WRITE(*,'(I2,A)') iscal_offset+ 7,'. StrainTensor'
    WRITE(*,'(I2,A)') iscal_offset+ 8,'. Strain 2S_ijS_ij (Ln)'
    WRITE(*,'(I2,A)') iscal_offset+ 9,'. StrainEquation'
    WRITE(*,'(I2,A)') iscal_offset+10,'. VelocityGradientInvariants'
    WRITE(*,'(I2,A)') iscal_offset+11,'. Space partition'
    WRITE(*,'(I2,A)') iscal_offset+12,'. Buoyancy'
    WRITE(*,'(I2,A)') iscal_offset+14,'. HorizontalDivergence'
    WRITE(*,'(I2,A)') iscal_offset+15,'. Turbulent quantities'
    WRITE(*,'(I2,A)') iscal_offset+16,'. Radiative forcing'
    WRITE(*,'(I2,A)') iscal_offset+17,'. Relative humidity'
    WRITE(*,'(I2,A)') iscal_offset+18,'. Particle Density'
    WRITE(*,'(I2,A)') iscal_offset+19,'. Thermodynamic quantities'
    WRITE(*,'(I2,A)') iscal_offset+20,'. Analysis of B and V'
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
  iread_flow = 0
  iread_scal = 0
  iread_part = 0
  inb_txc    = 0

  DO iv = 1,iopt_size
    IF ( opt_vec(iv) .EQ. 1              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
    IF ( opt_vec(iv) .EQ. 2              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
    IF ( opt_vec(iv) .EQ. 3              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
    IF ( opt_vec(iv) .EQ. 4              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,3); ENDIF
    IF ( opt_vec(iv) .EQ. 5              ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,1); ENDIF
    IF ( opt_vec(iv) .EQ. 6              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. 7              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
    IF ( opt_vec(iv) .EQ. 8              ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,7); ENDIF
    IF ( opt_vec(iv) .EQ. 9              ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,1); ENDIF
    IF ( opt_vec(iv) .GT. 9 .AND. opt_vec(iv) .LE. iscal_offset   ) THEN
                                                                 iread_scal = 1; inb_txc=MAX(inb_txc,4); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+1 ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+2 ) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,3); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+3 ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+4 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,4); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+5 ) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,7); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+6 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+7 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+8 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,3); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+9 ) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+10) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+12) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,4); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+14) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+15) THEN; iread_flow = 1;                 inb_txc=MAX(inb_txc,6); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+16) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+17) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+18) THEN; iread_part = 1;                 inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+19) THEN;                 iread_scal = 1; inb_txc=MAX(inb_txc,2); ENDIF
    IF ( opt_vec(iv) .EQ. iscal_offset+20) THEN; iread_flow = 1; iread_scal = 1; inb_txc=MAX(inb_txc,7); ENDIF
  ENDDO

  ! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Subdomain', '-1', sRes)

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
  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Format', '-1', sRes)

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
  opt_cond_relative = 0
  igate_size    = 0

  DO iv = 1,iopt_size
    IF ( opt_vec(iv) .EQ. iscal_offset+11 ) THEN
#include "dns_read_partition.h"
      IF ( opt_cond .GT. 1 ) inb_txc = MAX(inb_txc,5)
      EXIT
    ENDIF
  ENDDO

  ! -------------------------------------------------------------------
  ! Allocating memory space
  ! -------------------------------------------------------------------
  ALLOCATE(gate(isize_field))

  isize_wrk3d = isize_txc_field
#ifdef USE_MPI
  isize_wrk3d = isize_wrk3d + isize_field ! more space in wrk3d array needed in IO_WRITE_VISUALS
#endif
  IF ( icalc_part .eq. 1) THEN
    isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  END IF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  IF ( iread_part .EQ. 1 ) THEN ! Particle variables
    inb_part_txc = MAX(inb_part_txc,1)
    CALL PARTICLE_ALLOCATE(C_FILE_LOC)
  ENDIF

  ! -------------------------------------------------------------------
  ! Initialize
  ! -------------------------------------------------------------------
  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  IF ( ifourier .EQ. 1 .AND. inb_txc .GE. 1 ) THEN ! For Poisson solver
    CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( iread_flow .EQ. 1 .AND. inb_txc .GE. 3 ) THEN ! We need array space
    CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

  CALL FI_PROFILES_INITIALIZE(wrk1d) ! Initialize thermodynamic quantities

#ifdef USE_MPI
  CALL VISUALS_MPIO_AUX(opt_format, subdomain)
#else
  io_aux(:)%offset = 0
#endif

  MaskSize    = 6

  ! ###################################################################
  ! Grid
  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
    DO iv = 1,iopt_size
      IF ( opt_vec(iv) .EQ. 0 ) THEN
        CALL ENSIGHT_GRID('grid.ensight', g(1)%size, g(2)%size, g(3)%size, subdomain, g(1)%nodes,g(2)%nodes,g(3)%nodes)
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
    ENDIF

    IF ( iread_flow .EQ. 1 ) THEN ! Flow variables
      WRITE(flow_file,*) itime; flow_file = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(flow_file))
      CALL DNS_READ_FIELDS(flow_file, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)
    ENDIF

    CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

    IF ( iread_part .EQ. 1 ) THEN ! Particle variables
      WRITE(part_file,*) itime; part_file = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(part_file))
    ENDIF

    WRITE(sRes,100) rtime; sRes = 'Physical time '//TRIM(ADJUSTL(sRes))
    CALL IO_WRITE_ASCII(lfile, sRes)

    ! -------------------------------------------------------------------
    ! Calculate intermittency
    ! -------------------------------------------------------------------
    IF      ( opt_cond .EQ. 1 ) THEN ! Read external file
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
        txc(1:isize_field,1) = SQRT(q(1:isize_field,1)**2 + q(1:isize_field,2)**2 + q(1:isize_field,3)**2)
        plot_file = 'VelocityMagnitude'//time_str(1:MaskSize)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

      ENDIF

      ! ###################################################################
      ! Thermodynamic state
      ! ###################################################################
      IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density
          plot_file = 'Density'//time_str(1:MaskSize)
          IF    ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
            CALL THERMO_ANELASTIC_DENSITY(imax,jmax,kmax, s, epbackground,pbackground, txc(1,1))
          ELSE
            wrk1d(1:jmax,1) = C_0_R
            CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
            dummy = C_1_R /froude
            txc(1:isize_field,1) = txc(1:isize_field,1)*dummy + C_1_R
          ENDIF
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 7 .AND. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! temperature
          plot_file = 'Temperature'//time_str(1:MaskSize)
          CALL THERMO_ANELASTIC_TEMPERATURE(imax,jmax,kmax, s, epbackground, txc(1,1))
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

          IF ( damkohler(1) .GT. C_0_R ) THEN ! Supersaturated liquid; this is wrong
            plot_file = 'Supsat'//time_str(1:MaskSize)
            txc(1:isize_field,1:2) = s(1:isize_field,1:2)
            CALL THERMO_AIRWATER_PH(imax,jmax,kmax, txc(1,2),txc(1,1), epbackground,pbackground)
            txc(1:isize_field,3) = ( s(1:isize_field,3) -txc(1:isize_field,3) ) /s(1,3)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,3), wrk3d)
          ENDIF

        ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
          plot_file = 'Pressure'//time_str(1:MaskSize)
          CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

          plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
          CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc(1,1),txc(1,2), wrk3d, wrk2d,wrk3d)
          CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,3), wrk3d, wrk2d,wrk3d)
          CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc(1,1),txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(1:isize_field,2) =-( txc(1:isize_field,2)*q(1:isize_field,1) &
                                 + txc(1:isize_field,3)*q(1:isize_field,2) &
                                 + txc(1:isize_field,4)*q(1:isize_field,3) )
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

          CALL IO_WRITE_ASCII(lfile,'Computing pressure-strain correlation...')
          txc(1:isize_field,2) = txc(1:isize_field,1); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(1,2))

          plot_file = 'PressureStrainX'//time_str(1:MaskSize)
          txc(1:isize_field,3) = q(1:isize_field,1); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,3))
          CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc(1,3), txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(1:isize_field,3) = txc(1:isize_field,2)*txc(1:isize_field,4)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,3), wrk3d)

          plot_file = 'PressureStrainY'//time_str(1:MaskSize)
          txc(1:isize_field,3) = q(1:isize_field,2); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,3))
          CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,3), txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(1:isize_field,3) = txc(1:isize_field,2)*txc(1:isize_field,4)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,3), wrk3d)

          plot_file = 'PressureStrainZ'//time_str(1:MaskSize)
          txc(1:isize_field,3) = q(1:isize_field,3); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,3))
          CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc(1,3), txc(1,4), wrk3d, wrk2d,wrk3d)
          txc(1:isize_field,3) = txc(1:isize_field,2)*txc(1:isize_field,4)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,3), wrk3d)

          plot_file = 'PressureHydrostatic'//time_str(1:MaskSize)
          q = C_0_R
          CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,2), txc(1,3),txc(1,4), txc(1,5), wrk1d,wrk2d,wrk3d)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

          plot_file = 'PressureHydrodynamic'//time_str(1:MaskSize)
          txc(1:isize_field,1) = txc(1:isize_field,1) -txc(1:isize_field,2)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ENDIF

      ELSE ! compressible
        IF      ( opt_vec(iv) .EQ. 6 ) THEN ! density
          plot_file = 'Density'//time_str(1:MaskSize)
          txc(1:isize_field,1) = q(1:isize_field,5)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 7 ) THEN ! temperature
          plot_file = 'Temperature'//time_str(1:MaskSize)
          txc(1:isize_field,1) = q(1:isize_field,7)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ELSE IF ( opt_vec(iv) .EQ. 8 ) THEN ! pressure
          plot_file = 'Pressure'//time_str(1:MaskSize)
          txc(1:isize_field,1) = q(1:isize_field,6)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ENDIF

      ENDIF

      ! ###################################################################
      ! Scalars
      ! ###################################################################
      IF      ( opt_vec(iv) .EQ. 9 ) THEN ! All prognostic scalars
        DO is = 1,inb_scal_array
          WRITE(str,*) is; plot_file = 'Scalar'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
          txc(1:isize_field,1) = s(1:isize_field,is)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
        ENDDO

      ELSE IF ( opt_vec(iv) .GT. 9 .AND. opt_vec(iv) .LE. iscal_offset ) THEN ! Individual and diagnostic scalars
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! s(1,inb_scal+1) contains liquid mass fraction
          IF      ( opt_vec(iv) .EQ. 10 ) THEN ! vapor water mass fraction
            plot_file = TRIM(ADJUSTL(THERMO_SPNAME(1)))//time_str(1:MaskSize)
            s(:,1) = s(:,inb_scal) - s(:,inb_scal+1)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s, wrk3d)

          ELSE IF ( opt_vec(iv) .EQ. 11 ) THEN ! air mass fraction
            plot_file = TRIM(ADJUSTL(THERMO_SPNAME(2)))//time_str(1:MaskSize)
            s(:,1) = C_1_R - s(:,inb_scal)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s, wrk3d)

          ELSE IF ( opt_vec(iv) .EQ. 12 ) THEN ! liquid mass fraction
            plot_file = TRIM(ADJUSTL(THERMO_SPNAME(3)))//time_str(1:MaskSize)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, s(1,inb_scal+1), wrk3d)

          ENDIF

        ELSE ! Plot the chosen species
          is = opt_vec(iv) -9; plot_file = TRIM(ADJUSTL(THERMO_SPNAME(is)))//time_str(1:MaskSize)
          txc(1:isize_field,1) = s(1:isize_field,is)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        ENDIF

      ENDIF

      ! ###################################################################
      ! Scalar Derivatives
      ! ###################################################################
      IF ( opt_vec(iv) .GE. iscal_offset+1 .AND. opt_vec(iv) .LE. iscal_offset+3 ) THEN
        DO is = 1,inb_scal_array
          WRITE(str,*) is; str = 'Scalar'//TRIM(ADJUSTL(str))

          IF ( opt_vec(iv) .EQ. iscal_offset+1 ) THEN
            plot_file = TRIM(ADJUSTL(str))//'GradientVector'//time_str(1:MaskSize)
            CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s(1,is), txc(1,1), wrk3d, wrk2d,wrk3d)
            CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,is), txc(1,2), wrk3d, wrk2d,wrk3d)
            CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s(1,is), txc(1,3), wrk3d, wrk2d,wrk3d)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i3, subdomain, txc(1,1), wrk3d)
          ENDIF

          IF ( opt_vec(iv) .EQ. iscal_offset+2 .OR. opt_vec(iv) .EQ. iscal_offset+3 ) THEN ! Gradient magnitude
            plot_file = TRIM(ADJUSTL(str))//'Gradient'//time_str(1:MaskSize)
            CALL FI_GRADIENT(imax,jmax,kmax, s(1,is), txc(1,1),txc(1,2), wrk2d,wrk3d)
            IF ( opt_vec(iv) .EQ. iscal_offset+2 ) THEN
              plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
              txc(1:isize_field,1) = LOG(txc(1:isize_field,1)+C_SMALL_R)
            ENDIF
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
          ENDIF

          IF ( opt_vec(iv) .EQ. iscal_offset+3 .AND. is .LE. inb_scal ) THEN ! Scalar gradient equation
            IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
            ELSE;                                  diff = visc/schmidt(is)
            ENDIF

            CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
            plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
            CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s(1,is), q(1,1),q(1,2),q(1,3), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

            CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
            plot_file = TRIM(ADJUSTL(str))//'GradientDiffusion'//time_str(1:MaskSize)
            CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s(1,is), &
              txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
            txc(1:isize_field,1) = diff *txc(1:isize_field,1)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

          ENDIF

        ENDDO

      ENDIF

      ! ###################################################################
      ! Velocity Derivatives
      ! ###################################################################
      IF ( opt_vec(iv) .EQ. iscal_offset+4 ) THEN ! VorticityVector
        plot_file = 'VorticityVector'//time_str(1:MaskSize)
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i3, subdomain, txc(1,1), wrk3d)
      ENDIF

      IF ( opt_vec(iv) .EQ. iscal_offset+5 .OR. opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! Enstrophy
        plot_file = 'Enstrophy'//time_str(1:MaskSize)
        CALL FI_VORTICITY(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3), wrk2d,wrk3d)
        IF ( opt_vec(iv) .EQ. iscal_offset+5 ) THEN ! Natural log
          plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
          txc(1:isize_field,1)=LOG(txc(1:isize_field,1)+C_SMALL_R)
        ENDIF
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'LnPotentialEnstrophy'//time_str(1:MaskSize)
        IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
          CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,4))
        ELSE
          wrk1d(1:jmax,1) = C_0_R
          CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,4), wrk1d)
        ENDIF
        dummy =  C_1_R /froude
        txc(1:isize_field,4) = txc(1:isize_field,4) *dummy
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc(1,4),txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,4),txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc(1,4),txc(1,3), wrk3d, wrk2d,wrk3d)
        CALL FI_CURL(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1)*txc(1:isize_field,4) &
                             + txc(1:isize_field,2)*txc(1:isize_field,5) &
                             + txc(1:isize_field,3)*txc(1:isize_field,6)
        txc(1:isize_field,1) = LOG( txc(1:isize_field,1)*txc(1:isize_field,1) +C_SMALL_R )
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
      ENDIF

      IF ( opt_vec(iv) .EQ. iscal_offset+6 ) THEN ! EnstrophyEquation
        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
        plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
        CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
        plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
        CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),&
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,1) = visc *txc(1:isize_field,1)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

      ENDIF

      ! -------------------------------------------------------------------
      IF ( opt_vec(iv) .EQ. iscal_offset+7 ) THEN ! Strain Tensor
        plot_file = 'StrainTensor'//time_str(1:MaskSize)
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i6, subdomain, txc(1,1), wrk3d)
      ENDIF

      IF ( opt_vec(iv) .EQ. iscal_offset+8 .OR. opt_vec(iv) .EQ. iscal_offset+9 ) THEN ! Strain
        plot_file = 'Strain'//time_str(1:MaskSize)
        CALL FI_STRAIN(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)
        txc(1:isize_field,1) = C_2_R *txc(1:isize_field,1)
        IF ( opt_vec(iv) .EQ. iscal_offset+8 ) THEN ! Natural log
          plot_file = 'Ln'//TRIM(ADJUSTL(plot_file))
          txc(1:isize_field,1)=LOG( txc(1:isize_field,1) +C_SMALL_R )
        ENDIF
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
      ENDIF

      IF ( opt_vec(iv) .EQ. iscal_offset+9 ) THEN ! StrainEquation (I need the pressure)
        CALL IO_WRITE_ASCII(lfile,'Computing strain pressure...')
        plot_file = 'StrainPressure'//time_str(1:MaskSize)
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
          CALL IO_WRITE_ASCII(efile,'VISUALS. Strain eqn for incompressible undeveloped.')
          CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ELSE
          txc(:,6) = q(:,6)
        ENDIF
        CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, q(1,1),q(1,2),q(1,3),txc(1,6), &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        txc(1:isize_field,1) = C_2_R *txc(1:isize_field,1)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
        plot_file = 'StrainProduction'//time_str(1:MaskSize)
        CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
        txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,1)=C_2_R*txc(1:isize_field,1)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
        plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
        CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,1) = C_2_R *visc *txc(1:isize_field,1)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

      ENDIF

      ! -------------------------------------------------------------------
      IF ( opt_vec(iv) .EQ. iscal_offset+10 ) THEN ! Invariants
        plot_file = 'InvariantP'//time_str(1:MaskSize)
        CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'InvariantQ'//time_str(1:MaskSize)
        CALL FI_INVARIANT_Q(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'InvariantR'//time_str(1:MaskSize)
        CALL FI_INVARIANT_R(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
          txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
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
        plot_file = 'Buoyancy'//time_str(1:MaskSize)
        IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
          CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,1))
        ELSE
          wrk1d(1:jmax,1) = C_0_R
          CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
        ENDIF
        dummy =  C_1_R /froude
        txc(1:isize_field,1) = txc(1:isize_field,1) *dummy
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'Fvb'//time_str(1:MaskSize)     ! buoyancy flux along Oy
        txc(1:isize_field,2) = txc(1:isize_field,1) *q(1:isize_field,2)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

        plot_file = 'bPrime'//time_str(1:MaskSize)  ! buoyancy fluctuation
        CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(1,1))
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'Cvb'//time_str(1:MaskSize)     ! Covariance between b and v
        txc(1:isize_field,2) = q(1:isize_field,2); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,2))
        txc(1:isize_field,2) = txc(1:isize_field,1) *txc(1:isize_field,2)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)


        plot_file = 'LnBuoyancySource'//time_str(1:MaskSize)
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
          CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, txc(1,1),txc(1,2),txc(1,3))
          CALL FI_GRADIENT(imax,jmax,kmax, txc(1,1),txc(1,2), txc(1,4), wrk2d,wrk3d)
          dummy = buoyancy%parameters(inb_scal_array)
          txc(1:isize_field,2) = txc(1:isize_field,2) *txc(1:isize_field,3) *dummy
        ELSE
          CALL FI_GRADIENT(imax,jmax,kmax, s, txc(1,1),txc(1,2), wrk2d,wrk3d)
          CALL FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, txc(1,1), txc(1,2))
        ENDIF
        dummy =  visc /schmidt(1) /froude
        txc(1:isize_field,1) = txc(1:isize_field,2) *dummy
        txc(1:isize_field,1) = LOG( ABS(txc(1:isize_field,1)) +C_SMALL_R )
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

      ENDIF

      ! ###################################################################
      IF ( opt_vec(iv) .EQ. iscal_offset+14 ) THEN
        plot_file = 'HorizontalDivergence'//time_str(1:MaskSize)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), q(1,1), txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), q(1,3), txc(1,1), wrk3d, wrk2d,wrk3d)
        txc(1:isize_field,1) = txc(1:isize_field,1) + txc(1:isize_field,2)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
      ENDIF

      ! ###################################################################
      IF ( opt_vec(iv) .EQ. iscal_offset+15 ) THEN ! Turbulent quantities
        plot_file = 'LnDissipation'//time_str(1:MaskSize)
        CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), &
          txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
        txc(1:isize_field,1) = LOG( txc(1:isize_field,1) +C_SMALL_R )
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'Tke'//time_str(1:MaskSize)
        txc(1:isize_field,1) = q(1:isize_field,1); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,1))
        txc(1:isize_field,2) = q(1:isize_field,2); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,2))
        txc(1:isize_field,3) = q(1:isize_field,3); CALL FI_FLUCTUATION_INPLACE(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, txc(:,3))
        txc(1:isize_field,4) = C_05_R *( txc(1:isize_field,1)**2 +txc(1:isize_field,2)**2 +txc(1:isize_field,3)**2 )
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,4), wrk3d)

        plot_file = 'ReynoldsTensor'//time_str(1:MaskSize)
        txc(1:isize_field,4) = txc(1:isize_field,1) *txc(1:isize_field,2)
        txc(1:isize_field,5) = txc(1:isize_field,1) *txc(1:isize_field,3)
        txc(1:isize_field,6) = txc(1:isize_field,2) *txc(1:isize_field,3)
        txc(1:isize_field,1) = txc(1:isize_field,1) *txc(1:isize_field,1)
        txc(1:isize_field,2) = txc(1:isize_field,2) *txc(1:isize_field,2)
        txc(1:isize_field,3) = txc(1:isize_field,3) *txc(1:isize_field,3)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i6, subdomain, txc(1,1), wrk3d)

      ENDIF

      ! ###################################################################
      IF ( opt_vec(iv) .EQ. iscal_offset+16 ) THEN
        DO is = 1,inb_scal

          IF ( radiation%active(is) ) THEN
            WRITE(str,*) is; plot_file = 'Radiation'//TRIM(ADJUSTL(str))//time_str(1:MaskSize)
            IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
              CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, s(1,radiation%scalar(is)), txc(1,2))
              CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), txc(1,2),                 txc(1,1), wrk1d,wrk3d)
              CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, ribackground, txc(1,1))
            ELSE
              CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(1)), txc(1,1), wrk1d,wrk3d)
            ENDIF
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
          ENDIF

        ENDDO
      ENDIF

      ! ###################################################################
      IF (  opt_vec(iv) .EQ. iscal_offset+17) THEN
        plot_file = 'RelativeHumidity'//time_str(1:MaskSize)
        CALL THERMO_ANELASTIC_RELATIVEHUMIDITY(imax,jmax,kmax, s, epbackground,pbackground, wrk3d, txc(1,1))
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)
      ENDIF

      ! ###################################################################
      IF ( opt_vec(iv) .EQ. iscal_offset+18 ) THEN ! Particle information
        plot_file = 'ParticleDensity'//time_str(1:MaskSize)
        CALL IO_READ_PARTICLE(part_file, l_g, l_q)
        l_txc = C_1_R; ! We want density
        CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,1), wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 )  THEN
          txc(:,1) = txc(:,1) + 0.00000001
          DO is=1,2
            plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(is)))//time_str(1:MaskSize)
            CALL PARTICLE_TO_FIELD(l_q, l_q(1,3+is), txc(1,2), wrk2d,wrk3d)
            txc(:,2) = txc(:,2) /txc(:,1)
            CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
          END DO
        END IF

        IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN
          !inb_part_array is the last component -> residence times in bil_cloud_4
          plot_file = TRIM(ADJUSTL(LAGRANGE_SPNAME(3)))//time_str(1:MaskSize)
          CALL PARTICLE_TO_FIELD(l_q, l_q(1,inb_part_array-1), txc(1,2), wrk2d,wrk3d)
          CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
        ENDIF

      ENDIF

      ! ###################################################################
      IF (  opt_vec(iv) .EQ. iscal_offset+19 ) THEN ! Thermodynamic quantities
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

      ! ###################################################################
      IF (  opt_vec(iv) .EQ. iscal_offset+20 ) THEN
        plot_file = 'LaplacianV'//time_str(1:MaskSize)
        CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), q(1,2),txc(1,4), txc(1,5), wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), q(1,2),txc(1,3), txc(1,5), wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), q(1,2),txc(1,2), txc(1,5), wrk2d,wrk3d)
        txc(1:isize_field,2) = txc(1:isize_field,2) +txc(1:isize_field,3) +txc(1:isize_field,4)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

        plot_file = 'Buoyancy'//time_str(1:MaskSize)
        IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
          CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,1))
        ELSE
          wrk1d(1:jmax,1) = C_0_R
          CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,1), wrk1d)
        ENDIF
        dummy =  C_1_R /froude
        txc(1:isize_field,1) = txc(1:isize_field,1) *dummy
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'LaplacianB'//time_str(1:MaskSize)
        CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), txc(1,1),txc(1,4), txc(1,5), wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,3), txc(1,5), wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), txc(1,1),txc(1,2), txc(1,5), wrk2d,wrk3d)
        txc(1:isize_field,2) = txc(1:isize_field,2) +txc(1:isize_field,3) +txc(1:isize_field,4)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)

        plot_file = 'Pressure'//time_str(1:MaskSize)
        bbackground = C_0_R
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,1), wrk3d)

        plot_file = 'PressureGradientY'//time_str(1:MaskSize)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,1),txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL IO_WRITE_VISUALS(plot_file, opt_format, imax,jmax,kmax, i1, subdomain, txc(1,2), wrk3d)
      ENDIF

    ENDDO

  ENDDO

  100 FORMAT(G_FORMAT_R)
  CALL DNS_STOP(0)
END PROGRAM VISUALS

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

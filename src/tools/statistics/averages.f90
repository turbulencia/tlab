#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

#define C_FILE_LOC "AVERAGES"

PROGRAM AVERAGES

  USE DNS_TYPES,     ONLY : pointers_dt
  USE DNS_CONSTANTS, ONLY : efile,lfile,gfile, tag_flow,tag_scal,tag_part
  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
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
  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE         :: x,y,z
  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE, TARGET :: q, s, txc
  TREAL, DIMENSION(:),      ALLOCATABLE, SAVE         :: wrk1d,wrk3d
  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE         :: wrk2d

  TREAL, DIMENSION(:),      ALLOCATABLE, SAVE         :: mean, y_aux
  INTEGER(1), DIMENSION(:), ALLOCATABLE, SAVE         :: gate
  TREAL, DIMENSION(:,:,:),  ALLOCATABLE, SAVE         :: surface

  TYPE(pointers_dt), DIMENSION(16) :: vars

  ! Particle data
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_txc

  ! -------------------------------------------------------------------
  ! Local variables
  ! -------------------------------------------------------------------
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 varname(16)
  CHARACTER*64 str, line

  TINTEGER opt_main, opt_block, opt_order
  TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
  TINTEGER nfield, isize_wrk3d, is, ij, k, n, bcs(2,2)
  TREAL eloc1, eloc2, eloc3, cos1, cos2, cos3, dummy
  TINTEGER jmax_aux, iread_flow, iread_scal, ierr, idummy

  ! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  TREAL gate_threshold(igate_size_max)
  INTEGER(1) gate_level

  ! Reading variables
  CHARACTER*512 sRes

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

  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL(inifile)
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL AVG_DEFS_TEMPORAL
#ifdef USE_MPI
  ENDIF
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

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamAverages', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
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

     IF ( opt_main .GT. 2 ) THEN
        WRITE(*,*) 'Gate level to be used ?'
        READ(*,*) gate_level
        WRITE(*,*) 'Number of moments ?'
        READ(*,*) opt_order
     ENDIF

#endif
  ELSE
     opt_main = INT(opt_vec(1))
     IF ( iopt_size .GE. 2 ) opt_block = INT(opt_vec(2))
     IF ( iopt_size .GE. 3 ) gate_level= INT(opt_vec(3),KIND=1)
     IF ( iopt_size .GE. 3 ) opt_order = INT(opt_vec(4))

  ENDIF

  IF ( opt_main .LT. 0 ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'AVERAGES. Missing input [ParamAverages] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  IF ( opt_block .LT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile, 'AVERAGES. Invalid value of opt_block.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

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
     nfield = 10
     iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,8)
  CASE ( 4 )
     nfield = 6 +inb_scal
     iread_flow = 1; iread_scal = 1; inb_txc = MAX(inb_txc,3)
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) inb_txc = MAX(inb_txc,6)
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

  IF ( opt_main .EQ. 2 .OR. gate_level .NE. 0 ) THEN
#include "dns_read_partition.h"
     IF ( opt_cond .GT. 1 ) inb_txc = MAX(inb_txc,5)
  ENDIF

  ! -------------------------------------------------------------------
  ! Allocating memory space
  ! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(gate(isize_field))

  ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
  jmax_aux = g(2)%size/opt_block
  ALLOCATE(y_aux(g(2)%size)) ! Reduced vertical grid

  ! Subarray information to read and write envelopes
  IF ( opt_main .EQ. 2 .AND. opt_cond .GT. 1 .AND. igate_size .NE. 0 ) THEN
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

  ENDIF

  IF ( opt_main .EQ. 1 ) THEN
     ALLOCATE(mean(jmax*MAX_AVG_TEMPORAL))
  ELSE
     ALLOCATE(mean(2*opt_order*nfield))
  ENDIF

  isize_wrk3d = MAX(isize_field,2*opt_order*nfield*jmax)
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)
  IF ( icalc_part .EQ. 1) THEN
     isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  END IF
#include "dns_alloc_arrays.h"

  IF ( icalc_part .EQ. 1 ) THEN
#include "dns_alloc_larrays.h"
  ENDIF

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

  IF ( iread_flow .EQ. 1 ) THEN
     u => q(:,1)
     v => q(:,2)
     w => q(:,3)
  ENDIF

  ! ###################################################################
  ! Postprocess given list of files
  ! ###################################################################
  DO it = 1,itime_size
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

        IF ( opt_main .EQ. 2 ) rtime = params(1)

     ELSE IF ( opt_cond .GT. 1 ) THEN
        opt_cond_scal = 1
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
           opt_cond_scal = inb_scal_array
        ENDIF

        CALL IO_WRITE_ASCII(lfile,'Calculating partition...')
        CALL FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
             imax,jmax,kmax, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)
     ENDIF

     ! -------------------------------------------------------------------
     ! Type of averages
     ! -------------------------------------------------------------------
     SELECT CASE ( opt_main )

        ! ###################################################################
        ! Conventional statistics
        ! ###################################################################
     CASE ( 1 )
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,3), txc(1,1),txc(1,2), txc(1,4), wrk1d,wrk2d,wrk3d)

        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, q(1,4), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, q(1,5), wrk3d)! density
           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), q(1,7), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), q(1,7), q(1,6))
           IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, q(1,7), q(1,8))

        ENDIF

        IF ( icalc_flow .EQ. 1 ) THEN
           CALL AVG_FLOW_XZ(q,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
                txc(1,7),txc(1,8),txc(1,9), mean, wrk1d,wrk2d,wrk3d)
        ENDIF

        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
              CALL AVG_SCAL_XZ(is, q,s, s(1,is), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
           ENDDO

           ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
           IF ( buoyancy%TYPE .EQ. EQNS_BOD_QUADRATIC   .OR. &
                buoyancy%TYPE .EQ. EQNS_BOD_BILINEAR    .OR. &
                imixture .EQ. MIXT_TYPE_AIRWATER        .OR. &
                imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
              dummy = C_1_R/froude
              IF ( buoyancy%TYPE .EQ. EQNS_EXPLICIT ) THEN
                 CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, txc(1,7))
              ELSE
                 wrk1d(1:jmax) = C_0_R
                 CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, txc(1,7), wrk1d)
              ENDIF
              txc(1:isize_field,7) = txc(1:isize_field,7) *dummy

              CALL AVG_SCAL_XZ(is, q,s, txc(1,7), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)

           ENDIF

           IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
                 is = is + 1
                 CALL THERMO_ANELASTIC_THETA_L(imax,jmax,kmax, s, epbackground,pbackground, txc(1,7))
                 !                 CALL THERMO_ANELASTIC_STATIC_CONSTANTCP(imax,jmax,kmax, s, epbackground, txc(1,7))
                 CALL AVG_SCAL_XZ(is, q,s, txc(1,7), &
                      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
              ENDIF
           ENDIF

        ENDIF

        ! Lagrange Liquid and Liquid without diffusion
        IF ( icalc_part .EQ. 1 ) THEN
           IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN
              WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
              CALL IO_READ_PARTICLE(fname, l_g, l_q)

              l_txc = C_1_R; ! We want density
              CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,7), wrk2d,wrk3d)

              txc(:,7) = txc(:,7) + 0.00000001
              idummy = inb_part - 3 ! # scalar properties solved in the lagrangian
              DO is = inb_scal_array +1 +1, inb_scal_array+1 +idummy
                 sbg(is)%mean = C_1_R; sbg(is)%delta = C_0_R; sbg(is)%ymean = sbg(1)%ymean; schmidt(is) = schmidt(1)
                 l_txc(:,1)=l_q(:,3+is-inb_scal_array-1) !!! DO WE WANT l_txc(:,is) ???
                 CALL PARTICLE_TO_FIELD(l_q, l_txc, txc(1,8), wrk2d,wrk3d)
                 txc(:,8) = txc(:,8) /txc(:,7)
                 sbg(is)%mean  = sbg(inb_scal_array)%mean
                 sbg(is)%delta = sbg(inb_scal_array)%delta
                 sbg(is)%ymean = sbg(inb_scal_array)%ymean
                 CALL AVG_SCAL_XZ(is, q,s, txc(1,8), &
                      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), mean, wrk1d,wrk2d,wrk3d)
              ENDDO
           ENDIF
        ENDIF
        ! DO from is inb_scal_array+1 to inb_paticle+inb_scal_array
        ! PARTICLE TO FIELD to txc(5)
        ! CALL AVG_SCAL_XZ on new field (substutute by s(1,is)
        !ENDDO

        ! ###################################################################
        ! Partition of field
        ! ###################################################################
     CASE ( 2 )
        DO n = 1,igate_size
           WRITE(varname(n),*) n; varname(n) = 'Partition'//TRIM(ADJUSTL(varname(n)))
        ENDDO
        WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
        CALL INTER2D_N(fname, varname, rtime, imax,jmax,kmax, igate_size, y, gate)

        IF ( opt_cond .GT. 1 ) THEN ! write only if the gate information has not been read
           WRITE(fname,*) itime; fname = 'gate.'//TRIM(ADJUSTL(fname))
           params(1) = rtime; params(2) = M_REAL(igate_size); params_size = 2
           CALL IO_WRITE_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)

           DO n = 1,igate_size
              gate_level = INT(n,KIND=1)
              CALL BOUNDARY_LOWER_INT1(imax,jmax,kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1,2))
              DO k = 1,kmax ! rearranging
                 ij = (k-1) *imax +1
                 surface(1:imax,n           ,k) = wrk2d(ij:ij+imax-1,1)
              ENDDO
              CALL BOUNDARY_UPPER_INT1(imax,jmax,kmax, gate_level, y, gate, wrk3d, wrk2d, wrk2d(1,2))
              DO k = 1,kmax ! rearranging
                 ij = (k-1) *imax +1
                 surface(1:imax,n+igate_size,k) = wrk2d(ij:ij+imax-1,1)
              ENDDO
           ENDDO
           varname = ''
           WRITE(fname,*) itime; fname = 'envelopesJ.'//TRIM(ADJUSTL(fname))
           CALL IO_WRITE_SUBARRAY4(IO_SUBARRAY_ENVELOPES, fname, varname, surface, io_sizes, wrk3d)

        ENDIF

        ! ###################################################################
        ! Momentum equation
        ! ###################################################################
     CASE ( 3 )
       nfield = 0
       nfield = nfield+1; vars(nfield)%field => q(:,1);   varname(nfield) = 'U'
       nfield = nfield+1; vars(nfield)%field => q(:,3);   varname(nfield) = 'W'

       nfield = nfield+1; vars(nfield)%field => txc(:,1); varname(nfield) = 'Uy'
       nfield = nfield+1; vars(nfield)%field => txc(:,2); varname(nfield) = 'Uyy'
       CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,1), txc(1,2), txc(1,1), wrk2d,wrk3d)
       nfield = nfield+1; vars(nfield)%field => txc(:,3); varname(nfield) = 'Wy'
       nfield = nfield+1; vars(nfield)%field => txc(:,4); varname(nfield) = 'Wyy'
       CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,3), txc(1,4), txc(1,3), wrk2d,wrk3d)

       nfield = nfield+1; vars(nfield)%field => txc(:,5); varname(nfield) = '(UV)y'
       txc(1:isize_field,6) = q(1:isize_field,2) *q(1:isize_field,1)
       CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,6), txc(1,5), wrk3d, wrk2d,wrk3d)
       nfield = nfield+1; vars(nfield)%field => txc(:,6); varname(nfield) = '(WV)y'
       txc(1:isize_field,7) = q(1:isize_field,2) *q(1:isize_field,3)
       CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc(1,7), txc(1,6), wrk3d, wrk2d,wrk3d)

       nfield = nfield+1; vars(nfield)%field => txc(:,7); varname(nfield) = 'VUy'
       txc(1:isize_field,7) = q(1:isize_field,2) *txc(1:isize_field,1)
       nfield = nfield+1; vars(nfield)%field => txc(:,8); varname(nfield) = 'VWy'
       txc(1:isize_field,8) = q(1:isize_field,2) *txc(1:isize_field,3)

       IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
          DO is = 1,nfield
             CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
          ENDDO
       ENDIF

       WRITE(fname,*) itime; fname='cavg'//TRIM(ADJUSTL(fname))
       CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
            nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Main variables
        ! ###################################################################
     CASE ( 4 )
        nfield = 0
        nfield = nfield+1; vars(nfield)%field => u(:); varname(nfield) = 'U'
        nfield = nfield+1; vars(nfield)%field => v(:); varname(nfield) = 'V'
        nfield = nfield+1; vars(nfield)%field => w(:); varname(nfield) = 'W'

        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
           nfield = nfield+1; vars(nfield)%field => txc(:,1); varname(nfield) = 'P'

        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density
           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1))
           nfield = nfield+1; vars(nfield)%field => txc(:,2); varname(nfield) = 'R'
           nfield = nfield+1; vars(nfield)%field => txc(:,1); varname(nfield) = 'P'
           nfield = nfield+1; vars(nfield)%field => txc(:,3); varname(nfield) = 'T'

        ENDIF

        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
              nfield = nfield+1; vars(nfield)%field => s(:,is); WRITE(varname(nfield),*) is; varname(nfield) = 'Scalar'//TRIM(ADJUSTL(varname(nfield)))
           ENDDO
        ENDIF

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='cavg'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Enstrophy equation
        ! ###################################################################
     CASE ( 5 )
        CALL IO_WRITE_ASCII(lfile,'Computing baroclinic term...')
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           IF ( buoyancy%TYPE .EQ. EQNS_NONE ) THEN
              txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
           ELSE
              ! calculate buoyancy vector along Oy
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

           CALL THERMO_CALORIC_TEMPERATURE&
                (imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE&
                (imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1)) ! pressure in txc1
           ! result vector in txc4, txc5, txc6
           CALL FI_VORTICITY_BAROCLINIC(imax,jmax,kmax, txc(1,2),txc(1,1), txc(1,4), txc(1,3),txc(1,7), wrk2d,wrk3d)
        ENDIF
        ! result vector in txc1, txc2, txc3
        CALL FI_CURL(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3), txc(1,7), wrk2d,wrk3d)
        ! scalar product, store in txc8
        DO ij = 1,isize_field
           txc(ij,8) = txc(ij,1)*txc(ij,4) + txc(ij,2)*txc(ij,5) + txc(ij,3)*txc(ij,6)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
        CALL FI_VORTICITY_PRODUCTION(imax,jmax,kmax, u,v,w, txc(1,1),&
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
        CALL FI_VORTICITY_DIFFUSION(imax,jmax,kmax, u,v,w, txc(1,2),&
             txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,2) = visc*txc(ij,2)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy...')
        CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4),txc(1,5), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing dilatation term...')
        CALL FI_INVARIANT_P(imax,jmax,kmax, u,v,w, txc(1,4),txc(1,5), wrk2d,wrk3d)

        DO ij = 1,isize_field
           txc(ij,6) = txc(ij,4)*txc(ij,3) ! -w^2 div(u)
           txc(ij,5) = txc(ij,1)/txc(ij,3) ! production rate
           txc(ij,4) = LOG(txc(ij,3))      ! ln(w^2)
        ENDDO

        vars(1)%field => txc(:,3); varname(1) = 'EnstrophyW_iW_i'
        vars(2)%field => txc(:,4); varname(2) = 'LnEnstrophyW_iW_i'
        vars(3)%field => txc(:,1); varname(3) = 'ProductionW_iW_jS_ij'
        vars(4)%field => txc(:,2); varname(4) = 'DiffusionNuW_iLapW_i'
        vars(5)%field => txc(:,6); varname(5) = 'DilatationMsW_iW_iDivU'
        vars(6)%field => txc(:,8); varname(6) = 'Baroclinic'
        vars(7)%field => txc(:,5); varname(7) = 'RateAN_iN_jS_ij'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgW2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Strain equation
        ! ###################################################################
     CASE ( 6 )
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
        CALL FI_STRAIN_PRESSURE(imax,jmax,kmax, u,v,w, txc(1,1), &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,1) = C_2_R *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
        CALL FI_STRAIN_PRODUCTION(imax,jmax,kmax, u,v,w, &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        txc(1:isize_field,2) = C_2_R *txc(1:isize_field,2)

        CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
        CALL FI_STRAIN_DIFFUSION(imax,jmax,kmax, u,v,w, &
             txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk2d,wrk3d)
        txc(1:isize_field,3) = C_2_R *visc *txc(1:isize_field,3)

        CALL IO_WRITE_ASCII(lfile,'Computing strain...')
        CALL FI_STRAIN(imax,jmax,kmax, u,v,w, txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        txc(1:isize_field,4) = C_2_R *txc(1:isize_field,4)
        txc(1:isize_field,5) = LOG( txc(1:isize_field,4) )

        vars(1)%field => txc(:,4); varname(1) = 'Strain2S_ijS_i'
        vars(2)%field => txc(:,5); varname(2) = 'LnStrain2S_ijS_i'
        vars(3)%field => txc(:,2); varname(3) = 'ProductionMs2S_ijS_jkS_ki'
        vars(4)%field => txc(:,3); varname(4) = 'DiffusionNuS_ijLapS_ij'
        vars(5)%field => txc(:,1); varname(5) = 'Pressure2S_ijP_ij'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgS2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Scalar gradient equation
        ! ###################################################################
     CASE ( 7 )
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
        CALL FI_GRADIENT_PRODUCTION(imax,jmax,kmax, s, u,v,w, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        ! array u used as auxiliar
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
        CALL FI_GRADIENT_DIFFUSION(imax,jmax,kmax, s, &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),u, wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,2) = visc /schmidt(inb_scal) *txc(ij,2)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
        CALL FI_GRADIENT(imax,jmax,kmax, s,txc(1,3), txc(1,4), wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,5) = txc(ij,1)/txc(ij,3)
           txc(ij,4) = LOG(txc(ij,3))
        ENDDO

        vars(1)%field => txc(:,3); varname(1) = 'GradientG_iG_i'
        vars(2)%field => txc(:,4); varname(2) = 'LnGradientG_iG_i'
        vars(3)%field => txc(:,1); varname(3) = 'ProductionMsG_iG_jS_ij'
        vars(4)%field => txc(:,2); varname(4) = 'DiffusionNuG_iLapG_i'
        vars(5)%field => txc(:,5); varname(5) = 'StrainAMsN_iN_jS_ij'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgG2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Velocity gradient invariants
        ! ###################################################################
     CASE ( 8 )
        CALL IO_WRITE_ASCII(lfile,'Computing third invariant R...')
        CALL FI_INVARIANT_R(imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing second invariant Q...')
        CALL FI_INVARIANT_Q(imax,jmax,kmax, u,v,w, txc(1,2), txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing first invariant P...')
        CALL FI_INVARIANT_P(imax,jmax,kmax, u,v,w, txc(1,3), txc(1,4), wrk2d,wrk3d)

        vars(1)%field => txc(:,3); varname(1) = 'InvariantP'
        vars(2)%field => txc(:,2); varname(2) = 'InvariantQ'
        vars(3)%field => txc(:,1); varname(3) = 'InvariantR'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgInv'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Scalar gradient components
        ! ###################################################################
     CASE ( 9 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,3), wrk3d, wrk2d,wrk3d)
        ! Angles; s array is overwritten to save space
        DO ij = 1,isize_field
           dummy = txc(ij,2)/SQRT(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
           txc(ij,4) = ASIN(dummy)                 ! with Oy
           s(ij,1)  = ATAN2(txc(ij,3),txc(ij,1))  ! with Ox in plane xOz
        ENDDO

        vars(1)%field => txc(:,1); varname(1) = 'GradientX'
        vars(2)%field => txc(:,2); varname(2) = 'GradientY'
        vars(3)%field => txc(:,3); varname(3) = 'GradientZ'
        vars(4)%field => s(:,1);   varname(4) = 'Theta'
        vars(5)%field => txc(:,4); varname(5) = 'Phi'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgGi'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! eigenvalues of rate-of-strain tensor
        ! ###################################################################
     CASE ( 10 )
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...') ! txc6-txc9
        CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))

        vars(1)%field => txc(:,7); varname(1) = 'Lambda1'
        vars(2)%field => txc(:,8); varname(2) = 'Lambda2'
        vars(3)%field => txc(:,9); varname(3) = 'Lambda3'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgEig'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! eigenframe of rate-of-strain tensor
        ! ###################################################################
     CASE ( 11 )
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imax,jmax,kmax, u,v,w, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...')           ! txc7-txc9
        CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))

        CALL IO_WRITE_ASCII(lfile,'Computing eigenframe...')            ! txc1-txc6
        CALL FI_TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1,1), txc(1,7))

        ! local direction cosines of vorticity vector
        CALL IO_WRITE_ASCII(lfile,'Computing vorticity vector...')      ! txc7-txc9
        CALL FI_CURL(imax,jmax,kmax, u,v,w, txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk2d,wrk3d)

        DO ij = 1,isize_field
           dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           u(ij) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           v(ij) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           w(ij) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        ENDDO

        vars(1)%field => u; varname(1) = 'cos(w,lambda1)'
        vars(2)%field => v; varname(2) = 'cos(w,lambda2)'
        vars(3)%field => w; varname(3) = 'cos(w,lambda3)'

        ! local direction cosines of scalar gradient vector
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...') ! txc7-txc9
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s, txc(1,7), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s, txc(1,8), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s, txc(1,9), wrk3d, wrk2d,wrk3d)

        DO ij = 1,isize_field
           dummy = SQRT(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
           txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3
        ENDDO

        vars(4)%field => txc(:,7); varname(4) = 'cos(G,lambda1)'
        vars(5)%field => txc(:,8); varname(5) = 'cos(G,lambda2)'
        vars(6)%field => txc(:,9); varname(6) = 'cos(G,lambda3)'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgCos'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! longitudinal velocity derivatives
        ! ###################################################################
     CASE ( 12 )
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, txc(1,2), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, txc(1,3), wrk3d, wrk2d,wrk3d)

        vars(1)%field => txc(:,1); varname(1) = 'dudx'
        vars(2)%field => txc(:,2); varname(2) = 'dvdy'
        vars(3)%field => txc(:,3); varname(3) = 'dwdz'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgDer'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Vertical fluxes
        ! ###################################################################
     CASE ( 13 )
        nfield = 0

        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, txc(:,1), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, txc(:,2), wrk3d, wrk2d,wrk3d)
        txc(:,1) = ( txc(:,1) + txc(:,2) ) *visc
        nfield = nfield+1; vars(nfield)%field => txc(:,1); varname(nfield) = 'tauyx'

        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, txc(:,2), wrk3d, wrk2d,wrk3d)
        txc(:,2) =   txc(:,2) *C_2_R       *visc
        nfield = nfield+1; vars(nfield)%field => txc(:,2); varname(nfield) = 'tauyy'

        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, txc(:,3), wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, txc(:,4), wrk3d, wrk2d,wrk3d)
        txc(:,3) = ( txc(:,3) + txc(:,4) ) *visc
        nfield = nfield+1; vars(nfield)%field => txc(:,3); varname(nfield) = 'tauyz'

        DO is = 1,inb_scal_array
           CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(:,is), txc(:,3+is), wrk3d, wrk2d,wrk3d)
           txc(:,3+is) =   txc(:,3+is) *visc /schmidt(inb_scal)
           nfield = nfield+1; vars(nfield)%field => txc(:,3+is); WRITE(varname(nfield),*) is; varname(nfield) = 'tauy'//TRIM(ADJUSTL(varname(nfield)))
        ENDDO

        u = u*v
        nfield = nfield+1; vars(nfield)%field => u; varname(nfield) = 'vu'
        ! I need v below
        nfield = nfield+1; vars(nfield)%field => v; varname(nfield) = 'vv'
        w = w*v
        nfield = nfield+1; vars(nfield)%field => w; varname(nfield) = 'vw'
        DO is = 1,inb_scal_array
           s(:,is) = s(:,is)*v
           nfield = nfield+1; vars(nfield)%field => s(:,is); WRITE(varname(nfield),*) is; varname(nfield) = 'v'//TRIM(ADJUSTL(varname(nfield)))
        ENDDO
        v = v*v ! I need v' above for the scalar fluxes

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgMom'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Hydrostatic pressure
        ! ###################################################################
     CASE ( 14 )
        is = 0

        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,1), txc(1,2),txc(1,3), txc(1,4), wrk1d,wrk2d,wrk3d)
        is = is+1; vars(is)%field => txc(:,1); varname(is) = 'P'

        q = C_0_R
        CALL FI_PRESSURE_BOUSSINESQ(q,s, txc(1,2), txc(1,3),txc(1,4), txc(1,5), wrk1d,wrk2d,wrk3d)
        is = is+1; vars(is)%field => txc(:,2); varname(is) = 'Phydro'

        txc(:,3) = txc(:,1) - txc(:,2)
        is = is+1; vars(is)%field => txc(:,3); varname(is) = 'Padvec'

        IF ( nfield .NE. is ) THEN ! Check
           CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
           CALL DNS_STOP(DNS_ERROR_WRKSIZE)
        ENDIF

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgPre'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Dissipation
        ! ###################################################################
     CASE ( 15 )
        CALL FI_DISSIPATION(i1,imax,jmax,kmax, u,v,w, txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)

        vars(1)%field => txc(:,1); varname(1) = 'Eps'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgEps'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Covariances among scalars
        ! ###################################################################
     CASE ( 16 )
        CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, s(1,1))
        CALL REYFLUCT2D(imax,jmax,kmax, g(1)%jac,g(3)%jac, area, s(1,2))

        txc(1:isize_field,1) = s(1:isize_field,1)   *s(1:isize_field,2)
        txc(1:isize_field,2) = txc(1:isize_field,1) *s(1:isize_field,1)
        txc(1:isize_field,3) = txc(1:isize_field,1) *s(1:isize_field,2)

        is = 0
        is = is+1; vars(is)%field => txc(:,1); varname(is) = 's1s2'
        is = is+1; vars(is)%field => txc(:,2); varname(is) = 's1s2s1'
        is = is+1; vars(is)%field => txc(:,3); varname(is) = 's1s2s2'

        IF ( nfield .NE. is ) THEN ! Check
           CALL IO_WRITE_ASCII(efile, 'AVERAGES. Array space nfield incorrect.')
           CALL DNS_STOP(DNS_ERROR_WRKSIZE)
        ENDIF

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='cov'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)

        ! ###################################################################
        ! Potential vorticity
        ! ###################################################################
     CASE ( 17 )
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

        vars(1)%field => txc(:,1); varname(1) = 'PV'
        vars(2)%field => txc(:,2); varname(2) = 'Cos'

        IF (  jmax_aux*opt_block .NE. g(2)%size ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, i1,i1,i1, imax,jmax_aux*opt_block,kmax, vars(is)%field, wrk1d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgPV'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, gate_level, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, vars, mean)


     END SELECT
  ENDDO

  CALL DNS_END(0)

  STOP

END PROGRAM AVERAGES

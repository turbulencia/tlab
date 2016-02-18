#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "AVERAGES"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/05/16 - J.P. Mellado
!#              General input files sequence from dns.ini or stdin.
!#              Absoft argument readingfrom command line removed (v>=4.4.0)
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!# 2013/10/01 - J.P. Mellado
!#              Finalizing the management of conditional statistics
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM AVERAGES

  USE DNS_TYPES
  USE DNS_CONSTANTS
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
  TREAL, DIMENSION(:),      ALLOCATABLE, SAVE         :: x,y,z, dx,dy,dz
  TYPE(grid_structure), DIMENSION(3) :: g

  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE, TARGET :: q, s, txc
  TREAL, DIMENSION(:),      ALLOCATABLE, SAVE         :: wrk1d,wrk2d,wrk3d

  TREAL, DIMENSION(:),      ALLOCATABLE, SAVE         :: mean, y_aux
  INTEGER(1), DIMENSION(:), ALLOCATABLE, SAVE         :: gate
  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE         :: surface

  TYPE(pointers_structure), DIMENSION(16) :: data

! Particle data
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q
  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_txc
  INTEGER(8), DIMENSION(:),   ALLOCATABLE, SAVE :: l_tags

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 fname, inifile, bakfile
  CHARACTER*32 varname(16)
  CHARACTER*64 str, line

  TINTEGER opt_main, opt_block, opt_order, opt_bins, opt_bcs, opt_format, flag_buoyancy
  TINTEGER opt_cond, opt_threshold
  TINTEGER nfield, isize_wrk3d, is, ij, n
  TREAL diff, umin, umax, dummy, s_aux(MAX_NSP)
  TREAL eloc1, eloc2, eloc3, cos1, cos2, cos3
  TINTEGER jmax_aux, iread_flow, iread_scal, ierr, idummy

! Gates for the definition of the intermittency function (partition of the fields)
  TINTEGER igate_size
  INTEGER(1) opt_gate, igate_vec(igate_size_max)
  TREAL gate_threshold(igate_size_max)

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TREAL opt_vec(iopt_size_max)
  TREAL opt_vec2(iopt_size_max)

  TINTEGER params_size
  TREAL params(params_size_max)

#ifdef USE_MPI
  INTEGER icount
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w

!########################################################################
!########################################################################
  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
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
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))
  ALLOCATE(gate(isize_field))

  ALLOCATE(y_aux(jmax_total)) ! Reduced vertical grid

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Additional options
! -------------------------------------------------------------------
  opt_main  =-1 ! default values
  opt_block = 1
  opt_gate  = 1
  opt_order = 1
  opt_bins  = 0
  opt_bcs   = 1
  opt_format= 0

  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamAverages', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'Option ?'
     WRITE(*,*) ' 1. Conventional averages'
     WRITE(*,*) ' 2. Intermittency or gate function'
     WRITE(*,*) ' 3. Zonal average of scalar gradient G_iG_i conditioned on scalar'
     WRITE(*,*) ' 4. Zonal avarages of main variables'
     WRITE(*,*) ' 5. Zonal avarages of enstrophy W_iW_i/2 equation'
     WRITE(*,*) ' 6. Zonal avarages of strain 2S_ijS_ij/2 equation'
     WRITE(*,*) ' 7. Zonal avarages of scalar gradient G_iG_i/2 equation'
     WRITE(*,*) ' 8. Zonal avarages of velocity gradient invariants'
     WRITE(*,*) ' 9. Zonal avarages of scalar gradient components'
     WRITE(*,*) '10. Zonal avarages of eigenvalues of rate-of-strain tensor'
     WRITE(*,*) '11. Zonal avarages of eigenframe of rate-of-strain tensor'
     WRITE(*,*) '12. Zonal avarages of longitudinal velocity derivatives'
     WRITE(*,*) '13. Zonal avarages of momentum vertical transport'
     WRITE(*,*) '14. Zonal avarages of pressure partition'
     WRITE(*,*) '15. Zonal avarages of dissipation'
     READ(*,*) opt_main

     WRITE(*,*) 'Planes block size ?'
     READ(*,*) opt_block  

     IF ( opt_main .GT. 2 ) THEN
        WRITE(*,*) 'Gate level to be used ?'
        READ(*,*) opt_gate
        WRITE(*,*) 'Number of moments ?'
        READ(*,*) opt_order
     ENDIF

     IF ( opt_main .EQ. 3 ) THEN
        WRITE(*,*) 'Number of bins ?'
        READ(*,*) opt_bins
        WRITE(*,*) 'Local interval (1) or homogeneous (0) ?'
        READ(*,*) opt_bcs
     ENDIF

     IF ( opt_main .EQ. 2 ) THEN
        WRITE(*,*) 'File Format ?'
        WRITE(*,*) ' 0. Just averages'
        WRITE(*,*) ' 1. Binaries'
        READ(*,*) opt_format
     ENDIF

#endif
  ELSE
     opt_main = DINT(opt_vec(1))
     IF ( iopt_size .GE. 2 ) opt_block = DINT(opt_vec(2))
     IF ( opt_main .GT. 2 ) THEN
        IF ( iopt_size .GE. 3 ) opt_gate  = DINT(opt_vec(3))
        IF ( iopt_size .GE. 4 ) opt_order = DINT(opt_vec(4))
     ENDIF
     IF ( opt_main .EQ. 3 ) THEN
        opt_bins = DINT(opt_vec(5))
        opt_bcs  = DINT(opt_vec(6))
     ENDIF

     CALL SCANINIINT(bakfile, inifile, 'PostProcessing', 'Format', '0', opt_format)

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
! Defining gate levels for conditioning
! -------------------------------------------------------------------
  opt_cond      = 0 ! default values
  igate_size    = 0
  opt_threshold = 0

  IF ( opt_main .GT. 1 ) THEN

#include "dns_read_partition.h"

  ENDIF

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
! in case jmax_total is not divisible by opt_block, drop the upper most planes
  jmax_aux = jmax_total/opt_block

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
  inb_txc    = 0
  nfield     = 2

  IF      ( opt_cond .EQ. 2 ) THEN
     iread_scal = 1
     inb_txc    = 1
  ELSE IF ( opt_cond .EQ. 3 ) THEN
     inb_txc = MAX(inb_txc,3)
     iread_flow = 1
  ELSE IF ( opt_cond .EQ. 4 ) THEN
     inb_txc = MAX(inb_txc,3)
     iread_scal = 1
  ENDIF

  IF      ( opt_main .EQ. 1 ) THEN
     inb_txc = MAX(inb_txc,9)
     iread_flow = icalc_flow
     iread_scal = icalc_scal
  ELSE IF ( opt_main .EQ. 2 ) THEN
     ifourier = 0
  ELSE IF ( opt_main .EQ. 3 ) THEN
     nfield = 2
     inb_txc = MAX(inb_txc,3)
     iread_scal = 1
  ELSE IF ( opt_main .EQ. 4 ) THEN
     nfield = 6+inb_scal
     inb_txc = MAX(inb_txc,3)
     iread_scal = 1
     iread_flow = 1
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) inb_txc = MAX(inb_txc,5)
  ELSE IF ( opt_main .EQ. 5 ) THEN ! enstrophy
     nfield = 7
     inb_txc = MAX(inb_txc,8)
     iread_flow = 1
     iread_scal = 1
  ELSE IF ( opt_main .EQ. 6 ) THEN
     nfield = 4
     inb_txc = MAX(inb_txc,8)
     iread_flow = 1
     iread_scal = 1
  ELSE IF ( opt_main .EQ. 7 ) THEN ! scalar gradient
     nfield = 5
     inb_txc = MAX(inb_txc,6)
     iread_scal = 1
     iread_flow = 1
  ELSE IF ( opt_main .EQ. 8 ) THEN
     nfield = 3
     inb_txc = MAX(inb_txc,6)
     iread_flow = 1
  ELSE IF ( opt_main .EQ. 9 ) THEN
     nfield = 5
     inb_txc = MAX(inb_txc,4)
     iread_scal = 1
  ELSE IF ( opt_main .EQ.10 ) THEN ! eigenvalues
     nfield = 3
     inb_txc = MAX(inb_txc,9)
     iread_flow = 1
  ELSE IF ( opt_main .EQ.11 ) THEN ! eigenframe
     nfield = 6
     inb_txc = MAX(inb_txc,10)
     iread_flow = 1
     iread_scal = 1
  ELSE IF ( opt_main .EQ.12 ) THEN ! longitudinal velocity derivatives
     nfield = 3
     inb_txc = MAX(inb_txc,3)
     iread_flow = 1
     iread_scal = 0
  ELSE IF ( opt_main .EQ.13 ) THEN ! Momentum vertical flux
     nfield = 8
     inb_txc = MAX(inb_txc,4)
     iread_flow = 1
     iread_scal = 1
  ELSE IF ( opt_main .EQ.14 ) THEN ! pressure partition
     nfield = 3
     inb_txc = MAX(inb_txc,6)
     iread_flow = 1
     iread_scal = 1
  ELSE IF ( opt_main .EQ.15 ) THEN ! dissipation partition
     nfield = 1
     inb_txc = MAX(inb_txc,6)
     iread_flow = 1
     iread_scal = 0
  ENDIF

  IF ( opt_main .EQ. 1 ) THEN
     ALLOCATE(mean(jmax*MAX_AVG_TEMPORAL))
  ELSE
     IF ( opt_main .EQ. 15 ) THEN 
        ALLOCATE(mean(MAX(jmax*5,(MAX(opt_bins,2)*opt_order*nfield))))
     ELSE 
        ALLOCATE(mean(MAX(opt_bins,2)*opt_order*nfield))
     ENDIF
  ENDIF

  IF ( opt_main .EQ. 2 ) THEN
     IF ( opt_cond .EQ. 1 ) THEN; ALLOCATE(surface(imax*kmax,2*igate_size_max))
     ELSE;                        ALLOCATE(surface(imax*kmax,2*igate_size))
     ENDIF
  ENDIF

  isize_txc   = isize_txc_field*inb_txc
  isize_wrk3d = MAX(isize_field,opt_bins*opt_order*nfield*jmax)
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)

#include "dns_alloc_arrays.h"

  IF ( icalc_particle .EQ. 1 ) THEN
#include "dns_alloc_larrays.h"
  ENDIF
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
     y_aux(is) = y_aux(is) + y(ij)/M_REAL(opt_block)
  enddo

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( iread_flow .EQ. 1 ) THEN ! We need array space
     CALL DNS_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Define pointers
! ###################################################################
  IF ( iread_flow .EQ. 1 ) THEN
     u => q(:,1)
     v => q(:,2)
     w => q(:,3)
  ENDIF

! ###################################################################
! Calculating statistics
! ###################################################################
  DO it = 1,itime_size
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
#include "dns_calc_partition.h"
     IF ( opt_main .EQ. 2 .AND. opt_cond .EQ. 1 ) rtime = params(1)

! ###################################################################
! Conventional statistics
! ###################################################################
     IF      ( opt_main .EQ. 1 ) THEN
        IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN 
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_AIRWATER_PHAL(i1,i1,i1, mean_i(2), p_init, mean_i(1))
              CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)
              CALL THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(1,2), p_init, s(1,1))
           ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
              CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(:,inb_scal_array), wrk3d)
           ENDIF
           CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, u,v,w, s, txc(1,3), &
                txc(1,1),txc(1,2),txc(1,4), wrk1d,wrk2d,wrk3d)

        ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           txc(:,1) = C_1_R ! to be developed
           
        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, q(1,4), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, q(1,5), wrk3d)! density
           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, q(1,4), q(1,5), q(1,7), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, q(1,5), q(1,7), q(1,6))
           IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, q(1,7), q(1,8))

        ENDIF

        IF ( icalc_flow .EQ. 1 ) THEN
           CALL AVG_FLOW_TEMPORAL_LAYER(y,dx,dy,dz, q, s, &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
                txc(1,7),txc(1,8),txc(1,9), mean, wrk1d,wrk2d,wrk3d)
        ENDIF
        
        IF ( icalc_scal .EQ. 1 ) THEN
           IF ( imixture .EQ. MIXT_TYPE_BILAIRWATER .OR. imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT ) THEN
              CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,1),s(:,inb_scal_array)) ! Update the liquid function
           ENDIF
           DO is = inb_scal+1,inb_scal_array ! Add diagnostic fields, if any
              mean_i(is) = C_1_R; delta_i(is) = C_0_R; ycoor_i(is) = ycoor_i(1); schmidt(is) = schmidt(1)
           ENDDO
           DO is = 1,inb_scal_array          ! All, prognostic and diagnostic fields in array s
              CALL AVG_SCAL_TEMPORAL_LAYER(is, y, dx,dy,dz, q,s, s(1,is), &
                   txc(1,1),txc(1,2),txc(1,3),txc(1,4), mean, wrk1d,wrk2d,wrk3d)
           ENDDO

! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
           IF ( flag_buoyancy .EQ. 1 ) THEN
              wrk1d(1:jmax) = C_0_R 
              CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, txc(:,1), wrk1d)
              dummy = C_1_R/froude
              txc(1:isize_field,1) = txc(1:isize_field,1)*dummy
! mean values
              dummy = C_0_R
              s_aux(1:inb_scal) = mean_i(1:inb_scal) - C_05_R*delta_i(1:inb_scal)
              CALL FI_BUOYANCY(ibodyforce, i1,i1,i1, body_param, s_aux, umin, dummy)
              s_aux(1:inb_scal) = mean_i(1:inb_scal) + C_05_R*delta_i(1:inb_scal)
              CALL FI_BUOYANCY(ibodyforce, i1,i1,i1, body_param, s_aux, umax, dummy)
              mean_i(is) = (umax+umin)/froude; delta_i(is) = ABS(umax-umin)/froude; ycoor_i(is) = ycoor_i(1); schmidt(is) = schmidt(1)
              CALL AVG_SCAL_TEMPORAL_LAYER(is, y, dx,dy,dz, q,s, txc(1,1), &
                   txc(1,2),txc(1,3),txc(1,4),txc(1,5), mean, wrk1d,wrk2d,wrk3d)
              
           ENDIF

        ENDIF
        
! Lagrange Liquid and Liquid without diffusion
        IF ( icalc_particle .EQ. 1 ) THEN
           IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2 & 
                .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
              WRITE(fname,*) itime; fname = "particle."//TRIM(ADJUSTL(fname))
              CALL DNS_READ_PARTICLE(fname,l_q) 
                 
              l_txc = C_1_R; ! We want density
              CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,5))
              
              txc(:,5) = txc(:,5) + 0.00000001
              DO is = inb_scal_array+2,inb_scal_particle+inb_scal_array+1
                 l_txc(:,1)=l_q(:,3+is-inb_scal_array-1) !!! DO WE WANT l_txc(:,is) ???
                 CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,6))   
                 txc(:,6) = txc(:,6)/txc(:,5)
                 CALL AVG_SCAL_TEMPORAL_LAYER(is, y,dx,dy,dz, q,s, txc(1,6), &
                      txc(1,1),txc(1,2),txc(1,3),txc(1,4), mean, wrk1d,wrk2d,wrk3d)
              ENDDO
           ENDIF
        ENDIF
! DO from is inb_scal_array+1 to inb_paticle+inb_scal_array
! PARTICLE TO FIELD to txc(5)
! CALL AVG_SCAL_TEMPORAL_LAYER on new field (substutute by s(1,is)
!ENDDO 

! ###################################################################
! Partition of field
! ###################################################################
     ELSE IF ( opt_main .EQ. 2 ) THEN
        DO n = 1,igate_size
           WRITE(varname(n),*) n; varname(n) = 'Partition'//TRIM(ADJUSTL(varname(n)))
        ENDDO
        WRITE(fname,*) itime; fname='int'//TRIM(ADJUSTL(fname))
        CALL INTER2D_N(fname, varname, rtime, imax,jmax,kmax, igate_size, y, gate, igate_vec)
        
        IF ( opt_cond .GT. 1 ) THEN ! write only if the gate information has not been read
           WRITE(fname,*) itime; fname = 'par'//TRIM(ADJUSTL(fname))
           params(1) = rtime; params(2) = M_REAL(igate_size); params_size = 2
           CALL IO_WRITE_INT1(fname, i1, imax,jmax,kmax,itime, params_size,params, gate)
        ENDIF

! Envelopes
        DO n = 1,igate_size
           WRITE(str,*) n; str = 'Envelope'//TRIM(ADJUSTL(str))
           line = 'Calculating '//TRIM(ADJUSTL(str)); CALL IO_WRITE_ASCII(lfile,line)
           CALL BOUNDARY_LOWER_INT1(imax,jmax,kmax, igate_vec(n), y, gate, wrk3d, surface(:,n           ), wrk2d)
           CALL BOUNDARY_UPPER_INT1(imax,jmax,kmax, igate_vec(n), y, gate, wrk3d, surface(:,n+igate_size), wrk2d)
        ENDDO

        WRITE(fname,*) itime; fname = 'lower.'//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i1, imax,i1,kmax, igate_size, isize_wrk3d, surface,                 wrk3d)
        WRITE(fname,*) itime; fname = 'upper.'//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i1, imax,i1,kmax, igate_size, isize_wrk3d, surface(:,1+igate_size), wrk3d)

! ###################################################################
! Scalar gradient conditioned on Z
! ###################################################################
     ELSE IF ( opt_main .EQ. 3 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
        CALL FI_GRADIENT(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, s, txc(1,1), txc(1,2), txc(1,3), wrk1d, wrk2d, wrk3d)
        DO ij = 1,isize_field
           txc(ij,2) = log(txc(ij,1))
        ENDDO

        data(1)%field => txc(:,1); varname(1) = 'GradientG_iG_i'
        data(2)%field => txc(:,2); varname(2) = 'LnGradientG_iG_i'

        IF ( opt_bcs .EQ. 0 ) THEN
           CALL MINMAX(imax,jmax,kmax, s, umin,umax)
        ENDIF

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='cavgZ'//TRIM(ADJUSTL(fname))
        CALL CAVG2D_N(fname, opt_bcs, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_bins, opt_order, umin, umax, y_aux, gate, s, data, mean, wrk1d)

! ###################################################################
! Main variables
! ###################################################################
     ELSE IF ( opt_main .EQ. 4 ) THEN
        nfield = 0
        nfield = nfield+1; data(nfield)%field => u(:); varname(nfield) = 'U'
        nfield = nfield+1; data(nfield)%field => v(:); varname(nfield) = 'V'
        nfield = nfield+1; data(nfield)%field => w(:); varname(nfield) = 'W'

        IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN 
           CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, u,v,w, s, txc(1,1), &
                txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d)
           nfield = nfield+1; data(nfield)%field => txc(:,1); varname(nfield) = 'P'

        ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           txc(:,1) = C_1_R ! to be developed
           nfield = nfield+1; data(nfield)%field => txc(:,1); varname(nfield) = 'P'
           
        ELSE
           WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
           CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density
           CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
           CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, txc(1,2), txc(1,3), txc(1,1))
           nfield = nfield+1; data(nfield)%field => txc(:,2); varname(nfield) = 'R'
           nfield = nfield+1; data(nfield)%field => txc(:,1); varname(nfield) = 'P'
           nfield = nfield+1; data(nfield)%field => txc(:,3); varname(nfield) = 'T'

        ENDIF
     
        IF ( icalc_scal .EQ. 1 ) THEN
           nfield = nfield+1; data(nfield)%field => s(:,1);  varname(nfield) = 'Scalar1'
        ENDIF

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='cavg'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Enstrophy equation
! ###################################################################
     ELSE IF ( opt_main .EQ. 5 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing baroclinic term...')
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR.&
             imode_eqns .EQ. DNS_EQNS_ANELASTIC     )THEN
           IF ( ibodyforce .EQ. EQNS_NONE ) THEN
              txc(:,4) = C_0_R; txc(:,5) = C_0_R; txc(:,6) = C_0_R
           ELSE
! calculate buoyancy vector along Oy
           wrk1d(1:jmax) = C_0_R 
           CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
           DO ij = 1,isize_field
              s(ij,1) = wrk3d(ij)*body_vector(2)
           ENDDO

           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, s, txc(1,4), i0,i0, wrk1d,wrk2d,wrk3d)
           txc(:,4) =-txc(:,4)
           txc(:,5) = C_0_R
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, s, txc(1,6), i0,i0, wrk1d,wrk2d,wrk3d)
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
           CALL FI__BAROCLINIC(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
                dx,dy,dz, txc(1,2),txc(1,1), txc(1,4), txc(1,3),txc(1,7), wrk1d,wrk2d,wrk3d)
        ENDIF
! result vector in txc1, txc2, txc3
        CALL FI_CURL(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,1),txc(1,2),txc(1,3), txc(1,7), wrk1d,wrk2d,wrk3d)
! scalar product, store in txc8
        DO ij = 1,isize_field
           txc(ij,8) = txc(ij,1)*txc(ij,4) + txc(ij,2)*txc(ij,5) + txc(ij,3)*txc(ij,6)
        ENDDO
        
        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy production...')
        CALL FI_VORTICITY_PRODUCTION(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), &
             wrk1d,wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy diffusion...')
        CALL FI_VORTICITY_DIFFUSION&
             (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), &
             wrk1d,wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,2) = visc*txc(ij,2)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing enstrophy...')
        CALL FI_VORTICITY(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing dilatation term...')
        CALL FI_INVARIANT_P(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)

        DO ij = 1,isize_field
           txc(ij,6) = txc(ij,4)*txc(ij,3) ! -w^2 div(u)
           txc(ij,5) = txc(ij,1)/txc(ij,3) ! production rate 
           txc(ij,4) = log(txc(ij,3))      ! ln(w^2)
        ENDDO

        data(1)%field => txc(:,3); varname(1) = 'EnstrophyW_iW_i'
        data(2)%field => txc(:,4); varname(2) = 'LnEnstrophyW_iW_i'
        data(3)%field => txc(:,1); varname(3) = 'ProductionW_iW_jS_ij'
        data(4)%field => txc(:,2); varname(4) = 'DiffusionNuW_iLapW_i'
        data(5)%field => txc(:,6); varname(5) = 'DilatationMsW_iW_iDivU'
        data(6)%field => txc(:,8); varname(6) = 'Baroclinic'
        data(7)%field => txc(:,5); varname(7) = 'RateAN_iN_jS_ij'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgW2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Strain equation
! ###################################################################
     ELSE IF ( opt_main .EQ. 6 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname)) !need to read again thermo data
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i4, isize_wrk3d, txc(1,1), wrk3d)! energy
        CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i5,i5, isize_wrk3d, txc(1,2), wrk3d)! density

        CALL IO_WRITE_ASCII(lfile,'Computing strain pressure...')
        CALL THERMO_CALORIC_TEMPERATURE&
             (imax, jmax, kmax, s, txc(1,1), txc(1,2), txc(1,3), wrk3d)
        CALL THERMO_THERMAL_PRESSURE&
             (imax, jmax, kmax, s, txc(1,2), txc(1,3), txc(1,1)) ! pressure in txc1
        CALL FI_STRAIN_PRESSURE&
             (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx, dy, dz, u, v, w, txc(1,1), &
             txc(1,2), txc(1,3), txc(1,4), txc(1,5), txc(1,6), wrk1d,wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,1)=C_2_R*txc(ij,2)
        ENDDO
        
        CALL IO_WRITE_ASCII(lfile,'Computing strain production...')
        CALL FI_STRAIN_PRODUCTION(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk1d,wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,2)=C_2_R*txc(ij,2)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing strain diffusion...')
        CALL FI_STRAIN_DIFFUSION&
             (iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, &
             txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk1d,wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,3)=C_2_R*visc*txc(ij,3)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing strain...')
        CALL FI_STRAIN(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
        DO ij = 1,isize_field
           txc(ij,4)=C_2_R*txc(ij,4)
        ENDDO
              
        data(1)%field => txc(:,4); varname(1) = 'Strain2S_ijS_i'
        data(2)%field => txc(:,2); varname(2) = 'ProductionMs2S_ijS_jkS_ki'
        data(3)%field => txc(:,3); varname(3) = 'DiffusionNuS_ijLapS_ij'
        data(4)%field => txc(:,1); varname(4) = 'Pressure2S_ijP_ij'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgS2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Scalar gradient equation
! ###################################################################
     ELSE IF ( opt_main .EQ. 7 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient production...')
        CALL FI_GRADIENT_PRODUCTION(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, s, u, v, w, txc(1,1), txc(1,2), txc(1,3), txc(1,4), txc(1,5), txc(1,6),&
             wrk1d, wrk2d, wrk3d)

! array u used as auxiliar
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient diffusion...')
        CALL FI_GRADIENT_DIFFUSION&
             (iunifx, iunify, iunifz, imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, s, txc(1,2), txc(1,3), txc(1,4), txc(1,5), txc(1,6), u, &
             wrk1d, wrk2d, wrk3d)
        DO ij = 1,isize_field
           txc(ij,2) = diff*txc(ij,2)
        ENDDO

        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient...')
        CALL FI_GRADIENT(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, s, txc(1,3), txc(1,4), txc(1,5), wrk1d, wrk2d, wrk3d)
        DO ij = 1,isize_field
           txc(ij,5) = txc(ij,1)/txc(ij,3)
           txc(ij,4) = log(txc(ij,3))
        ENDDO

        data(1)%field => txc(:,3); varname(1) = 'GradientG_iG_i'
        data(2)%field => txc(:,4); varname(2) = 'LnGradientG_iG_i'
        data(3)%field => txc(:,1); varname(3) = 'ProductionMsG_iG_jS_ij'
        data(4)%field => txc(:,2); varname(4) = 'DiffusionNuG_iLapG_i'
        data(5)%field => txc(:,5); varname(5) = 'StrainAMsN_iN_jS_ij'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgG2'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Velocity gradient invariants
! ###################################################################
     ELSE IF ( opt_main .EQ. 8 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing third invariant R...')
        CALL FI_INVARIANT_R(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, u, v, w, txc(1,1), &
             txc(1,2), txc(1,3), txc(1,4), txc(1,5), txc(1,6), wrk1d, wrk2d, wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing second invariant Q...')
        CALL FI_INVARIANT_Q(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, u, v, w, txc(1,2), &
             txc(1,3), txc(1,4), txc(1,5), wrk1d, wrk2d, wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing first invariant P...')
        CALL FI_INVARIANT_P(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, u, v, w, txc(1,3), &
             txc(1,4), wrk1d, wrk2d, wrk3d)

        data(1)%field => txc(:,3); varname(1) = 'InvariantP'
        data(2)%field => txc(:,2); varname(2) = 'InvariantQ'
        data(3)%field => txc(:,1); varname(3) = 'InvariantR'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgInv'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Scalar gradient components
! ###################################################################
     ELSE IF ( opt_main .EQ. 9 ) THEN
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, s, txc(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s, txc(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, s, txc(1,3), i0,i0, wrk1d,wrk2d,wrk3d)
! Angles; s array is overwritten to save space
        DO ij = 1,isize_field
           dummy = txc(ij,2)/sqrt(txc(ij,1)*txc(ij,1)+txc(ij,2)*txc(ij,2)+txc(ij,3)*txc(ij,3))
           txc(ij,4) = asin(dummy)                 ! with Oy
           s(ij,1)  = atan2(txc(ij,3),txc(ij,1))  ! with Ox in plane xOz
        ENDDO

        data(1)%field => txc(:,1); varname(1) = 'GradientX'
        data(2)%field => txc(:,2); varname(2) = 'GradientY'
        data(3)%field => txc(:,3); varname(3) = 'GradientZ'
        data(4)%field => s(:,1);   varname(4) = 'Theta'
        data(5)%field => txc(:,4); varname(5) = 'Phi'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgGi'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! eigenvalues of rate-of-strain tensor
! ###################################################################
     ELSE IF ( opt_main .EQ. 10 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
             dx, dy, dz, u, v, w, txc(1,1), wrk1d, wrk2d, wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...') ! txc6-txc9
        CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))

        data(1)%field => txc(:,7); varname(1) = 'Lambda1'
        data(2)%field => txc(:,8); varname(2) = 'Lambda2'
        data(3)%field => txc(:,9); varname(3) = 'Lambda3'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgEig'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! eigenframe of rate-of-strain tensor
! ###################################################################
     ELSE IF ( opt_main .EQ. 11 ) THEN
        CALL IO_WRITE_ASCII(lfile,'Computing rate-of-strain tensor...') ! txc1-txc6
        CALL FI_STRAIN_TENSOR(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,1), wrk1d,wrk2d,wrk3d)

        CALL IO_WRITE_ASCII(lfile,'Computing eigenvalues...')           ! txc7-txc9
        CALL FI_TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1,1), txc(1,7))

        CALL IO_WRITE_ASCII(lfile,'Computing eigenframe...')            ! txc1-txc6
        CALL FI_TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1,1), txc(1,7))

! local direction cosines of vorticity vector
        CALL IO_WRITE_ASCII(lfile,'Computing vorticity vector...')      ! txc7-txc9
        CALL FI_CURL(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
             dx,dy,dz, u,v,w, txc(1,7),txc(1,8),txc(1,9), txc(1,10), wrk1d,wrk2d,wrk3d)

        DO ij = 1,isize_field
           dummy = sqrt(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           u(ij) = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           v(ij) = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           w(ij) = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
        ENDDO

        data(1)%field => u; varname(1) = 'cos(w,lambda1)'
        data(2)%field => v; varname(2) = 'cos(w,lambda2)'
        data(3)%field => w; varname(3) = 'cos(w,lambda3)'

! local direction cosines of scalar gradient vector
        CALL IO_WRITE_ASCII(lfile,'Computing scalar gradient vector...') ! txc7-txc9
        CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
             dx, s, txc(1,7), i0, i0, wrk1d, wrk2d, wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
             dy, s, txc(1,8), i0, i0, wrk1d, wrk2d, wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
             dz, s, txc(1,9), i0, i0, wrk1d, wrk2d, wrk3d)

        DO ij = 1,isize_field
           dummy = sqrt(txc(ij,7)*txc(ij,7)+txc(ij,8)*txc(ij,8)+txc(ij,9)*txc(ij,9))
           cos1  = (txc(ij,7)*txc(ij,1) + txc(ij,8)*txc(ij,2) + txc(ij,9)*txc(ij,3))/dummy
           cos2  = (txc(ij,7)*txc(ij,4) + txc(ij,8)*txc(ij,5) + txc(ij,9)*txc(ij,6))/dummy
           eloc1 = txc(ij,2)*txc(ij,6)-txc(ij,5)*txc(ij,3)
           eloc2 = txc(ij,3)*txc(ij,4)-txc(ij,6)*txc(ij,1)
           eloc3 = txc(ij,1)*txc(ij,5)-txc(ij,4)*txc(ij,2)
           cos3  = (txc(ij,7)*eloc1 + txc(ij,8)*eloc2 + txc(ij,9)*eloc3)/dummy
           txc(ij,7) = cos1; txc(ij,8) = cos2; txc(ij,9) = cos3 
        ENDDO

        data(4)%field => txc(:,7); varname(4) = 'cos(G,lambda1)'
        data(5)%field => txc(:,8); varname(5) = 'cos(G,lambda2)'
        data(6)%field => txc(:,9); varname(6) = 'cos(G,lambda3)'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgCos'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! longitudinal velocity derivatives
! ###################################################################
     ELSE IF ( opt_main .EQ. 12 ) THEN
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, txc(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, txc(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, txc(1,3), i0,i0, wrk1d,wrk2d,wrk3d)

        data(1)%field => txc(:,1); varname(1) = 'dudx'
        data(2)%field => txc(:,2); varname(2) = 'dvdy'
        data(3)%field => txc(:,3); varname(3) = 'dwdz'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgDer'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Momentum vertical transport
! ###################################################################
     ELSE IF ( opt_main .EQ. 13 ) THEN
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, txc(:,1), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, txc(:,2), i0,i0, wrk1d,wrk2d,wrk3d)
        txc(:,1) = ( txc(:,1) + txc(:,2) ) *visc

        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, txc(:,2), i0,i0, wrk1d,wrk2d,wrk3d)
        txc(:,2) =   txc(:,2) *C_2_R       *visc

        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, txc(:,3), i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, txc(:,4), i0,i0, wrk1d,wrk2d,wrk3d)
        txc(:,3) = ( txc(:,3) + txc(:,4) ) *visc

        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s(:,1), txc(:,4), i0,i0, wrk1d,wrk2d,wrk3d)
        txc(:,4) =   txc(:,4) *diff

        is = 0
        is = is+1; data(is)%field => txc(:,1); varname(is) = 'tauyx'
        is = is+1; data(is)%field => txc(:,2); varname(is) = 'tauyy'
        is = is+1; data(is)%field => txc(:,3); varname(is) = 'tauyz'
        is = is+1; data(is)%field => txc(:,4); varname(is) = 'diffz'

        CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, u)
        CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, v)
        CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, w)
        CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, s(:,1))
        u = u*v
        w = w*v
        s(:,1) = s(:,1)*v
        v = v*v

        is = is+1; data(is)%field => u; varname(is) = 'vu'
        is = is+1; data(is)%field => v; varname(is) = 'vv'
        is = is+1; data(is)%field => w; varname(is) = 'vw'
        is = is+1; data(is)%field => s(:,1); varname(is) = 'vs'

        IF ( nfield .NE. is ) THEN ! Check
           CALL IO_WRITE_ASCII(efile, 'AVERAGES. Array space nfield incorrect.')
           CALL DNS_STOP(DNS_ERROR_WRKSIZE)
        ENDIF

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgMom'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Hydrostatic pressure
! ###################################################################
     ELSE IF ( opt_main .EQ. 14 ) THEN
        is = 0

        CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, u,v,w, s, txc(1,1), &
             txc(1,2),txc(1,3),txc(1,4), wrk1d,wrk2d,wrk3d)
        is = is+1; data(is)%field => txc(:,1); varname(is) = 'P'
        
        txc(:,3) = C_0_R
        CALL FI_PRESSURE_BOUSSINESQ(y,dx,dy,dz, txc(1,3),txc(1,3),txc(1,3), s, txc(1,2), &
             txc(1,4),txc(1,5),txc(1,6), wrk1d,wrk2d,wrk3d)
        is = is+1; data(is)%field => txc(:,2); varname(is) = 'Phydro'

        txc(:,3) = txc(:,1) - txc(:,2)
        is = is+1; data(is)%field => txc(:,3); varname(is) = 'Padvec'

        IF ( nfield .NE. is ) THEN ! Check
           CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
           CALL DNS_STOP(DNS_ERROR_WRKSIZE)
        ENDIF

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgPre'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)

! ###################################################################
! Dissipation
! ###################################################################
     ELSE IF ( opt_main .EQ. 15 ) THEN 

        txc(1:isize_field,5) = C_1_R

        CALL FI_DISSIPATION(i1,imode_fdm,imax,jmax,kmax,i1bc,j1bc,k1bc, & 
             area,visc,dx,dy,dz,txc(:,5),u,v,w,txc(:,1), & 
             txc(:,2),txc(:,3),txc(:,4),mean,wrk1d,wrk2d,wrk3d)

        data(1)%field => txc(:,1); varname(1) = 'Eps'

        IF (  jmax_aux*opt_block .NE. jmax_total ) THEN
           DO is = 1,nfield
              CALL REDUCE_BLOCK_INPLACE(imax,jmax,kmax, imax,jmax_aux*opt_block,kmax, data(is)%field, wrk3d)
           ENDDO
        ENDIF

        WRITE(fname,*) itime; fname='avgEps'//TRIM(ADJUSTL(fname))
        CALL AVG2D_N(fname, varname, opt_gate, rtime, imax*opt_block, jmax_aux, kmax, &
             nfield, opt_order, y_aux, gate, data, mean)
     ENDIF
  ENDDO

  CALL DNS_END(0)

  STOP

100 FORMAT(G_FORMAT_R)

END PROGRAM AVERAGES

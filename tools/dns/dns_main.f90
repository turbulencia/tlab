#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "DNS"

PROGRAM DNS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE LAGRANGE_GLOBAL
  USE DNS_LOCAL 
  USE DNS_TOWER
  USE BOUNDARY_INFLOW
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS
  USE PARTICLE_TRAJECTORIES
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
! Grid and associated arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z

! Flow/Scalar variables and RHS space
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q,s, h_q,h_s, txc

! Particle data
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_hq, l_txc
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: l_comm

! Work arrays
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk1d,wrk2d,wrk3d

! Auxiliar memory space
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: vaux

! Inflow arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x_inf, y_inf, z_inf
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q_inf, s_inf

  TARGET q

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: e, rho, p, T, vis
  
  CHARACTER*32 fname, inifile
  CHARACTER*128 str, line
  TINTEGER idummy, ig
  TINTEGER ierr, isize_wrk3d, isize_vaux, isize_loc
  TREAL dummy

! ###################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL(inifile)
  ENDIF
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(inifile)
#endif
#ifdef LES
  CALL LES_READ_INI(inifile)
#endif
  CALL DNS_READ_LOCAL(inifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
  IF ( imode_rhs .EQ. EQNS_RHS_NONBLOCKING ) CALL DNS_NB3DFFT_INITIALIZE
#endif

  itime        = nitera_first
  logs_data(1) = 0 ! Status

! #######################################################################
! Memory management
! #######################################################################
  isize_loc = MAX(g_inf(1)%size,MAX(g_inf(2)%size,g_inf(3)%size))
  isize_wrk1d = MAX(isize_wrk1d,isize_loc)

  isize_loc = MAX(g_inf(1)%size*g_inf(2)%size,MAX(g_inf(1)%size*g_inf(3)%size,g_inf(2)%size*g_inf(3)%size))
  isize_wrk2d = MAX(isize_wrk2d, isize_loc)

! txc
  inb_txc = 9
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN 
                                             inb_txc = 6
     IF ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) inb_txc = inb_txc+1
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL       .AND. &
            iadvection .EQ. EQNS_SKEWSYMMETRIC      .AND. &
            iviscous   .EQ. EQNS_EXPLICIT                 ) THEN
                                             inb_txc = 6
  ENDIF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .GT. C_0_R ) inb_txc = inb_txc + 1

  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN ! because of the statistics
     inb_txc = MAX(inb_txc,7)

     IF ( fstavg .EQ. 1 ) THEN
        idummy =  nstatavg*jmax*MAX_AVG_SPATIAL / isize_txc_field
        IF ( MOD( nstatavg*jmax*MAX_AVG_SPATIAL , isize_txc_field ) .GT. 0 ) idummy = idummy+1
        inb_txc = MAX(inb_txc,idummy)
     ENDIF

  ENDIF

#ifdef USE_PSFFT
  IF ( imode_rhs .EQ. EQNS_RHS_NONBLOCKING ) inb_txc = MAX(inb_txc,15)
#endif 

#ifdef LES
  IF ( iles .EQ. 1 ) THEN ! this number needs to be revised
     isize_loc = 13

     IF ( iles_type_regu .EQ. LES_REGU_SMGDYN .OR. iles_type_regu .EQ. LES_REGU_SMGDYNRMS ) THEN
        IF ( iles_type_tran .EQ. LES_TRAN_NONE ) isize_loc = 10
        isize_loc = isize_loc + 4 ! space for aux_sg in dynamic smagorinsky
     ENDIF
     IF ( iles_type_chem .EQ. LES_CHEM_QUASIBS ) isize_loc = isize_loc + 4 ! space for chi

     idummy =  isize_wrk1d*isize_wrk1d*3 / isize_txc_field
     IF ( MOD( isize_wrk1d*isize_wrk1d*3 , isize_txc_field ) .GT. 0 ) idummy = idummy+1
     isize_loc = MAX(isize_loc,idummy)

     inb_txc = MAX(isize_loc, inb_txc)

  ENDIF
#endif

  isize_txc = inb_txc*isize_txc_field

! wkr3d
  isize_wrk3d = MAX(imax,g_inf(1)%size)*MAX(jmax,g_inf(2)%size)*kmax
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)
  IF ( icalc_part .eq. 1) THEN
     isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
     isize_wrk3d = MAX(isize_wrk3d,(jmax*(kmax+1)*inb_particle_interp*2))
     isize_wrk3d = MAX(isize_wrk3d,(jmax*(imax+1)*inb_particle_interp*2))
  END IF
  IF ( tower_mode .EQ. 1 ) THEN 
     isize_wrk3d = MAX(isize_wrk3d,nitera_save*(g(2)%size+2))
  ENDIF

#ifdef LES
#ifdef USE_MPI
  IF ( ims_npro .GT. 1 ) THEN ! wrk3d for OZ filter in PARALLEL mode may require more space in LES
     isize_wrk3d = MAX(isize_wrk3d,imax*jmax*(isgs_f0size + isgs_f1size + kmax))
  ENDIF
#endif
#endif

! -------------------------------------------------------------------
! Allocating basic memory space
! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

#include "dns_alloc_arrays.h"

  ALLOCATE(x_inf(g_inf(1)%size,g_inf(1)%inb_grid)) ! Inflow fields for spatial simulations
  ALLOCATE(y_inf(g_inf(2)%size,g_inf(2)%inb_grid))
  ALLOCATE(z_inf(g_inf(3)%size,g_inf(3)%inb_grid))
  ALLOCATE(q_inf(g_inf(1)%size*g_inf(2)%size*kmax,inb_flow_array))
  ALLOCATE(s_inf(g_inf(1)%size*g_inf(2)%size*kmax,inb_scal_array))

! -------------------------------------------------------------------
! Allocating additional memory space
! -------------------------------------------------------------------
! Rhs
  WRITE(str,*) inb_flow; line = 'Allocating array rhs flow of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_q(isize_field,inb_flow),    stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_q.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) inb_scal; line = 'Allocating array rhs scal of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_s(isize_field,inb_scal),    stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_s.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

! Array vaux
  CALL DNS_VAUX(isize_vaux)

  WRITE(str,*) isize_vaux; line = 'Allocating array vaux of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(vaux(isize_vaux),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for vaux.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

! Lagrangian part
  IF ( icalc_part .EQ. 1 ) THEN
#include "dns_alloc_larrays.h"
     
     WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(l_comm(isize_l_comm), stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
     
     WRITE(str,*) isize_particle; line = 'Allocating array l_hq of size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) inb_part; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(l_hq(isize_particle,inb_part),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_hq.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  
  ENDIF

! ###################################################################
! Subarray for writing planes
! ###################################################################
#ifdef USE_MPI
  CALL DNS_MPIO_AUX
#else
  io_aux(:)%offset = 0
#endif

! #######################################################################
! Initializing tower stuff 
! #######################################################################
  IF ( tower_mode .EQ. 1 ) THEN 
     CALL DNS_TOWER_INITIALIZE(tower_stride)  
  ENDIF

! #######################################################################
! Log files
! #######################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN; CALL AVG_DEFS_TEMPORAL
     ELSE;                                         CALL AVG_DEFS_SPATIAL; ENDIF
#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Initialize arrays to zero
! ###################################################################
  q(:,:) = C_0_R; h_q(:,:) = C_0_R ! Basic fields and rhs
  s(:,:) = C_0_R; h_s(:,:) = C_0_R

  q_inf(:,:) = C_0_R               ! Inflow fields for spatial simulations
  s_inf(:,:) = C_0_R

  vaux(:) = C_0_R                  ! Auxiliary information

  IF ( icalc_part .EQ. 1 ) THEN ! Lagrangian
     l_q = C_0_R; l_hq = C_0_R
  ENDIF

! ###################################################################
! Read the grid 
! ###################################################################
#include "dns_read_grid.h"

  IF ( g_inf(1)%size .GT. 1 ) THEN ! Inflow fields for spatial simulations
     CALL IO_READ_GRID('grid.inf', g_inf(1)%size, g_inf(2)%size, g_inf(3)%size,  &
          g_inf(1)%scale,g_inf(2)%scale,g_inf(3)%scale, x_inf,y_inf,z_inf)
     CALL FDM_INITIALIZE(x_inf, g_inf(1), wrk1d)
     g_inf(2)%nodes => y_inf(:,1)
     g_inf(3)%nodes => z_inf(:,1)
  ENDIF
  
! ###################################################################
! Initialize filters
! ###################################################################  
  DO ig = 1,3
     CALL OPR_FILTER_INITIALIZE( g(ig), FilterDomain(ig), wrk1d )
  END DO

! ###################################################################
! Initialize spectral Poisson solver
! ###################################################################
  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Check
! ###################################################################
  CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

! ###################################################################
! Read fields
! ###################################################################
! Flow fields
  iviscchg  = 0 ! Default is no continuous change in viscosity
  viscstop  = visc
  WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
  CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)
  viscstart = visc

  IF ( viscstart .NE. viscstop ) THEN ! check if viscosity has been changed
     WRITE(str,*) viscstart
     str = 'Changing original viscosity '//TRIM(ADJUSTL(str))//' to new value.'
     CALL IO_WRITE_ASCII(lfile,str)
     IF ( visctime .GT. C_0_R ) THEN ! Continuous change
        iviscchg = 1; visctime = (viscstart - viscstop)/visctime
     ELSE
        visc = viscstop
     ENDIF
  ENDIF

! Scalar fields
  IF ( icalc_scal .EQ. 1 ) THEN
     dummy = rtime ! flow data controls global data
     WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
     rtime = dummy

     IF ( itime .NE. nitera_first ) THEN ! files may contain different values
        CALL IO_WRITE_ASCII(wfile, 'DNS_MAIN: Scalar ItNumber size mismatch; enforced.')
        itime = nitera_first
     ENDIF

  ENDIF

! Particle fields
  IF ( icalc_part .EQ. 1 ) THEN
     WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
     CALL IO_READ_PARTICLE(fname, l_g, l_q)
     
! set boundarys for residence time pdf 
     IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN
        l_y_lambda =  (g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(1)%ymean - C_2_R
        l_y_base =   ((g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(1)%ymean -(g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(3)%ymean )/C_2_R &
             +  (g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(3)%ymean
        IF (residence_reset .EQ. 1) THEN
           l_q(:,6:7) = C_0_R
        ENDIF
     ENDIF
     
     IF ( itrajectory .NE. LAG_TRAJECTORY_NONE ) THEN
        CALL PARTICLE_TRAJECTORIES_INITIALIZE(nitera_save, nitera_last)
     END IF
     
  END IF
  
! Running average field
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. nitera_stats_spa .GT. 0 ) THEN
     WRITE(fname,*) nitera_first; fname = 'st'//TRIM(ADJUSTL(fname))
     CALL DNS_READ_AVGIJ(fname, vaux(vindex(VA_MEAN_WRK))) 
  ENDIF

! ###################################################################
! Initialize thermodynamic quantities
! ###################################################################
  CALL FI_PROFILES(wrk1d)
     
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Calculate q_l
        CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2), s(1,1), epbackground,pbackground)         

     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))

     ENDIF

  ELSE
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)

     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)

     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)

#ifdef CHEMISTRY
     IF ( ireactive .NE. CHEM_NONE ) THEN ! Calculate TGFM if reactive case
        CALL THERMO_GAMMA(imax,jmax,kmax, s, T, wrk3d)
        CALL CHEM_BURKESCHUMANN(imax,jmax,kmax, s, T, wrk3d, h_q(1,1))
     ENDIF
#endif

  ENDIF

! ###################################################################
! Check
! ###################################################################
  CALL DNS_CONTROL(i0, q,s, txc, wrk2d,wrk3d)
     
! ###################################################################
! Initialize data for boundary conditions
! ###################################################################
  CALL BOUNDARY_BUFFER_INITIALIZE(q,s, txc, wrk3d)

  CALL BOUNDARY_BCS_INITIALIZE(wrk3d)

  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL BOUNDARY_INFLOW_INITIALIZE(rtime, q_inf,s_inf, txc, wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Initialize LES 
! ###################################################################
#ifdef LES
  IF ( iles .EQ. 1 ) CALL LES_INI(q,s,h_q,h_s, txc, vaux, wrk1d,wrk2d,wrk3d)
#endif

! ###################################################################
! Initialize time step dt
! ###################################################################
  CALL TIME_COURANT(q,s, wrk3d)

! ###################################################################
! Initialize logfiles
! ###################################################################
  CALL DNS_LOGS(i1) ! headers
  CALL DNS_LOGS(i2) ! first line
  IF ( INT(logs_data(1)) .NE. 0 ) CALL DNS_STOP(INT(logs_data(1)))

! ###################################################################
! Do simulation: Integrate equations
! ###################################################################
  CALL TIME_INTEGRATION(q,h_q, s,h_s, q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
       l_q, l_hq, l_txc, l_comm)

! ###################################################################
#ifdef USE_FFTW
  IF ( ifourier .EQ. 1 ) THEN
     CALL dfftw_destroy_plan(fft_plan_fx)
     CALL dfftw_destroy_plan(fft_plan_bx)
     IF ( g(3)%size .GT. 1 ) THEN
        CALL dfftw_destroy_plan(fft_plan_fz)
        CALL dfftw_destroy_plan(fft_plan_bz)
     ENDIF
  ENDIF
#endif

  CALL DNS_END(0)

  STOP
END PROGRAM DNS

! ###################################################################
! ###################################################################
#ifdef USE_MPI

SUBROUTINE DNS_MPIO_AUX()

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : inb_flow_array, inb_scal_array
  USE DNS_LOCAL,  ONLY : nplanes_i,nplanes_j,nplanes_k, planes_i,planes_k
  USE DNS_MPI
  
  IMPLICIT NONE

#include "mpif.h" 

! -----------------------------------------------------------------------
  TINTEGER                :: ndims, idummy, id
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

! #######################################################################
  mpio_aux(:)%active = .FALSE. ! defaults
  mpio_aux(:)%offset = 0

! ###################################################################
! Subarray information to write planes
! ###################################################################
  idummy = inb_flow_array +inb_scal_array

  id = MPIO_SUBARRAY_PLANES_XOY
  IF ( nplanes_k .GT. 0 ) THEN ! Saving full vertical xOy planes; writing only info of PE containing the first plane
     IF ( ims_pro_k .EQ. ( planes_k(1) /kmax) ) mpio_aux(id)%active = .TRUE.
     mpio_aux(id)%communicator = ims_comm_x
     
     ndims = 2
     sizes(1)   = imax *ims_npro_i; sizes(2)   = jmax *nplanes_k *idummy
     locsize(1) = imax;             locsize(2) = jmax *nplanes_k *idummy
     offset(1)  = ims_offset_i;     offset(2)  = 0
     
     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(id)%subarray, ims_err)
     CALL MPI_Type_commit(mpio_aux(id)%subarray, ims_err)
     
  ENDIF

  id = MPIO_SUBARRAY_PLANES_ZOY
  IF ( nplanes_i .GT. 0 ) THEN ! Saving full vertical zOy planes; writing only info of PE containing the first plane
     IF ( ims_pro_i .EQ.  ( planes_i(1) /imax) ) mpio_aux(id)%active = .TRUE.
     mpio_aux(id)%communicator = ims_comm_z
     
     ndims = 2
     sizes(1)   = jmax *nplanes_i *idummy; sizes(2)   = kmax *ims_npro_k 
     locsize(1) = jmax *nplanes_i *idummy; locsize(2) = kmax 
     offset(1)  = 0;                       offset(2)  = ims_offset_k
     
     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(id)%subarray, ims_err)
     CALL MPI_Type_commit(mpio_aux(id)%subarray, ims_err)

  ENDIF

  id = MPIO_SUBARRAY_PLANES_XOZ
  IF ( nplanes_j .GT. 0 ) THEN ! Saving full blocks xOz planes
     mpio_aux(id)%active = .TRUE.
     mpio_aux(id)%communicator = MPI_COMM_WORLD
     
     ndims = 3 ! Subarray for the output of the 2D data
     sizes(1)  =imax *ims_npro_i; sizes(2)   = nplanes_j*idummy; sizes(3)   = kmax *ims_npro_k
     locsize(1)=imax;             locsize(2) = nplanes_j*idummy; locsize(3) = kmax
     offset(1) =ims_offset_i;     offset(2)  = 0;                offset(3)  = ims_offset_k
     
     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
          MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(id)%subarray, ims_err)
     CALL MPI_Type_commit(mpio_aux(id)%subarray, ims_err)

  ENDIF

  RETURN
END SUBROUTINE DNS_MPIO_AUX

#endif

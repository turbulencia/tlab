#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "DNS"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/06/12 - J.P. Mellado
!#              Formulation in terms of the energy instead of pressure
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM DNS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE LAGRANGE_GLOBAL
  USE DNS_LOCAL 
  USE DNS_TOWER
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
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: x,y,z, dx,dy,dz

! Flow/Scalar variables and RHS space
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q,s, h_q,h_s, txc  

! Particle data
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_hq

! Auxiliar memory space
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: vaux

  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_txc
  TREAL,      DIMENSION(:,:,:), ALLOCATABLE, SAVE :: l_trajectories
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: l_comm
  INTEGER(8), DIMENSION(:),   ALLOCATABLE, SAVE :: l_tags, l_trajectories_tags

! Work arrays
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk1d,wrk2d,wrk3d

! Inflow arrays
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: x_inf, y_inf, z_inf
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q_inf, s_inf

  TARGET q

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER           :: u, v, w, e, rho, p, T, vis

  TINTEGER iread_flow, iread_scal, idummy
  TINTEGER ierr, isize_wrk3d, isize_vaux, isize_loc
#ifdef USE_MPI
  TINTEGER id
#endif
  TREAL dummy
  CHARACTER*32 fname
  CHARACTER*64 str, line

#ifdef USE_MPI
  TINTEGER ndims_l, sizes_l(3), locsize_l(3), offset_l(3)
#endif

! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL('dns.ini')
#endif
#ifdef LES
  CALL LES_READ_INI('dns.ini')
#endif
  CALL DNS_READ_LOCAL('dns.ini')

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
  IF ( imode_rhs .EQ. EQNS_RHS_NONBLOCKING ) CALL DNS_NB3DFFT_INITIALIZE
#endif

  itime        = nitera_first
  logs_data(1) = 0 ! Status

! #######################################################################
! Definining types for parallel mode
! #######################################################################
#ifdef USE_MPI
! -------------------------------------------------------------------
! Filters at boundaries
! -------------------------------------------------------------------
  IF ( buff_nps_imax .GT. 1 ) THEN ! Required for outflow explicit filter in Ox
     CALL IO_WRITE_ASCII(lfile,'Initialize MPI types for Ox BCs explicit filter.')
     id    = DNS_MPI_K_OUTBCS
     isize_loc = buff_nps_imax*jmax
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

  IF ( buff_nps_jmin .GT. 1 ) THEN ! Required for outflow explicit filter in Oy
     CALL IO_WRITE_ASCII(lfile,'Initialize MPI types for Oy BCs explicit filter.')
     id    = DNS_MPI_K_TOPBCS
     isize_loc = imax*buff_nps_jmin
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

  IF ( ifilt_inflow .EQ. 1 ) THEN !  Required for inflow explicit filter
     CALL IO_WRITE_ASCII(lfile,'Initialize MPI types for inflow filter.')
     id    = DNS_MPI_K_INFLOW
     isize_loc = ifilt_inflow_iwidth*ifilt_inflow_jwidth
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF


! -------------------------------------------------------------------
! Characteristic BCs in compressible mode
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN

  IF ( i1bc .NE. DNS_BCS_PERIODIC ) THEN ! Required for NRBCs in Ox
     id    = DNS_MPI_K_NRBCX
     isize_loc = MOD(jmax,ims_npro_k)
     ims_bcs_imax = 2*(inb_flow+inb_scal_array)
     DO WHILE ( MOD(isize_loc*ims_bcs_imax,ims_npro_k) .GT. 0 ) 
        ims_bcs_imax = ims_bcs_imax + 1
     ENDDO
     WRITE(str,*) ims_bcs_imax
     str = 'Initialize MPI types for Ox BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
     CALL IO_WRITE_ASCII(lfile,str)
     isize_loc = ims_bcs_imax*jmax
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

  IF ( j1bc .NE. DNS_BCS_PERIODIC ) THEN ! Required for NRBCs in Oy
     id    = DNS_MPI_K_NRBCY
     isize_loc = MOD(imax,ims_npro_k)
     ims_bcs_jmax = 2*(inb_flow+inb_scal_array)
     DO WHILE ( MOD(isize_loc*ims_bcs_jmax,ims_npro_k) .GT. 0 ) 
        ims_bcs_jmax = ims_bcs_jmax + 1
     ENDDO
     WRITE(str,*) ims_bcs_jmax
     str = 'Initialize MPI types for Oy BCs transverse terms. '//TRIM(ADJUSTL(str))//' planes.'
     CALL IO_WRITE_ASCII(lfile,str)
     isize_loc = imax*ims_bcs_jmax
     CALL DNS_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, i1, i1, i1, i1, &
          ims_size_k(id), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF

  ENDIF

#endif

! #######################################################################
! Memory management
! #######################################################################
  isize_loc = MAX(imax_inf,MAX(jmax_inf,kmax_inf))
  isize_wrk1d = MAX(isize_wrk1d,isize_loc)

  isize_loc = MAX(imax_inf*jmax_inf,MAX(imax_inf*kmax_inf,jmax_inf*kmax_inf))
  isize_wrk2d = MAX(isize_wrk2d, isize_loc)
  IF ( icalc_particle .eq. 1) THEN
    isize_wrk2d = MAX(isize_wrk2d, jmax*inb_lag_total_interp)
  END IF

! txc
  inb_txc = 9
  IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR.  &
            imode_eqns .EQ. DNS_EQNS_ANELASTIC            ) THEN 
                                             inb_txc = 6
     IF ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) inb_txc = inb_txc+1
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL       .AND. &
            iadvection .EQ. EQNS_SKEWSYMMETRIC      .AND. &
            iviscous   .EQ. EQNS_EXPLICIT                 ) THEN
                                             inb_txc = 6
  ENDIF
  IF ( imixture .EQ. MIXT_TYPE_SUPSAT   ) inb_txc = inb_txc + 1

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
  isize_wrk3d = MAX(imax,imax_inf)*MAX(jmax,jmax_inf)*kmax
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)
  IF ( icalc_particle .eq. 1) THEN
     isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
     isize_wrk3d = MAX(isize_wrk3d,(jmax*(kmax+1)*inb_lag_total_interp*2))
     isize_wrk3d = MAX(isize_wrk3d,(jmax*(imax+1)*inb_lag_total_interp*2))
  END IF
  IF ( tower_mode .EQ. 1 ) THEN 
     isize_wrk3d = MAX(isize_wrk3d,nitera_save*(jmax_total+2))
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
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

! inflow
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     ALLOCATE(x_inf(imax_inf))
     ALLOCATE(y_inf(jmax_inf))
     ALLOCATE(z_inf(kmax_total))
     ALLOCATE(q_inf(imax_inf*jmax_inf*kmax_inf,inb_flow))
     ALLOCATE(s_inf(imax_inf*jmax_inf*kmax_inf,inb_scal_array))
  ELSE
     ALLOCATE(x_inf(1),y_inf(1),z_inf(1),q_inf(1,1),s_inf(1,1))
  ENDIF

! The case wiht icalc_scal = 0 in Ekman flow gives an error...
  iread_flow = icalc_flow; iread_scal = 1 !icalc_scal
#include "dns_alloc_arrays.h"

! -------------------------------------------------------------------
! Allocating additional memory space
! -------------------------------------------------------------------
! Rhs
  WRITE(str,*) inb_flow; line = 'Allocating array rhs flow. Size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_q(isize_field,inb_flow),    stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_q.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) inb_scal; line = 'Allocating array rhs scal. Size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_s(isize_field,inb_scal),    stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_s.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

! Array vaux
  CALL DNS_VAUX(isize_vaux)

  WRITE(str,*) isize_vaux; line = 'Allocating array vaux. Size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(vaux(isize_vaux),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for vaux.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

! Lagrangian part
  IF ( icalc_particle .EQ. 1 ) THEN
#include "dns_alloc_larrays.h"

     ALLOCATE(l_trajectories(3,num_trajectories,nitera_save),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF

     ALLOCATE(l_trajectories_tags(num_trajectories),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories_tags.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
     
     ALLOCATE(l_comm(isize_l_comm), stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
     
     ALLOCATE(l_hq(isize_particle,inb_particle),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_hq.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  
  ENDIF

! ###################################################################
! SUBARRAY FOR WRITE OF PLANES
! ###################################################################
#ifdef USE_MPI
  ndims_l = 3 ! Subarray for the output of the 2D data
  sizes_l(1)  =imax_total;   sizes_l(2)   = npln_j*(inb_flow+inb_scal); sizes_l(3)   = kmax_total
  locsize_l(1)=imax;         locsize_l(2) = npln_j*(inb_flow+inb_scal); locsize_l(3) = kmax
  offset_l(1) =ims_offset_i; offset_l(2)  = 0;                          offset_l(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims_l, sizes_l, locsize_l, offset_l, &
       MPI_ORDER_FORTRAN, MPI_REAL4, ims_subarray_pln_j, ims_err)
  CALL MPI_Type_commit(ims_subarray_pln_j, ims_err)
#endif

! #######################################################################
! Initializing tower stuff 
! #######################################################################
  IF ( tower_mode .EQ. 1 ) THEN 
!     CALL DNS_TOWER_INITIALIZE(x,y,z,tower_stride)  
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
  vaux(:) = C_0_R

  q(:,:)  = C_0_R; h_q(:,:) = C_0_R
  s(:,:)  = C_0_R; h_s(:,:) = C_0_R

  IF ( icalc_particle .EQ. 1 ) THEN ! Lagrangian
     l_q = C_0_R; l_hq = C_0_R
     l_trajectories = C_0_R; l_trajectories_tags = C_0_R
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN ! inflow domain for spatial case
     q_inf(:,:) = C_0_R
     s_inf(:,:) = C_0_R
  ENDIF

! ###################################################################
! Read the grid 
! ###################################################################
#include "dns_read_grid.h"

! ###################################################################
! Initialize constant vectors for filters
! ###################################################################
  IF      ( ifilt_domain .EQ. DNS_FILTER_4E  .OR. &
            ifilt_domain .EQ. DNS_FILTER_ADM      ) THEN  
     CALL FILT4E_INI(imax_total, i1bc, scalex, x, vaux(vindex(VA_FLT_CX)))
     CALL FILT4E_INI(jmax_total, j1bc, scaley, y, vaux(vindex(VA_FLT_CY)))
     CALL FILT4E_INI(kmax_total, k1bc, scalez, z, vaux(vindex(VA_FLT_CZ)))

  ELSE IF ( ifilt_domain .EQ. DNS_FILTER_COMPACT  ) THEN
     CALL FILT4C_INI(imax_total, i1bc, ifilt_alpha, dx, vaux(vindex(VA_FLT_CX)))
     CALL FILT4C_INI(jmax_total, j1bc, ifilt_alpha, dy, vaux(vindex(VA_FLT_CY)))
     CALL FILT4C_INI(kmax_total, k1bc, ifilt_alpha, dz, vaux(vindex(VA_FLT_CZ)))

  ENDIF

! ####################################################################
! Initializing position of l_q(particles)
! ####################################################################
  IF ( icalc_particle .EQ. 1 ) THEN
    WRITE(fname,*) nitera_first; fname = "particle_id."//TRIM(ADJUSTL(fname))
    CALL DNS_READ_PARTICLE_TAGS(fname,l_tags)
    
    WRITE(fname,*) nitera_first; fname = "particle."//TRIM(ADJUSTL(fname))
    CALL DNS_READ_PARTICLE(fname,l_q) ! h_particle only as dummy
    ! set boundarys for residence time pdf 
    IF (inb_particle_aux .EQ. 1) THEN
       l_y_lambda = ((y(jmax)-y(1))*ycoor_i(1)-2)
       l_y_base = ( ((y(jmax)-y(1))*ycoor_i(1)-(y(jmax)-y(1))*ycoor_i(3) )/2 + (y(jmax)-y(1))*ycoor_i(3))
       IF (residence_reset .EQ. 1) THEN
          l_q(:,6) = C_0_R
       ENDIF
    ENDIF
  END IF

! ###################################################################
! Initialize spectral Poisson solver
! ###################################################################
  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Check
! ###################################################################
  CALL DNS_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

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
     str = 'Changing original viscosity '//TRIM(ADJUSTL(str))//'to new value.'
     CALL IO_WRITE_ASCII(lfile,str)
     IF ( visctime .GT. C_0_R ) THEN ! Continuous change
        iviscchg = 1; visctime = (viscstart - viscstop)/visctime
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

! Running average field
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. frunstat .EQ. 1 ) THEN
     WRITE(fname,*) nitera_first; fname = 'st'//TRIM(ADJUSTL(fname))
     CALL DNS_READ_AVGIJ(fname, vaux(vindex(VA_MEAN_WRK))) 
  ENDIF

! ###################################################################
! Define pointers
! ###################################################################
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)

     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)

  ENDIF

! ###################################################################
! Initialize thermodynamic quantities
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL DNS_PROFILES(y, vaux(vindex(VA_BCS_VI)), wrk1d)

     IF      ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_AIRWATER_PHAL(i1, i1,i1, mean_i(2), p_init, mean_i(1)) !Calculate the mean water in eq
        CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)  !Calculate density (using the mean values)
        CALL THERMO_AIRWATER_PHAL(imax, jmax, kmax, s(:,2), p_init, s(:,1)) !Calculate the water in eq in the domain

     ELSE IF ( imixture .EQ. MIXT_TYPE_SUPSAT   ) THEN
        IF (damkohler(1) .GT. C_0_R)  THEN !supersaturation condition
	    CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)  !Calculate density (using the mean values)	    
	ELSE
	    CALL THERMO_AIRWATER_PHAL(i1, i1,i1, mean_i(2), p_init, mean_i(1)) !Calculate the mean water in eq
            CALL THERMO_THERMAL_DENSITY_HP_ALWATER(i1,i1,i1, mean_i(2),mean_i(1),p_init,mean_rho)  !Calculate density (using the mean values)
            CALL THERMO_AIRWATER_PHAL(imax, jmax, kmax, s(:,2), p_init, s(:,1)) !Calculate the water in eq in the domain
        ENDIF

     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))

     ENDIF

  ELSE
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

! check max/min values
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL DNS_CONTROL_FLOW(q,s, wrk3d)
  ENDIF
  CALL DNS_CONTROL_SCAL(s)
     
! ###################################################################
! Initialize data for boundary conditions
! ###################################################################
! buffer zone and reference pressures
  CALL BOUNDARY_INIT(vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
                     vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), &
                     vaux(vindex(VA_BCS_HT)),  vaux(vindex(VA_BCS_HB)),  &
                     vaux(vindex(VA_BCS_VI)),  vaux(vindex(VA_BCS_VO)),  &
                     x,y,z, dx,dy,dz, q,s, txc, wrk3d)

! inflow
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL BOUNDARY_INFLOW_INIT(rtime, y, x_inf,y_inf,z_inf, q_inf,s_inf, txc, wrk1d,wrk2d,wrk3d)
  ENDIF

! ###################################################################
! Initialize LES 
! ###################################################################
#ifdef LES
  IF ( iles .EQ. 1 ) CALL LES_INI(x,y,z,dx,dy,dz, q,s,h_q,h_s, txc, vaux, wrk1d,wrk2d,wrk3d)
#endif

! ###################################################################
! Initialize time step dt
! ###################################################################
  CALL TIME_COURANT(dx,dy,dz, q,s, wrk2d,wrk3d)

! ###################################################################
! Initialize logfiles
! ###################################################################
  CALL DNS_LOGS(i1) ! headers

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL FI_INVARIANT_P(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, q(1,1),q(1,2),q(1,3), h_q(1,1),h_q(1,2), wrk1d,wrk2d,wrk3d)
     CALL MINMAX(imax,jmax,kmax, h_q(1,1), logs_data(11),logs_data(10))
     logs_data(10)=-logs_data(10); logs_data(11)=-logs_data(11)
     IF ( MAX(ABS(logs_data(10)),ABS(logs_data(11))) .GT. d_bound_max ) THEN
        logs_data(1) = 1
     ENDIF
  ENDIF

  CALL DNS_LOGS(i2) ! first line

  IF ( logs_data(1) .NE. 0 ) CALL DNS_STOP(DNS_ERROR_DILATATION)

! ###################################################################
! Do simulation: Integrate equations
! ###################################################################
  CALL TIME_INTEGRATION(x,y,z,dx,dy,dz, q,h_q,s,h_s, &
       x_inf,y_inf,z_inf,q_inf,s_inf, txc, vaux, wrk1d,wrk2d,wrk3d, &
       l_q, l_hq, l_txc, l_tags, l_comm, l_trajectories, l_trajectories_tags)

! ###################################################################
#ifdef USE_FFTW
  IF ( ifourier .EQ. 1 ) THEN
     CALL dfftw_destroy_plan(fft_plan_fx)
     CALL dfftw_destroy_plan(fft_plan_bx)
     IF ( kmax_total .GT. 1 ) THEN
        CALL dfftw_destroy_plan(fft_plan_fz)
        CALL dfftw_destroy_plan(fft_plan_bz)
     ENDIF
  ENDIF
#endif

  CALL DNS_END(0)

  STOP
END PROGRAM DNS

#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE PARTICLE_READ_GLOBAL(inifile)
    
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_npro
#endif

  IMPLICIT NONE

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*64 lstr
  CHARACTER*32 bakfile
  TINTEGER idummy

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL IO_WRITE_ASCII(lfile, 'Reading particle input data.')

! -------------------------------------------------------------------
  CALL IO_WRITE_ASCII(bakfile,  '#')
  CALL IO_WRITE_ASCII(bakfile,  '#[Lagrange]')
  CALL IO_WRITE_ASCII(bakfile,  '#Type=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_number=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_rnd_mode=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Y_Particle_Pos=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Y_Particle_Width=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_bumper=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Jmax_part=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Jmin_part=<value>')

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'Type', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'              ) THEN; ilagrange = LAG_TYPE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tracer'            ) THEN; ilagrange = LAG_TYPE_TRACER
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'simplesettling'    ) THEN; ilagrange = LAG_TYPE_SIMPLE_SETT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinearcloudthree') THEN; ilagrange = LAG_TYPE_BIL_CLOUD_3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinearcloudfour' ) THEN; ilagrange = LAG_TYPE_BIL_CLOUD_4
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'PARTICLE_READ_GLOBAL. Wrong lagrangian model.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINILONGINT(bakfile, inifile, 'Lagrange', 'Particle_number', '0', particle_number  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Particle_bumper', '2.0', particle_bumper  )
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Jmax_part', '1', jmax_part  )
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Jmin_part', '1', jmin_part  )

! -------------------------------------------------------------------
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Particle_rnd_mode', '1', particle_rnd_mode  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_Pos', '0.5', y_particle_pos  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_Width', '1.0', y_particle_width  )

! -------------------------------------------------------------------
  inb_trajectory = 3 ! Default, just position
  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'TrajectoryType', 'first', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'first'     ) THEN; itrajectory = LAG_TRAJECTORY_FIRST
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'vorticity' ) THEN; itrajectory = LAG_TRAJECTORY_VORTICITY
     inb_trajectory = 3 + 3 + 1 ! position + vorticity + buoyancy
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'largest'   ) THEN; itrajectory = LAG_TRAJECTORY_LARGEST
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; itrajectory = LAG_TRAJECTORY_NONE
  ELSE
     CALL IO_WRITE_ASCII(efile,'PARTICLE_READ_GLOBAL. Invalid option in TrajectoryType')
     CALL DNS_STOP(DNS_ERROR_CALCTRAJECTORIES)
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'TrajectoryNumber', '0', isize_trajectory)
  IF ( isize_trajectory .GT. particle_number ) THEN
     CALL IO_WRITE_ASCII(efile,'PARTICLE_READ_GLOBAL. Number of trajectories must be less or equal than number of particles.')
     CALL DNS_STOP(DNS_ERROR_CALCTRAJECTORIES)
  ENDIF
  IF ( isize_trajectory .LE. 0 ) itrajectory = LAG_TRAJECTORY_NONE
  IF ( icalc_part .EQ. 0   ) itrajectory = LAG_TRAJECTORY_NONE
     
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'CalculateParticlePDF', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_part_pdf = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_part_pdf = 0
  ENDIF
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_PDF_Pos', '0.5', y_particle_pdf_pos  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_PDF_Width', '1.0', y_particle_pdf_width  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'X_Particle_PDF_Pos', '0.0', x_particle_pdf_pos  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'X_Particle_PDF_Width', '0.0', x_particle_pdf_width  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Z_Particle_PDF_Pos', '0.0', z_particle_pdf_pos  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Z_Particle_PDF_Width', '0.0', z_particle_pdf_width  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Particle_PDF_Max', '10', particle_pdf_max  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Particle_PDF_Interval', '0.5', particle_pdf_interval  )

  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'ResidenceReset', 'yes', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; residence_reset = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; residence_reset = 0
  ELSE
     CALL IO_WRITE_ASCII(efile,'PARTICLE_READ_GLOBAL. ResidenceReset must be yes or no')
     CALL DNS_STOP(DNS_ERROR_RESIDENCERESET)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'Parameters', '0.0', sRes)
  idummy = MAX_LAGPARAM
  CALL LIST_REAL(sRes, idummy, lagrange_param)

! ###################################################################
! Initializing size of Lagrangian arrays
! ###################################################################
  CALL PARTICLE_TYPE_INITIALIZE

  WRITE(lstr,*) inb_particle 
  CALL IO_WRITE_ASCII(lfile, 'Initialize inb_particle = '//TRIM(ADJUSTL(lstr)))
  WRITE(lstr,*) inb_lag_total_interp 
  CALL IO_WRITE_ASCII(lfile, 'Initialize inb_lag_total_interp = '//TRIM(ADJUSTL(lstr)))
  WRITE(lstr,*) inb_lag_aux_field
  CALL IO_WRITE_ASCII(lfile, 'Initialize inb_lag_aux_field = '//TRIM(ADJUSTL(lstr)))

#ifdef USE_MPI
!     isize_particle=INT(particle_number/INT(ims_npro,KIND=8)*INT(particle_bumper,KIND=8))
  isize_particle=INT(particle_number/INT(ims_npro,KIND=8)*INT(particle_bumper*100,KIND=8)/INT(100,KIND=8))
#else
  isize_particle=INT(particle_number)
#endif

  isize_hf_1 = 2     *jmax*kmax   *inb_lag_total_interp 
  isize_hf_2 =   imax*jmax     *2 *inb_lag_total_interp 
  isize_hf_3 = 2     *jmax     *2 *inb_lag_total_interp 
  isize_max_hf  = isize_hf_1+isize_hf_2+isize_hf_3
  isize_pbuffer = int(isize_particle/4*(inb_particle*2+1) ) !same size for both buffers
  isize_l_comm  = isize_hf_1+isize_hf_2+isize_hf_3+2*isize_pbuffer

  idummy = MAX((imax+1)*jmax, MAX((imax+1)*kmax,jmax*(kmax+1)))
  isize_wrk2d = MAX(isize_wrk2d,idummy)

  RETURN  
END SUBROUTINE PARTICLE_READ_GLOBAL

! ###################################################################
! ###################################################################
SUBROUTINE PARTICLE_TYPE_INITIALIZE
  
  USE LAGRANGE_GLOBAL

  USE DNS_GLOBAL, ONLY : inb_particle, inb_particle_txc
  IMPLICIT NONE

! -------------------------------------------------------------------
  inb_particle_txc = 0

  IF   (ilagrange .EQ. LAG_TYPE_TRACER) THEN
     inb_particle_evolution = 3
     inb_particle_aux = 0          
     inb_particle_txc = 0
     inb_lag_aux_field = 0
     inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
     
  ELSEIF  (ilagrange .EQ. LAG_TYPE_SIMPLE_SETT) THEN
     inb_particle_evolution = 3
     inb_particle_aux = 0          
     inb_particle_txc = 0
     inb_lag_aux_field = 0
     inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
     
  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3) THEN
     inb_particle_evolution = 5    !amount of particle properties 
     inb_particle_aux = 0          !amount of particle properties without runge kutta (only sent and sorted)
     inb_particle_txc = 1          !l_txc properties
     inb_lag_aux_field = 4         !field data on txc
     inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
     LAGRANGE_SPNAME(1) = 'droplet_diff_3'
     LAGRANGE_SPNAME(2) = 'droplet_nodiff_3'
     
  ELSEIF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
     inb_particle_evolution = 5    !amount of particle properties with runge kutta
     inb_particle_aux = 1          !amount of particle properties without runge kutta (only sent and sorted)
     inb_particle_txc = 1          !l_txc properties
     inb_lag_aux_field = 4         !field data on txc
     inb_particle = inb_particle_evolution + inb_particle_aux    !amount of particle properties which are sent
     LAGRANGE_SPNAME(1) = 'droplet_diff_3'
     LAGRANGE_SPNAME(2) = 'droplet_nodiff_3'
     LAGRANGE_SPNAME(3) = 'residence_part'
     
  END IF
  
  inb_scal_particle = inb_particle_evolution - 3          !Number of scalar properties solved in the lagrangian
  inb_lag_total_interp = inb_lag_aux_field + 3  !Number of fields needed by lagrangian (no extra memory, usually txc fields)

  RETURN
END SUBROUTINE PARTICLE_TYPE_INITIALIZE

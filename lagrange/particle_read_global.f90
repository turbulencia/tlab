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

! ###################################################################
! Main block
! ###################################################################
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Lagrange', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'              ) THEN; ilagrange = LAG_TYPE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tracer'            ) THEN; ilagrange = LAG_TYPE_TRACER
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'simplesettling'    ) THEN; ilagrange = LAG_TYPE_SIMPLE_SETT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinearcloudthree') THEN; ilagrange = LAG_TYPE_BIL_CLOUD_3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinearcloudfour' ) THEN; ilagrange = LAG_TYPE_BIL_CLOUD_4
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong lagrangian model.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! ###################################################################
! Lagrange block
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,  '#')
  CALL IO_WRITE_ASCII(bakfile,  '#[Lagrange]')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_number=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_rnd_mode=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#inb_particle=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Y_Particle_Pos=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Y_Particle_Width=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Particle_bumper=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Jmax_part=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Jmin_part=<value>')

  CALL SCANINILONGINT(bakfile, inifile, 'Lagrange', 'Particle_number', '0', particle_number  )
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Particle_rnd_mode', '1', particle_rnd_mode  )
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Jmax_part', '1', jmax_part  )
  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Jmin_part', '1', jmin_part  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_Pos', '0.5', y_particle_pos  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Y_Particle_Width', '1.0', y_particle_width  )
  CALL SCANINIREAL(bakfile, inifile, 'Lagrange', 'Particle_bumper', '2.0', particle_bumper  )

  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'CalculateTrajectories', 'none', sRes)
  !Plot all droplets. This version is NOT optimized and writes to disc every time step
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'all' ) THEN; icalc_trajectories = LAG_TRAJECTORY_ALL 
  !plot the 'Num_trajectories' largest droplets at the simulation last time step from a file created lagrange_traject.x. This function is optimized and does not write
  ! every tinme step to disc 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'largest') THEN; icalc_trajectories = LAG_TRAJECTORY_LARGEST 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no' .OR. TRIM(ADJUSTL(sRes)) .eq. 'none' ) THEN; icalc_trajectories = LAG_TRAJECTORY_NONE !Alberto
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Invalid option in CalculateTrajectories')
     CALL DNS_STOP(DNS_ERROR_CALCTRAJECTORIES)
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Lagrange', 'Num_trajectories', '50', num_trajectories  )

  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'CalculateParticlePDF', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_particle_pdf = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_particle_pdf = 0
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
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. ResidenceReset must be yes or no')
     CALL DNS_STOP(DNS_ERROR_RESIDENCERESET)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Lagrange', 'Parameters', '0.0', sRes)
  idummy = MAX_LAGPARAM
  CALL LIST_REAL(sRes, idummy, lagrange_param)

! ###################################################################
! Initializing size of Lagrangian arrays
! ###################################################################
  CALL LAGRANGE_TYPE_INITIALIZE

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

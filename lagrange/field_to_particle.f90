#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
SUBROUTINE  FIELD_TO_PARTICLE &
    (nvar, data_in, data_out, l_q,l_hq,l_tags,l_comm, wrk2d,wrk3d)

  USE DNS_CONSTANTS,  ONLY : efile, lfile
  USE DNS_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax, isize_particle
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI,        ONLY:  ims_err, ims_pro, ims_size_p
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar)                 :: data_in
  TYPE(pointers_dt),   DIMENSION(nvar)                 :: data_out
  TREAL,               DIMENSION(isize_particle,*)     :: l_q, l_hq
  INTEGER(8),          DIMENSION(isize_particle)       :: l_tags
  TREAL,               DIMENSION(isize_l_comm), TARGET :: l_comm
  TREAL,               DIMENSION(*)                    :: wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal 
  TINTEGER npar_start, npar
  TINTEGER ip1,ip2,ip3, np1,np2,np3, iv

  TYPE(pointers3d_dt), DIMENSION(nvar) :: data_halo1, data_halo2, data_halo3

! Check
  CHARACTER*(32) str
  CHARACTER*(128) line
  TINTEGER particle_number_local
  
!#######################################################################
  IF ( nvar .GT. inb_particle_interp ) THEN
     CALL IO_WRITE_ASCII(efile,'FIELD_TO_PARTICLE. Not enough memory.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  np1 = 2*jmax*kmax; ip1 = 1
  np2 = imax*jmax*2; ip2 = 1 +isize_hf_1
  np3 = 2   *jmax*2; ip3 = 1 +isize_hf_1 +isize_hf_2
  DO iv = 1,nvar
     data_halo1(iv)%field(1:2,1:jmax,1:kmax) => l_comm(ip1:ip1+np1-1); ip1 = ip1 +np1
     data_halo2(iv)%field(1:imax,1:jmax,1:2) => l_comm(ip2:ip2+np2-1); ip2 = ip2 +np2
     data_halo3(iv)%field(1:2,1:jmax,1:2)    => l_comm(ip3:ip3+np3-1); ip3 = ip3 +np3
  ENDDO

! #######################################################################
! Setting fields in halo regions
! ######################################################################
#ifdef USE_MPI
  CALL PARTICLE_HALO_K(nvar, data_in, data_halo2(1)%field, wrk3d(1), wrk3d(imax*jmax*nvar+1), &
       wrk2d(1), wrk2d(2*(jmax*nvar+1)))

  CALL PARTICLE_HALO_I(nvar, data_in, data_halo1(1)%field, data_halo3(1)%field, wrk3d(1), wrk3d(jmax*(kmax+1)*nvar+1),&
       wrk2d(1), wrk2d(jmax*nvar+1), wrk2d(2*(jmax*nvar+1)))

#else
  CALL PARTICLE_HALO_SERIAL(nvar, data_in, data_halo1(1)%field, data_halo2(1)%field, data_halo3(1)%field)

#endif
  
!#######################################################################
! Sorting and counting particles for each zone
!#######################################################################
  CALL PARTICLE_SORT_HALO(grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal,&
       l_hq, l_tags, l_q)

! Check
#ifdef USE_MPI
   particle_number_local = ims_size_p(ims_pro+1)
#else
   particle_number_local = INT(particle_number)
#endif

  npar = grid_zone+ halo_zone_x+ halo_zone_z+ halo_zone_diagonal
  IF ( particle_number_local .NE. npar ) THEN
     WRITE(str,*) npar
     line = 'FIELD_TO_PARTICLE. Npar is '//TRIM(ADJUSTL(str))
     WRITE(str,*) particle_number_local
     line = TRIM(ADJUSTL(line))//' whereas np_local is '//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(efile,line)
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF
! End of check

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
#endif
  
!#######################################################################
! Interpolating
!#######################################################################
  npar_start = 1
  npar       = grid_zone
  CALL PARTICLE_INTERPOLATION(i0, nvar, data_in, data_out, l_q, npar_start, npar)

  IF ( halo_zone_x .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_x
     CALL PARTICLE_INTERPOLATION(i1, nvar, data_halo1, data_out, l_q, npar_start,npar)
  END IF
  
  IF ( halo_zone_z .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_z
     CALL PARTICLE_INTERPOLATION(i2, nvar, data_halo2, data_out, l_q, npar_start,npar)
 END IF

  IF ( halo_zone_diagonal .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_diagonal
     CALL PARTICLE_INTERPOLATION(i3, nvar, data_halo3, data_out, l_q, npar_start,npar)
  END IF

  DO iv = 1,nvar
     NULLIFY(data_halo1(iv)%field,data_halo2(iv)%field,data_halo3(iv)%field)
  ENDDO
  
  RETURN
END SUBROUTINE FIELD_TO_PARTICLE

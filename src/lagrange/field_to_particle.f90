#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
SUBROUTINE  FIELD_TO_PARTICLE &
    (nvar, data_in, data_out, l_g,l_q,l_comm, wrk3d)

  USE TLAB_CONSTANTS,  ONLY : efile, lfile
  USE TLAB_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE TLAB_VARS,     ONLY : imax,jmax,kmax, isize_particle
  USE TLAB_PROCS
  USE LAGRANGE_VARS,ONLY : particle_dt, isize_l_comm, inb_particle_interp
#ifdef USE_MPI
  USE TLAB_MPI_VARS,        ONLY:  ims_err
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar)                 :: data_in
  TYPE(pointers_dt),   DIMENSION(nvar)                 :: data_out
  TYPE(particle_dt)                                    :: l_g
  TREAL,               DIMENSION(isize_particle,*)     :: l_q
  TREAL,               DIMENSION(isize_l_comm), TARGET :: l_comm
  TREAL,               DIMENSION(*)                    :: wrk3d

! -------------------------------------------------------------------
  TINTEGER grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal
  TINTEGER npar_start, npar
  TINTEGER ip1,ip2,ip3, np1,np2,np3, iv

  TYPE(pointers3d_dt), DIMENSION(nvar) :: data_halo_i, data_halo_k, data_halo_ik

!#######################################################################
  IF ( nvar .GT. inb_particle_interp ) THEN
     CALL TLAB_WRITE_ASCII(efile,'FIELD_TO_PARTICLE. Not enough memory.')
     CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  np1 = 2*jmax*kmax; ip1 = 1
  np2 = imax*jmax*2; ip2 = ip1 +np1 *nvar !isize_hf_1
  np3 = 2   *jmax*2; ip3 = ip2 +np2 *nvar !isize_hf_2
  DO iv = 1,nvar
     data_halo_i(iv)%field(1:2,1:jmax,1:kmax) => l_comm(ip1:ip1+np1-1); ip1 = ip1 +np1
     data_halo_k(iv)%field(1:imax,1:jmax,1:2) => l_comm(ip2:ip2+np2-1); ip2 = ip2 +np2
     data_halo_ik(iv)%field(1:2,1:jmax,1:2)   => l_comm(ip3:ip3+np3-1); ip3 = ip3 +np3
  ENDDO

! -------------------------------------------------------------------
! Setting fields in halo regions
  CALL PARTICLE_HALO_K(nvar, data_in, data_halo_k(1)%field, wrk3d(1), wrk3d(imax*jmax*nvar+1))
  CALL PARTICLE_HALO_I(nvar, data_in, data_halo_i(1)%field, data_halo_k(1)%field, data_halo_ik(1)%field, wrk3d(1), wrk3d(jmax*(kmax+1)*nvar+1))

! -------------------------------------------------------------------
! Sorting and counting particles for each zone
  CALL PARTICLE_SORT_HALO(l_g,l_q, nvar,data_out, grid_zone,halo_zone_x,halo_zone_z,halo_zone_diagonal)

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
#endif

! -------------------------------------------------------------------
! Interpolating
  npar_start = 1
  npar       = grid_zone
  CALL FIELD_TO_PARTICLE_INTERPOLATE(i0, nvar, data_in, data_out, l_q, npar_start, npar)

  IF ( halo_zone_x .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_x
     CALL FIELD_TO_PARTICLE_INTERPOLATE(i1, nvar, data_halo_i, data_out, l_q, npar_start,npar)
  END IF

  IF ( halo_zone_z .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_z
     CALL FIELD_TO_PARTICLE_INTERPOLATE(i2, nvar, data_halo_k, data_out, l_q, npar_start,npar)
 END IF

  IF ( halo_zone_diagonal .NE. 0 ) THEN
     npar_start = npar +1
     npar       = npar +halo_zone_diagonal
     CALL FIELD_TO_PARTICLE_INTERPOLATE(i3, nvar, data_halo_ik, data_out, l_q, npar_start,npar)
  END IF

! -------------------------------------------------------------------
  DO iv = 1,nvar
     NULLIFY(data_halo_i(iv)%field,data_halo_k(iv)%field,data_halo_ik(iv)%field)
  ENDDO

  RETURN
END SUBROUTINE FIELD_TO_PARTICLE

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_K(nvar, data, halo_field_k, buffer_send, buffer_recv)

  USE TLAB_TYPES,  ONLY : pointers3d_dt
  USE TLAB_VARS, ONLY : imax,jmax,kmax

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro, ims_pro_k, ims_npro_k, ims_map_k
  USE TLAB_MPI_VARS, ONLY : ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar):: data
  TREAL, DIMENSION(imax,jmax,2,nvar)  :: halo_field_k
  TREAL, DIMENSION(imax,jmax,1,nvar)  :: buffer_send, buffer_recv

! -------------------------------------------------------------------
  TINTEGER i

#ifdef USE_MPI
  integer source, dest, l, size
  integer mpireq(ims_npro*2+2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)
#endif

! ######################################################################
#ifdef USE_MPI
  IF ( ims_npro_k .EQ. 1 ) THEN
#endif
     DO i = 1,nvar
        halo_field_k(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        halo_field_k(1:imax,1:jmax,2,i) = data(i)%field(1:imax,1:jmax,1   )
     ENDDO

#ifdef USE_MPI
  ELSE
     DO i = 1,nvar
        halo_field_k(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        buffer_send(1:imax,1:jmax,1,i)  = data(i)%field(1:imax,1:jmax,1   ) ! data to be transfered
     ENDDO

     mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
     l      = 2*ims_pro +1
     dest   = ims_map_k(MOD(ims_pro_k-1 +ims_npro_k,ims_npro_k) +1)
     source = ims_map_k(MOD(ims_pro_k+1 +ims_npro_k,ims_npro_k) +1)
     size   = imax*jmax*nvar
     CALL MPI_ISEND(buffer_send,size,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,size,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

     halo_field_k(1:imax,1:jmax,2,1:nvar) = buffer_recv(1:imax,1:jmax,1,1:nvar)

  END IF
#endif

  RETURN
END SUBROUTINE PARTICLE_HALO_K

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_I(nvar, data, halo_field_i, halo_field_k, halo_field_ik, buffer_send, buffer_recv)

  USE TLAB_TYPES,      ONLY : pointers3d_dt
  USE TLAB_VARS,     ONLY : imax,jmax,kmax

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro, ims_pro_i, ims_npro_i, ims_map_i
  USE TLAB_MPI_VARS, ONLY : ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar):: data
  TREAL, DIMENSION(2,jmax,kmax,nvar)  :: halo_field_i
  TREAL, DIMENSION(imax,jmax,2,nvar)  :: halo_field_k
  TREAL, DIMENSION(2,jmax,2,nvar)     :: halo_field_ik
  TREAL, DIMENSION(1,jmax,kmax+1,nvar):: buffer_send, buffer_recv

! -------------------------------------------------------------------
  TINTEGER i

#ifdef USE_MPI
  integer source, dest, l, size
  integer mpireq(ims_npro*2+2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)
#endif

! ######################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .EQ. 1 ) THEN
#endif
     DO i = 1,nvar
        halo_field_i (1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax)
        halo_field_i (2,1:jmax,1:kmax,i) = data(i)%field(1,   1:jmax,1:kmax)
        halo_field_ik(2,1:jmax,2     ,i) = halo_field_k (1,   1:jmax,2     ,i) ! top-right corner
     END DO

#ifdef USE_MPI
  ELSE
     DO i = 1,nvar
        halo_field_i(1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax)
        buffer_send(1,1:jmax,1:kmax,i)  = data(i)%field(1,   1:jmax,1:kmax)   ! data to be transfered
        buffer_send(1,1:jmax,kmax+1,i)  = halo_field_k (1,   1:jmax,2     ,i)
     ENDDO

     mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
     l      = 2*ims_pro +1
     dest   = ims_map_i(MOD(ims_pro_i-1 +ims_npro_i,ims_npro_i) +1)
     source = ims_map_i(MOD(ims_pro_i+1 +ims_npro_i,ims_npro_i) +1)
     size   = jmax*(kmax+1)*nvar
     CALL MPI_ISEND(buffer_send,size,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,size,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

     halo_field_i (2,1:jmax,1:kmax,1:nvar) = buffer_recv(1,1:jmax,1:kmax,1:nvar)
     halo_field_ik(2,1:jmax,2     ,1:nvar) = buffer_recv(1,1:jmax,kmax+1,1:nvar) ! top-right corner

  END IF
#endif

  halo_field_ik(1,1:jmax,1,1:nvar) = halo_field_i(1,   1:jmax,kmax,1:nvar)
  halo_field_ik(2,1:jmax,1,1:nvar) = halo_field_i(2,   1:jmax,kmax,1:nvar)
  halo_field_ik(1,1:jmax,2,1:nvar) = halo_field_k(imax,1:jmax,2   ,1:nvar)

  RETURN
END SUBROUTINE PARTICLE_HALO_I

!########################################################################
!########################################################################
SUBROUTINE FIELD_TO_PARTICLE_INTERPOLATE &
     (iflag, nvar, data_in, data_out, l_q, grid_start, grid_end)

  USE TLAB_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE TLAB_VARS,     ONLY : isize_particle
  USE TLAB_VARS,     ONLY : g
  USE LAGRANGE_VARS,ONLY : l_g
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY: ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"

  TINTEGER iflag, nvar, grid_start, grid_end
  TYPE(pointers3d_dt), DIMENSION(nvar)             :: data_in
  TYPE(pointers_dt),   DIMENSION(nvar)             :: data_out
  TREAL,               DIMENSION(isize_particle,3) :: l_q

! -------------------------------------------------------------------
  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  g_p(10), g1loc, g2loc, g5loc, g6loc
  TINTEGER i, j
  TREAL dx_loc_inv, dz_loc_inv

! ######################################################################
  dx_loc_inv = M_REAL( g(1)%size ) /g(1)%scale
  dz_loc_inv = M_REAL( g(3)%size ) /g(3)%scale

! Managing the iflag option outside the loop
  g_p(7) =1
  g_p(8) =2
  g_p(9) =1
  g_p(10)=2

  g1loc = 1
  g2loc = 2
  g5loc = 5
  g6loc = 6
  IF      ( iflag .EQ. 1 ) THEN ! halo_zone_x
     g1loc = 7
     g2loc = 8
  ELSE IF ( iflag .EQ. 2 ) THEN ! halo_zone_z
     g5loc = 9
     g6loc = 10
  ELSE IF ( iflag .EQ. 3 ) THEN ! halo_zone_diagonal
     g1loc = 7
     g2loc = 8
     g5loc = 9
     g6loc = 10
  ENDIF

! ######################################################################
  IF  ( g(3)%size .NE. 1 ) THEN

     DO i = grid_start,grid_end ! loop over all particles

        length_g_p(1) = l_q(i,1) *dx_loc_inv            ! Local X position
        g_p(1)        = FLOOR( length_g_p(1) )
        length_g_p(1) = length_g_p(1) -M_REAL( g_p(1) )
#ifdef USE_MPI
        g_p(1)        = g_p(1) +1 -ims_offset_i
#else
        g_p(1)        = g_p(1) +1
#endif
        g_p(2)        = g_p(1) +1
        length_g_p(2) = C_1_R -length_g_p(1)

        length_g_p(5) = l_q(i,3) *dz_loc_inv            ! Local Z position
        g_p(5)        = FLOOR( length_g_p(5) )
        length_g_p(5) = length_g_p(5) -M_REAL( g_p(5) )
#ifdef USE_MPI
        g_p(5)        = g_p(5) +1 -ims_offset_k
#else
        g_p(5)        = g_p(5) +1
#endif
        g_p(6)        = g_p(5) +1
        length_g_p(6) = C_1_R -length_g_p(5)

        g_p(3)        = l_g%nodes(i)                    ! Local Y position
        g_p(4)        = g_p(3) +1
        length_g_p(3) =(l_q(i,2) - g(2)%nodes(g_p(3))) /(g(2)%nodes(g_p(4))-g(2)%nodes(g_p(3)))
        length_g_p(4) = C_1_R -length_g_p(3)

        cube_g_p(1) = length_g_p(1) *length_g_p(3) ! bilear cubes for X and Y
        cube_g_p(2) = length_g_p(1) *length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

! -------------------------------------------------------------------
! Trilinear interpolation
! Two bilinear interpolations for each k plane (g_p(5) and g_p(6)
! Then multipled by (1-length) for trilinear aspect
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                ((cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),g_p(g5loc)) &
                 +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),g_p(g5loc)) &
                 +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),g_p(g5loc)) &
                 +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),g_p(g5loc)))*length_g_p(6)) &
               +((cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),g_p(g6loc)) &
                 +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),g_p(g6loc)) &
                 +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),g_p(g6loc)) &
                 +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),g_p(g6loc)))*length_g_p(5))
        ENDDO

     END DO

! ######################################################################
  ELSE !2D case

     DO i = grid_start,grid_end

        length_g_p(1) = l_q(i,1) *dx_loc_inv
        g_p(1)        = FLOOR( length_g_p(1) )
        length_g_p(1) = length_g_p(1) -M_REAL( g_p(1) )
#ifdef USE_MPI
        g_p(1)        = g_p(1) +1 -ims_offset_i
#else
        g_p(1)        = g_p(1) +1
#endif
        g_p(2)        = g_p(1) +1
        length_g_p(2) = C_1_R -length_g_p(1)

        g_p(3)        = l_g%nodes(i)
        g_p(4)        = g_p(3) +1
        length_g_p(3) =(l_q(i,2) - g(2)%nodes(g_p(3))) /(g(2)%nodes(g_p(4))-g(2)%nodes(g_p(3)))
        length_g_p(4) = C_1_R -length_g_p(3)

        cube_g_p(1) = length_g_p(1) *length_g_p(3)
        cube_g_p(2) = length_g_p(1) *length_g_p(4)
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

! -------------------------------------------------------------------
! Bilinear interpolation
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                (cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),1) &
                +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),1) &
                +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),1) &
                +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),1))
        ENDDO

     END DO

  ENDIF

  RETURN
END SUBROUTINE FIELD_TO_PARTICLE_INTERPOLATE

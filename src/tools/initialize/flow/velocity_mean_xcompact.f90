#include "types.h"
#include "dns_const.h"

subroutine VELOCITY_MEAN_XCOMPACT(u,v,w)

  use TLAB_VARS,     only: qbg, g, imax,jmax,kmax

#ifdef USE_MPI 
  use TLAB_MPI_VARS, only: ims_pro, ims_offset_j
#endif 

  implicit none

#include "integers.h"

#ifdef USE_MPI
#include "mpif.h"
#else
  TINTEGER, parameter                           :: ims_pro=0, ims_offset_j=0
#endif

  TREAL, dimension(imax,jmax,kmax), intent(out) :: u, v, w

  ! -------------------------------------------------------------------
  TINTEGER                                      :: i, j, k, ii, code, jstart
  TREAL                                         :: um, y
  TREAL, parameter                              :: init_noise = 0.125

  !########################################################################

  ! global array indicies for each mpi task (indices start with 0)
  jstart = ims_offset_j ! jend = ims_offset_j + jmax - 1

  ! ini velocities with zeros
  u = C_0_R; v = C_0_R; w = C_0_R

  ! random
  call system_clock(count = code)
  call random_seed( size  = ii)
  call random_seed( put   = code + ims_pro*(/ (i - 1, i = 1, ii) /))
  !
  call random_number(u)
  call random_number(v)
  call random_number(w)

  !modulation of the random noise + initial velocity profile (requires ly=2*delta)
  do k=1,kmax
     do j=1,jmax

        y  = g(2)%nodes(j+jstart) - g(2)%nodes(g(2)%size) / C_2_R
        
        um = exp(-0.2 * y**C_2_R)

        do i=1,imax
           u(i,j,k) = (init_noise * um * (C_2_R * u(i,j,k) - C_1_R) + C_1_R - y**C_2_R) * qbg(1)%delta
           v(i,j,k) = (init_noise * um * (C_2_R * v(i,j,k) - C_1_R)                   ) * qbg(1)%delta
           w(i,j,k) = (init_noise * um * (C_2_R * w(i,j,k) - C_1_R)                   ) * qbg(1)%delta
        end do
     end do
  end do

  ! -------------------------------------------------------------------
  if ( g(3)%size .EQ. 1 ) w = C_0_R

  return
end subroutine VELOCITY_MEAN_XCOMPACT
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2021/XX/XX - J.K. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!#
!#
!#
!#
!########################################################################
subroutine FI_CHANNEL_CFR_FORCING(u, v, h1, h2, wrk1d, wrk3d)

  use TLAB_VARS, only : jmax, isize_field
  use TLAB_VARS, only : ubulk, ubulk_parabolic, delta_ubulk

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u     
  TREAL, dimension(isize_field), intent(inout) :: v     
  TREAL, dimension(isize_field), intent(inout) :: h1    
  TREAL, dimension(isize_field), intent(inout) :: h2      
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

! -----------------------------------------------------------------------

  TINTEGER                                     :: ij

! #######################################################################

  ! current ubulk of the flow
  call FI_CHANNEL_UBULK(u, wrk1d, wrk3d) 

  delta_ubulk = ubulk_parabolic - ubulk

  ! constant flow rate forcing
  do ij = 1, isize_field
    h1(ij) = h1(ij) + delta_ubulk
  end do

  ! rotate turbulent channelflow (reducing spinup time)
  call FI_CHANNEL_SPINUP(u, v, h1, h2)
  
  return
end subroutine FI_CHANNEL_CFR_FORCING
!########################################################################
subroutine FI_CHANNEL_CPG_FORCING(u, v, h1, h2, wrk1d, wrk3d)

  use TLAB_VARS, only : jmax, isize_field
  use TLAB_VARS, only : fcpg

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u     
  TREAL, dimension(isize_field), intent(inout) :: v     
  TREAL, dimension(isize_field), intent(inout) :: h1    
  TREAL, dimension(isize_field), intent(inout) :: h2      
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d
! -----------------------------------------------------------------------

  TINTEGER                                     :: ij

! #######################################################################

  ! current ubulk of the flow
  call FI_CHANNEL_UBULK(u, wrk1d, wrk3d) 

  ! add constant streamwise pressure gradient
  do ij = 1, isize_field
    h1(ij) = h1(ij) + fcpg
  end do

  ! rotate turbulent channel flow (reducing spinup time)
  call FI_CHANNEL_SPINUP(u, v, h1, h2)
  
  return
end subroutine FI_CHANNEL_CPG_FORCING
!########################################################################
subroutine FI_CHANNEL_SPINUP(u, v, h1, h2)

  use TLAB_VARS,      only : isize_field, itime
  use TLAB_VARS,      only : channel_rot, nitera_spinup, spinuptime
  use TLAB_PROCS,     only : TLAB_WRITE_ASCII
  use TLAB_CONSTANTS, only : wfile

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u     
  TREAL, dimension(isize_field), intent(inout) :: v     
  TREAL, dimension(isize_field), intent(inout) :: h1    
  TREAL, dimension(isize_field), intent(inout) :: h2    
! -----------------------------------------------------------------------

  TINTEGER                                     :: ij, itime_sub
  character*32                                 :: str
  character*128                                :: line

! #######################################################################

  ! rotate
  if (itime <= spinuptime + nitera_spinup) then
    if (itime > itime_sub) then
      write(str,*) itime; line = 'Rotating turbulent channel flow at iteration step '//trim(adjustl(str))//'.'
      call TLAB_WRITE_ASCII(wfile,line)
      itime_sub = itime
    end if
    !
    do ij = 1, isize_field
      h1(ij) = h1(ij) - channel_rot * v(ij)
      h2(ij) = h2(ij) + channel_rot * u(ij)
    end do
  end if

  ! other possible option, not tested!!
  ! ip   = 1
  ! um   = exp(-0.5 * (g(2)%nodes - 1)**C_2_R) - exp(-0.5)
  ! DO ij = 1, isize_field
  !   h1(ip:ip+jmax-1) = h1(ip:ip+jmax-1) - spinup_rot * um(:)
  !   ! h1(ij) = h1(ij) - 0.12 * v(ij)
  !   h2(ip:ip+jmax-1) = h2(ip:ip+jmax-1) + spinup_rot * um(:)
  !   ! h2(ij) = h2(ij) + 0.12 * u(ij)
  !   ip = ip + jmax
  ! ENDDO

  return
end subroutine FI_CHANNEL_SPINUP
!########################################################################
subroutine FI_CHANNEL_UBULK(u,wrk1d,wrk3d)

  use TLAB_VARS, only : imax,jmax,kmax, isize_field
  use TLAB_VARS, only : g, area
  use TLAB_VARS, only : ubulk

  implicit none

#include "integers.h"

  TREAL, dimension(isize_field), intent(inout) :: u
  TREAL, dimension(jmax),        intent(inout) :: wrk1d
  TREAL, dimension(isize_field), intent(inout) :: wrk3d

! -----------------------------------------------------------------------
  TREAL, dimension(jmax)                       :: uxz
  TREAL                                        :: SIMPSON_NU

! #######################################################################

  call AVG_IK_V(imax,jmax,kmax, jmax, u, g(1)%jac,g(3)%jac, uxz, wrk1d, area)

  ubulk = (C_1_R / g(2)%nodes(g(2)%size)) * SIMPSON_NU(jmax, uxz, g(2)%nodes)

  return
end subroutine FI_CHANNEL_UBULK
!########################################################################
subroutine FI_CHANNEL_INITIALIZE()

  use TLAB_VARS,      only : g, qbg, ubulk_parabolic
  use TLAB_VARS,      only : itime, spinuptime, channel_rot, nitera_spinup
  use TLAB_CONSTANTS, only : efile, wfile, lfile
  use TLAB_PROCS

  implicit none

#include "integers.h"

  character*32  :: str
  character*128 :: line
  character*10  :: clock(2)

! -----------------------------------------------------------------------
  ! initialize bulk velocity
  ! laminar streamwise bulk velocity, make sure that the initial parabolic 
  ! velocity profile is correct (u(y=0)=u(y=Ly)=0)

  if(qbg(1)%type == PROFILE_PARABOLIC .or. qbg(1)%type == PROFILE_PARABOLIC_XCOMPACT3D) then
    ubulk_parabolic = (C_1_R / g(2)%nodes(g(2)%size)) * (C_4_R/C_3_R) * qbg(1)%delta 
  else 
    call TLAB_WRITE_ASCII(efile, 'Analytical ubulk cannot be computed. Check initial velocity profile.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
  end if

  ! initialize turbulent channel flow rotation
  if (nitera_spinup > i0) then
    CALL TLAB_WRITE_ASCII(wfile, '############################################################################')
    CALL DATE_AND_TIME(clock(1),clock(2))
    line='Initializing channel spinup rotation on '//TRIM(ADJUSTL(clock(1)(1:8)))//' at '//TRIM(ADJUSTL(clock(2)))
    CALL TLAB_WRITE_ASCII(wfile,line)
    call TLAB_WRITE_ASCII(wfile, 'Turbulent channel flow rotation is turned on in order to reduce spinup time.')
    write(str,*) channel_rot; line = 'Channel rotation is set to '//trim(adjustl(str))//'.'
    call TLAB_WRITE_ASCII(wfile, line)
    spinuptime = itime
  end if

  return
end subroutine FI_CHANNEL_INITIALIZE
!########################################################################
subroutine CHANNEL_RESCALE_GRID(imax,jmax,kmax, scalex,scaley,scalez, x,y,z, area)

  use TLAB_VARS, ONLY : isize_wrk1d, g

  implicit none

  TINTEGER,                   intent(inout) :: imax, jmax, kmax
  TREAL,                      intent(inout) :: scalex, scaley, scalez
  TREAL,     dimension(imax), intent(inout) :: x
  TREAL,     dimension(jmax), intent(inout) :: y
  TREAL,     dimension(kmax), intent(inout) :: z
  TREAL,     optional,        intent(inout) :: area

  ! -----------------------------------------------------------------------
  
  TREAL,     dimension(:),    allocatable   :: work1,work2
  TREAL                                     :: scaley_old, scaley_new
  TREAL                                     :: dxmx, dxmn, axmx, axmn
  TINTEGER                                  :: idir, n
  character*32                              :: sfile, name
  ! #######################################################################

  ! rescale grid with delta==1 (channel half height, channel height == scaley_new === 2)  
  scaley_new = C_2_R
  scaley_old = scaley

  ! x - nodes and scale
  x      = x      / scaley_old
  scalex = scalex / scaley_old

  ! y - nodes and scale
  y      = (y / scaley_old) * scaley_new
  scaley =                    scaley_new  ! udate scaley

  ! z - nodes and scale
  z      = z      / scaley_old
  scalez = scalez / scaley_old
  
  ! new area
  if ( present(area) ) then
    area = scalex
    if ( kmax > 1 ) area = area * scalez ! 3d case
  end if

  ! #######################################################################
  name = 'grid_rescale'

  ! write out rescaled grid_new file
  call IO_WRITE_GRID(name, imax,jmax,kmax, scalex,scaley,scalez, x,y,z)
  
  ! write out new grid_new.sts file
    ! g%nodes pointer is not asigned yet (in fdm_initialize)
    ! g%size and g%name is the same (old and new)
    ! g%scale is updated (in this subroutine)
  sfile    = trim(adjustl(name))//'.sts'
  !
  allocate(work1(isize_wrk1d))
  allocate(work2(isize_wrk1d))
  !
  open(20,file=sfile)
  do idir = 1,3
    write(20,3000) '['//trim(adjustl(g(idir)%name))//'-direction]'
    if ( g(idir)%size > 1 .and. g(idir)%name == 'x' ) then
      work1(2) = x(2) - x(1)
      do n = 3,g(idir)%size
        work1(n) = x(n) - x(n-1)
        work2(n) = work1(n)/work1(n-1)
      end do
      dxmx = MAXVAL(work1(2:g(idir)%size)); dxmn = MINVAL(work1(2:g(idir)%size))
      axmx = MAXVAL(work2(3:g(idir)%size)); axmn = MINVAL(work2(3:g(idir)%size))
      write(20,2000) 'number of points .......: ',g(idir)%size
      write(20,1000) 'origin .................: ',x(1)
      write(20,1000) 'end point ..............: ',x(g(idir)%size)
      write(20,1000) 'scale ..................: ',g(idir)%scale
      write(20,1000) 'minimum step ...........: ',dxmn
      write(20,1000) 'maximum step ...........: ',dxmx
      write(20,1000) 'minimum stretching .....: ',axmn
      write(20,1000) 'maximum stretching .....: ',axmx
    else if ( g(idir)%size > 1 .and. g(idir)%name == 'y' ) then
      work1(2) = y(2) - y(1)
      do n = 3,g(idir)%size
        work1(n) = y(n) - y(n-1)
        work2(n) = work1(n)/work1(n-1)
      end do
      dxmx = MAXVAL(work1(2:g(idir)%size)); dxmn = MINVAL(work1(2:g(idir)%size))
      axmx = MAXVAL(work2(3:g(idir)%size)); axmn = MINVAL(work2(3:g(idir)%size))
      write(20,2000) 'number of points .......: ',g(idir)%size
      write(20,1000) 'origin .................: ',y(1)
      write(20,1000) 'end point ..............: ',y(g(idir)%size)
      write(20,1000) 'scale ..................: ',g(idir)%scale
      write(20,1000) 'minimum step ...........: ',dxmn
      write(20,1000) 'maximum step ...........: ',dxmx
      write(20,1000) 'minimum stretching .....: ',axmn
      write(20,1000) 'maximum stretching .....: ',axmx
    else if ( g(idir)%size > 1 .and. g(idir)%name == 'z') then
      work1(2) = z(2) - z(1)
      do n = 3,g(idir)%size
        work1(n) = z(n) - z(n-1)
        work2(n) = work1(n)/work1(n-1)
      end do
      dxmx = MAXVAL(work1(2:g(idir)%size)); dxmn = MINVAL(work1(2:g(idir)%size))
      axmx = MAXVAL(work2(3:g(idir)%size)); axmn = MINVAL(work2(3:g(idir)%size))
      write(20,2000) 'number of points .......: ',g(idir)%size
      write(20,1000) 'origin .................: ',z(1)
      write(20,1000) 'end point ..............: ',z(g(idir)%size)
      write(20,1000) 'scale ..................: ',g(idir)%scale
      write(20,1000) 'minimum step ...........: ',dxmn
      write(20,1000) 'maximum step ...........: ',dxmx
      write(20,1000) 'minimum stretching .....: ',axmn
      write(20,1000) 'maximum stretching .....: ',axmx
    else
      write(20,'(a7)') '2D case'
    end if
  end do
  close(20)
  1000 format(a25,e12.5)
  2000 format(a25,i5)
  3000 format(a13)
  
  return
end subroutine CHANNEL_RESCALE_GRID
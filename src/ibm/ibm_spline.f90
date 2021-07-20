#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  
!#  
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################

subroutine IBM_SPLINE_X(u, u_mod, nlines, g)
  
  use DNS_IBM
  use DNS_GLOBAL, only: isize_field
  use DNS_TYPES,  only: grid_dt

  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
#endif

type(grid_dt),                  intent(in) :: g

TREAL, dimension(isize_field), intent(in)      :: u 
TREAL, dimension(isize_field), intent(out)     :: u_mod 

TINTEGER, intent(in) :: nlines

  
  ! ================================================================== !

  ! nothing implemented yet
  u_mod(:) = u(:)

  return
end subroutine IBM_SPLINE_X

!########################################################################

subroutine IBM_SPLINE_Y(u, u_mod, nlines, g)
  
  use DNS_IBM
  use DNS_GLOBAL, only: isize_field
  use DNS_TYPES,  only: grid_dt

  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
#endif

  type(grid_dt),                  intent(in) :: g

  TREAL, dimension(isize_field), intent(in)      :: u 
  TREAL, dimension(isize_field), intent(out)     :: u_mod 
  
  TINTEGER, intent(in) :: nlines
  ! ================================================================== !

  ! nothing implemented yet
  u_mod(:) = u(:)

  return
end subroutine IBM_SPLINE_Y

!########################################################################

subroutine IBM_SPLINE_Z(u, u_mod, nlines, g)
  
  use DNS_IBM,       only: nobk, nobk_b, nobk_e
  use DNS_IBM,       only: nflu, wrk_ibm, iwrk_ibm
  use DNS_IBM,       only: xa, xb, ya, yb
  use DNS_GLOBAL,    only: isize_field
  use DNS_CONSTANTS, only: efile
  use DNS_TYPES,     only: grid_dt

! MPI just for debugging
#ifdef USE_MPI
  use DNS_MPI,       only: ims_pro
#endif  
   
  implicit none
  
#include "integers.h"

! MPI just for debugging
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
#else
  TINTEGER, parameter :: ims_pro=0, ims_npro=1
#endif
  
  TREAL, dimension(isize_field), intent(in)  :: u 
  TREAL, dimension(isize_field), intent(out) :: u_mod 
  type(grid_dt),                 intent(in)  :: g
  TINTEGER,                      intent(in)  :: nlines 

  TINTEGER                                   :: l, ij, ip, ib, iob, ip_mod
  logical                                    :: splines
  character, dimension(128)                  :: line

  ! ================================================================== !
  ! debug
  if (ims_pro == 0) write(*,*) '========================================================='
  ! ================================================================== !

  ! modify u field with splines in solid region (here: equally spaced grid in z) 

  ! index convention on k-axis
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  splines = .true. ! 1. case doesn't need splines
  u_mod   = u

  do ij = 1, nlines            ! index of ij-plane, loop over plane and check for objects in each line
    if(nobk(ij) /= i0) then ! if k-line contains immersed object(s) --yes-->  spline interpolation
      ip = i0
      do iob = 1, nobk(ij)  ! loop over immersed object(s)

        ! select different cases [1...4] of immersed objects
        if(nobk_b(ip+ij) == i1) then
        ! ================================================================== !
          if(nobk_e(ip+ij) == g%size) then
            ! 1. case: object over full extend of z-line
            
            ! debug
            if (ims_pro == 0) write(*,*) 'IBM_SPLINE_Z Case 1'
            
            ! do nothing
            splines = .false.

          else if((nobk_e(ip+ij) <= (g%size - nflu)) .eqv. g%periodic) then
            ! 2. case: object is semi-immersed (z-direction needs to be periodic)
            
            ! debug
            if (ims_pro == 0) write(*,*) 'IBM_SPLINE_Z Case 2'
            
            ! build vectors for spline generation
            splines = .false. ! not implemented yet

          else
            ! not implemented yet, break
            ! either not enough fluidpoints or no periodic BC in z
            write(line, *) 'IBM_SPLINE_Z case not implemented yet'
            call IO_WRITE_ASCII(efile, line)
            call DNS_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        ! ================================================================== !
        else if(nobk_b(ip+ij) >= (nflu+1)) then
          if(nobk_e(ip+ij) <= (g%size - nflu)) then 
            ! 3. case: object is fully immersed

            ! debug
            if (ims_pro == 0) write(*,*) 'IBM_SPLINE_Z Case 3'

            ! build vectors for spline generation
            call IBM_SPLINE_VECTOR_3XZ(u, xa, ya, xb, ib, nobk_b(ip+ij), nobk_e(ip+ij), nlines, ij) 
            
          else if((nobk_e(ip+ij) == g%size) .eqv. g%periodic) then
            ! 4. case: object is semi-immersed (z-direction needs to be periodic)
            
            ! debug
            if (ims_pro == 0) write(*,*) 'IBM_SPLINE_Z Case 4'

            ! build vectors for spline generation
            splines = .false. ! not implemented yet
            
          else
            ! not implemented yet, break
            ! either not enough fluidpoints or no periodic BC in z
            write(line, *) 'IBM_SPLINE_Z case not implemented yet'
            call IO_WRITE_ASCII(efile, line)
            call DNS_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        else
          ! not implemented yet, break
          write(line, *) 'IBM_SPLINE_Z case not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_IBM_SPLINE)
        end if
        ! ================================================================== !

        ! spline interpolation and fill gap in u_ibm
        if (splines) then
          ip_mod = i0 
          call IBM_SPLINE(xa, ya, ib, xb(1:ib), wrk_ibm, iwrk_ibm, yb(1:ib))     
          do l = 1, ib
            u_mod(nobk_b(ip+ij) + l * nlines)  = yb(l)
          if (ims_pro == 0) write(*,*) 'u_mod', u_mod(nobk_b(ip+ij) + l * nlines)
          end do
        end if

        ! ================================================================== !
        ! debug
        if (ims_pro == 0) write(*,*) 'xa', xa
        if (ims_pro == 0) write(*,*) 'ya', ya
        if (ims_pro == 0) write(*,*) 'xb', xb(1:ib)
        if (ims_pro == 0) write(*,*) 'yb', yb(1:ib)
        if (ims_pro == 0) write(*,*) '=================='
        ! ================================================================== !

        ip = ip + nlines
      end do
    end if
  end do

  if (ims_pro == 0) write(*,*) '========================================================='
  
  return
end subroutine IBM_SPLINE_Z

!########################################################################

subroutine IBM_SPLINE(xa, ya, ib, xb, wrk_ibm, iwrk_ibm, yb)
  
  use DNS_IBM,       only: wrk_ibm_size, iwrk_ibm_size, nest, nsp, kspl
  use DNS_CONSTANTS, only: efile
   
  implicit none
  
#include "integers.h"

  TREAL,    dimension(nsp),            intent(in)    :: xa, ya 
  TINTEGER,                            intent(in)    :: ib
  TREAL,    dimension(ib),             intent(in)    :: xb
  TREAL,    dimension(wrk_ibm_size),   intent(inout) :: wrk_ibm
  TINTEGER, dimension(iwrk_ibm_size),  intent(inout) :: iwrk_ibm
  TREAL,    dimension(ib),             intent(out)   :: yb 

  TREAL                                              :: xstart, xend, s, fp
  TINTEGER                                           :: iopt, n, l, ier
  TINTEGER                                           :: ip1, ip2, ip3, ip4

  character, dimension(128)                          :: line
  
  ! ================================================================== !
  ! spline function parameter
  iopt = i0    ! (iopt=0 or 1) smoothing spline, weighted least-squares spline (iopt=-1)
  s    = C_0_R ! control the tradeoff between closeness of fit and smoothness

  ! set interval for spline approximation
  xstart = xa(1);  xend = xa(nsp)
 
  ! define working arrays and their relative positions
  ip1  = 1          ! w(nsp)
  ip2  = ip1 + nsp  ! t(nest)
  ip3  = ip2 + nest ! c(nest)
  ip4  = ip3 + nest ! wrk(nsp*(kspl+1)+nest*(7+3*kspl))

  ! weights of data points w(nsp)
  do l = 1, nsp
    wrk_ibm(l) = C_1_R ! here: all weights are equal
  end do

  ! ================================================================== !
    ! evaluation of spline function [curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)]
    !s    : (in case iopt>=0) s must specify the smoothing factor. 
    !       if s=0 then interpolation spline, data points coincident with spline points
    !t    : array,length n, which contains the position of the knots.
    !n    : integer, giving the total number of knots of s(x). [output]
    !c    : array,length n, which contains the b-spline coefficients.
    !k    : integer, giving the degree of s(x).
    !x    : array,length m, which contains the points where s(x) must
    !fp   : contains the weighted sum of squared residuals of the spline approximation [output]
    !ier  : ier contains a non-positive value on exit [-2,-1,0], if error ier=[1,2,3,10]

  !    curfit(iopt, m,   x,  y,  w,            xb,     xe,   k,    s, nest, n, &
  call curfit(iopt, nsp, xa, ya, wrk_ibm(ip1), xstart, xend, kspl, s, nest, n, & 
  !           t,            c,            fp, wrk,          lwrk,          iwrk,     ier)
              wrk_ibm(ip2), wrk_ibm(ip3), fp, wrk_ibm(ip4), iwrk_ibm_size, iwrk_ibm, ier)
  
  if ( (ier /= 0) .and. (ier /= -1) ) then
    write(line, *) 'INTERPOLATE_1D. Curfit error code = ', ier
    call IO_WRITE_ASCII(efile, line)
    call DNS_STOP(DNS_ERROR_CURFIT)
  end if

  ! ================================================================== !
    ! evaluation of the spline [call splev(t,n,c,k,x,y,m,ier)] function to evaluate a B-spline or its derivatives
    !###### input parameters:
    !t    : array,length n, which contains the position of the knots.
    !n    : integer, giving the total number of knots of s(x).
    !c    : array,length n, which contains the b-spline coefficients.
    !k    : integer, giving the degree of s(x).
    !x    : array,length m, which contains the points where s(x) must be evaluated.
    !m    : integer, giving the number of points where s(x) must be evaluated.
    !###### output parameter:
    !y    : array,length m, giving the value of s(x) at the different points.
    !ier  : error flag: ier = 0 : normal return, ier =10 : invalid input data (see restrictions)

  !    splev(t,            n, c,            k,    x,  y,  m,  ier)
  call splev(wrk_ibm(ip2), n, wrk_ibm(ip3), kspl, xb, yb, ib, ier)
  
  if (ier /= 0) then
    write(line, *) 'INTERPOLATE_1D. Splev error code = ', ier
    call IO_WRITE_ASCII(efile, line)
    call DNS_STOP(DNS_ERROR_CURFIT)
  end if

  return
end subroutine IBM_SPLINE

!########################################################################

subroutine IBM_SPLINE_VECTOR_3XZ(u, xa, ya, xb, ib, ip_il, ip_ir, nlines, plane) 

  use DNS_IBM,    only: nflu, wrk1d_ibm_size, nsp 
  use DNS_GLOBAL, only: isize_field
   
  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field),    intent(in)  :: u
  TREAL, dimension(nsp),            intent(out) :: xa
  TREAL, dimension(nsp),            intent(out) :: ya
  TREAL, dimension(wrk1d_ibm_size), intent(out) :: xb
  TINTEGER,                         intent(out) :: ib
  TINTEGER,                         intent(in)  :: ip_il, ip_ir, nlines, plane

  TINTEGER                                      :: ia, kflu, kint
  TINTEGER                                      :: ip_fl, iu_fl, iu_ir

  ! ================================================================== !

  ! index to remember current position in vectors
  ia = i1
  ib = i0

  ! needed indices
  ip_fl = ip_il - nflu                ! k-axis index of most left fluid point
  iu_fl = ip_fl * nlines + plane      ! u-index of most left fluid point
  iu_ir = ip_ir * nlines + plane      ! u-index of right interface  point

  ! build left half of xa, ya (from left to right)
  do kflu = 1, nflu
    xa(ia) = dble(ip_fl + (kflu - 1)) ! other option: grid node positions      
    ya(ia) =    u(iu_fl + (kflu - 1) * nlines)
    ia     = ia + 1
  end do

  ! set interfaces (left and right)
  xa(ia) = dble(ip_il)   
  ya(ia) = C_0_R       ! velocity is zero 
  ia     = ia + 1
  !
  xa(ia) = dble(ip_ir) 
  ya(ia) = C_0_R       ! velocity is zero 
  ia     = ia + 1

  ! build right half of xa, ya (from left to right)
  do kflu = 1, nflu
    xa(ia) = dble(ip_ir + kflu)
    ya(ia) =    u(iu_ir + kflu * nlines)
    ia     = ia + 1
  end do

  ! build gap vector where splines are evaluated (here: with interface points)
  do kint = ip_il, ip_ir
    ib = ib + 1
    xb(ib) = dble(kint)
  end do

  ! not needed indices --> just for debugging
  ! ip_fr = nobk_e(ip+plane) + nflu   ! k-axis index of most right fluid point
  ! iu_il = ip_il * nlines + plane)   ! u-index of left interface point
  ! iu_fr = ip_fr * nlines + plane    ! u-index of most right fluid point

  return
end subroutine IBM_SPLINE_VECTOR_3XZ
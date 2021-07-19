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
  use DNS_IBM,       only: nflu
  use DNS_IBM,       only: xa, xb, ya, yb
  use DNS_IBM,       only: wrk_ibm, iwrk_ibm
  use DNS_GLOBAL,    only: imax, jmax, kmax 
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

  type(grid_dt),                  intent(in) :: g

  TINTEGER,                       intent(in) :: nlines 
  TINTEGER                                   :: k, ij, ip, ia, ib, iob, kflu, kint
  TINTEGER                                   :: ip_fl, iu_fl, ip_il, ip_ir, iu_ir
  TINTEGER                                   :: nxy    ! nyz, nxz
  character, dimension(64)                   :: line

  ! ================================================================== !
  if (ims_pro == 0) write(*,*) '========================================================='

  ! modify u field with splines in solid region
  ! equally spaced grid in z 

  ! index convention on k-axis
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  nxy = nlines

  do ij = 1, 1!nxy            ! index of ij-plane, loop over plane and check for objects in each line
    if(nobk(ij) /= i0) then ! if k-line contains immersed object(s) --yes-->  spline interpolation
      ip = i0
      do iob = 1, nobk(ij)  ! loop over immersed object(s)
        ! build vectors for spline interpolation and  each single object 
        ! ================================================================== !
        ! index to remember current position in vectors
        ia = i1

        ! first interface node (ip_il = nobk_b(ip+ij) ) (build left half of vectors)     
        if(nobk_b(ip+ij) >= (nflu+1)) then ! obj is immersed
          
          ! k-axis indices
          ip_il = nobk_b(ip+ij)          ! k-axis index of left interface point
          ip_fl = ip_il - nflu   ! k-axis index of most left fluid point

          ! coresponding indices in u-field (interface iu_il = ip_il * nxy + ij)
          iu_fl = ip_fl * nxy + ij       ! u-index of most left fluid point

          ! loop over nflu-fluid points (from left to right)
          do kflu = 1, nflu
            xa(ia) = dble(ip_fl + (kflu - 1))     ! could be also done with actual grid nod positions      
            ya(ia) =    u(iu_fl + (kflu - 1) * nxy)
            ia     = ia + 1
          end do

          ! set left interface
          xa(ia) = dble(nobk_b(ip+ij))
          ya(ia) = C_0_R
          ia     = ia + 1
         
        ! not implemented yet  
        else if((nobk_b(ip+ij) < (nflu+1)) .and. ((nobk_b(ip+ij) > i0))) then
          ! periodic conditions for xa and ya needed!
          ! not implemented yet
          write(line, *) 'IBM_spline Z not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_IBM_SPLINE) 

        else ! if (nobk_b(ip+ij) == i0)
          ! object is half immersed or object exents over full length
          ! not implemented yet
          write(line, *) 'IBM_spline Z not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_IBM_SPLINE) 
        end if 

        ! ================================================================== !

        ! find second interface
        if(nobk_e(ip+ij) <= (g%size - nflu)) then ! obj is immersed
          
          ! k-axis indices
          ip_ir = nobk_e(ip+ij)        ! k-axis index of right interface  point
          ! ip_fr = nobk_e(ip+ij) + nflu   ! k-axis index of most right fluid point

          ! coresponding indices in u-field
          iu_ir = ip_ir * nxy + ij     ! u-index of right interface  point
          ! iu_fr = ip_fr * nxy + ij       ! u-index of most right fluid point

          ! set right interface
          xa(ia) = dble(ip_ir) 
          ya(ia) = C_0_R
          ia     = ia + 1

          ! loop over nflu fluid points, from left to right
          do kflu = 1, nflu
            xa(ia) = dble(ip_ir + kflu)
            ya(ia) =    u(iu_ir + kflu * nxy)
            ia     = ia + 1
          end do

          ! build gap vector where splines are evaluated, with interface points
          ib = i0
          do kint = ip_il, ip_ir
            ib = ib + 1
            xb(ib) = dble(kint)
          end do
          
        ! not implemented yet  
        else if((nobk_e(ip+ij) > (g%size - nflu)) .and. ((nobk_e(ip+ij) < g%size))) then
          ! periodic conditions for xa and ya needed!
          ! not implemented yet
          write(line, *) 'IBM_spline Z not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_IBM_SPLINE) 

        else ! if (nobk_e(ip+ij) == g%size)
          ! object is half immersed or object exents over full length
          ! not implemented yet
          write(line, *) 'IBM_spline Z not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_IBM_SPLINE) 
        end if 

        ! debug
        if (ims_pro == 0) write(*,*) 'xa', xa
        if (ims_pro == 0) write(*,*) 'ya', ya
        if (ims_pro == 0) write(*,*) 'xb', xb(1:ib)
        ! if (ims_pro == 0) write(*,*) 'ib', ib

        ! ================================================================== !
        ! spline interpolation
        call IBM_SPLINE(xa, ya, ib, xb(1:ib), wrk_ibm, iwrk_ibm, yb(1:ib))

        ! ================================================================== !
        ! fill gap in u_ibm 


        ! ================================================================== !
        ip = ip + nxy
      end do
    end if
  end do

  if (ims_pro == 0) write(*,*) '========================================================='
  
  return
end subroutine IBM_SPLINE_Z

!########################################################################

subroutine IBM_SPLINE(xa, ya, ib, xb, wrk_ibm, iwrk_ibm, yb)
  
  use DNS_IBM,       only: wrk_ibm_size, iwrk_ibm_size, wrk1d_ibm_size
  use DNS_IBM,       only: nest, nsp
  use DNS_IBM,       only: nflu, kspl
  use DNS_CONSTANTS, only: efile
   
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

  TREAL,    dimension(nsp),            intent(in)    :: xa, ya 
  TINTEGER,                            intent(in)    :: ib
  TREAL,    dimension(ib),             intent(in)    :: xb
  TREAL,    dimension(wrk_ibm_size),   intent(inout) :: wrk_ibm
  TINTEGER, dimension(iwrk_ibm_size),  intent(inout) :: iwrk_ibm

  TREAL,    dimension(ib),             intent(out)   :: yb 

  TREAL                                              :: xstart, xend, s, fp

  TINTEGER                                           :: iopt, n, l, ier
  TINTEGER                                           :: ip1, ip2, ip3, ip4

  character, dimension(64)                           :: line
  
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
  
  ! ================================================================== !
  ! debug
  if (ims_pro == 0) write(*,*) 'yb', yb
  if (ims_pro == 0) write(*,*) '=================='

  ! ================================================================== !

  return
end subroutine IBM_SPLINE

!########################################################################
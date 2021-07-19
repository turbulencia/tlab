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
  use DNS_GLOBAL,    only: imax, jmax, kmax 
  use DNS_GLOBAL,    only: isize_field
  use DNS_CONSTANTS, only: efile
  use DNS_TYPES,     only: grid_dt


#ifdef USE_MPI
  use DNS_MPI,       only: ims_pro ! for debugging
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

  ! vectors for spline interpolation
  TREAL, dimension(2*nflu+2) :: xa, ya
  TREAL, dimension(2*nflu+2) :: xb, yb

  type(grid_dt),                  intent(in) :: g

  TINTEGER, intent(in) :: nlines ! nlines == nxy
  TINTEGER                                   :: k, ij, ip, ia, iob, kflu
  TINTEGER                                   :: ip_fl, iu_fl, ip_ir, iu_ir
  TINTEGER                                   :: nyz, nxz, nxy
  character, dimension(64)                   :: line

  ! ================================================================== !
  if (ims_pro == 0) write(*,*) '========================================================='

  ! modify u field with splines in solid region
  ! equally spaced grid in z 

  ! index convention on k-axis
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  nxy = nlines

  do ij = 1, nxy            ! index of ij-plane, loop over plane and check for objects in each line
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
          ip_fl = nobk_b(ip+ij) - nflu   ! k-axis index of most left fluid point

          ! coresponding indices in u-field (interface iu_il = ip_il * nxy + ij)
          iu_fl = ip_fl * nxy + ij       ! u-index of most left fluid point

          ! loop over nflu-fluid points (from left to right)
          do kflu = 1, nflu
            xa(ia) =   ip_fl + (kflu - 1)     ! could be also done with actual grid nod positions      
            ya(ia) = u(iu_fl + (kflu - 1) * nxy)
            ia     = ia + 1
          end do

          ! set left interface
          xa(ia) = nobk_b(ip+ij)
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
          xa(ia) = ip_ir 
          ya(ia) = C_0_R
          ia     = ia + 1

          ! loop over nflu fluid points, from left to right
          do kflu = 1, nflu
            xa(ia) =   ip_ir + kflu
            ya(ia) = u(iu_ir + kflu * nxy)
            ia     = ia + 1
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
        if (ims_pro == 0) write(*,*) 'xa', xa
        if (ims_pro == 0) write(*,*) 'ya', ya

        ! ================================================================== !
        ! spline interpolation
        call IBM_SPLINE(xa, ya, xb, yb)

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

subroutine IBM_SPLINE(xa, ya, xb, yb) 
  
  use DNS_IBM, only: nflu, kspl
   
  implicit none
  
#include "integers.h"

  TREAL, dimension(2*nflu+2), intent(in)  :: xa, ya
  TREAL, dimension(2*nflu+2), intent(out) :: xb, yb 
  
  ! ================================================================== !


!   ! wrk_ibm
!         ! contains the following arrays for for curfit():
!           ! w(m)
!           ! t(nest)
!           ! c(nest)
!           ! wrk(m*(k+1)+nest*(7+3*k))
!           ! iwrk(lwrk=m*(k+1)+nest*(7+3*k))    
!             ! m            = 2 * (nflu + 1) ! number of data points m > k (k==kspl)
!             ! nest         = m + kspl + 1
!             ! nest         = 2 * (nflu + 1) + kspl + 1
!             ! wkk_ibm_size = m + 2*nest + 2*(m*(k+1)+nest*(7+3*k))
!     wrk_ibm_size = (2 * (nflu + 1)) + 2*(2 * (nflu + 1) + kspl + 1) + &
!     2*((2 * (nflu + 1))*(kspl+1)+(2 * (nflu + 1) + kspl + 1)*(7+3*kspl))

  

!   ! #######################################################################
!   iopt = 0
!   xb = x_org(1); xe = x_org(imax)                  ! just used in the non-periodic case
!   imax1 = imax+1; x_org(imax1) = x_org(1) + scalex ! the periodic case
!   kx = 5                                           ! order of the splines interpolation
!   s = C_0_R
!   nest = imax1 + 2*kx                              ! the periodic case requires maximum sizes
!   lwrk = imax1*(kx+1)+nest*(8+5*kx)

! ! Working array relative possitions
!   ip1 = 1           ! w
!   ip2 = ip1 + imax1 ! t
!   ip3 = ip2 + nest  ! c
!   ip4 = ip3 + nest  ! wrk
!   ip5 = ip4 + lwrk  ! iwkr
!   ip6 = ip5 + nest  ! returned value

!   IF ( isize_wrk .LT. ip6 ) THEN
!      CALL IO_WRITE_ASCII(efile, 'INTERPOLATE_1D. Temporary Array not large enough')
!      CALL DNS_STOP(DNS_ERROR_CURFIT)
!   ENDIF





! ! ###################################################################
! ! spline function parameter
!   iopt = 0     ! (iopt=0 or 1) smoothing spline, weighted least-squares spline (iopt=-1)
!   s    = C_0_R ! control the tradeoff between closeness of fit and smoothness


!     ! evaluation of spline function [curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)]
!     !s    : (in case iopt>=0) s must specify the smoothing factor. 
!     !       if s=0 then interpolation spline, data points coincident with spline points
!     !t    : array,length n, which contains the position of the knots.
!     !n    : integer, giving the total number of knots of s(x).
!     !c    : array,length n, which contains the b-spline coefficients.
!     !k    : integer, giving the degree of s(x).
!     !x    : array,length m, which contains the points where s(x) must
!     !fp   : contains the weighted sum of squared residuals of the spline approximation
!     !ier  : ier contains a non-positive value on exit [-2,-1,0], if error ier=[1,2,3,10]
!   call curfit(iopt,imax,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  
!   if ( (ier /= 0) .and. (ier /= -1) ) then
!     write(line, *) 'INTERPOLATE_1D. Curfit error code = ', ier
!     call IO_WRITE_ASCII(efile, line)
!     call DNS_STOP(DNS_ERROR_CURFIT)
!   end if


!     ! evaluation of the spline [call splev(t,n,c,k,x,y,m,ier)] function to evaluate a B-spline or its derivatives
!     !###### input parameters:
!     !t    : array,length n, which contains the position of the knots.
!     !n    : integer, giving the total number of knots of s(x).
!     !c    : array,length n, which contains the b-spline coefficients.
!     !k    : integer, giving the degree of s(x).
!     !x    : array,length m, which contains the points where s(x) must be evaluated.
!     !m    : integer, giving the number of points where s(x) must be evaluated.
!     !###### output parameter:
!     !y    : array,length m, giving the value of s(x) at the different points.
!     !ier  : error flag: ier = 0 : normal return, ier =10 : invalid input data (see restrictions)
!   call splev(t,n,c,k,x_new,y_sp,imax_new,ier)
  
!   if (ier /= 0) then
!     write(line, *) 'INTERPOLATE_1D. Splev error code = ', ier
!     call IO_WRITE_ASCII(efile, line)
!     call DNS_STOP(DNS_ERROR_CURFIT)
!   end if
  





















  return
end subroutine IBM_SPLINE

!########################################################################
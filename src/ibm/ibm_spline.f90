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

subroutine IBM_SPLINE_Y(fld, fld_mod, nlines, g)
  
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

  TREAL, dimension(isize_field), intent(in)      :: fld 
  TREAL, dimension(isize_field), intent(out)     :: fld_mod 
  
  TINTEGER, intent(in) :: nlines
  ! ================================================================== !

  ! nothing implemented yet
  fld_mod(:) = fld(:)

  return
end subroutine IBM_SPLINE_Y

!########################################################################

subroutine IBM_SPLINE_XZ(fld, fld_mod, g, nlines, isize_nob, isize_nob_be, nob, nob_b, nob_e)

  use DNS_IBM,       only: xa, xb, ya, yb, nflu
  use DNS_GLOBAL,    only: isize_field
  use DNS_CONSTANTS, only: efile
  use DNS_TYPES,     only: grid_dt

  ! MPI just for debugging
  use DNS_GLOBAL,    only: imax, jmax, kmax ! debug
#ifdef USE_MPI
  use DNS_MPI,       only: ims_pro
  use DNS_MPI,       only: ims_size_i, ims_size_j, ims_size_k    
  use DNS_MPI,       only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i 
  use DNS_MPI,       only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use DNS_MPI,       only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
#endif  
   
  implicit none
  
#include "integers.h"

! MPI just for debugging
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
  TINTEGER, parameter                            :: idi = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter                            :: idj = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter                            :: idk = DNS_MPI_K_PARTIAL 
#else
  TINTEGER, parameter                            :: ims_pro=0, ims_npro=1
#endif
  
  TREAL,    dimension(isize_field),  intent(in)  :: fld 
  TREAL,    dimension(isize_field),  intent(out) :: fld_mod 
  type(grid_dt),                     intent(in)  :: g
  TINTEGER,                          intent(in)  :: nlines, isize_nob, isize_nob_be
  TINTEGER, dimension(isize_nob),    intent(in)  :: nob
  TINTEGER, dimension(isize_nob_be), intent(in)  :: nob_b, nob_e

  TINTEGER                                       :: l, ii, ip, ib, iob, iu_il
  logical                                        :: splines
  character, dimension(128)                      :: line

  ! debug
  TREAL, dimension(isize_field)                  :: wrk3d      ! debug
  TREAL, dimension(isize_field)                  :: fld_mod_tr ! debug

  ! ================================================================== !
  ! debug
  if (ims_pro == 0) write(*,*) '========================================================='
  ! ================================================================== !

  ! modify field with splines in solid region 

  ! index convention on contiguous lines
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  splines = .true. ! 1. case doesn't need splines
  fld_mod = fld

  ! index ii (dummy index; for x: ii == jk, for y: ii == ik, for z: ii == ij)

  do ii = 1, nlines        ! index of ii-plane, loop over plane and check for objects in each line
    if(nob(ii) /= i0) then ! if line contains immersed object(s) --yes-->  spline interpolation
      ip = i0
      do iob = 1, nob(ii)  ! loop over immersed object(s)

        ! select different cases [1...4] of immersed objects
        if(nob_b(ip+ii) == i1) then
        ! ================================================================== !
          if(nob_e(ip+ii) == g%size) then
            ! 1. case: object over full extend of line
            
            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 1'
            
            ! do nothing
            splines = .false.

          else if((nob_e(ip+ii) <= (g%size - nflu)) .eqv. g%periodic) then
            ! 2. case: object is semi-immersed
            
            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 2'
            
            ! build vectors for spline generation
            splines = .false. ! not implemented yet

          else if((nob_e(ip+ii) <= (g%size - nflu)) .neqv. g%periodic) then
            ! 3. case: object is semi-immersed
            
            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 3'
            
            ! build vectors for spline generation
            splines = .false. ! not implemented yet

          else
            write(line, *) 'IBM_SPLINE not enough fluid points right of the right interface'
            call IO_WRITE_ASCII(efile, line)
            call DNS_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        ! ================================================================== !
        else if(nob_b(ip+ii) >= (nflu+1)) then
          if(nob_e(ip+ii) <= (g%size - nflu)) then 
            ! 4. case: object is fully immersed

            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 4'

            ! build vectors for spline generation
            call IBM_SPLINE_VECTOR_4(fld, xa, ya, xb, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            
          else if((nob_e(ip+ii) == g%size) .eqv. g%periodic) then
            ! 5. case: object is semi-immersed
            
            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 5'

            ! build vectors for spline generation
            splines = .false. ! not implemented yet
            
          else if((nob_e(ip+ii) == g%size) .neqv. g%periodic) then
            ! 6. case: object is semi-immersed
            
            ! debug
            ! if (ims_pro == 0) write(*,*) 'IBM_SPLINE_XZ Case 6'

            ! build vectors for spline generation
            splines = .false. ! not implemented yet

          else
            write(line, *) 'IBM_SPLINE not enough fluid points left of the left interface'
            call IO_WRITE_ASCII(efile, line)
            call DNS_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        else
          write(line, *) 'IBM_SPLINE this case is not implemented yet'
          call IO_WRITE_ASCII(efile, line)
          call DNS_STOP(DNS_ERROR_NOTIMPL)
        end if

        ! ================================================================== !

        ! spline interpolation and fill gap in fld_ibm
        if (splines) then

          ! generate splines
          call IBM_SPLINE(xa, ya, ib, xb(1:ib), yb(1:ib)) 
          
          ! fld index of left interface
          iu_il = (nob_b(ip+ii) - 1) * nlines + ii     
          
          ! replace splines in solid gaps
          do l = 1, ib
            fld_mod(iu_il + (l-1) * nlines) = yb(l)
          end do
          
          ! ! ================================================================== !
          ! ! debug
          ! if (ims_pro == 0) write(*,*) 'xa', xa
          ! if (ims_pro == 0) write(*,*) 'ya', ya
          ! if (ims_pro == 0) write(*,*) 'xb', xb(1:ib)
          ! if (ims_pro == 0) write(*,*) 'yb', yb(1:ib)
          ! if (ims_pro == 0) write(*,*) '=================='
          ! ! ================================================================== !
        
        end if
        ip = ip + nlines
      end do
    end if
  end do

  if (g%name == 'z') then
    ! ================================================================== !
    ! debug
    if (ims_pro == 0) write(*,*) '========================================================='

    ! write out fld_mod for debugging
#ifdef USE_MPI
    if ( ims_npro_k > 1 ) then
      call DNS_MPI_TRPB_K(fld_mod, fld_mod_tr, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    endif
    call DNS_WRITE_FIELDS('fld_mod', i2, imax,jmax,kmax, i1, imax*jmax*kmax, fld_mod_tr, wrk3d)
#else
    fld_mod_tr = fld_mod
    call DNS_WRITE_FIELDS('fld_mod', i2, imax,jmax,kmax, i1, imax*jmax*kmax, fld_mod_tr, wrk3d)
#endif

    ! stop after writing field
    call DNS_STOP(DNS_ERROR_IBM_SPLINE)
  ! ================================================================== !
  end if

  return
end subroutine IBM_SPLINE_XZ

!########################################################################

subroutine IBM_SPLINE(xa, ya, ib, xb, yb)
  
  use DNS_IBM,       only: isize_iwrk_ibm, nest, nsp, kspl
  use DNS_IBM,       only: wrk_ibm, iwrk_ibm
  use DNS_CONSTANTS, only: efile
   
  implicit none
  
#include "integers.h"

  TREAL,    dimension(nsp), intent(in)  :: xa, ya 
  TINTEGER,                 intent(in)  :: ib
  TREAL,    dimension(ib),  intent(in)  :: xb
  TREAL,    dimension(ib),  intent(out) :: yb 

  TREAL                                 :: xstart, xend, s, fp
  TINTEGER                              :: iopt, n, l, ier
  TINTEGER                              :: ip1, ip2, ip3, ip4

  character, dimension(128)             :: line
  
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
  !           t,            c,            fp, wrk,          lwrk,           iwrk,     ier)
              wrk_ibm(ip2), wrk_ibm(ip3), fp, wrk_ibm(ip4), isize_iwrk_ibm, iwrk_ibm, ier)
  
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

  ! force yb at interface to zero (without this yb approx 10^-16 at interface)
  yb(1)  = C_0_R
  yb(ib) = C_0_R
  
  if (ier /= 0) then
    write(line, *) 'INTERPOLATE_1D. Splev error code = ', ier
    call IO_WRITE_ASCII(efile, line)
    call DNS_STOP(DNS_ERROR_CURFIT)
  end if

  return
end subroutine IBM_SPLINE

!########################################################################

subroutine IBM_SPLINE_VECTOR_4(fld, xa, ya, xb, ib, ip_il, ip_ir, nlines, plane) 

  use DNS_IBM,    only: nflu, isize_wrk1d_ibm, nsp 
  use DNS_GLOBAL, only: isize_field
   
  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field),     intent(in)  :: fld
  TREAL, dimension(nsp),             intent(out) :: xa
  TREAL, dimension(nsp),             intent(out) :: ya
  TREAL, dimension(isize_wrk1d_ibm), intent(out) :: xb
  TINTEGER,                          intent(out) :: ib
  TINTEGER,                          intent(in)  :: ip_il, ip_ir, nlines, plane

  TINTEGER                                       :: ia, kflu, kint
  TINTEGER                                       :: ip_fl, iu_fl, iu_ir

  ! ================================================================== !

  ! index to remember current position in vectors
  ia = i1
  ib = i0

  ! needed indices
  ip_fl = ip_il - nflu                      ! k-axis index of most left fluid point
  iu_fl = (ip_fl - 1) * nlines + plane      ! fld-index of most left fluid point
  iu_ir = (ip_ir - 1) * nlines + plane      ! fld-index of right interface  point

  ! build left half of xa, ya (from left to right)
  do kflu = 1, nflu
    xa(ia) = dble(ip_fl + (kflu - 1)) ! other option: grid node positions      
    ya(ia) =  fld(iu_fl + (kflu - 1) * nlines)
    ia     = ia + 1
  end do

  ! set interfaces (left and right)
  xa(ia) = dble(ip_il)   
  ya(ia) = C_0_R       
  ia     = ia + 1
  !
  xa(ia) = dble(ip_ir) 
  ya(ia) = C_0_R       
  ia     = ia + 1

  ! build right half of xa, ya (from left to right)
  do kflu = 1, nflu
    xa(ia) = dble(ip_ir + kflu)
    ya(ia) =  fld(iu_ir + kflu * nlines)
    ia     = ia + 1
  end do

  ! build gap vector where splines are evaluated (here: with interface points)
  do kint = ip_il, ip_ir
    ib = ib + 1
    xb(ib) = dble(kint)
  end do

  return
end subroutine IBM_SPLINE_VECTOR_4
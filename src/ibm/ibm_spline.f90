#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI 
#include "dns_const_mpi.h"  
#endif

!########################################################################
!# HISTORY / AUTHORS
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

subroutine IBM_SPLINE_XYZ(is, fld, fld_mod, g, nlines, isize_nob, isize_nob_be, nob, nob_b, nob_e, wrk3d)

  use DNS_IBM,        only: xa, xb, ya, yb, nflu, ibmscaljmin
  use TLAB_VARS,      only: isize_field
  use TLAB_CONSTANTS, only: efile
  use TLAB_TYPES,     only: grid_dt
  use TLAB_PROCS

  ! MPI just for debugging
#ifdef IBM_DEBUG
  use IO_FIELDS
  use TLAB_VARS,      only: imax, jmax, kmax
#ifdef USE_MPI 
  use MPI
  use TLAB_MPI_VARS,  only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use TLAB_MPI_VARS,  only: ims_npro_k, ims_pro
  use TLAB_MPI_PROCS
#endif
#endif

  implicit none
  
#include "integers.h"

  TINTEGER,                          intent(in)    :: is     ! scalar index; if 0, then velocity
  TREAL,    dimension(isize_field),  intent(in)    :: fld 
  TREAL,    dimension(isize_field),  intent(out)   :: fld_mod 
  type(grid_dt),                     intent(in)    :: g
  TINTEGER,                          intent(in)    :: nlines, isize_nob, isize_nob_be
  TINTEGER, dimension(isize_nob),    intent(in)    :: nob
  TINTEGER, dimension(isize_nob_be), intent(in)    :: nob_b, nob_e
  TREAL,    dimension(isize_field),  intent(inout) :: wrk3d  ! for cubic splines subroutine

  TINTEGER                                         :: l, ii, ip, ia, ib, iob, iu_il
  logical                                          :: splines

! MPI just for debugging
#ifdef IBM_DEBUG
#ifdef USE_MPI 
  TINTEGER, parameter                              :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter                              :: idj = TLAB_MPI_J_PARTIAL 
  TINTEGER, parameter                              :: idk = TLAB_MPI_K_PARTIAL 
#else
  TINTEGER, parameter                              :: ims_pro=0, ims_npro=1
#endif
  ! debug
  TREAL, dimension(isize_field)                    :: fld_mod_tr ! debug
  TINTEGER                                         :: m
#endif
  ! ================================================================== !  
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) '========================================================='
#endif
  ! ================================================================== !

  ! modify field with splines in solid region 

  ! index convention on contiguous lines
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  splines = .true.  ! 1. case doesn't need splines
  fld_mod = fld

  ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)

  do ii = 1, nlines        ! index of ii-plane, loop over plane and check for objects in each line
    if(nob(ii) /= i0) then ! if line contains immersed object(s) --yes-->  spline interpolation
      ip = i0
      do iob = 1, nob(ii)  ! loop over immersed object(s)
        ! select different cases of immersed objects
        if(nob_b(ip+ii) == i1) then
        ! ================================================================== !
          if(nob_e(ip+ii) == g%size) then
            ! 1. case: object over full extend of line
            splines = .false. ! do nothing
            ! .............................................................. !
          else if((nob_e(ip+ii) <= (g%size - nflu)) .eqv. g%periodic) then
            ! 2. case: object is semi-immersed - periodic case (not implemented yet)
            call IBM_SPLINE_VECTOR(is, i2, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if((nob_e(ip+ii) <= (g%size - nflu)) .neqv. g%periodic) then ! in j-direction
            ! 3. case: object is semi-immersed - non-periodic case - lower boundary
            call IBM_SPLINE_VECTOR(is, i3, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else
            call TLAB_WRITE_ASCII(efile, 'IBM_SPLINE not enough fluid points right of the right interface.')
            call TLAB_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        ! ================================================================== !
        else if(nob_b(ip+ii) >= (nflu+1)) then
          if(nob_e(ip+ii) <= (g%size - nflu)) then 
            ! 4. case: object is fully immersed
            call IBM_SPLINE_VECTOR(is, i4, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if((nob_e(ip+ii) == g%size) .eqv. g%periodic) then
            ! 5. case: object is semi-immersed - periodic case (not implemented yet)
            call IBM_SPLINE_VECTOR(is, i5, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if((nob_e(ip+ii) == g%size) .neqv. g%periodic) then  ! in j-direction
            ! 6. case: object is semi-immersed - non-periodic case - upper boundary
            call IBM_SPLINE_VECTOR(is, i6, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else
            call TLAB_WRITE_ASCII(efile, 'IBM_SPLINE not enough fluid points left of the left interface.')
            call TLAB_STOP(DNS_ERROR_IBM_SPLINE)
          end if
          ! .............................................................. !
        else
          call TLAB_WRITE_ASCII(efile, 'IBM_SPLINE this case is not implemented yet.')
          call TLAB_STOP(DNS_ERROR_NOTIMPL)
        end if

        ! ================================================================== !

        ! spline interpolation and fill gap in fld_ibm
        if ( splines ) then

          ! generate splines
          call IBM_SPLINE( ia, ib, xa(1:ia), ya(1:ia), xb(1:ib), yb(1:ib)) ! old fitpack routines
          ! call CUBIC_SPLINE(ia, ib, xa(1:ia), ya(1:ia), xb(1:ib), yb(1:ib), wrk3d) 

          ! force yb at interface to physical BCs (without this yb approx delta 10^-16 wrong at interface)
          if ( is /= i0 ) then
            yb(1)  = ibmscaljmin(is)
            yb(ib) = ibmscaljmin(is)
          else
            yb(1)  = C_0_R
            yb(ib) = C_0_R
          endif
                  
          ! fld index of left interface
          iu_il = (nob_b(ip+ii) - 1) * nlines + ii     
          
          ! replace splines in solid gaps
          do l = 1, ib
            fld_mod(iu_il + (l - 1) * nlines) = yb(l)
          end do
          ! ================================================================== !
          ! debug
          ! if (ims_pro == 0) write(*,*) 'size(xa)', size(xa(1:ia))
          ! if (ims_pro == 0) write(*,*) 'size(ya)', size(ya(1:ia))
          ! if (ims_pro == 0) write(*,*) 'size(xb)', size(xb(1:ib))
          ! if (ims_pro == 0) write(*,*) 'size(yb)', size(yb(1:ib))
          ! if (ims_pro == 0) write(*,*) 'xa', xa(1:ia)
          ! if (ims_pro == 0) write(*,*) 'ya', ya(1:ia)
          ! if (ims_pro == 0) write(*,*) 'xb', xb(1:ib)
          ! if (ims_pro == 0) write(*,*) 'yb', yb(1:ib)
          ! if (ims_pro == 0) write(*,*) '=================='
          ! call TLAB_STOP(DNS_ERROR_IBM_SPLINE)
          ! ================================================================== !
        end if
        ip = ip + nlines
      end do
    end if
  end do

#if 0
  ! ================================================================== !
  ! debug -- Block comment with #if 0 and #endif
  ! ================================================================== !
  if (g%name == 'y') then
    ! ================================================================== !
    ! debug
    if (ims_pro == 0) write(*,*) '========================================================='

    ! write out fld_mod for debugging
    call DNS_TRANSPOSE(fld_mod, kmax, imax*jmax, kmax, fld_mod_tr, imax*jmax)
    call IO_WRITE_FIELDS('fld_mod',  IO_SCAL, imax,jmax,kmax, i1, fld_mod_tr, wrk3d)

    ! stop after writing field
    call TLAB_STOP(i0)
    ! ================================================================== !
  end if

  if (g%name == 'z') then
    ! ================================================================== !
    ! debug
    if (ims_pro == 0) write(*,*) '========================================================='

    ! write out fld_mod for debugging
#ifdef USE_MPI
    if ( ims_npro_k > 1 ) then
      call TLAB_MPI_TRPB_K(fld_mod, fld_mod_tr, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    endif
    call IO_WRITE_FIELDS('fld_mod',  IO_FLOW, imax,jmax,kmax, i1, fld_mod_tr, wrk3d)
#else
    fld_mod_tr = fld_mod
    call IO_WRITE_FIELDS('fld_mod',  IO_FLOW, imax,jmax,kmax, i1, fld_mod_tr, wrk3d)
#endif

    ! stop after writing field
    call TLAB_STOP(i0)
    ! ================================================================== !
  end if
    ! ================================================================== !
#endif
  return
end subroutine IBM_SPLINE_XYZ

!########################################################################

subroutine IBM_SPLINE(ia, ib, xa, ya, xb, yb)
  
  use DNS_IBM,        only: isize_iwrk_ibm, nest, nsp, kspl
  use DNS_IBM,        only: wrk_ibm, iwrk_ibm
  use TLAB_CONSTANTS, only: efile
  use TLAB_PROCS
   
  implicit none
  
#include "integers.h"

  TINTEGER,                intent(in)  :: ia
  TINTEGER,                intent(in)  :: ib
  TREAL,    dimension(ia), intent(in)  :: xa 
  TREAL,    dimension(ia), intent(in)  :: ya 
  TREAL,    dimension(ib), intent(in)  :: xb
  TREAL,    dimension(ib), intent(out) :: yb 

  TREAL                                :: xstart, xend, s, fp
  TINTEGER                             :: iopt, n, l, ier
  TINTEGER                             :: ip1, ip2, ip3, ip4

  character(len=128)                   :: line
  
  ! ================================================================== !
  ! spline function parameter
  iopt = i0    ! (iopt=0 or 1) smoothing spline, weighted least-squares spline (iopt=-1)
  s    = C_0_R ! control the tradeoff between closeness of fit and smoothness

  ! set interval for spline approximation
  xstart = xa(1);  xend = xa(ia)
 
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
  call curfit(iopt, ia,  xa, ya, wrk_ibm(ip1), xstart, xend, kspl, s, nest, n, & 
  !           t,            c,            fp, wrk,          lwrk,           iwrk,     ier)
              wrk_ibm(ip2), wrk_ibm(ip3), fp, wrk_ibm(ip4), isize_iwrk_ibm, iwrk_ibm, ier)
  
  if ( (ier /= 0) .and. (ier /= -1) ) then
    write(line, *) 'IBM_SPLINE. Curfit error code = ', ier
    call TLAB_WRITE_ASCII(efile, line)
    call TLAB_STOP(DNS_ERROR_CURFIT)
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
    write(line, *) 'IBM_SPLINE. Splev error code = ', ier
    call TLAB_WRITE_ASCII(efile, line)
    call TLAB_STOP(DNS_ERROR_CURFIT)
  end if

  return
end subroutine IBM_SPLINE

!########################################################################

subroutine IBM_SPLINE_VECTOR(is, case, fld, g, xa, ya, xb, ia, ib, ip_il, ip_ir, nlines, plane) 

  use DNS_IBM,        only: nflu, isize_wrk1d_ibm, nsp, ibmscaljmin
  use TLAB_VARS,      only: isize_field
  use TLAB_TYPES,     only: grid_dt
  use TLAB_PROCS
  use TLAB_CONSTANTS, only: efile

   
  implicit none
  
#include "integers.h"

  TINTEGER,                          intent(in)  :: is     ! scalar index; if 0, then velocity
  TINTEGER,                          intent(in)  :: case
  TREAL, dimension(isize_field),     intent(in)  :: fld 
  type(grid_dt),                     intent(in)  :: g   
  TREAL, dimension(nsp),             intent(out) :: xa ! max size (not always needed)
  TREAL, dimension(nsp),             intent(out) :: ya  
  TREAL, dimension(isize_wrk1d_ibm), intent(out) :: xb  
  TINTEGER,                          intent(out) :: ia
  TINTEGER,                          intent(out) :: ib
  TINTEGER,                          intent(in)  :: ip_il, ip_ir, nlines, plane

  TINTEGER                                       :: kflu, gap
  TINTEGER                                       :: ip_fl, iu_fl, iu_ir

  character(len=128)                             :: str, line

  ! ================================================================== !
  ! index to remember current position in vectors
  ia = i0
  ib = i0

  ! needed indices
  ip_fl =  ip_il - nflu                     ! k-axis index of most left fluid point
  iu_fl = (ip_fl - 1) * nlines + plane      ! fld-index of most left fluid point
  iu_ir = (ip_ir - 1) * nlines + plane      ! fld-index of right interface  point
  ! ================================================================== !
  ! vectors are built from left to right
  ! -----------------------------------------------------------------
  ! build left half of xa, ya
  select case (case)
  case(i2, i5)
    write(str,*) case; line = 'IBM_SPLINE this case is not implemented yet. Current spline case '//trim(adjustl(str))
    call TLAB_WRITE_ASCII(efile,line)
    call TLAB_STOP(DNS_ERROR_NOTIMPL)
  case(i3)
    ! mirror nflu points on the ground + zeros at left interface
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = - g%nodes((nflu + 2) - kflu) 
      ! ya(ia) =   C_0_R
      if ( is /= i0 ) then
         ya(ia) = ibmscaljmin(is)
      else
        ya(ia) =   C_0_R
      endif
    end do
    ! add zero at right interface (needed if mirroring is not used)
    ! ia     = ia + 1
    ! ya(ia) = C_0_R  
  case(i4, i6)
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(ip_fl + (kflu - 1)) ! dble(ip_fl + (kflu - 1))   
      ya(ia) =     fld(iu_fl + (kflu - 1) * nlines)
    end do
  end select
  ! -----------------------------------------------------------------
  ! set interfaces (left and right)
  ia     = ia + 1
  xa(ia) = g%nodes(ip_il) ! dble(ip_il)   
  ! ya(ia) =   C_0_R
  if ( is /= i0 ) then
    ya(ia) = ibmscaljmin(is)
  else
    ya(ia) =   C_0_R
  endif 
  !
  ia     = ia + 1
  xa(ia) = g%nodes(ip_ir) 
  ! ya(ia) =   C_0_R
  if ( is /= i0 ) then
    ya(ia) = ibmscaljmin(is)
  else
    ya(ia) =   C_0_R
  endif 
 ! -----------------------------------------------------------------
  ! build right half of xa, ya
  select case (case)
  case(i2, i5)
    write(str,*) case; line = 'IBM_SPLINE this case is not implemented yet. Current spline case '//trim(adjustl(str))
    call TLAB_WRITE_ASCII(efile,line)
    call TLAB_STOP(DNS_ERROR_NOTIMPL)
  case(i3, i4)
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(ip_ir + kflu)
      ya(ia) =     fld(iu_ir + kflu * nlines)
    end do
  case(i6)
    ! mirror nflu points on the top + zeros at right interface
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(g%size) + ( g%nodes(g%size) - g%nodes(g%size - kflu) )
      ! ya(ia) =   C_0_R
      if ( is /= i0 ) then
        ya(ia) = ibmscaljmin(is)
      else
        ya(ia) =   C_0_R
      endif 
    end do
    ! add zero at right interface (needed if mirroring is not used)
    ! ia     = ia + 1
    ! ya(ia) = C_0_R  
  end select
! ================================================================== !
! build gap vector where splines are evaluated (here: with interface points)
  do gap = ip_il, ip_ir
    ib      = ib + 1
    xb(ib)  = g%nodes(gap) ! dble(ip_int)
  end do
  ! ================================================================== !
  ! select case (case)
  ! case(i3) ! one interface on the right, extrapolation in gap
  !   ! mirror nflu points on the ground + zeros at left interface
  !   do kflu  = 1, nflu
  !     ia     =   ia + 1
  !     xa(ia) = - g%nodes((nflu + 2) - kflu) 
  !     ya(ia) =   C_0_R
  !   end do 
  !   ! add zero at left interface (needed if mirroring is not used)
  !   ! ia     = ia + 1
  !   ! ya(ia) = C_0_R       
  !   ! set interfaces (left and right)
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_il) ! dble(ip_il)   
  !   ya(ia) = C_0_R   
  !   !
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_ir) 
  !   ya(ia) = C_0_R       
  !   ! build right half of xa, ya (from left to right)
  !   do kflu  = 1, nflu
  !     ia     = ia + 1
  !     xa(ia) = g%nodes(ip_ir + kflu)
  !     ya(ia) =     fld(iu_ir + kflu * nlines)
  !   end do
  ! ! -----------------------------------------------------------------
  ! case(i4)
  !   ! build left half of xa, ya (from left to right)
  !   do kflu  = 1, nflu
  !     ia     = ia + 1
  !     xa(ia) = g%nodes(ip_fl + (kflu - 1)) ! dble(ip_fl + (kflu - 1))   
  !     ya(ia) =     fld(iu_fl + (kflu - 1) * nlines)
  !   end do
  !   ! set interfaces (left and right)
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_il) ! dble(ip_il)   
  !   ya(ia) = C_0_R       
  !   !
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_ir) ! dble(ip_ir) 
  !   ya(ia) = C_0_R       
  !   ! build right half of xa, ya (from left to right)
  !   do kflu  = 1, nflu
  !     ia     = ia + 1
  !     xa(ia) = g%nodes(ip_ir + kflu) ! dble(ip_ir + kflu)
  !     ya(ia) =     fld(iu_ir + kflu * nlines)
  !   end do
  ! ! -----------------------------------------------------------------
  ! case(i6) ! one interface on the left, extrapolation in gap
  !   ! build left half of xa, ya (from left to right)
  !   do kflu  = 1, nflu
  !     ia     = ia + 1
  !     xa(ia) = g%nodes(ip_fl + (kflu - 1)) ! dble(ip_fl + (kflu - 1))   
  !     ya(ia) =     fld(iu_fl + (kflu - 1) * nlines)
  !   end do
  !   ! set interfaces (left and right)
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_il) ! dble(ip_il)   
  !   ya(ia) = C_0_R   
  !   !
  !   ia     = ia + 1
  !   xa(ia) = g%nodes(ip_ir) 
  !   ya(ia) = C_0_R       
  !   ! mirror nflu points on the top + zeros at right interface
  !   do kflu  = 1, nflu
  !     ia     = ia + 1
  !     xa(ia) = g%nodes(g%size) + ( g%nodes(g%size) - g%nodes(g%size - kflu) )
  !     ya(ia) = C_0_R
  !   end do
  !   ! add zero at right interface (needed if mirroring is not used)
  !   ! ia     = ia + 1
  !   ! ya(ia) = C_0_R       
  ! end select
  return
end subroutine IBM_SPLINE_VECTOR
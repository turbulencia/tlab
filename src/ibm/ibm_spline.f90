#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef IBM_DEBUG
#ifdef USE_MPI 
#include "dns_const_mpi.h"  
#endif
#endif
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   cubic spline reconstruction in solid regions
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_SPLINE_XYZ(is, fld, fld_mod, g, nlines, isize_nob, isize_nob_be, nob, nob_b, nob_e, wrk3d)

  use DNS_IBM,        only : xa, xb, ya, yb, nflu, ibmscaljmin
  use TLAB_VARS,      only : isize_field
  use TLAB_CONSTANTS, only : efile
  use TLAB_TYPES,     only : grid_dt
  use TLAB_PROCS

  implicit none
  
#include "integers.h"

  TINTEGER,                          intent(in   ) :: is     ! scalar index; if 0, then velocity
  TREAL,    dimension(isize_field),  intent(in   ) :: fld 
  TREAL,    dimension(isize_field),  intent(  out) :: fld_mod 
  type(grid_dt),                     intent(in   ) :: g
  TINTEGER,                          intent(in   ) :: nlines, isize_nob, isize_nob_be
  TINTEGER, dimension(isize_nob),    intent(in   ) :: nob
  TINTEGER, dimension(isize_nob_be), intent(in   ) :: nob_b, nob_e
  TREAL,    dimension(isize_field),  intent(inout) :: wrk3d  ! for cubic splines subroutine

  TINTEGER                                         :: l, ii, ip, ia, ib, iob, iu_il
  logical                                          :: splines
  TINTEGER, dimension(2)                           :: bc      
  TREAL,    dimension(2)                           :: bcval 
  TREAL                                            :: m1, m2

  ! ================================================================== !
  ! index convention on contiguous lines
  ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

  splines = .true.  ! 1. case doesn't need splines
  fld_mod = fld     ! never modify u,v,w,s directly !

  ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
  do ii = 1, nlines        ! index of ii-plane, loop over plane and check for objects in each line
    if ( nob(ii) /= 0 ) then ! if line contains immersed object(s) --yes-->  spline interpolation
      ip = i0
      do iob = 1, nob(ii)  ! loop over immersed object(s)
        ! select different cases of immersed objects
        if ( nob_b(ip+ii ) == 1) then
        ! ================================================================== !
          if ( nob_e(ip+ii) == g%size ) then
            ! 1. case: object over full extend of line
            splines = .false. ! do nothing
            ! .............................................................. !
          else if ( ( nob_e(ip+ii) <= (g%size - nflu) ) .and. ( g%periodic .eqv. .true. ) ) then
            ! 2. case: object is semi-immersed - periodic case
            call IBM_SPLINE_VECTOR(is, i2, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if ( ( nob_e(ip+ii) <= (g%size - nflu) ) .and. ( g%periodic .neqv. .true. ) ) then ! e.g. in vertical direction
            ! 3. case: object is semi-immersed - non-periodic case - lower boundary
            call IBM_SPLINE_VECTOR(is, i3, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else
            call TLAB_WRITE_ASCII(efile, 'IBM_SPLINE not enough fluid points right of the right interface.')
            call TLAB_STOP(DNS_ERROR_IBM_SPLINE)
          end if
        ! ================================================================== !
        else if ( nob_b(ip+ii) >= (nflu+1) ) then
          if ( nob_e(ip+ii) <= (g%size - nflu) ) then 
            ! 4. case: object is fully immersed
            call IBM_SPLINE_VECTOR(is, i4, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if ( (nob_e(ip+ii) == g%size) .eqv. g%periodic ) then
            ! 5. case: object is semi-immersed - periodic case
            call IBM_SPLINE_VECTOR(is, i5, fld, g, xa, ya, xb, ia, ib, nob_b(ip+ii), nob_e(ip+ii), nlines, ii) 
            ! .............................................................. !
          else if ( (nob_e(ip+ii) == g%size) .neqv. g%periodic ) then  ! e.g. in vertical direction
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
          ! generate splines (other possibility: natural boundary conditions)
          bc(:) = i2 ! fixed first derivative at endpoints
          m1 = (ya(2)  - ya(1)   ) / (xa(2)  - xa(1)   ); bcval(1) = m1 
          m2 = (ya(ia) - ya(ia-1)) / (xa(ia) - xa(ia-1)); bcval(2) = m2
          call CUBIC_SPLINE(bc, bcval, ia, ib, xa(1:ia), ya(1:ia), xb(1:ib), yb(1:ib), wrk3d)
          ! force yb at interface to physical BCs again, to get exact boundary values here
          if ( is /= 0 ) then
            yb(1)  = ibmscaljmin(is)
            yb(ib) = ibmscaljmin(is)
          else
            yb(1)  = C_0_R
            yb(ib) = C_0_R
          end if
          ! fld index of left interface
          iu_il = (nob_b(ip+ii) - 1) * nlines + ii     
          ! replace splines in solid gaps
          do l = 1, ib
            fld_mod(iu_il + (l - 1) * nlines) = yb(l)
          end do
        end if
        ip = ip + nlines
      end do
    end if
  end do

  return
end subroutine IBM_SPLINE_XYZ

!########################################################################

subroutine IBM_SPLINE_VECTOR(is, case, fld, g, xa, ya, xb, ia, ib, ip_il, ip_ir, nlines, plane) 

  use DNS_IBM,        only : nflu, isize_wrk1d_ibm, nspl, ibmscaljmin
  use TLAB_VARS,      only : isize_field
  use TLAB_TYPES,     only : grid_dt
   
  implicit none
  
#include "integers.h"

  TINTEGER,                             intent(in ) :: is
  TINTEGER,                             intent(in ) :: case
  TREAL,    dimension(isize_field),     intent(in ) :: fld 
  type(grid_dt),                        intent(in ) :: g   
  TREAL,    dimension(nspl),            intent(out) :: xa
  TREAL,    dimension(nspl),            intent(out) :: ya  
  TREAL,    dimension(isize_wrk1d_ibm), intent(out) :: xb  
  TINTEGER,                             intent(out) :: ia
  TINTEGER,                             intent(out) :: ib
  TINTEGER,                             intent(in ) :: ip_il, ip_ir, nlines, plane

  TINTEGER                                          :: kflu, gap
  TINTEGER                                          :: ip_fl, iu_fl, iu_ir

  ! ================================================================== !
  ! indices
  ia = i0; ib = i0 ! remember current position in vectors
  select case (case)
  case(i2) ! semi-immersed + periodic
    ip_fl = g%size - nflu                 ! k-axis index of most left fluid point
    iu_fl =  ip_fl      * nlines + plane  ! fld-index of most left fluid point
    iu_ir = (ip_ir - 1) * nlines + plane  ! fld-index of right interface  point
  case(i5) ! semi-immersed + periodic
    ip_fl = ip_il  - nflu                 
    iu_fl = (ip_fl - 1) * nlines + plane  
    iu_ir =                        plane  
  case default
    ip_fl = ip_il  - nflu                 
    iu_fl = (ip_fl - 1) * nlines + plane  
    iu_ir = (ip_ir - 1) * nlines + plane  
  end select

  ! ================================================================== !
  ! vectors are built from left to right
  ! -----------------------------------------------------------------
  ! build left half of xa, ya
  select case (case)
  case(i2) ! semi-immersed + periodic
    ! mirror nflu points
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = - (g%scale - g%nodes(g%size - nflu + kflu))
      ya(ia) = fld(iu_fl + (kflu - 1) * nlines)
    end do
  case(i3) ! semi-immersed + non-periodic
    ! mirror nflu points on the ground + zeros at left interface
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = - g%nodes((nflu + 2) - kflu) 
      if ( is /= i0 ) then
        ya(ia) = ibmscaljmin(is)
      else
        ya(ia) = C_0_R
      endif
    end do
  case(i4, i5, i6) 
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(ip_fl + (kflu - 1))
      ya(ia) =     fld(iu_fl + (kflu - 1) * nlines)
    end do
  end select
  ! -----------------------------------------------------------------
  ! set interfaces (left and right)
  ia     = ia + 1
  xa(ia) = g%nodes(ip_il)  
  if ( is /= i0 ) then
    ya(ia) = ibmscaljmin(is)
  else
    ya(ia) = C_0_R
  endif 
  !
  ia     = ia + 1
  xa(ia) = g%nodes(ip_ir) 
  if ( is /= i0 ) then
    ya(ia) = ibmscaljmin(is)
  else
    ya(ia) = C_0_R
  endif 
  ! -----------------------------------------------------------------
  ! build right half of xa, ya
  select case (case)
  case(i2,i3, i4)
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(ip_ir + kflu)
      ya(ia) =     fld(iu_ir + kflu * nlines)
    end do
  case(i5) ! semi-immersed + periodic
    ! mirror nflu points
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(g%size) + g%nodes(kflu + 1)
      ya(ia) =     fld(iu_ir + (kflu-1) * nlines)
    end do
  case(i6) ! semi-immersed + non-periodic
    ! mirror nflu points on the top + zeros at right interface
    do kflu  = 1, nflu
      ia     = ia + 1
      xa(ia) = g%nodes(g%size) + ( g%nodes(g%size) - g%nodes(g%size - kflu) )
      if ( is /= i0 ) then
        ya(ia) = ibmscaljmin(is)
      else
        ya(ia) = C_0_R
      endif 
    end do
  end select
  ! -----------------------------------------------------------------
  ! build gap vector where splines are evaluated (here: with interface points)
  do gap = ip_il, ip_ir
    ib      = ib + 1
    xb(ib)  = g%nodes(gap)
  end do

  return
end subroutine IBM_SPLINE_VECTOR

!########################################################################
#include "types.h"
#include "dns_const.h"

! #include "dns_error.h"        ! needed if DNS_STOP is used, in this case: define specific error 
                                ! (e.g. number of objects do not fit with grid)
! #ifdef USE_MPI
! #include "dns_const_mpi.h"    ! most likely not needed, if no specific MPI calls are present
! #endif

!########################################################################
!# HISTORY
!#
!# 2021/05/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#
!# Calculate the boundary values for immersed object(s) in the flow  
!# domain. The epsilon field epsi can be seen as an indicator field:
!#    epsi(i,j,k) = 0   for  flow
!#    epsi(i,j,k) = 1   for  solid
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
module DNS_IBM

  ! use 

  implicit none

  ! TINTEGER
  ! TREAL 

  ! all functions/subroutines are private by default, make neeeded ones public 
  private 
  public :: BOUNDARY_BCS_IBM_FLOW

contains
  !########################################################################
  subroutine INITIALIZE_IBM()
    !
    implicit none
    !
    return
  end subroutine INITIALIZE_IBM
  !########################################################################
  subroutine BOUNDARY_BCS_IBM_FLOW(nx, ny, nz, g1, g2, g3, epsi)  ! (ibc, nx,ny,nz, g, u, bcs_hb,bcs_ht, wrk1d,tmp1,tmp2)

    use DNS_TYPES, only     : grid_dt
    use DNS_CONSTANTS, only : lfile

    use DNS_GLOBAL,only     : imax,jmax,kmax

#ifdef USE_MPI 
    use DNS_MPI, only       : ims_offset_i, ims_offset_j, ims_offset_k
    use DNS_MPI, only       : ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k          ! each number of each proc
    use DNS_MPI, only       : ims_npro, ims_npro_i, ims_npro_j, ims_npro_k         ! total numbers of proc
#endif 

    implicit none

#include "integers.h"

#ifdef USE_MPI
#include "mpif.h"
#else
    TINTEGER, parameter  :: ims_offset_i=0, ims_offset_j=0, ims_offset_k=0
    TINTEGER, parameter  :: ims_pro_i=0,    ims_pro_j=0,    ims_pro_k=0,    ims_pro=0  
    TINTEGER, parameter  :: ims_npro_i=1,   ims_npro_j=1,   ims_npro_k=1,   ims_npro=0 
#endif

    ! later in dns.ini file
    TINTEGER, parameter        :: nbars = 4                ! number of bars in
    TINTEGER, parameter        :: hbar  = 10, wbar = 10    ! width, height of bars in gridpoints
    TREAL                      :: zcenter_bar
    TINTEGER, dimension(nbars) :: zstart_bar, zend_bar

    ! ini
    TINTEGER,intent(in)                           :: nx,ny,nz!, ibc
    TYPE(grid_dt),                    intent(in)  :: g1, g2, g3  ! grid informations
    TREAL, dimension(imax,jmax,kmax), intent(out) :: epsi
    ! -------------------------------------------------------------------
    TINTEGER                                      :: istart, iend, jstart, jend, kstart, kend
    TINTEGER                                      :: i,j,k,l

    
    CHARACTER*250 line1


    !###################################################################

    ! global array indicies for each mpi task (indices starting with 1)
    ! istart = ims_offset_i+1; iend = ims_offset_i + imax    
    ! jstart = ims_offset_j+1; jend = ims_offset_j + jmax 
    ! kstart = ims_offset_k+1; kend = ims_offset_k + kmax 
    istart = ims_offset_i; iend = ims_offset_i + imax - 1    
    jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
    kstart = ims_offset_k; kend = ims_offset_k + kmax - 1
    ! debugging
    if ( ims_pro .EQ. 0 ) then 
      write(*,*) '======== Initialization of Grid and Decomposition ======='
      write(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax,' x ', kmax*ims_npro_k 
      write(*,*) 'DECOMP: ranks', ims_npro_i,' x ',1,' x ',ims_npro_k 
      write(*,*) '        grid ', imax,      ' x ',jmax,      ' x ',kmax 
      ! WRITE(*,*) '======== Global Grid Information of each Task     ======='
    end if
    ! IF ( ims_pro .EQ. 0 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 1 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 2 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 3 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF

    ! if square bars in x-/streamwise direction equally distributed
        ! -------------------------------------------------------------------
        ! Requirenments (interfaces of bars have to be on gridpoints, center doesn't matter), if:
        !   kmax_total/(2*nbars) = integer         :: wbar = even   (center is on gridpoints)
        !   kmax_total/(2*nbars) = x * 1/2         :: wbar = uneven (center is between grid points)
        !
        ! 
        ! -------------------------------------------------------------------


    ! global z-positions of bars
    do l = 1, nbars
      zcenter_bar   = g3%size / nbars * (l - 0.5)
      zstart_bar(l) = int(zcenter_bar - 0.5 * wbar)
      zend_bar(l)   = int(zcenter_bar + 0.5 * wbar)
    end do

    if ( ims_pro .EQ. 0 ) then
      write(*,*) '======== Z - Positions of Bars                    ======='
      do l = 1, nbars
        write(*,*) zstart_bar(l), ' ---- ', zend_bar(l)
      end do
    end if
    
    ! ini epsi
    ! write(*,*) 'everthing fine 1'
    epsi(:,:,:) = C_0_R
    ! write(*,*) 'everthing fine 2'
    
    do j = 1, hbar   
      do k = 1, kmax 
        do l = 1, nbars 
          if( ((k+kstart).ge.zstart_bar(l)) .and. ((k+kstart).le.zend_bar(l)) ) then 
            do i = 1, imax
              epsi(i,j,k) = C_1_R
            end do
          end if
        end do 
      end do
    end do 

    ! WRITE(*,*) '=== array positions of bars==='
    ! write(*,*) 'epsilon size x ', size(epsi,1)
    ! write(*,*) 'epsilon size y ', size(epsi,2)
    ! write(*,*) 'epsilon size z ', size(epsi,3)
    ! write(*,*) 'iend-istart    ', iend - istart
    ! write(*,*) 'jend-jstart    ', jend - jstart
    ! write(*,*) 'kend-kstart    ', kend - kstart

    return
  end subroutine BOUNDARY_BCS_IBM_FLOW
  !########################################################################
  subroutine WRITE_GEOMETRY()
    !
    implicit none
    !
    return
  end subroutine WRITE_GEOMETRY
  !########################################################################
  subroutine BOUNDARY_BCS_IBM_SCAL()
    !
    implicit none
    !
    return
  end subroutine BOUNDARY_BCS_IBM_SCAL
  !########################################################################
  subroutine IBM_FINALIZE()
    !
    implicit none
    !
    return
  end subroutine IBM_FINALIZE
  !########################################################################
!########################################################################
end module DNS_IBM
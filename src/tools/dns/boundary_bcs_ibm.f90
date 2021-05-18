#include "types.h"
#include "dns_const.h"

!  #include "dns_error.h"  ! needed if DNS_STOP is used, in this case: define specific error 

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/05/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#
!# epsilon field 'epsi' is an indicator field:
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
!# SUBROUTINES IN MODULE
!#
!# 
!#                            
!#                           
!#
!########################################################################
module DNS_IBM

  use DNS_GLOBAL, only: imax,jmax,kmax    
  use DNS_LOCAL,  only: xbars_geo

  implicit none

  TREAL, dimension(:,:,:), allocatable :: xepsi, zepsi
  TREAL, dimension(:,:),   allocatable :: nobjx, nobjy, nobjz ! number of objects in x/y/z 
  TREAL, dimension(:,:,:), allocatable :: xobji, yobji, zobji ! location of starting interfaces of object (fluid --> solid)
  TREAL, dimension(:,:,:), allocatable :: xobjf, yobjf, zobjf ! location of ending   interfaces of object (solid --> fluid)
  TINTEGER                             :: ximax, xjmax, xkmax, zimax, zjmax, zkmax ! new array dimensions after transpositions

  ! all functions/subroutines are private by default, make needed subroutines public 
  private 
  public  :: INITIALIZE_GEOMETRY

contains
  !########################################################################
  subroutine INITIALIZE_IBM()       ! called once in dns_main.f90 before time integration starts, initialize IBM procedure
    !
    implicit none
    !     
! #ifdef USE_MPI 
! #include "mpif.h"
! #endif 
    !
    return
  end subroutine INITIALIZE_IBM
  !########################################################################
  subroutine INITIALIZE_GEOMETRY(epsi)  ! should be called once in dns_main.f90 before time integration starts
                                        ! after calling, all relevant geometry information needed for the IBM are available
    !
    implicit none
    !
#ifdef USE_MPI 
#include "mpif.h"
#endif 
    
    TREAL, dimension(imax,jmax,kmax):: epsi

    ! ================================================================== !

    ! generate geometry (epsi) of immersed objects (define own routine)
    call GEOMETRY_XBARS(epsi) 

    ! transpose geometry (epsi) and allocate memory
    call TRANSPOSE_GEOMETRY(epsi)  ! xepsi, zepsi

    ! generate relevant geometry fields             
    call GENERATE_GEOMETRY(epsi)   ! xepsi, zepsi, nobjx, nobjy, nobjz, xobji, xobjf, yobji, yobjf, zobji, zobjf
    
    ! write/read geometry fields to/from disk
    ! call WRITE_GEOMETRY()           
    
    return
  end subroutine INITIALIZE_GEOMETRY
  !########################################################################
  subroutine GEOMETRY_XBARS(epsi)
    ! -------------------------------------------------------------------
    ! if square bars in x-/streamwise direction equally distributed
    ! -------------------------------------------------------------------
    ! Requirenments (interfaces of bars have to be on gridpoints, center doesn't matter), if:
    !   kmax_total/(2*nbars) = integer         :: wbar = even   (center is on gridpoints)
    !   kmax_total/(2*nbars) = x * 1/2         :: wbar = uneven (center is between grid points)
    ! -------------------------------------------------------------------

    use DNS_GLOBAL,    only: g

#ifdef USE_MPI 
    use DNS_MPI,       only: ims_offset_i, ims_offset_j, ims_offset_k
    use DNS_MPI,       only: ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k          ! each number of each proc
    use DNS_MPI,       only: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k         ! total numbers of proc
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

    TINTEGER                                         :: nbars, hbar, wbar
    TREAL,    dimension(imax,jmax,kmax), intent(out) :: epsi ! or: intent(inout)?
    TREAL                                            :: zcenter_bar
    TINTEGER, dimension(xbars_geo(1))                :: zstart_bar, zend_bar
    TINTEGER                                         :: istart, iend, jstart, jend, kstart, kend
    TINTEGER                                         :: i,j,k,l

    ! debugging 
    CHARACTER(32)                                    :: fname
    TREAL, dimension(imax*jmax*kmax)                 :: wrk3d

    ! ================================================================== !

    ! global array indicies for each mpi task (indices starting with 0)
    istart = ims_offset_i; iend = ims_offset_i + imax - 1    
    jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
    kstart = ims_offset_k; kend = ims_offset_k + kmax - 1

    ! debugging
    if ( ims_pro .eq. 0 ) then 
      write(*,*) '======== Initialization of Grid and Decomposition ======='
      write(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax,' x ', kmax*ims_npro_k 
      write(*,*) 'DECOMP: ranks', ims_npro_i,' x ',1,' x ',ims_npro_k 
      write(*,*) '        grid ', imax,      ' x ',jmax,      ' x ',kmax 
    end if
    ! IF ( ims_pro .EQ. 0 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 1 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 2 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
    ! IF ( ims_pro .EQ. 3 ) THEN; write(*,*) 'Task ', ims_pro, 'with indices: ', istart,'x',iend,' # ',jstart,'x',jend,' # ',kstart,'x',kend;ENDIF
   
    ! geometry informtion from dns.ini
    nbars=xbars_geo(1); hbar=xbars_geo(2); wbar=xbars_geo(3)  
   
    ! global z-positions of bars, equally distributed on gridpoints
    do l = 1, nbars
      zcenter_bar   = g(3)%size / nbars * (l - 0.5)
      zstart_bar(l) = int(zcenter_bar - 0.5 * wbar)
      zend_bar(l)   = int(zcenter_bar + 0.5 * wbar)
    end do

    ! debugging
    if ( ims_pro .eq. 0 ) then
      write(*,*) '======== Z - Positions of Bars =========================='
      do l = 1, nbars
        write(*,*)'bar', l, ' start:', zstart_bar(l), '  end:', zend_bar(l)
      end do
    end if
    
    ! ini epsi
    epsi(:,:,:) = C_0_R
    
    ! streamwise aligned square bars, equaly spaced, interfaces on gridpoints (define own geometry)
    do j = 1, hbar   
      do k = 1, kmax 
        do l = 1, nbars 
          if( ((k+kstart).gt.zstart_bar(l)) .and. ((k+kstart).le.zend_bar(l)) ) then 
            do i = 1, imax
              epsi(i,j,k) = C_1_R
            end do
          end if
        end do 
      end do
    end do 

    ! write epsi field
    write(fname,*) i0; 
    fname = TRIM(ADJUSTL('epsi'))//TRIM(ADJUSTL(fname))
    if (ims_pro .eq. 0) then
      write(*,*) '========================================================='  
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, 1, imax*jmax*kmax, epsi, wrk3d)

    return
  end subroutine GEOMETRY_XBARS
  !########################################################################
  subroutine GENERATE_GEOMETRY(epsi)!, nobjx, nobjy, nobjz, xobji, xobjf, yobji, yobjf, zobji, zobj)  
    !
#ifdef USE_MPI
    use DNS_MPI, only: ims_pro
#endif    
    !
    implicit none
    !
#include "integers.h"
    !
#ifdef USE_MPI 
#include "mpif.h"
#endif 

    TREAL, dimension(imax,jmax,kmax)  :: epsi   ! intent() ?
    TINTEGER                          :: nobjx_max, nobjy_max, nobjz_max, num
    TINTEGER                          :: i,j,k,l

    CHARACTER(32)                     :: fname
    TREAL, DIMENSION(imax*kmax)       :: wrk3d

    ! ================================================================== !

    ! ini nobj arrays
    nobjx(:,:) = C_0_R; nobjy(:,:) = C_0_R; nobjz(:,:) = C_0_R
    nobjx_max  = C_0_R; nobjy_max  = C_0_R; nobjz_max  = C_0_R

    ! number of objects in y-direction
    do i = 1, imax
      do k = 1, kmax
        num = 0
        if(epsi(i,1,k) .eq. 1) then
          nobjy(i,k) = 1
          num = 1
        end if
        do j = 1, jmax-1       
          if((epsi(i,j,k) .eq. 0) .and. (epsi(i,j+1,k) .eq. 1)) then
            nobjy(i,k) = nobjy(i,k) + 1
            num = num + 1
          end if
        end do
        if (num .gt. nobjy_max) then
          nobjy_max = num    
        end if
      end do    
    end do

    ! debugging
    if (ims_pro .eq. 0) then
      write(*,*) '========================================================='
      write(*,*) 'nobjy_max   =', nobjy_max
      ! write nobjy field
    end if
    write(fname,*) i0; 
    fname = TRIM(ADJUSTL('nobjy'))//TRIM(ADJUSTL(fname))
    if (ims_pro .eq. 0) then  
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, imax,1,kmax, 1, imax*1*kmax, nobjy, wrk3d)

    ! number of objects in x-direction

    return
  end subroutine GENERATE_GEOMETRY
    !########################################################################
  subroutine TRANSPOSE_GEOMETRY(epsi) !,xepsi,yepsi,zepsi 
    !
#ifdef USE_MPI
    use DNS_MPI, only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use DNS_MPI, only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use DNS_MPI, only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
#endif    
    !
    implicit none
    !
#include "integers.h"
    !
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
#endif 

    TREAL, dimension(imax,jmax,kmax)  :: epsi

#ifdef USE_MPI
    TINTEGER, parameter               :: idi = DNS_MPI_I_PARTIAL
    TINTEGER, parameter               :: idk = DNS_MPI_K_PARTIAL
#endif

    ! ================================================================== !

    ! new array indices
#ifdef USE_MPI
    ximax = ims_npro_i * imax ! imax_total
    xjmax = jmax / ims_npro_i 
    xkmax = kmax
    !
    zimax = imax
    zjmax = jmax / ims_npro_k
    zkmax = ims_npro_k * kmax ! kmax_total
#else
    ximax = imax; xjmax = jmax; xkmax = kmax
    zimax = imax; zjmax = jmax; zkmax = kmax
#endif

    ! allocate memory for arrays
    allocate(xepsi(ximax,xjmax,xkmax), zepsi(zimax,zjmax,zkmax))
    allocate(nobjx(xjmax,xkmax), nobjy(imax,kmax), nobjz(zimax,zjmax))
    allocate(xobji(xbars_geo(1),xjmax,xkmax), yobji(xbars_geo(1),imax,kmax), zobji(xbars_geo(1),zimax,zjmax))
    allocate(xobjf(xbars_geo(1),xjmax,xkmax), yobjf(xbars_geo(1),imax,kmax), zobjf(xbars_geo(1),zimax,zjmax))

    ! debugging
    if ( ims_pro .eq. 0 ) then 
      write(*,*) '======== TRANSPOSING GEOMETRY ==========================='
      write(*,*) 'Array xepsi: ', size(xepsi,1), ' x ', size(xepsi,2), ' x ', size(xepsi,3)
      write(*,*) 'Array yepsi: ', size(epsi,1), ' x ', size(epsi,2), ' x ', size(epsi,3)
      write(*,*) 'Array zepsi: ', size(zepsi,1), ' x ', size(zepsi,2), ' x ', size(zepsi,3)    
    end if

    ! transpose in x-direction
#ifdef USE_MPI
    if ( ims_npro_i .gt. 1 ) then
      call DNS_MPI_TRPF_I(epsi, xepsi, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
    else
#endif
      xepsi = epsi
#ifdef USE_MPI
    end if
#endif

    ! transpose in z-direction
#ifdef USE_MPI
    if ( ims_npro_k .gt. 1 ) then
      call DNS_MPI_TRPF_K(epsi, zepsi, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    else
#endif
      zepsi = epsi
#ifdef USE_MPI
    end if
#endif

    if ( ims_pro .eq. 0 ) then 
      write(*,*) '========================================================='
      write(*,*) 'Transposing and Allocating arrays: done'
      write(*,*) '========================================================='    
      write(*,*) 'Generating geometry'     
    end if

    return
  end subroutine TRANSPOSE_GEOMETRY
  !########################################################################
  subroutine BOUNDARY_BCS_IBM_FLOW()

    implicit none
    
#ifdef USE_MPI 
#include "mpif.h"
#endif 

    return
  end subroutine BOUNDARY_BCS_IBM_FLOW
  !########################################################################
  subroutine WRITE_GEOMETRY()
    !
    implicit none
  
#ifdef USE_MPI 
#include "mpif.h"
#endif 

    ! write epsi field
    ! m = m + 1
    ! WRITE(fname,*) m; 
    ! fname = TRIM(ADJUSTL('epsi'))//TRIM(ADJUSTL(fname))
    ! write(*,*) fname ! debugging
    ! CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, 1, imax*jmax*kmax, epsi, wrk3d)
    
    return
  end subroutine WRITE_GEOMETRY
  !########################################################################
  subroutine READ_GEOMETRY() ! if needed: restart/run with already generated geometry
    !
    implicit none
  
#ifdef USE_MPI 
#include "mpif.h"
#endif 
    !
    return
  end subroutine READ_GEOMETRY
  !########################################################################  
  subroutine BOUNDARY_BCS_IBM_SCAL()
    !
    implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
    !
    return
  end subroutine BOUNDARY_BCS_IBM_SCAL
  !########################################################################
  subroutine IBM_FINALIZE() ! dealloc of arrays? maybe not needed...
    !
    implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
    !
    return
  end subroutine IBM_FINALIZE
  !########################################################################
!########################################################################
end module DNS_IBM
#include "types.h"
#include "dns_const.h"

!  #include "dns_error.h"  ! needed if DNS_STOP is used, in this case: define specific error 

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/06/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#
!# epsilon field 'eps' is an indicator field:
!#    eps(i,j,k) = 0   for  flow
!#    eps(i,j,k) = 1   for  solid
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

  use DNS_GLOBAL, only: imax, jmax, kmax, isize_field, isize_txc_field, inb_txc
  use DNS_LOCAL,  only: xbars_geo

  implicit none

  TREAL, dimension(:,:,:), allocatable :: eps_aux                ! eps_aux field for debugging and geometry generation
  TREAL, dimension(:),     allocatable :: epsi, epsj, epsk, eps  ! eps transposed in i/j/k
  TREAL, dimension(:),     allocatable :: nobi, nobj, nobk       ! number of objects        in i/j/k 
  ! TREAL, dimension(:,:,:), allocatable :: ! location of starting interfaces of object (fluid --> solid)
  ! TREAL, dimension(:,:,:), allocatable :: ! location of ending   interfaces of object (solid --> fluid)


  ! all functions/subroutines are private by default
  ! puplish only needed subroutines here 
  private 
  public  :: ALLOCATE_IBM, INITIALIZE_GEOMETRY, eps, epsi, epsj, epsk 

contains
  !########################################################################
  subroutine ALLOCATE_IBM(allocated)  ! called once in dns_main.f90 l.210 to allocate needed memory
    !
    use DNS_CONSTANTS, only: lfile, efile
    !
#ifdef USE_MPI
    use DNS_MPI, only: ims_size_i, ims_size_j, ims_size_k 
#endif    
    !  
    implicit none

#include "dns_error.h"
#include "integers.h"
    !
    !
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
    TINTEGER, parameter       :: idi = DNS_MPI_I_PARTIAL ! = 1
    TINTEGER, parameter       :: idj = DNS_MPI_J_PARTIAL ! = 1
    TINTEGER, parameter       :: idk = DNS_MPI_K_PARTIAL ! = 1 --> idi = idj = idk
#endif
    ! 
    logical, intent(inout)    :: allocated       ! flag, just allocate memory space once
    TINTEGER                  :: ierr, inb_ibm
    TINTEGER                  :: nyz, nxz, nxy
    character(128)            :: str, line
   
    ! ================================================================== !
    inb_ibm = i1 ! can be also defined in dns_read_local.f90 and module dns_global.f90 (cf. inb_flow ...)

    ! npages (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI 
    nyz = ims_size_i(idi) ! local
    nxz = ims_size_j(idj) 
    nxy = ims_size_k(idk) 
#else
    nyz = jmax * kmax     ! global
    nxz = imax * kmax     
    nxy = imax * jmax     
#endif      

    ! ================================================================== !

    ! allocate here all ibm related arrays
    if ( .not. allocated ) then 

      ! eps_aux (3D-array), for debugging and easy generating geometries
      write(str,*) inb_ibm; line = 'Allocating array IBM eps_aux of size '//trim(adjustl(str))//'x'
      ! write(str,*) inb_ibm; line = 'Allocating array IBM eps_aux  of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(eps_aux(imax,jmax,kmax), stat=ierr)
      ! allocate(eps_aux(isize_field,inb_ibm), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps_aux.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! eps
      write(str,*) inb_ibm; line = 'Allocating array IBM eps of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(eps(isize_field), stat=ierr)
      ! allocate(eps(isize_field,inb_ibm), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! ================================================================== !

      ! epsi
      write(str,*) inb_ibm; line = 'Allocating array IBM epsi of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(epsi(isize_field), stat=ierr)
      ! allocate(epsi(isize_field,inb_ibm), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsi.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! epsj
      write(str,*) inb_ibm; line = 'Allocating array IBM epsj of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(epsj(isize_field), stat=ierr)
      ! allocate(epsj(isize_field,inb_ibm), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsj.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! epsk
      write(str,*) inb_ibm; line = 'Allocating array IBM epsk of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(epsk(isize_field), stat=ierr)
      ! allocate(epsk(isize_field,inb_ibm), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsk.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! ================================================================== !

      ! nobi
      write(str,*) inb_ibm; line = 'Allocating array IBM nobi of size '//trim(adjustl(str))//'x'
      write(str,*) nyz; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(nobi(nyz), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! nobj
      write(str,*) inb_ibm; line = 'Allocating array IBM nobj of size '//trim(adjustl(str))//'x'
      write(str,*) nxz; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(nobj(nxz), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if
      
      ! nobk
      write(str,*) inb_ibm; line = 'Allocating array IBM nobk of size '//trim(adjustl(str))//'x'
      write(str,*) nxy; line = trim(adjustl(line))//trim(adjustl(str))
      call IO_WRITE_ASCII(lfile,line)
      allocate(nobk(nxy), stat=ierr)
      if ( ierr .ne. 0 ) then
        call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk.')
        call DNS_STOP(DNS_ERROR_ALLOC)
      end if

      ! ================================================================== !

      ! set flag
      allocated = .true.
    end if
    !
  return
  end subroutine ALLOCATE_IBM
  !########################################################################
  subroutine INITIALIZE_IBM()       ! called once in dns_main.f90 before time integration starts
                                    ! initialize IBM procedure
    !
    implicit none
    !     
! #ifdef USE_MPI 
! #include "mpif.h"
! #endif 
    !
    return
  end subroutine 
  !########################################################################
  subroutine INITIALIZE_GEOMETRY(txc, wrk3d)  ! should be called once in dns_main.f90 before time integration starts
                                    ! after calling, all relevant geometry information needed for the IBM are available and written to disk
    !
    implicit none
    !
#ifdef USE_MPI 
#include "mpif.h"
#endif 
    
    TREAL, dimension(*), intent(inout) :: txc, wrk3d

    ! ================================================================== !

    ! generate native 3d-geometry field (eps_aux) of immersed objects (define your own geomtry here)
    call GEOMETRY_XBARS(wrk3d) 

    ! transpose eps in epsi, epsj, epsk and allocate neccessary memory
    call TRANSPOSE_GEOMETRY(txc)

    ! generate relevant geometry fields for IBM routines (nobi, nobj, nobk)
    call GENERATE_GEOMETRY(wrk3d)
    
    ! read/write geometry fields from/to disk
    ! call READ_GEOMETRY() ! call WRITE_GEOMETRY()           
    
    return
  end subroutine INITIALIZE_GEOMETRY
  !########################################################################
  subroutine GEOMETRY_XBARS(wrk3d)
    ! -------------------------------------------------------------------
    ! here: square bars in x-/streamwise direction, equally distributed + spaced in z-/spanwise direction
    ! -------------------------------------------------------------------
    !   Requirenments (interfaces of bars have to be on gridpoints, center doesn't matter), if:
    !     kmax_total/(2*nbars) = integer         :: wbar = even   (center is on gridpoints)
    !     kmax_total/(2*nbars) = x * 1/2         :: wbar = uneven (center is between grid points)
    ! -------------------------------------------------------------------

    use DNS_GLOBAL, only: g ! type grid, all grid related informations  

#ifdef USE_MPI 
    use DNS_MPI,    only: ims_offset_i, ims_offset_j, ims_offset_k
    use DNS_MPI,    only: ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k          ! each number of each proc
    use DNS_MPI,    only: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k         ! total numbers of proc
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

    TINTEGER                                     :: nbars, hbar, wbar
    TREAL                                        :: zcenter_bar
    TINTEGER, dimension(xbars_geo(1))            :: zstart_bar, zend_bar
    TINTEGER                                     :: istart, iend, jstart, jend, kstart, kend
    TINTEGER                                     :: i,j,k,l

    ! debugging 
    character(32)                                :: fname
    TREAL, dimension(isize_field), intent(inout) :: wrk3d

    ! ================================================================== !

    ! global array indicies for each mpi task (indices start with 0)
    istart = ims_offset_i; iend = ims_offset_i + imax - 1    
    jstart = ims_offset_j; jend = ims_offset_j + jmax - 1
    kstart = ims_offset_k; kend = ims_offset_k + kmax - 1

    ! debugging
    if ( ims_pro .eq. 0 ) then 
      write(*,*) '======== Initialization of Grid and Decomposition ======='
      write(*,*) 'GRID:        ', imax*ims_npro_i,' x ', jmax, ' x ', kmax*ims_npro_k
      ! write(*,*) 'GRID g(...): ', g(1)%size,      ' x ', g(2)%size,' x ', g(3)%size 
      write(*,*) 'DECOMP: ranks', ims_npro_i,     ' x ', 1,    ' x ', ims_npro_k 
      write(*,*) '        grid ', imax,           ' x ', jmax, ' x ', kmax 
    end if
     
    ! geometry (from dns.ini)
    nbars=xbars_geo(1); hbar=xbars_geo(2); wbar=xbars_geo(3)  
   
    ! global z-positions of bars, equally distributed on gridpoints with equal spacing
    do l = 1, nbars
      zcenter_bar   = g(3)%size / nbars * (l - 0.5)
      zstart_bar(l) = int(zcenter_bar - 0.5 * wbar)
      zend_bar(l)   = int(zcenter_bar + 0.5 * wbar)
    end do

    ! debugging
    if ( ims_pro .eq. 0 ) then
      write(*,*) '======== Z - Positions of streamwise aligned Bars ======='
      do l = 1, nbars
        write(*,*)'bar nr.', l, ' start:', zstart_bar(l), ' end:', zend_bar(l)
      end do
    end if
    
    ! ini eps_aux
    eps_aux(:,:,:) = C_0_R
    
    ! streamwise aligned square bars, equaly spaced, interfaces on gridpoints (define own geometry)
    do j = 1, hbar   
      do k = 1, kmax 
        do l = 1, nbars 
          if( ((k+kstart).gt.zstart_bar(l)) .and. ((k+kstart).le.zend_bar(l)) ) then 
            do i = 1, imax
              eps_aux(i,j,k) = C_1_R
            end do
          end if
        end do 
      end do
    end do 

    ! reshape 3D-eps_aux field into 1D-eps
    eps = reshape(eps_aux,(/isize_field/))

    ! write eps field
    write(fname,*) i0; 
    fname = trim(adjustl('eps'))//trim(adjustl(fname))
    if (ims_pro .eq. 0) then
      write(*,*) '========================================================='  
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, 1, imax*jmax*kmax, eps, wrk3d)

    return
  end subroutine GEOMETRY_XBARS
!########################################################################
  subroutine TRANSPOSE_GEOMETRY(txc)
    !
    use DNS_GLOBAL, only: g ! type grid, all grid related informations  
    !
#ifdef USE_MPI
    use DNS_MPI, only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use DNS_MPI, only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use DNS_MPI, only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
    use DNS_MPI, only: ims_size_i, ims_size_j, ims_size_k 
#endif    
    
    implicit none
    
#include "integers.h"
    
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
    TINTEGER, parameter                                  :: idi = DNS_MPI_I_PARTIAL ! = 1
    TINTEGER, parameter                                  :: idj = DNS_MPI_J_PARTIAL ! = 1
    TINTEGER, parameter                                  :: idk = DNS_MPI_K_PARTIAL ! = 1 --> idi = idj = idk
#else
    TINTEGER, parameter                                  :: ims_pro = 0             ! task in serial mode
#endif

    TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc 
    TREAL, dimension(isize_field)                        :: tmp1                    !... tmp6    

    TINTEGER                                             :: nyz, nxz, nxy
    
    ! ================================================================== !
    
    ! tmp aux 
    tmp1(:) = txc(:,1)

    ! npages (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI
    nyz = ims_size_i(idi) ! npage = jmax*kmax
    nxz = ims_size_j(idj) ! npage = imax*kmax
    nxy = ims_size_k(idk) ! npage = imax*jmax
#else
    nyz = jmax * kmax     ! global
    nxz = imax * kmax     ! global
    nxy = imax * jmax     ! global
#endif

    ! debugging
    if ( ims_pro .eq. 0 ) then 
      write(*,*) '======== DEBUG TRANSPOSING GEOMETRY ====================='
      write(*,*) 'ims_size_i(1)    ', nyz
      write(*,*) 'ims_size_j(1)    ', nxz
      write(*,*) 'ims_size_k(1)    ', nxy
      write(*,*) 'isize_field      ', isize_field
      write(*,*) 'isize_txc_field  ', isize_txc_field
      write(*,*) 'inb_txc          ', inb_txc
    end if

    ! MPI transposition and local transposition (make x-direction the last one) in x 
#ifdef USE_MPI
    if ( ims_npro_i .GT. 1 ) then
      call DNS_MPI_TRPF_I(eps, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
    else
#endif
    tmp1 = eps
#ifdef USE_MPI
    end if
#endif
  
#ifdef USE_ESSL
    call DGETMO       (tmp1, g(1)%size, g(1)%size, nyz,       epsi, nyz)
#else
    call DNS_TRANSPOSE(tmp1, g(1)%size, nyz,       g(1)%size, epsi, nyz)
#endif

    ! local transposition (make y-direction the last one) in y 
    ! (no MPI transposition needed, native standing in y)
#ifdef USE_ESSL
    call DGETMO       (eps, nxy, nxy,  kmax, epsj, kmax)
#else

    call DNS_TRANSPOSE(eps, nxy, kmax, nxy,  epsj, kmax)
#endif

    ! MPI transposition in z 
    ! (no local transposition needed, already in right order)
#ifdef USE_MPI
    if ( ims_npro_k .GT. 1 ) then
      call DNS_MPI_TRPF_K(eps, epsk, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    else
#endif
    epsk = eps
#ifdef USE_MPI
    end if
#endif

    ! debugging
    if ( ims_pro .eq. 0 ) then 
      write(*,*) '========================================================='
      write(*,*) 'Transposing and Allocating arrays: done'
      write(*,*) '========================================================='    
      write(*,*) 'Generating geometry'     
    end if

    return
  end subroutine TRANSPOSE_GEOMETRY
!############# DEBUGGING, BLOCK COMMENT (until end of module) ###########
! #if 0
!########################################################################
  !########################################################################
  subroutine GENERATE_GEOMETRY(wrk3d)
    !
    use DNS_GLOBAL, only: g ! type grid, all grid related informations  
    !
#ifdef USE_MPI
    use DNS_MPI, only: ims_pro
    use DNS_MPI, only: ims_size_i, ims_size_j, ims_size_k 
#endif    
    !
    implicit none
    !
#include "integers.h"
    !
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
    TINTEGER, parameter                          :: idi = DNS_MPI_I_PARTIAL ! = 1
    TINTEGER, parameter                          :: idj = DNS_MPI_J_PARTIAL ! = 1
    TINTEGER, parameter                          :: idk = DNS_MPI_K_PARTIAL ! = 1 --> idi = idj = idk
#else
    TINTEGER, parameter                          :: ims_pro=0  
#endif

    TREAL, dimension(isize_field), intent(inout) :: wrk3d 

    TINTEGER                                     :: nobi_max, nobj_max, nobk_max
    TINTEGER                                     :: i, j, k, ij, ik, jk, ip
    TINTEGER                                     :: nyz, nxz, nxy

    CHARACTER(32)                                :: fname

    ! DEBUG
    TREAL, dimension(imax*kmax):: nobjj 


    ! ================================================================== !

    ! npages (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI
    nyz = ims_size_i(idi) ! local
    nxz = ims_size_j(idj) 
    nxy = ims_size_k(idk) 
#else
    nyz = jmax * kmax     ! global
    nxz = imax * kmax     
    nxy = imax * jmax     
#endif

    ! initialize 
    nobi(:)  = C_0_R; nobj(:)  = C_0_R; nobk(:)  = C_0_R
    nobi_max = i0;    nobj_max = i0;    nobk_max = i0

    ! ================================================================== !
    ! number of objects in x-direction
    ip = 1
    do i = 1, g(1)%size - 1     ! contiguous i-lines
      do jk = 1, nyz            ! pages of   i-lines
        if((ip .eq. 1) .and. (epsi(jk) .eq. 1)) then ! exception: check first plane for objects
          nobi(jk) = C_1_R
        end if 
        if((epsi(ip+jk-1) .eq. 0) .and. (epsi(ip+jk-1+nyz) .eq. 1)) then ! check for interface 
          nobi(jk) = nobi(jk) + C_1_R
        end if
      end do
      ip = ip + nyz
    end do

    ! ! nobi_max = int(maxval(nobi)) 
    ! ! CALL MPI_REDUCE(idummy,t_dif,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD) 
 
    ! ================================================================== !
    ! number of objects in y-direction
    ip = 1
    do j = 1, g(2)%size - 1     ! contiguous j-lines
      do ik = 1, nxz            ! pages of   j-lines
        if((ip .eq. 1) .and. (epsj(ik) .eq. 1)) then ! exception: check first plane for objects
          nobj(ik) = C_1_R
        end if 
        if((epsj(ip+ik-1) .eq. 0) .and. (epsj(ip+ik-1+nxz) .eq. 1)) then ! check for interface 
          nobj(ik) = nobj(ik) + C_1_R
        end if
      end do
      ip = ip + nxz
    end do
    ! nobj_max = int(maxval(nobj))   
    ! ================================================================== !
    ! number of objects in z-direction
    ip = 1
    do k = 1, g(3)%size - 1     ! contiguous k-lines
      do ij = 1, nxy            ! pages of   k-lines
        if((ip .eq. 1) .and. (epsk(ij) .eq. 1)) then ! exception: check first plane for objects
          nobk(ij) = C_1_R
        end if 
        if((epsk(ip+ij-1) .eq. 0) .and. (epsk(ip+ij-1+nxy) .eq. 1)) then ! check for interface 
          nobk(ij) = nobk(ij) + C_1_R
        end if
      end do
      ip = ip + nxy
    end do
    ! ! nobk_max = int(maxval(nobk))    
    ! ! ================================================================== !

    ! debugging
    if (ims_pro .eq. 0) then
      write(*,*) '========================================================='
      write(*,*) 'nobi_max   =', nobi_max
      write(*,*) 'nobj_max   =', nobj_max
      write(*,*) 'nobk_max   =', nobk_max
      write(*,*) '========================================================='
    end if

    ! write nobi fields
    write(fname,*) i0; 
    fname = trim(adjustl('nobi'))//trim(adjustl(fname))
    if (ims_pro .eq. 0) then  
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, i1,   jmax, kmax, i1, jmax*kmax,   nobi, wrk3d)
    

    ! write nobj fields
    write(fname,*) i0; 
    fname = trim(adjustl('nobj'))//trim(adjustl(fname))
    if (ims_pro .eq. 0) then  
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, imax,i1,   kmax, i1, imax*kmax,   nobj, wrk3d)
    
    ! write nobk fields
    write(fname,*) i0; 
    fname = trim(adjustl('nobk'))//trim(adjustl(fname))
    if (ims_pro .eq. 0) then
      write(*,*) fname ! debugging
    end if
    call DNS_WRITE_FIELDS(fname, i2, imax,jmax,i1, i1, imax*jmax, nobk, wrk3d)

    return
  end subroutine GENERATE_GEOMETRY
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

    ! write eps field
    ! m = m + 1
    ! WRITE(fname,*) m; 
    ! fname = trim(adjustl('eps'))//trim(adjustl(fname))
    ! write(*,*) fname ! debugging
    ! CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, 1, imax*jmax*kmax, eps, wrk3d)
    
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
! #endif 
end module DNS_IBM
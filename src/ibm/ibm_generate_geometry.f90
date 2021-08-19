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

subroutine IBM_GENERATE_GEOMETRY(wrk3d,txc)
  
  use DNS_IBM
  use TLAB_VARS,        only: g   
  use TLAB_VARS,        only: imax, jmax, kmax 
  use TLAB_VARS,        only: isize_field, inb_txc
  
#ifdef USE_MPI
  use TLAB_MPI_VARS,    only: ims_pro, ims_npro
  use TLAB_MPI_VARS,    only: ims_size_i, ims_size_j, ims_size_k    
  use TLAB_MPI_VARS,    only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i 
  use TLAB_MPI_VARS,    only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use TLAB_MPI_VARS,    only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
  use TLAB_MPI_PROCS
#endif    
  
  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
  TINTEGER, parameter                                :: idi        = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter                                :: idj        = TLAB_MPI_J_PARTIAL 
  TINTEGER, parameter                                :: idk        = TLAB_MPI_K_PARTIAL 
  TINTEGER, parameter                                :: idi_nob    = TLAB_MPI_I_IBM_NOB 
  TINTEGER, parameter                                :: idk_nob    = TLAB_MPI_K_IBM_NOB 
  TINTEGER, parameter                                :: idi_nob_be = TLAB_MPI_I_IBM_NOB_BE 
  TINTEGER, parameter                                :: idk_nob_be = TLAB_MPI_K_IBM_NOB_BE
  TINTEGER                                           :: ims_err, max_nob
  ! debugging arrays    
  TREAL, dimension(ims_size_i(idi))                  :: nobi_out     ! DEBUG
  TREAL, dimension(ims_size_j(idj))                  :: nobj_out     ! DEBUG
  TREAL, dimension(ims_size_k(idk))                  :: nobk_out     ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo%number) :: nobi_b_out   ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo%number) :: nobj_b_out   ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo%number) :: nobk_b_out   ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo%number) :: nobi_e_out   ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo%number) :: nobj_e_out   ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo%number) :: nobk_e_out   ! DEBUG
  ! for type casting (tlab subroutines just with real)
  TREAL, dimension(ims_size_i(idi))                  :: nobi_real    ! DEBUG
  TREAL, dimension(ims_size_j(idj))                  :: nobj_real    ! DEBUG
  TREAL, dimension(ims_size_k(idk))                  :: nobk_real    ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo%number) :: nobi_b_real  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo%number) :: nobj_b_real  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo%number) :: nobk_b_real  ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo%number) :: nobi_e_real  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo%number) :: nobj_e_real  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo%number) :: nobk_e_real  ! DEBUG
#else
  TINTEGER, parameter                                :: ims_pro=0, ims_npro=1
  ! debugging arrays
  TREAL, dimension(jmax * kmax)                      :: nobi_out     ! DEBUG
  TREAL, dimension(imax * kmax)                      :: nobj_out     ! DEBUG
  TREAL, dimension(imax * jmax)                      :: nobk_out     ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo%number)   :: nobi_b_out   ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo%number)   :: nobj_b_out   ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo%number)   :: nobk_b_out   ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo%number)   :: nobi_e_out   ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo%number)   :: nobj_e_out   ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo%number)   :: nobk_e_out   ! DEBUG
  ! for type casting (tlab subroutines just with real)
  TREAL, dimension(jmax * kmax)                      :: nobi_real    ! DEBUG
  TREAL, dimension(imax * kmax)                      :: nobj_real    ! DEBUG
  TREAL, dimension(imax * jmax)                      :: nobk_real    ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo%number)   :: nobi_b_real  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo%number)   :: nobj_b_real  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo%number)   :: nobk_b_real  ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo%number)   :: nobi_e_real  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo%number)   :: nobj_e_real  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo%number)   :: nobk_e_real  ! DEBUG
#endif

  TREAL, dimension(isize_field), intent(inout)       :: wrk3d 
    
  TREAL                                              :: nobi_max, nobj_max, nobk_max
  TINTEGER                                           :: i, j, k, ij, ik, jk, ip, inum
  TINTEGER                                           :: nyz, nxz, nxy
  TINTEGER                                           :: nob_max
    
  CHARACTER(len=32)                                  :: fname
    
  ! DEBUG
  TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc   ! DEBUG
  TREAL, dimension(isize_field)                        :: tmp1, tmp2, tmp3, tmp4   ! DEBUG

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

  ! ================================================================== !
 
  ! max(nobi_max,nobj_max,nobk_max) from dns.ini file
  nob_max = xbars_geo%number

  ! initialize 
  nobi(:)   = i0; nobj(:)   = i0; nobk(:)   = i0
  nobi_b(:) = i0; nobj_b(:) = i0; nobk_b(:) = i0
  nobi_e(:) = i0; nobj_e(:) = i0; nobk_e(:) = i0
  nobi_max  = i0; nobj_max  = i0; nobk_max  = i0

  ! ================================================================== !
  ! number of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((ip == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for objects
        nobi(jk) = i1
      end if 
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! check for interface 
        nobi(jk) = nobi(jk) + i1
      end if
    end do
    ip = ip + nyz
  end do

  nobi_max = maxval(nobi)
#ifdef USE_MPI
  max_nob = nobi_max
  call MPI_ALLREDUCE(max_nob, nobi_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG

  nobi_real = dble(nobi)
  nobi_out  = nobi_real
  
  call DNS_TRANSPOSE(nobi_out, nyz, i1, nyz, nobi_real, i1)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(nobi_real, nobi_out, ims_ds_i(1,idi_nob), ims_dr_i(1,idi_nob), ims_ts_i(1,idi_nob), ims_tr_i(1,idi_nob))
  endif
#else 
  nobi_out = nobi_real
#endif

  ! ================================================================== !
  ! begin and end of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((i == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for interface
        nobi_b(jk) = i ! nobi_b
      end if
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then     ! nobi_b check for interface 
        inum = i0
        do while (nobi_b(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_b(inum+jk) = i + i1
      elseif((epsi(ip+jk-1) == C_1_R) .and. (epsi(ip+jk-1+nyz) == C_0_R)) then ! nobi_e check for interface 
        inum = i0
        do while (nobi_e(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = i        
      end if
      if((i == (g(1)%size - 1)) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobi_e(inum+jk) /= i0)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = g(1)%size    
      end if
    end do
    ip = ip + nyz
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG
  nobi_b_real = dble(nobi_b)
  nobi_e_real = dble(nobi_e)
  nobi_b_out  = nobi_b_real
  nobi_e_out  = nobi_e_real
  call DNS_TRANSPOSE(nobi_b_out, nyz, nob_max, nyz, nobi_b_real, nob_max)
  call DNS_TRANSPOSE(nobi_e_out, nyz, nob_max, nyz, nobi_e_real, nob_max)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(nobi_b_real, nobi_b_out, ims_ds_i(1,idi_nob_be), ims_dr_i(1,idi_nob_be), ims_ts_i(1,idi_nob_be), ims_tr_i(1,idi_nob_be))
    call TLAB_MPI_TRPB_I(nobi_e_real, nobi_e_out, ims_ds_i(1,idi_nob_be), ims_dr_i(1,idi_nob_be), ims_ts_i(1,idi_nob_be), ims_tr_i(1,idi_nob_be))
  endif
#else
  nobi_b_out = nobi_b_real
  nobi_e_out = nobi_e_real
#endif

! ================================================================== !
  ! number of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((ip == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for objects
        nobj(ik) = i1
      end if 
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! check for interface 
        nobj(ik) = nobj(ik) + i1
      end if
    end do
    ip = ip + nxz
  end do

  nobj_max = maxval(nobj)
#ifdef USE_MPI
  max_nob = nobj_max
  call MPI_ALLREDUCE(max_nob, nobj_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG
  nobj_real = dble(nobj)
  call DNS_TRANSPOSE(nobj_real, kmax, imax, kmax, nobj_out, imax)

  ! ================================================================== !
  ! begin and end of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((j == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for interface
        nobj_b(ik) = j ! nobj_b
      end if
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then     ! nobj_b check for interface 
        inum = i0
        do while (nobj_b(inum+ik) /= i0)
          inum = inum + nxz            
        end do 
        nobj_b(inum+ik) = j + i1
      elseif((epsj(ip+ik-1) == C_1_R) .and. (epsj(ip+ik-1+nxz) == C_0_R)) then ! nobj_e check for interface 
        inum = i0
        do while (nobj_e(inum+ik) /= i0)
          inum = inum + nxz            
        end do 
        nobj_e(inum+ik) = j        
      end if
      if((j == (g(2)%size - 1)) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobj_e(inum+ik) /= i0)
          inum = inum + nxz           
        end do 
        nobj_e(inum+ik) = g(2)%size    
      end if
    end do
    ip = ip + nxz
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG
  nobj_b_real = dble(nobj_b)
  nobj_e_real = dble(nobj_e)
  call DNS_TRANSPOSE(nobj_b_real, kmax, imax * nob_max, kmax, nobj_b_out, imax * nob_max)
  call DNS_TRANSPOSE(nobj_e_real, kmax, imax * nob_max, kmax, nobj_e_out, imax * nob_max)
  
  ! ================================================================== !
  ! number of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
        nobk(ij) = i1
      end if 
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
        nobk(ij) = nobk(ij) + i1
      end if
    end do
    ip = ip + nxy
  end do

  nobk_max = maxval(nobk)
#ifdef USE_MPI
  max_nob = nobk_max
  call MPI_ALLREDUCE(max_nob, nobk_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG
  nobk_real = dble(nobk)
#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(nobk_real, nobk_out, ims_ds_k(1,idk_nob), ims_dr_k(1,idk_nob), ims_ts_k(1,idk_nob), ims_tr_k(1,idk_nob))
  endif
#else 
  nobk_out = nobk_real
#endif

#ifdef IBM_DEBUG
  if (ims_pro == 0) then
    write(*,*) 'nobk_b     =', size(nobk_b)
    write(*,*) 'nobk_b_out =', size(nobk_b_out)
    write(*,*) 'nobk_e     =', size(nobk_e)
    write(*,*) 'nobk_e_out =', size(nobk_e_out)
  end if
#endif

  ! ! ================================================================== !
  ! ! begin and end of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
        nobk_b(ij) = k ! nobj_b
      end if
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
        inum = i0
        do while (nobk_b(inum+ij) /= i0)
          inum = inum + nxy            
        end do 
        nobk_b(inum+ij) = k + i1
      elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
        inum = i0
        do while (nobk_e(inum+ij) /= i0)
          inum = inum + nxy            
        end do 
        nobk_e(inum+ij) = k        
      end if
      if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobk_e(inum+ij) /= i0)
          inum = inum + nxy           
        end do 
        nobk_e(inum+ij) = g(3)%size    
      end if
    end do
    ip = ip + nxy
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG
  nobk_b_real = dble(nobk_b)
  nobk_e_real = dble(nobk_e)
#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(nobk_b_real, nobk_b_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
    call TLAB_MPI_TRPB_K(nobk_e_real, nobk_e_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
  endif
#else 
  nobk_b_out = nobk_b_real
  nobk_e_out = nobk_e_real
#endif

  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !
  ! DEBUG nob fields with  2D fields
  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !
#ifdef IBM_DEBUG
  if (ims_pro == 0) then
    write(*,*) '========================================================='
    write(*,*) 'nbars      =', xbars_geo%number
    write(*,*) 'nobi_max   =', int(nobi_max)
    write(*,*) 'nobj_max   =', int(nobj_max)
    write(*,*) 'nobk_max   =', int(nobk_max)
    write(*,*) '======== Writing geometry ==============================='
  end if
#endif

  ! ------------------------------------------------------------------ !

  ! write nobi fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, i1, g(2)%size/ims_npro, g(3)%size, i1, nyz, nobi_out, wrk3d)

  ! write nobi_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi_b'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, nob_max, g(2)%size/ims_npro, g(3)%size, i1, nyz*nob_max, nobi_b_out, wrk3d)     

  ! write nobi_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi_e'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, nob_max, g(2)%size/ims_npro, g(3)%size, i1, nyz*nob_max, nobi_e_out, wrk3d)     
  
  ! ------------------------------------------------------------------ !

  ! write nobj fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, i1, g(3)%size, i1, nxz, nobj_out, wrk3d)

  ! write nobj_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj_b'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, nob_max, g(3)%size, i1, nxz*nob_max, nobj_b_out, wrk3d)

  ! write nobj_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj_e'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, nob_max, g(3)%size, i1, nxz*nob_max, nobj_e_out, wrk3d)
  
  ! ------------------------------------------------------------------ !

  ! write nobk fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, i1, i1, nxy, nobk_out, wrk3d)

  ! write nobk_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk_b'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, nob_max, i1, nxy*nob_max, nobk_b_out, wrk3d)

  ! write nobk_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk_e'))//trim(adjustl(fname))
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) fname
#endif
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, nob_max, i1, nxy*nob_max, nobk_e_out, wrk3d)

  ! ================================================================== !
  ! ================================================================== !
  ! Block comment: begin
#if 0
  ! ================================================================== !
  ! ================================================================== !

  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !
  ! DEBUG nob fields with 3D fields
  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !

  if (ims_pro == 0) then
    write(*,*) '======== Writing geometry in 3D fields =================='
  end if

  txc(:,:) = C_0_R

  ! tmp aux 
  tmp1(:) = txc(:,1) 
  tmp2(:) = txc(:,2)
  tmp3(:) = txc(:,3) 
  tmp4(:) = txc(:,4) 

  ! ================================================================== !
  ! number of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((ip == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for objects
        tmp1(jk) = C_1_R
      end if 
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! check for interface 
        tmp1(jk) = tmp1(jk) + C_1_R
      end if
    end do
    ip = ip + nyz
  end do

  call DNS_TRANSPOSE(tmp1, nyz, g(1)%size, nyz,        tmp2, g(1)%size)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(tmp2, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  endif
  call DNS_WRITE_FIELDS('nobi3d', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp1, wrk3d)
#else
  call DNS_WRITE_FIELDS('nobi3d', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
#endif
  tmp1(:) = C_0_R
  tmp2(:) = C_0_R

  ! ================================================================== !
  ! number of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((ip == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for objects
        tmp1(ik) = C_1_R
      end if 
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! check for interface 
        tmp1(ik) = tmp1(ik) + C_1_R
      end if
    end do
    ip = ip + nxz
  end do

  call DNS_TRANSPOSE(tmp1, kmax, imax * jmax, kmax, tmp2, imax * jmax)
  call DNS_WRITE_FIELDS('nobj3d', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
  tmp1(:) = C_0_R
  tmp2(:) = C_0_R

  ! ================================================================== !
  ! number of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
        tmp1(ij) = C_1_R
      end if 
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
        tmp1(ij) = tmp1(ij) + C_1_R
      end if
    end do
    ip = ip + nxy
  end do

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(tmp1, tmp2, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  endif
  call DNS_WRITE_FIELDS('nobk3d', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
#else
  call DNS_WRITE_FIELDS('nobk3d', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp1, wrk3d)
#endif

  tmp1(:) = C_0_R
  tmp2(:) = C_0_R

  if (ims_pro == 0) write(*,*) 'done writing nobi3d,   nobj3d,   nobk3d'
  ! ================================================================== !
  ! ================================================================== !
  ! begin and end of objects in x-direction

  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((i == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for interface
        tmp1(jk) = dble(i) ! nobi_b
      end if
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then     ! nobi_b check for interface 
        inum = i0
        do while (tmp1(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp1(inum+jk) = dble(i + 1)
      elseif((epsi(ip+jk-1) == C_1_R) .and. (epsi(ip+jk-1+nyz) == C_0_R)) then ! nobi_e check for interface 
        inum = i0
        do while (tmp2(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp2(inum+jk) = dble(i)        
      end if
      if((i == (g(1)%size - 1)) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        tmp2(inum+jk) = dble(g(1)%size)    
      end if
    end do
    ip = ip + nyz
  end do

  call DNS_TRANSPOSE(tmp1, nyz, g(1)%size, nyz,        tmp3, g(1)%size)
  call DNS_TRANSPOSE(tmp2, nyz, g(1)%size, nyz,        tmp4, g(1)%size)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLAB_MPI_TRPB_I(tmp3, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
    call TLAB_MPI_TRPB_I(tmp4, tmp2, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  endif
  call DNS_WRITE_FIELDS('nobi3d_b', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp1, wrk3d)
  call DNS_WRITE_FIELDS('nobi3d_e', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
#else
  call DNS_WRITE_FIELDS('nobi3d_b', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp3, wrk3d)
  call DNS_WRITE_FIELDS('nobi3d_e', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp4, wrk3d)
#endif

  tmp1(:) = C_0_R
  tmp2(:) = C_0_R
  tmp3(:) = C_0_R
  tmp4(:) = C_0_R
  if (ims_pro == 0) write(*,*) 'done writing nobi3d_b, nobi3d_e'
  ! ================================================================== !
  ! begin and end of objects in y-direction

  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((j == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for interface
        tmp1(ik) = dble(j) ! nobj_b
      end if
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then     ! nobj_b check for interface 
        inum = i0
        do while (tmp1(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        tmp1(inum+ik) = dble(j + 1)
      elseif((epsj(ip+ik-1) == C_1_R) .and. (epsj(ip+ik-1+nxz) == C_0_R)) then ! nobj_e check for interface 
        inum = i0
        do while (tmp2(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        tmp2(inum+ik) = dble(j)        
      end if
      if((j == (g(2)%size - 1)) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+ik) /= C_0_R)
          inum = inum + nxz           
        end do 
        tmp2(inum+ik) = dble(g(2)%size)    
      end if
    end do
    ip = ip + nxz
  end do

  call DNS_TRANSPOSE(tmp1, kmax, imax * jmax, kmax, tmp3, imax * jmax)
  call DNS_TRANSPOSE(tmp2, kmax, imax * jmax, kmax, tmp4, imax * jmax)

  call DNS_WRITE_FIELDS('nobj3d_b', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp3, wrk3d)
  call DNS_WRITE_FIELDS('nobj3d_e', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp4, wrk3d)

  tmp1(:) = C_0_R
  tmp2(:) = C_0_R
  tmp3(:) = C_0_R
  tmp4(:) = C_0_R  
  if (ims_pro == 0) write(*,*) 'done writing nobj3d_b, nobj3d_e'

  ! ================================================================== !
  ! begin and end of objects in z-direction

  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
        tmp1(ij) = dble(k)    ! nobk_b
      end if
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
        inum = i0
        do while (tmp1(ij+inum) /= C_0_R)
          inum = inum + nxy
        end do
        tmp1(ij+inum) = dble(k + 1)
      elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
        inum = i0
        do while (tmp2(ij+inum) /= C_0_R)
          inum = inum + nxy
        end do
        tmp2(ij+inum) = dble(k)
      end if
      if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (tmp2(inum+ij) /= C_0_R)
          inum = inum + nxy           
        end do 
        tmp2(inum+ij) = dble(g(3)%size)    
      end if
    end do
    ip = ip + nxy
  end do

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLAB_MPI_TRPB_K(tmp1, tmp3, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    call TLAB_MPI_TRPB_K(tmp2, tmp4, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  endif
  call DNS_WRITE_FIELDS('nobk3d_b', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp3, wrk3d)
  call DNS_WRITE_FIELDS('nobk3d_e', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp4, wrk3d)
#else
  call DNS_WRITE_FIELDS('nobk3d_b', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp1, wrk3d)
  call DNS_WRITE_FIELDS('nobk3d_e', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
#endif

  tmp1(:) = C_0_R
  tmp2(:) = C_0_R
  tmp3(:) = C_0_R
  tmp4(:) = C_0_R  

  if (ims_pro == 0) write(*,*) 'done writing nobk3d_b, nobk3d_e'
  if (ims_pro == 0) write(*,*) '========================================================='

  ! ================================================================== !
  ! ================================================================== !
  ! Block comment: end
#endif
  ! ================================================================== !
  ! ================================================================== !

  return
end subroutine IBM_GENERATE_GEOMETRY

!########################################################################

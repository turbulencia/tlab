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
  use DNS_GLOBAL, only: g   
  use DNS_GLOBAL, only: imax, jmax, kmax 
  use DNS_GLOBAL, only: isize_field, inb_txc !,isize_txc_field
  
#ifdef USE_MPI
  use DNS_MPI,    only: ims_pro, ims_npro
  use DNS_MPI,    only: ims_size_i, ims_size_j, ims_size_k    
  use DNS_MPI,    only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i 
  use DNS_MPI,    only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use DNS_MPI,    only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
#endif    
  
  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"  
  TINTEGER, parameter                            :: idi        = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter                            :: idj        = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter                            :: idk        = DNS_MPI_K_PARTIAL 
  TINTEGER, parameter                            :: idi_nob    = DNS_MPI_I_IBM_NOB 
  TINTEGER, parameter                            :: idk_nob    = DNS_MPI_K_IBM_NOB 
  TINTEGER, parameter                            :: idi_nob_be = DNS_MPI_I_IBM_NOB_BE 
  TINTEGER, parameter                            :: idk_nob_be = DNS_MPI_K_IBM_NOB_BE
  TINTEGER                                       :: ims_err, max_nob
  TREAL, dimension(ims_size_i(idi))              :: nobi_out    ! DEBUG
  TREAL, dimension(ims_size_j(idj))              :: nobj_out    ! DEBUG
  TREAL, dimension(ims_size_k(idk))              :: nobk_out    ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo(1)) :: nobi_b_out  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo(1)) :: nobj_b_out  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo(1)) :: nobk_b_out  ! DEBUG
  TREAL, dimension(ims_size_i(idi)*xbars_geo(1)) :: nobi_e_out  ! DEBUG
  TREAL, dimension(ims_size_j(idj)*xbars_geo(1)) :: nobj_e_out  ! DEBUG
  TREAL, dimension(ims_size_k(idk)*xbars_geo(1)) :: nobk_e_out  ! DEBUG
#else
  TINTEGER, parameter                            :: ims_pro=0, ims_npro=1
  TREAL, dimension(jmax * kmax)                  :: nobi_out    ! DEBUG
  TREAL, dimension(imax * kmax)                  :: nobj_out    ! DEBUG
  TREAL, dimension(imax * jmax)                  :: nobk_out    ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo(1))   :: nobi_b_out  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo(1))   :: nobj_b_out  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo(1))   :: nobk_b_out  ! DEBUG
  TREAL, dimension(jmax * kmax * xbars_geo(1))   :: nobi_e_out  ! DEBUG
  TREAL, dimension(imax * kmax * xbars_geo(1))   :: nobj_e_out  ! DEBUG
  TREAL, dimension(imax * jmax * xbars_geo(1))   :: nobk_e_out  ! DEBUG
#endif

  TREAL, dimension(isize_field), intent(inout)   :: wrk3d 

  TINTEGER                                       :: nobi_max, nobj_max, nobk_max
  TINTEGER                                       :: i, j, k, ij, ik, jk, ip, inum
  TINTEGER                                       :: nyz, nxz, nxy
  TINTEGER                                       :: nob_max

  CHARACTER(32)                                  :: fname

  ! DEBUG
  TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc 
  TREAL, dimension(isize_field)                        :: tmp1, tmp2, tmp3, tmp4 

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
  nob_max = xbars_geo(1)

  ! initialize 
  nobi(:)    = C_0_R; nobj(:)    = C_0_R; nobk(:)    = C_0_R
  nobi_b(:)  = C_0_R; nobj_b(:)  = C_0_R; nobk_b(:)  = C_0_R
  nobi_e(:)  = C_0_R; nobj_e(:)  = C_0_R; nobk_e(:)  = C_0_R
  nobi_max   = i0;    nobj_max   = i0;    nobk_max   = i0

  ! ================================================================== !
  ! number of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((ip == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for objects
        nobi(jk) = C_1_R
      end if 
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! check for interface 
        nobi(jk) = nobi(jk) + C_1_R
      end if
    end do
    ip = ip + nyz
  end do

  nobi_max = int(maxval(nobi))
#ifdef USE_MPI
  max_nob = nobi_max
  call MPI_ALLREDUCE(max_nob, nobi_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG

  nobi_out = nobi
  
  call DNS_TRANSPOSE(nobi_out, nyz, i1, nyz, nobi, i1)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call DNS_MPI_TRPB_I(nobi, nobi_out, ims_ds_i(1,idi_nob), ims_dr_i(1,idi_nob), ims_ts_i(1,idi_nob), ims_tr_i(1,idi_nob))
  endif
#else 

  nobi_out = nobi
#endif

  ! ================================================================== !
  ! begin and end of objects in x-direction
  ip = i1
  do i = 1, g(1)%size - 1     ! contiguous i-lines
    do jk = 1, nyz            ! pages of   i-lines
      if((i == 1) .and. (epsi(jk) == C_1_R)) then ! exception: check first plane for interface
        nobi_b(jk) = dble(i) ! nobi_b
      end if
      if((epsi(ip+jk-1) == C_0_R) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then     ! nobi_b check for interface 
        inum = i0
        do while (nobi_b(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        nobi_b(inum+jk) = dble(i)
      elseif((epsi(ip+jk-1) == C_1_R) .and. (epsi(ip+jk-1+nyz) == C_0_R)) then ! nobi_e check for interface 
        inum = i0
        do while (nobi_e(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = dble(i)        
      end if
      if((i == (g(1)%size - 1)) .and. (epsi(ip+jk-1+nyz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobi_e(inum+jk) /= C_0_R)
          inum = inum + nyz            
        end do 
        nobi_e(inum+jk) = dble(g(1)%size)    
      end if
    end do
    ip = ip + nyz
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG

  nobi_b_out = nobi_b
  nobi_e_out = nobi_e
  call DNS_TRANSPOSE(nobi_b_out, nyz, nob_max, nyz, nobi_b, nob_max)
  call DNS_TRANSPOSE(nobi_e_out, nyz, nob_max, nyz, nobi_e, nob_max)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call DNS_MPI_TRPB_I(nobi_b, nobi_b_out, ims_ds_i(1,idi_nob_be), ims_dr_i(1,idi_nob_be), ims_ts_i(1,idi_nob_be), ims_tr_i(1,idi_nob_be))
    call DNS_MPI_TRPB_I(nobi_e, nobi_e_out, ims_ds_i(1,idi_nob_be), ims_dr_i(1,idi_nob_be), ims_ts_i(1,idi_nob_be), ims_tr_i(1,idi_nob_be))
  endif
#else
  nobi_b_out = nobi_b
  nobi_e_out = nobi_e
#endif

! ================================================================== !
  ! number of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((ip == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for objects
        nobj(ik) = C_1_R
      end if 
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! check for interface 
        nobj(ik) = nobj(ik) + C_1_R
      end if
    end do
    ip = ip + nxz
  end do

  nobj_max = int(maxval(nobj))
#ifdef USE_MPI
  max_nob = nobj_max
  call MPI_ALLREDUCE(max_nob, nobj_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG

  call DNS_TRANSPOSE(nobj, kmax, imax, kmax, nobj_out, imax)

  ! ================================================================== !
  ! begin and end of objects in y-direction
  ip = i1
  do j = 1, g(2)%size - 1     ! contiguous j-lines
    do ik = 1, nxz            ! pages of   j-lines
      if((j == 1) .and. (epsj(ik) == C_1_R)) then ! exception: check first plane for interface
        nobj_b(ik) = dble(j) ! nobj_b
      end if
      if((epsj(ip+ik-1) == C_0_R) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then     ! nobj_b check for interface 
        inum = i0
        do while (nobj_b(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        nobj_b(inum+ik) = dble(j)
      elseif((epsj(ip+ik-1) == C_1_R) .and. (epsj(ip+ik-1+nxz) == C_0_R)) then ! nobj_e check for interface 
        inum = i0
        do while (nobj_e(inum+ik) /= C_0_R)
          inum = inum + nxz            
        end do 
        nobj_e(inum+ik) = dble(j)        
      end if
      if((j == (g(2)%size - 1)) .and. (epsj(ip+ik-1+nxz) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobj_e(inum+ik) /= C_0_R)
          inum = inum + nxz           
        end do 
        nobj_e(inum+ik) = dble(g(2)%size)    
      end if
    end do
    ip = ip + nxz
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG

  call DNS_TRANSPOSE(nobj_b, kmax, imax * nob_max, kmax, nobj_b_out, imax * nob_max)
  call DNS_TRANSPOSE(nobj_e, kmax, imax * nob_max, kmax, nobj_e_out, imax * nob_max)
  
  ! ================================================================== !
  ! number of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((ip == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for objects
        nobk(ij) = C_1_R
      end if 
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! check for interface 
        nobk(ij) = nobk(ij) + C_1_R
      end if
    end do
    ip = ip + nxy
  end do

  nobk_max = int(maxval(nobk))
#ifdef USE_MPI
  max_nob = nobk_max
  call MPI_ALLREDUCE(max_nob, nobk_max, i0, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

  ! ------------------------------------------------------------------ !
  ! DEBUG

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call DNS_MPI_TRPB_K(nobk, nobk_out, ims_ds_k(1,idk_nob), ims_dr_k(1,idk_nob), ims_ts_k(1,idk_nob), ims_tr_k(1,idk_nob))
  endif
#else 
  nobk_out = nobk
#endif

  if (ims_pro == 0) then
    write(*,*) 'nobk_b     =', size(nobk_b)
    write(*,*) 'nobk_b_out =', size(nobk_b_out)
    write(*,*) 'nobk_e     =', size(nobk_e)
    write(*,*) 'nobk_e_out =', size(nobk_e_out)
  end if

  ! ! ================================================================== !
  ! ! begin and end of objects in z-direction
  ip = i1
  do k = 1, g(3)%size - 1     ! contiguous k-lines
    do ij = 1, nxy            ! pages of   k-lines
      if((k == 1) .and. (epsk(ij) == C_1_R)) then ! exception: check first plane for interface
        nobk_b(ij) = dble(k) ! nobj_b
      end if
      if((epsk(ip+ij-1) == C_0_R) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then     ! nobk_b check for interface 
        inum = i0
        do while (nobk_b(inum+ij) /= C_0_R)
          inum = inum + nxy            
        end do 
        nobk_b(inum+ij) = dble(k)
      elseif((epsk(ip+ij-1) == C_1_R) .and. (epsk(ip+ij-1+nxy) == C_0_R)) then ! nobk_e check for interface 
        inum = i0
        do while (nobk_e(inum+ij) /= C_0_R)
          inum = inum + nxy            
        end do 
        nobk_e(inum+ij) = dble(k)        
      end if
      if((k == (g(3)%size - 1)) .and. (epsk(ip+ij-1+nxy) == C_1_R)) then ! exception: check last plane for interface
        inum = i0
        do while (nobk_e(inum+ij) /= C_0_R)
          inum = inum + nxy           
        end do 
        nobk_e(inum+ij) = dble(g(3)%size)    
      end if
    end do
    ip = ip + nxy
  end do

  ! ------------------------------------------------------------------ !
  ! DEBUG

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call DNS_MPI_TRPB_K(nobk_b, nobk_b_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
    call DNS_MPI_TRPB_K(nobk_e, nobk_e_out, ims_ds_k(1,idk_nob_be), ims_dr_k(1,idk_nob_be), ims_ts_k(1,idk_nob_be), ims_tr_k(1,idk_nob_be))
  endif
#else 
  nobk_b_out = nobk_b
  nobk_e_out = nobk_e
#endif

  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !
  ! DEBUG nob fields with  2D fields
  ! ================================================================== !
  ! ================================================================== !
  ! ================================================================== !
  if (ims_pro == 0) then
    write(*,*) '========================================================='
    write(*,*) 'nbars      =', xbars_geo(1)
    write(*,*) 'nobi_max   =', nobi_max
    write(*,*) 'nobj_max   =', nobj_max
    write(*,*) 'nobk_max   =', nobk_max
    write(*,*) '======== Writing geometry ==============================='
  end if

  ! ------------------------------------------------------------------ !

  ! write nobi fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, i1, g(2)%size/ims_npro, g(3)%size, i1, nyz, nobi_out, wrk3d)

  ! write nobi_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi_b'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, nob_max, g(2)%size/ims_npro, g(3)%size, i1, nyz*nob_max, nobi_b_out, wrk3d)     

  ! write nobi_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobi_e'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, nob_max, g(2)%size/ims_npro, g(3)%size, i1, nyz*nob_max, nobi_e_out, wrk3d)     
  
  ! ------------------------------------------------------------------ !

  ! write nobj fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, i1, g(3)%size, i1, nxz, nobj_out, wrk3d)

  ! write nobj_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj_b'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, nob_max, g(3)%size, i1, nxz*nob_max, nobj_b_out, wrk3d)

  ! write nobj_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobj_e'))//trim(adjustl(fname))
  if (ims_pro == 0) then  
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, nob_max, g(3)%size, i1, nxz*nob_max, nobj_e_out, wrk3d)
  
  ! ------------------------------------------------------------------ !

  ! write nobk fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk'))//trim(adjustl(fname))
  if (ims_pro == 0) then
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, i1, i1, nxy, nobk_out, wrk3d)

  ! write nobk_b fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk_b'))//trim(adjustl(fname))
  if (ims_pro == 0) then
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, nob_max, i1, nxy*nob_max, nobk_b_out, wrk3d)

  ! write nobk_e fields
  write(fname,*) i0; 
  fname = trim(adjustl('nobk_e'))//trim(adjustl(fname))
  if (ims_pro == 0) then
    write(*,*) fname ! DEBUG
  end if
  call DNS_WRITE_FIELDS(fname, i2, g(1)%size/ims_npro, g(2)%size, nob_max, i1, nxy*nob_max, nobk_e_out, wrk3d)

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
    call DNS_MPI_TRPB_I(tmp2, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
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
    call DNS_MPI_TRPB_K(tmp1, tmp2, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
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
        tmp1(inum+jk) = dble(i)
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
    call DNS_MPI_TRPB_I(tmp3, tmp1, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
    call DNS_MPI_TRPB_I(tmp4, tmp2, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
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
        tmp1(inum+ik) = dble(j)
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
        tmp1(ij+inum) = dble(k)
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
    call DNS_MPI_TRPB_K(tmp1, tmp3, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
    call DNS_MPI_TRPB_K(tmp2, tmp4, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
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

  return
end subroutine IBM_GENERATE_GEOMETRY

!########################################################################

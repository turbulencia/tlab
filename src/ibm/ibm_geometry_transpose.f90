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

subroutine IBM_GEOMETRY_TRANSPOSE(wrk3d,txc)
  
  use DNS_IBM
  use TLAB_VARS,        only: g
  use TLAB_VARS,        only: imax, jmax, kmax 
  use TLAB_VARS,        only: isize_field, isize_txc_field, inb_txc 
  
#ifdef USE_MPI
  use TLAB_MPI_VARS,    only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  use TLAB_MPI_VARS,    only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use TLAB_MPI_VARS,    only: ims_npro_i, ims_npro_j, ims_npro_k, ims_pro
  use TLAB_MPI_VARS,    only: ims_size_i, ims_size_j, ims_size_k 
  use TLAB_MPI_PROCS
#endif

  implicit none
  
#include "integers.h"
  
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
  TINTEGER, parameter                                  :: idi = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter                                  :: idj = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter                                  :: idk = DNS_MPI_K_PARTIAL 
#else
  TINTEGER, parameter                                  :: ims_pro = 0         
#endif

  TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc 
  TREAL, dimension(isize_field)                        :: tmp1 
  TREAL, dimension(isize_field),         intent(inout) :: wrk3d ! debug 

  TINTEGER                                             :: nyz, nxy

  ! DEBUG
  ! TREAL, dimension(isize_field)                        :: tmp2, tmp3, tmp4

  ! ================================================================== !
  
  ! tmp aux 
  tmp1(:) = txc(:,1)

  ! npages for dns_transpose function (cf. dns_mpi_initialize.f90)
#ifdef USE_MPI
  nyz = ims_size_i(idi) 
#else
  nyz = jmax * kmax   
#endif
  nxy = imax * jmax

  ! ================================================================== !
  ! DEBUG
#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) then 
    write(*,*) '======== Debug transposing geometry ====================='
#ifdef USE_MPI
    write(*,*) 'ims_size_i(1)    ', ims_size_i(idi)
    write(*,*) 'ims_size_j(1)    ', ims_size_j(idj)
    write(*,*) 'ims_size_k(1)    ', ims_size_j(idk)
#endif      
    write(*,*) 'nyz              ', jmax * kmax
    write(*,*) 'nxz              ', imax * kmax
    write(*,*) 'nxy              ', imax * jmax 
    write(*,*) 'isize_field      ', isize_field
    write(*,*) 'isize_txc_field  ', isize_txc_field
    write(*,*) 'inb_txc          ', inb_txc
  end if
#endif
  ! ================================================================== !

  ! MPI transposition and local transposition in x (make x-direction the last one)
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
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

  ! ------------------------------------------------------------------ !

  ! local transposition in y (make y-direction the last one, native standing in y, no MPI transp)
#ifdef USE_ESSL
  call DGETMO       (eps, nxy, nxy, kmax, epsj, kmax)
#else
  call DNS_TRANSPOSE(eps, nxy, kmax, nxy, epsj, kmax)
#endif

  ! ------------------------------------------------------------------ !

  ! MPI transposition in z (no local transposition needed, already in right order)
#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call DNS_MPI_TRPF_K(eps, epsk, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  else
#endif
  epsk = eps
#ifdef USE_MPI
  end if
#endif

  ! ================================================================== !
  ! ================================================================== !
  ! DEBUG transposed fields of eps
  ! ================================================================== !
  ! ================================================================== !
#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) then 
    write(*,*) '========================================================='
    write(*,*) 'Transposing arrays: done'
    write(*,*) '======== Generating geometry informations ==============='
  end if
  
  ! ! tmp
  ! tmp2(:) = txc(:,2);  tmp3(:) = txc(:,3); tmp4(:) = txc(:,4)

  ! ! check transpositions in x! 
  ! call DNS_TRANSPOSE(epsi, nyz, g(1)%size, nyz, tmp2, g(1)%size)
  ! call DNS_MPI_TRPB_I(tmp2, tmp3, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
  ! call DNS_WRITE_FIELDS('epsi0', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp3, wrk3d)
  
  ! ! check epsj
  ! call DNS_TRANSPOSE(epsj, kmax, nxy, kmax, tmp2, nxy)
  ! call DNS_WRITE_FIELDS('epsj0', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)

  ! ! check transpositions in z! 
  ! call DNS_MPI_TRPB_K(epsk, tmp2, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  ! call DNS_WRITE_FIELDS('epsk0', i2, imax,jmax,kmax, i1, imax*jmax*kmax, tmp2, wrk3d)
#endif
  ! ================================================================== !
  ! ================================================================== !
  
return
end subroutine IBM_GEOMETRY_TRANSPOSE

!########################################################################
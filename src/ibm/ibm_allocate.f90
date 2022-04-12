#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 l.210 to allocate needed memory for ibm
!#   module
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

subroutine IBM_ALLOCATE(C_FILE_LOC, allocated)

  use DNS_IBM
  use TLAB_CONSTANTS, only : lfile, efile
  use TLAB_VARS,      only : g, isize_field, istagger
  use TLAB_PROCS  
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_size_i, ims_size_j, ims_size_k 
#else
  use TLAB_VARS,      only : imax, jmax, kmax
#endif    

  implicit none

#include "integers.h"

  character(len=128), intent(in   ) :: C_FILE_LOC
  logical,            intent(inout) :: allocated

#ifdef USE_MPI 
  TINTEGER, parameter               :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter               :: idj = TLAB_MPI_J_PARTIAL 
  TINTEGER, parameter               :: idk = TLAB_MPI_K_PARTIAL 
#endif
  TINTEGER                          :: ierr, inb_ibm
  TINTEGER                          :: nyz, nxz, nxy    
  character(len=128)                :: str, line

  ! ================================================================== !
  ! npages
#ifdef USE_MPI 
  nyz = ims_size_i(idi)
  nxz = ims_size_j(idj) 
  nxy = ims_size_k(idk) 
#else
  nyz = jmax * kmax  
  nxz = imax * kmax     
  nxy = imax * jmax
#endif  

  ! ================================================================== !
  ! array sizes
  inb_ibm         = i1 
  isize_nobi      = nyz    
  isize_nobj      = nxz  
  isize_nobk      = nxy
  isize_nobi_be   = nyz * nob_max    
  isize_nobj_be   = nxz * nob_max  
  isize_nobk_be   = nxy * nob_max
  nspl            = 2 * nflu + 2    ! data points (incl. 2 interface points)
  isize_wrk1d_ibm = max(g(1)%size, max(g(2)%size, g(3)%size)) ! gap size unknown (max size assumed)

  ! ================================================================== !
  ! allocate here all ibm related arrays
  if ( .not. allocated ) then 
    ! eps
    write(str,*) inb_ibm; line = 'Allocating array IBM eps of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(eps(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  eps.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! epsp
    if ( istagger == 1 ) then
      write(str,*) inb_ibm; line = 'Allocating array IBM epsp of size '//trim(adjustl(str))//'x'
      write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
      call TLAB_WRITE_ASCII(lfile,line)
      allocate(epsp(isize_field), stat=ierr)
      if ( ierr /= 0 ) then
      call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  epsp.')
      call TLAB_STOP(DNS_ERROR_ALLOC)
      end if
    end if
    ! ------------------------------------------------------------------ !
    ! fld_ibm
    write(str,*) inb_ibm; line = 'Allocating array IBM fld_ibm of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(fld_ibm(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  fld_ibm.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    ! ------------------------------------------------------------------ !
    ! nobi
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobi(isize_nobi), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobi.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobj(isize_nobj), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobj.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobk(isize_nobk), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobk.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    ! ------------------------------------------------------------------ !
    ! nobi_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobi_b(isize_nobi_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobi_b.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobj_b(isize_nobj_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobj_b.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobk_b(isize_nobk_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobk_b.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    ! ------------------------------------------------------------------ !
    ! nobi_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobi_e(isize_nobi_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobi_e.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobj_e(isize_nobj_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobj_e.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk_be; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(nobk_e(isize_nobk_be), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  nobk_e.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    ! ------------------------------------------------------------------ !
    ! xa, ya spline arrays input
    write(str,*) inb_ibm; line = 'Allocating array IBM xa of size '//trim(adjustl(str))//'x'
    write(str,*) nspl; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(xa(nspl), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  xa.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    !
    write(str,*) inb_ibm; line = 'Allocating array IBM ya of size '//trim(adjustl(str))//'x'
    write(str,*) nspl; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(ya(nspl), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  ya.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    ! xb, yb spline arrays output
    write(str,*) inb_ibm; line = 'Allocating array IBM xb of size '//trim(adjustl(str))//'x'
    write(str,*) isize_wrk1d_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(xb(isize_wrk1d_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  xb.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    !
    write(str,*) inb_ibm; line = 'Allocating array IBM yb of size '//trim(adjustl(str))//'x'
    write(str,*) isize_wrk1d_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call TLAB_WRITE_ASCII(lfile,line)
    allocate(yb(isize_wrk1d_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call TLAB_WRITE_ASCII(efile,  C_FILE_LOC//'.  Error while allocating memory space for  yb.')
    call TLAB_STOP(DNS_ERROR_ALLOC)
    end if
    ! ------------------------------------------------------------------ !
    ! set alloc flag: done
    allocated = .true.
  end if

  return
end subroutine IBM_ALLOCATE

!########################################################################
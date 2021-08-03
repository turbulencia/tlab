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
!#   called once in dns_main.f90 l.210 to allocate needed memory
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

subroutine IBM_ALLOCATE(allocated)

  use DNS_IBM
  use DNS_CONSTANTS, only: lfile, efile
  use DNS_GLOBAL,    only: imax, jmax, kmax
  use DNS_GLOBAL,    only: g
  use DNS_GLOBAL,    only: isize_field

#ifdef USE_MPI
  use DNS_MPI,       only: ims_size_i, ims_size_j, ims_size_k 
#endif    

  implicit none

#include "integers.h"

#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"
  TINTEGER, parameter       :: idi = DNS_MPI_I_PARTIAL 
  TINTEGER, parameter       :: idj = DNS_MPI_J_PARTIAL 
  TINTEGER, parameter       :: idk = DNS_MPI_K_PARTIAL 
#endif

  logical, intent(inout)    :: allocated       ! flag, just allocate memory space once
  
  TINTEGER                  :: ierr, inb_ibm
  TINTEGER                  :: nyz, nxz, nxy
  TINTEGER                  :: nob_max
  
  character(128)            :: str, line

  ! ================================================================== !

  inb_ibm = i1 ! can be also defined in dns_read_local.f90 and module dns_global.f90 (cf. inb_flow ...)

  ! max(nobi_max,nobj_max,nobk_max) from dns.ini file
  nob_max = xbars_geo%number

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
    
  ! array sizes
  isize_nobi      = nyz    
  isize_nobj      = nxz  
  isize_nobk      = nxy
  isize_nobi_be   = nyz * nob_max    
  isize_nobj_be   = nxz * nob_max  
  isize_nobk_be   = nxy * nob_max
  !
  if (ibm_spline_global) then
    nsp           = max(g(1)%size, max(g(2)%size, g(3)%size)) 
  else
    nsp           = 2 * nflu + 2    ! number of data points (with 2 interface points) nsp > kspl
  endif
  ! cf. fitpack package (routines curfit and splev)
  nest            = nsp + kspl + 1
  isize_wrk_ibm   = nsp + 2 * nest + nsp * (kspl + 1) + nest * (7 + 3 * kspl)
  isize_iwrk_ibm  =                  nsp * (kspl + 1) + nest * (7 + 3 * kspl)
  isize_wrk1d_ibm = max(g(1)%size, max(g(2)%size, g(3)%size)) ! gap size unknown (max size assumed)

  ! ================================================================== !
    
  ! allocate here all ibm related arrays
  if ( .not. allocated ) then 

    ! eps_aux 3d-array
    write(str,*) inb_ibm; line = 'Allocating array IBM eps_aux of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(eps_aux(imax,jmax,kmax), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps_aux.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! eps
    write(str,*) inb_ibm; line = 'Allocating array IBM eps of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(eps(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for eps.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! epsi
    write(str,*) inb_ibm; line = 'Allocating array IBM epsi of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsi(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsi.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! epsj
    write(str,*) inb_ibm; line = 'Allocating array IBM epsj of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsj(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsj.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! epsk
    write(str,*) inb_ibm; line = 'Allocating array IBM epsk of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(epsk(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for epsk.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi(isize_nobi), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj(isize_nobj), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk(isize_nobk), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi_b(isize_nobi_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj_b(isize_nobj_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_b
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_b of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk_b(isize_nobk_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk_b.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! nobi_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobi_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobi_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobi_e(isize_nobi_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobi_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobj_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobj_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobj_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobj_e(isize_nobj_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobj_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! nobk_e
    write(str,*) inb_ibm; line = 'Allocating array IBM nobk_e of size '//trim(adjustl(str))//'x'
    write(str,*) isize_nobk_be; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(nobk_e(isize_nobk_be), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for nobk_e.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! fld_ibm
    write(str,*) inb_ibm; line = 'Allocating array IBM fld_ibm of size '//trim(adjustl(str))//'x'
    write(str,*) isize_field; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(fld_ibm(isize_field), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for fld_ibm.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! wrk_ibm [contains: w(nsp); t(nest); c(nest); wrk(nsp*(kspl+1)+nest*(7+3*kspl))], cf. fitpack
    write(str,*) inb_ibm; line = 'Allocating array IBM wrk_ibm of size '//trim(adjustl(str))//'x'
    write(str,*) isize_wrk_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(wrk_ibm(isize_wrk_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for wrk_ibm.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! iwrk_ibm [contains: iwrk(lwrk=m*(k+1)+nest*(7+3*k)) --> integer array]   
    write(str,*) inb_ibm; line = 'Allocating array IBM iwrk_ibm of size '//trim(adjustl(str))//'x'
    write(str,*) isize_iwrk_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(iwrk_ibm(isize_iwrk_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for iwrk_ibm.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! xa, ya spline arrays input
    write(str,*) inb_ibm; line = 'Allocating array IBM xa of size '//trim(adjustl(str))//'x'
    write(str,*) nsp; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(xa(nsp), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for xa.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if
    !
    write(str,*) inb_ibm; line = 'Allocating array IBM ya of size '//trim(adjustl(str))//'x'
    write(str,*) nsp; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(ya(nsp), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for ya.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! xb, yb spline arrays output
    write(str,*) inb_ibm; line = 'Allocating array IBM xb of size '//trim(adjustl(str))//'x'
    write(str,*) isize_wrk1d_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(xb(isize_wrk1d_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for xb.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if
    !
    write(str,*) inb_ibm; line = 'Allocating array IBM yb of size '//trim(adjustl(str))//'x'
    write(str,*) isize_wrk1d_ibm; line = trim(adjustl(line))//trim(adjustl(str))
    call IO_WRITE_ASCII(lfile,line)
    allocate(yb(isize_wrk1d_ibm), stat=ierr)
    if ( ierr /= 0 ) then
    call IO_WRITE_ASCII(efile,'DNS. Not enough memory for yb.')
    call DNS_STOP(DNS_ERROR_ALLOC)
    end if

    ! ------------------------------------------------------------------ !

    ! set alloc flag: done
    allocated = .true.

  end if

  return
end subroutine IBM_ALLOCATE

!########################################################################
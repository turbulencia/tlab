#include "types.h"
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
  use TLAB_VARS,      only : g, isize_field, istagger
  use TLAB_VARS,      only : imax, jmax, kmax
  use TLAB_PROCS  
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS,  only : ims_size_i, ims_size_k 
  use TLAB_MPI_VARS,  only : ims_npro_i, ims_npro_k 
#endif    

  implicit none

#include "integers.h"

  character(len=128), intent(in   ) :: C_FILE_LOC
  logical,            intent(inout) :: allocated

#ifdef USE_MPI 
  TINTEGER, parameter               :: idi = TLAB_MPI_I_PARTIAL 
  TINTEGER, parameter               :: idk = TLAB_MPI_K_PARTIAL 
#endif

  TINTEGER                          :: inb_ibm

  ! ================================================================== !
  ! npages
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    nyz = ims_size_i(idi)
  else
#endif
  nyz = jmax * kmax 
#ifdef USE_MPI
  end if
#endif

  nxz = imax * kmax     

#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    nxy = ims_size_k(idk)
  else
#endif
  nxy = imax * jmax
#ifdef USE_MPI
  end if
#endif

  ! array sizes
  isize_nobi      = nyz    
  isize_nobj      = nxz  
  isize_nobk      = nxy
  !
  isize_nobi_be   = nyz * nob_max    
  isize_nobj_be   = nxz * nob_max  
  isize_nobk_be   = nxy * nob_max
  !
  nspl            = 2 * nflu + 2    ! data points (incl. 2 interface points)
  isize_wrk1d_ibm = max(g(1)%size, max(g(2)%size, g(3)%size)) ! gap size unknown (max size assumed)

  ! allocate here all ibm related arrays
  if ( .not. allocated ) then
    ! eps          (geometry fields)
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC,   eps,      isize_field, 'eps'     )
    if ( istagger == 1 ) then
      call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, epsp,     isize_field, 'epsp'    )
    end if

    ! fld_ibm      (copying modified field)
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC,   fld_ibm,  isize_field, 'fld_ibm' )

    ! nob(i/j/k)   (number of objects)
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobi,   isize_nobi,    'nobi'  )
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobj,   isize_nobj,    'nobj'  )
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobk,   isize_nobk,    'nobk'  )
    
    ! nob(i/j/k)_b (beginnging objects)
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobi_b, isize_nobi_be, 'nobi_b')
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobj_b, isize_nobj_be, 'nobj_b')
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobk_b, isize_nobk_be, 'nobk_b')

    ! nob(i/j/k)_e (end of objects)
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobi_e, isize_nobi_be, 'nobi_e')
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobj_e, isize_nobj_be, 'nobj_e')
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, nobk_e, isize_nobk_be, 'nobk_e')
    
    ! xa, ya (spline arrays input)
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, xa, nspl,            'xa')
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, ya, nspl,            'ya')

    ! xb, yb (spline arrays output)
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, xb, isize_wrk1d_ibm, 'xb')
    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, yb, isize_wrk1d_ibm, 'yb')

    ! set alloc flag: done
    allocated = .true.
  end if

  return
end subroutine IBM_ALLOCATE

!########################################################################
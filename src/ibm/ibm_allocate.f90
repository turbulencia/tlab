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
!#   called once in dns_main.f90 to allocate needed memory for ibm module
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

subroutine IBM_ALLOCATE(C_FILE_LOC)

  use IBM_VARS
  use TLAB_VARS,      only : g, isize_field, inb_scal, stagger_on
  use TLAB_VARS,      only : imax, jmax, kmax
  use TLAB_CONSTANTS, only : wi
  use TLab_WorkFlow  
  use TLab_Memory
#ifdef USE_MPI
  use MPI
  use TLabMPI_VARS,  only : ims_size_i, ims_size_k 
  use TLabMPI_VARS,  only : ims_npro_i, ims_npro_k 
#endif    

  implicit none

  character(len=128), intent(in) :: C_FILE_LOC

#ifdef USE_MPI 
  integer(wi), parameter         :: idi = TLabMPI_I_PARTIAL 
  integer(wi), parameter         :: idk = TLabMPI_K_PARTIAL 
#endif
  integer(wi)                    :: nyz, nxz, nxy 

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

  ! eps          (geometry fields)
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   eps,     [isize_field], 'eps'    )

  if ( stagger_on ) then
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, epsp,    [isize_field], 'epsp'   )
  end if

  ! fld_ibm      (copy modified field)
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC,   fld_ibm, [isize_field], 'fld_ibm')

  ! nob(i/j/k)   (number of objects)
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobi,   [isize_nobi],    'nobi'  )
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobj,   [isize_nobj],    'nobj'  )
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobk,   [isize_nobk],    'nobk'  )
  
  ! nob(i/j/k)_b (beginning objects)
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobi_b, [isize_nobi_be], 'nobi_b')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobj_b, [isize_nobj_be], 'nobj_b')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobk_b, [isize_nobk_be], 'nobk_b')

  ! nob(i/j/k)_e (end of objects)
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobi_e, [isize_nobi_be], 'nobi_e')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobj_e, [isize_nobj_be], 'nobj_e')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, nobk_e, [isize_nobk_be], 'nobk_e')
  
  ! xa, ya (spline arrays input)
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, xa, [nspl],            'xa')
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, ya, [nspl],            'ya')

  ! xb, yb (spline arrays output)
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, xb, [isize_wrk1d_ibm], 'xb')
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, yb, [isize_wrk1d_ibm], 'yb')

  ! gammas for conditional averages
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, dy,       [jmax-1], 'dy'   )
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, facu,     [jmax-2], 'facu' )
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, facl,     [jmax-2], 'facl' )
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, gamma_0,  [jmax],   'eps_0')
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, gamma_1,  [jmax],   'eps_1')
  call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, scal_bcs, [jmax, inb_scal], 'scal_bcs')

  ! IBM case
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, ibm_case_x, [isize_nobi_be], 'ibm_case_x')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, ibm_case_y, [isize_nobj_be], 'ibm_case_y')
  call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, ibm_case_z, [isize_nobk_be], 'ibm_case_z')

  return
end subroutine IBM_ALLOCATE

!########################################################################
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
!#   transposes geometry field eps (epsi, epsj, epsk)
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

subroutine IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, tmp)

    use IBM_VARS
    use FDM, only: g
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLab_Constants, only: wi, wp
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_Transpose
#endif

    implicit none

    real(wp), dimension(isize_field), intent(out) :: epsi, epsj, epsk
    real(wp), dimension(isize_field), intent(inout) :: tmp

#ifdef USE_MPI
    integer(wi), parameter :: idi = TLAB_MPI_TRP_I_PARTIAL
    integer(wi), parameter :: idk = TLAB_MPI_TRP_K_PARTIAL
#endif
    integer(wi) :: nyz, nxy

    ! ================================================================== !
    ! MPI  and local transposition in x
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        call TLabMPI_TransposeI_Forward(eps, tmp, idi)
        ! nyz = ims_size_i(idi)
        nyz = ims_trp_plan_i(idi)%nlines
    else
#endif
        tmp = eps
        nyz = jmax*kmax
#ifdef USE_MPI
    end if
#endif

#ifdef USE_ESSL
    call DGETMO(tmp, g(1)%size, g(1)%size, nyz, epsi, nyz)
#else
    call TLab_Transpose(tmp, g(1)%size, nyz, g(1)%size, epsi, nyz)
#endif
    ! -------------------------------------------------------------------
    ! local transposition in y
    nxy = imax*jmax
#ifdef USE_ESSL
    call DGETMO(eps, nxy, nxy, kmax, epsj, kmax)
#else
    call TLab_Transpose(eps, nxy, kmax, nxy, epsj, kmax)
#endif
    ! -------------------------------------------------------------------
    ! MPI transposition in z
#ifdef USE_MPI
    if (ims_npro_k > 1) then
        call TLabMPI_TransposeK_Forward(eps, epsk, idk)
    else
#endif
        epsk = eps
#ifdef USE_MPI
    end if
#endif

    return
end subroutine IBM_GEOMETRY_TRANSPOSE

!########################################################################

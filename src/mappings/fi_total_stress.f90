#include "dns_const.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2023/03/28 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   Calculate the total stress tensor
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

module FI_TOTAL_STRESS
    
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS,      only: g, visc, imode_ibm
    use IBM_VARS,       only: ibm_partial
    use OPR_PARTIAL

    implicit none
    
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_TOTAL_STRESS_TENSOR

contains

!########################################################################

    subroutine FI_TOTAL_STRESS_TENSOR(nx, ny, nz, u, v, w, p, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
        
        integer(wi),                   intent(in   ) :: nx, ny, nz
        real(wp), dimension(nx,ny,nz), intent(in   ) :: u, v, w, p
        real(wp), dimension(nx,ny,nz), intent(  out) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 ! Sxx,Syy,Szz,Sxy,Sxz,Syz
        
        ! ================================================================== !
    
        ! OPR_PARTIAL_X/Y/Z uses modified fields for derivatives
        if ( imode_ibm == 1 ) ibm_partial = .true.

        ! Off diagonal terms

        ! Uy, Vx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp4)
        tmp4 = visc*(tmp4 + tmp1)

        ! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp5)
        tmp5 = visc*(tmp5 + tmp1)

        ! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp6)
        tmp6 = visc*(tmp6 + tmp1)

        ! diagonal terms
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        tmp1 =  -p + 2.0_wp*visc*tmp1
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        tmp2 =  -p + 2.0_wp*visc*tmp2
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3)
        tmp3 =  -p + 2.0_wp*visc*tmp3

        return

    end subroutine FI_TOTAL_STRESS_TENSOR
    
end module FI_TOTAL_STRESS

!     subroutine FI_WALL_SHEAR_IBM(nx, ny, nz, u, v, w, p, result, tmp1, tmp2, tmp3, tmp4, tmp5)
       
!         real(wp), dimension(nx,ny,nz), intent(IN   ) :: u, v, w, p
!         integer(wi),                   intent(IN   ) :: nx, ny, nz
!         real(wp), dimension(nx,ny,nz), intent(  OUT) :: result
!         real(wp), dimension(nx,ny,nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5 
        
!         target :: tmp5

!         real(wp), dimension(:), pointer :: tau
!         real(wp), dimension(:), pointer :: g_s, g_f, wrk1d

!         integer(wi) :: i, j, k, ip


!         ! ================================================================== !
!         ! assign pointers
!         tmp5(:, :, :) = 0.0_wp; i = 1_wi; k = 1_wi
!         tau   => tmp5(i, :, k); i = i + 1_wi
!         wrk1d => tmp5(i, :, k); i = i + 1_wi

!         ! OPR_PARTIAL_X/Y/Z uses modified fields for derivatives
!         ibm_partial = .true.
        
!         ! 1. Wall shear stress (tau_w[:,1 & ny,:] = nu*(dudy**2+dwdy**2)**0.5) on lower ground
!         call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
!         call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2)

!         result = visc*sqrt(tmp1**2 + tmp2**2)

!         call IBM_BCS_FIELD(result)

!         call WALL_SHEAR_UPSIDE(nx, ny, nz, result, nobj, nobj_b, nobj_e)
!         ! write top of elements in upper/lower

!         ! call IBM_BCS_FIELD(result)





!         ! call AVG_IK_V(nx, ny, nz, ny, tmp3, g(1)%jac, g(3)%jac, tau, wrk1d, area)


!         ! write(*,*)'tau(1)  = ', tau(1)
!         ! write(*,*)'tau(ny) = ', tau(ny)

        
!         ibm_partial = .false.

!         nullify (wrk1d,tau)
        
!         return

!     end subroutine FI_WALL_SHEAR_IBM

! !########################################################################

!     subroutine WALL_SHEAR_UPSIDE(nx, ny, nz, result, nobj, nobj_b, nobj_e)
       
!         integer(wi),                           intent(IN   ) :: nx, ny, nz
!         real(wp),    dimension(nx,ny,nz),      intent(INOUT) :: result
!         integer(wi), dimension(nz,nx),         intent(IN   ) :: nobj
!         integer(wi), dimension(nz,nx,nob_max), intent(IN   ) :: nobj_b, nobj_e
        
!         integer(wi) :: i,j,k, iob
!         ! ================================================================== !

!         do i = 1, nx
!             do k = 1, nz
!                 if ( nobj(k,i) /= 0 ) then   ! line contains immersed object(s)
!                     if ( nobj(k,i) > 2) then
!                         call TLAB_WRITE_ASCII(efile, 'IBM_WALL_SHEAR. Not implemented.')
!                         call TLAB_STOP(DNS_ERROR_IBM_SHEAR)
!                     end if 
!                     do iob = 1, nobj(k,i)    ! loop over immersed object(s)
!                         if ( nobj_b(k,i,iob) == 1 ) then
!                             result(i,1,k) = result(i,nobj_e(k,i,iob)+1,k)
!                         else if ( nobj_b(k,i,iob) /= 1 ) then
!                             result(i,2,k) = result(i,nobj_b(k,i,iob)-1,k)
!                         end if
!                     end do
!                 end if
!             end do
!         end do

!         return

!     end subroutine WALL_SHEAR_UPSIDE

! !########################################################################

!     subroutine FI_WALL_SHEAR_NOIBM(nx, ny, nz, u, w, result, tmp1, tmp2)
       
!         integer(wi),                   intent(IN   ) :: nx, ny, nz
!         real(wp), dimension(nx,ny,nz), intent(IN   ) :: u, w
!         real(wp), dimension(nx,ny,nz), intent(  OUT) :: result
!         real(wp), dimension(nx,ny,nz), intent(INOUT) :: tmp1, tmp2
        
!         ! ================================================================== !

!         ! 1. Wall shear stress (tau_w[:,1 & ny,:] = nu*(dudy**2+dwdy**2)**0.5) 
!         call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
!         call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2)

!         result = visc*sqrt(tmp1**2 + tmp2**2)
    
!         ! 2. only values on upper/lower boundary
!         ! result(:,2:ny-1,:) = 0.0_wp
!         result(:,2,:) = result(:,ny,:)
!         ! result(:,3:ny,:) = 0.0_wp

!         return

!     end subroutine FI_WALL_SHEAR_NOIBM
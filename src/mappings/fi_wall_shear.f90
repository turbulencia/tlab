#include "dns_const.h"
#include "dns_error.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2023/03/28 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   Compute Wall shear stress of fields
!#
!#
!#   This routine only works for if objects are fixed  on the lower domain boundary!
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

module FI_SHEAR
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS,      only: g, area, visc
    use IBM_VARS!,       only: eps, gamma_f, gamma_s, ibm_partial
    use AVGS,           only: AVG_IK_V
    use OPR_PARTIAL

    use TLAB_CONSTANTS, only : efile
    use TLAB_PROCS

    implicit none
    
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_WALL_SHEAR_IBM, FI_WALL_SHEAR_NOIBM

contains

!########################################################################
    
    subroutine FI_WALL_SHEAR_IBM(nx, ny, nz, u, v, w, p, result, tmp1, tmp2, tmp3, tmp4, tmp5)
       
        real(wp), dimension(nx,ny,nz), intent(IN   ) :: u, v, w, p
        integer(wi),                   intent(IN   ) :: nx, ny, nz
        real(wp), dimension(nx,ny,nz), intent(  OUT) :: result
        real(wp), dimension(nx,ny,nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5 
        
        target :: tmp5

        real(wp), dimension(:), pointer :: tau
        real(wp), dimension(:), pointer :: g_s, g_f, wrk1d

        integer(wi) :: i, j, k, ip


        ! ================================================================== !
        ! assign pointers
        tmp5(:, :, :) = 0.0_wp; i = 1_wi; k = 1_wi
        tau   => tmp5(i, :, k); i = i + 1_wi
        wrk1d => tmp5(i, :, k); i = i + 1_wi

        ! OPR_PARTIAL_X/Y/Z uses modified fields for derivatives
        ibm_partial = .true.
        
        ! 1. Wall shear stress (tau_w[:,1 & ny,:] = nu*(dudy**2+dwdy**2)**0.5) on lower ground
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2)

        result = visc*sqrt(tmp1**2 + tmp2**2)

        call IBM_BCS_FIELD(result)

        call WALL_SHEAR_UPSIDE(nx, ny, nz, result, nobj, nobj_b, nobj_e)
        ! write top of elements in upper/lower

        ! call IBM_BCS_FIELD(result)





        ! call AVG_IK_V(nx, ny, nz, ny, tmp3, g(1)%jac, g(3)%jac, tau, wrk1d, area)


        ! write(*,*)'tau(1)  = ', tau(1)
        ! write(*,*)'tau(ny) = ', tau(ny)

        
        ibm_partial = .false.

        nullify (wrk1d,tau)
        
        return

    end subroutine FI_WALL_SHEAR_IBM

!########################################################################

    subroutine WALL_SHEAR_UPSIDE(nx, ny, nz, result, nobj, nobj_b, nobj_e)
       
        integer(wi),                            intent(IN   ) :: nx, ny, nz
        real(wp),    dimension(nx,ny,nz),       intent(INOUT) :: result
        integer(wi), dimension(nz,nx),          intent(IN   ) :: nobj
        integer(wi), dimension(nz,nx, nob_max), intent(IN   ) :: nobj_b, nobj_e
        
        integer(wi) :: i,j,k, iob
        ! ================================================================== !



        do i = 1, nx
            do k = 1, nz
                if ( nobj(k,i) /= 0 ) then   ! line contains immersed object(s)
                    do iob = 1, nobj(k,i)    ! loop over immersed object(s)
                        ! write(*,*) 'nobj:1  ',nobj(:,1) 
                        ! write(*,*) 'nobj_b',nobj_b(:,1,2) 
                        ! stop 
                        if ( nobj_b(k,i,iob) == 1 ) then
                            ! result(i,1,k) = result(i,nobj_e(k,i,iob),k)
                            result(i,1,k) = result(i,nobj_e(k,i,iob)+1,k)
                        
                        else if ( nobj_b(k,i,iob) /= 1 ) then
                            ! result(i,ny,k) = result(i,nobj_b(k,i,iob),k)
                            result(i,ny,k) = result(i,nobj_b(k,i,iob)-1,k)

                        end if
                        
                        if ( iob > 2) then
                            call TLAB_WRITE_ASCII(efile, 'IBM_WALL_SHEAR. .....')
                            call TLAB_STOP(DNS_ERROR_IBM_SHEAR)
                        end if 

                    end do
                end if
            end do
        end do

        return

    end subroutine WALL_SHEAR_UPSIDE

!########################################################################

    subroutine FI_WALL_SHEAR_NOIBM(nx, ny, nz, u, w, result, tmp1, tmp2)
       
        integer(wi),                   intent(IN   ) :: nx, ny, nz
        real(wp), dimension(nx,ny,nz), intent(IN   ) :: u, w
        real(wp), dimension(nx,ny,nz), intent(  OUT) :: result
        real(wp), dimension(nx,ny,nz), intent(INOUT) :: tmp1, tmp2
        
        ! ================================================================== !

        ! 1. Wall shear stress (tau_w[:,1 & ny,:] = nu*(dudy**2+dwdy**2)**0.5) 
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2)

        result = visc*sqrt(tmp1**2 + tmp2**2)
    
        ! 2. only values on upper/lower boundary
        result(:,2:ny-1,:) = 0.0_wp

        return

    end subroutine FI_WALL_SHEAR_NOIBM

end module FI_SHEAR

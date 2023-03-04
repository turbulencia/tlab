#include "dns_const.h"

!########################################################################
!# Calculate the strain tensor and its magnitude as given by s_ij s_ij:
!# u,x^2 + v,y^2 + w,z^2 +1/2( (u,y+v,x)^2 + (u,z+w,x)^2 + (v,z+w,y)^2 )
!# and terms in its evolution equation.
!########################################################################
module FI_STRAIN_EQN
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: g
    use OPR_PARTIAL

    implicit none
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_STRAIN_TENSOR
    public :: FI_STRAIN
    public :: FI_STRAIN_PRODUCTION
    public :: FI_STRAIN_DIFFUSION
    public :: FI_STRAIN_PRESSURE
    public :: FI_STRAIN_A

contains
!########################################################################
! Calculate tensor components
!########################################################################
    subroutine FI_STRAIN_TENSOR(nx, ny, nz, u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 ! components Sxx, Syy, Szz, Sxy, Sxz, Syz

! ###################################################################
! Off diagonal terms; result1 used as auxiliar array
! Uy, Vx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp4)
        tmp4 = 0.5_wp*(tmp4 + tmp1)

! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp5)
        tmp5 = 0.5_wp*(tmp5 + tmp1)

! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp6)
        tmp6 = 0.5_wp*(tmp6 + tmp1)

! diagonal terms
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3)

        return
    end subroutine FI_STRAIN_TENSOR

!########################################################################
! Calculate tensor magnitude
!########################################################################
    subroutine FI_STRAIN(nx, ny, nz, u, v, w, result, tmp1, tmp2)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

! ###################################################################
! Ux, Vy, Wz
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        result = tmp1*tmp1 + tmp2*tmp2

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp2)
        result = result + tmp2*tmp2

! Uy, Vx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        result = result + 0.5_wp*(tmp1 + tmp2)*(tmp1 + tmp2)

! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp2)
        result = result + 0.5_wp*(tmp1 + tmp2)*(tmp1 + tmp2)

! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        result = result + 0.5_wp*(tmp1 + tmp2)*(tmp1 + tmp2)

        return
    end subroutine FI_STRAIN

!########################################################################
!# Calculate the strain production term as given by -s_ij s_jk s_ki
!# minus 1/4 of the vorticity production. The first term is:
!# s_11^3 + s_22^3 +  s_33^3 +
!# 3 s_11 (s_12^2+s_13^2) + 3 s_22 (s_21^2+s_23^2) + 3 s_33 (s_31^2+s_32^2) +
!# 2 s_12 s_23 s_31
!# or
!# 2 s_12 s_23 s_31 +
!# s_11 ( s_11^2 + 3 ( s_12^2 + s_13^2 ) ) +
!# s_22 ( s_22^2 + 3 ( s_12^2 + s_23^2 ) ) +
!# s_33 ( s_33^2 + 3 ( s_13^2 + s_23^2 ) )
!########################################################################
    subroutine FI_STRAIN_PRODUCTION(nx, ny, nz, u, v, w, result, s_12, s_13, s_23, tmp1, tmp2)
        use FI_VORTICITY_EQN, only: FI_VORTICITY_PRODUCTION
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: s_12, s_13, s_23
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

! ###################################################################
! Vorticity production part
! ###################################################################
        call FI_VORTICITY_PRODUCTION(nx, ny, nz, u, v, w, result, s_12, s_13, s_23, tmp1, tmp2)

        result = 0.25_wp*result

! ###################################################################
! Pure strain term
! ###################################################################
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        s_12 = 0.5_wp*(tmp1 + tmp2)

! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp2)
        s_13 = 0.5_wp*(tmp1 + tmp2)

! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        s_23 = 0.5_wp*(tmp1 + tmp2)

! -------------------------------------------------------------------
! Production term
! -------------------------------------------------------------------
        result = result + 2.0_wp*s_12*s_13*s_23

        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        result = result + tmp1*(tmp1*tmp1 + 3.0_wp*(s_12*s_12 + s_13*s_13))

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp1)
        result = result + tmp1*(tmp1*tmp1 + 3.0_wp*(s_12*s_12 + s_23*s_23))

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp1)
        result = result + tmp1*(tmp1*tmp1 + 3.0_wp*(s_13*s_13 + s_23*s_23))

! right sign
        result = -result

        return
    end subroutine FI_STRAIN_PRODUCTION

!########################################################################
!# Calculate the strain diffusion term as given by s_ij lap s_ij
!# The kinematic viscosity \nu is not multiplied here.
!########################################################################
    subroutine FI_STRAIN_DIFFUSION(nx, ny, nz, u, v, w, result, strain, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: strain
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! ###################################################################
! -------------------------------------------------------------------
! S_11 = U,x
! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, strain)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = strain*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! S_22 = V,y
! -------------------------------------------------------------------
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, strain)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = result + strain*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! S_33 = W,z
! -------------------------------------------------------------------
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, strain)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = result + strain*(tmp1 + tmp2 + tmp3)
! -------------------------------------------------------------------
! S_12 = (V,x + U,y)/2
! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        strain = tmp1 + tmp2

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = result + 0.5_wp*strain*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! S_13 = (U,z + W,x)/2
! -------------------------------------------------------------------
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)
        strain = tmp1 + tmp2

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = result + 0.5_wp*strain*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! S_23 = (W,y + V,z)/2
! -------------------------------------------------------------------
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        strain = tmp1 + tmp2

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), strain, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), strain, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), strain, tmp1, tmp4)

        result = result + 0.5_wp*strain*(tmp1 + tmp2 + tmp3)

        return
    end subroutine FI_STRAIN_DIFFUSION

!########################################################################
!# Calculate the strain-pressure term as given by -s_ij p,ij:
!########################################################################
    subroutine FI_STRAIN_PRESSURE(nx, ny, nz, u, v, w, p, result, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w, p
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), p, tmp2, tmp4)
        result = tmp1*tmp2

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), p, tmp2, tmp4)
        result = result + tmp1*tmp2

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), p, tmp2, tmp4)
        result = result + tmp1*tmp2

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! 2 s_12 p,xy
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), p, tmp2)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp2, tmp1)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp2)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp3)
        result = result + tmp1*(tmp2 + tmp3)

! 2 s_13 p,xz
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), p, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), tmp2, tmp1)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp3)
        result = result + tmp1*(tmp2 + tmp3)

! 2 s_23 p,yz
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), p, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), tmp2, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp3)
        result = result + tmp1*(tmp2 + tmp3)

! right sign
        result = -result

        return
    end subroutine FI_STRAIN_PRESSURE

!########################################################################
!# Compute the scalar quantity da/dx_i*du_j/dx_i*da/dx_j (strain1) and
!# that normalized by grad(a)*grad(a) (strain2)
!########################################################################
    subroutine FI_STRAIN_A(nx, ny, nz, a, u, v, w, strain1, strain2, normal1, normal2, normal3, tmp1)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: a, u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: strain1, strain2
        real(wp), dimension(nx*ny*nz), intent(OUT) :: normal1, normal2, normal3
        real(wp), dimension(nx*ny*nz), intent(OUT) :: tmp1 ! Returns grad(a)*grad(a)

        integer(wi) i
        
        ! ###################################################################
        ! Compute normals
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), a, normal1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), a, normal2)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), a, normal3)

        ! Compute gradient terms with u
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        strain1 = normal1*tmp1*normal1
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
        strain1 = strain1 + normal2*tmp1*normal1
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
        strain1 = strain1 + normal3*tmp1*normal1

        ! Compute gradient terms with v
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        strain1 = strain1 + normal1*tmp1*normal2
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp1)
        strain1 = strain1 + normal2*tmp1*normal2
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp1)
        strain1 = strain1 + normal3*tmp1*normal2

        ! Compute gradient terms with 3
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        strain1 = strain1 + normal1*tmp1*normal3
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        strain1 = strain1 + normal2*tmp1*normal3
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp1)
        strain1 = strain1 + normal3*tmp1*normal3

        ! norm of the gradient of a
        tmp1 = normal1**2 + normal2**2 + normal3**2

        ! normalize
        do i = 1, nx*ny*nz
            if (tmp1(i) > 0.0_wp) then
                strain2(i) = strain1(i)/tmp1(i)
            else
                strain2(i) = strain1(i)
            end if
        end do

        return
    end subroutine FI_STRAIN_A
end module FI_STRAIN_EQN

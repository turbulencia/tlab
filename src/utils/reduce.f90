subroutine REDUCE(imax, jmax, kmax, a, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, np                           ! np is the number of sampled planes
    integer(wi), dimension(np), intent(IN) :: p ! array with the i location of the planes
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: a ! input array (big)
    real(wp), dimension(np, jmax, kmax), intent(OUT) :: b ! output array (small)

    integer(wi) i, j, k, n

    do k = 1, kmax
        do j = 1, jmax
            do n = 1, np
                i = p(n)
                b(n, j, k) = a(i, j, k)
            end do
        end do
    end do

    return
end subroutine REDUCE

! #######################################################################
! #######################################################################
subroutine REDUCE_SUM(imax, jmax, kmax, a1, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, np
    integer(wi), dimension(np), intent(IN) :: p
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: a1, a2
    real(wp), dimension(np, jmax, kmax), intent(OUT) :: b

    integer(wi) i, j, k, n

    do k = 1, kmax
        do j = 1, jmax
            do n = 1, np
                i = p(n)
                b(n, j, k) = a1(i, j, k) + a2(i, j, k)
            end do
        end do
    end do

    return
end subroutine REDUCE_SUM

! #######################################################################
! #######################################################################
subroutine REDUCE_SUB(imax, jmax, kmax, a1, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, np
    integer(wi), dimension(np), intent(IN) :: p
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: a1, a2
    real(wp), dimension(np, jmax, kmax), intent(OUT) :: b

    integer(wi) i, j, k, n

    do k = 1, kmax
        do j = 1, jmax
            do n = 1, np
                i = p(n)
                b(n, j, k) = a1(i, j, k) - a2(i, j, k)
            end do
        end do
    end do

    return
end subroutine REDUCE_SUB

! #######################################################################
! #######################################################################
subroutine REDUCE_MUL(imax, jmax, kmax, a1, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, np
    integer(wi), dimension(np), intent(IN) :: p
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: a1, a2
    real(wp), dimension(np, jmax, kmax), intent(OUT) :: b

    integer(wi) i, j, k, n

    do k = 1, kmax
        do j = 1, jmax
            do n = 1, np
                i = p(n)
                b(n, j, k) = a1(i, j, k)*a2(i, j, k)
            end do
        end do
    end do

    return
end subroutine REDUCE_MUL

! #######################################################################
! #######################################################################
subroutine REDUCE_DIV(imax, jmax, kmax, a1, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, np
    integer(wi), dimension(np), intent(IN) :: p
    real(wp), dimension(imax, jmax, kmax), intent(IN) :: a1, a2
    real(wp), dimension(np, jmax, kmax), intent(OUT) :: b

    integer(wi) i, j, k, n

    do k = 1, kmax
        do j = 1, jmax
            do n = 1, np
                i = p(n)
                b(n, j, k) = a1(i, j, k)/a2(i, j, k)
            end do
        end do
    end do

    return
end subroutine REDUCE_DIV

! #######################################################################
! #######################################################################
subroutine REDUCE_BLOCK_INPLACE(nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst, a, wrk1d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst
    real(wp), dimension(nx*ny*nz), intent(INOUT) :: a
    real(wp), dimension(nx_dst), intent(INOUT) :: wrk1d

! -------------------------------------------------------------------
    integer(wi) j, k, nxy, nxy_dst, ip, ip_dst

! -------------------------------------------------------------------
    nxy = nx*ny
    nxy_dst = nx_dst*ny_dst

    do k = 1, nz_dst
        ip = (k - 1)*nxy + (nz1 - 1)*nxy + (ny1 - 1)*nx + (nx1 - 1) + 1
        ip_dst = (k - 1)*nxy_dst + 1
        do j = 1, ny_dst
            wrk1d(1:nx_dst) = a(ip:ip + nx_dst - 1)
            a(ip_dst:ip_dst + nx_dst - 1) = wrk1d(1:nx_dst)

            ip = ip + nx
            ip_dst = ip_dst + nx_dst
        end do
    end do

    return
end subroutine REDUCE_BLOCK_INPLACE

! #######################################################################
! #######################################################################
subroutine REDUCE_BLOCK_INPLACE_INT1(nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst, a, wrk1d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst
    integer(1), dimension(nx*ny*nz), intent(INOUT) :: a
    integer(1), dimension(nx_dst), intent(INOUT) :: wrk1d

! -------------------------------------------------------------------
    integer(wi) j, k, nxy, nxy_dst, ip, ip_dst

! -------------------------------------------------------------------
    nxy = nx*ny
    nxy_dst = nx_dst*ny_dst

    do k = 1, nz_dst
        ip = (k - 1)*nxy + (nz1 - 1)*nxy + (ny1 - 1)*nx + (nx1 - 1) + 1
        ip_dst = (k - 1)*nxy_dst + 1
        do j = 1, ny_dst
            wrk1d(1:nx_dst) = a(ip:ip + nx_dst - 1)
            a(ip_dst:ip_dst + nx_dst - 1) = wrk1d(1:nx_dst)

            ip = ip + nx
            ip_dst = ip_dst + nx_dst
        end do
    end do

    return
end subroutine REDUCE_BLOCK_INPLACE_INT1

! #######################################################################
! #######################################################################
subroutine REDUCE_BLOCK(nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst, a, wrk3d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst
    real(wp), dimension(nx*ny*nz), intent(INOUT) :: a
    real(wp), dimension(nx_dst*ny_dst*nz_dst), intent(INOUT) :: wrk3d

! -------------------------------------------------------------------
    integer(wi) j, k, nxy, nxy_dst, ip, ip_dst

! -------------------------------------------------------------------
    nxy = nx*ny
    nxy_dst = nx_dst*ny_dst

    do k = 1, nz_dst
        ip = (k - 1)*nxy + (nz1 - 1)*nxy + (ny1 - 1)*nx + (nx1 - 1) + 1
        ip_dst = (k - 1)*nxy_dst + 1
        do j = 1, ny_dst
            wrk3d(ip_dst:ip_dst + nx_dst - 1) = a(ip:ip + nx_dst - 1)

            ip = ip + nx
            ip_dst = ip_dst + nx_dst
        end do
    end do

    return
end subroutine REDUCE_BLOCK

! #######################################################################
! #######################################################################
subroutine REDUCE_Y_ALL(nx, ny, nz, nvar1, a1, nvar2, a2, aux, np, np_aux, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) nx, ny, nz, np, np_aux, nvar1, nvar2                           ! np is the number of sampled planes
    ! np_aux is one or zero, additional data
    integer(wi), dimension(np), intent(IN) :: p      ! array with the j location of the planes
    real(wp), dimension(nx, ny, nz, *), intent(IN) :: a1, a2 ! input array (big)
    real(wp), dimension(nx, nz, nvar1 + nvar2), intent(IN) :: aux    ! additional data with different structure
    real(wp), dimension(nx, np, nvar1 + nvar2, nz), intent(OUT) :: b      ! output array (small)

    integer(wi) j, j_loc, ivar, k

    do k = 1, nz
        do ivar = 1, nvar1
            do j = 1, np - np_aux
                j_loc = p(j)
                b(1:nx, j, ivar, k) = a1(1:nx, j_loc, k, ivar)
            end do
            if (np_aux > 0) b(1:nx, j, ivar, k) = aux(1:nx, k, ivar) ! Additional data, if needed
        end do

        do ivar = 1, nvar2 ! if nvar is 0, then array a2 is not used
            do j = 1, np - np_aux
                j_loc = p(j)
                b(1:nx, j, ivar + nvar1, k) = a2(1:nx, j_loc, k, ivar)
            end do
            if (np_aux > 0) b(1:nx, j, ivar + nvar1, k) = aux(1:nx, k, ivar + nvar1) ! Additional data, if needed
        end do

    end do

    return
end subroutine REDUCE_Y_ALL

! #######################################################################
! #######################################################################
subroutine REDUCE_Z_ALL(nx, ny, nz, nvar1, a1, nvar2, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) nx, ny, nz, np, nvar1, nvar2                               ! np is the number of sampled planes
    integer(wi), dimension(np), intent(IN) :: p      ! array with the j location of the planes
    real(wp), dimension(nx*ny, nz, *), intent(IN) :: a1, a2 ! input array (big)
    real(wp), dimension(nx*ny, nvar1 + nvar2, np), intent(OUT) :: b      ! output array (small)

    integer(wi) k, k_loc, ivar

    do k = 1, np
        k_loc = p(k)
        do ivar = 1, nvar1
            b(1:nx*ny, ivar, k) = a1(1:nx*ny, k_loc, ivar)
        end do

        do ivar = 1, nvar2 ! if nvar is 0, then array a2 is not used
            b(1:nx*ny, ivar + nvar1, k) = a2(1:nx*ny, k_loc, ivar)
        end do

    end do

    return
end subroutine REDUCE_Z_ALL

! #######################################################################
! #######################################################################
subroutine REDUCE_X_ALL(nx, ny, nz, nvar1, a1, nvar2, a2, np, p, b)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) nx, ny, nz, np, nvar1, nvar2                               ! np is the number of sampled planes
    integer(wi), dimension(np), intent(IN) :: p      ! array with the j location of the planes
    real(wp), dimension(nx, ny, nz, *), intent(IN) :: a1, a2 ! input array (big)
    real(wp), dimension(np, ny, nvar1 + nvar2, nz), intent(OUT) :: b      ! output array (small)

    integer(wi) i, i_loc, j, ivar, k

    do k = 1, nz
        do ivar = 1, nvar1
            do j = 1, ny
                do i = 1, np
                    i_loc = p(i)
                    b(i, j, ivar, k) = a1(i_loc, j, k, ivar)
                end do
            end do
        end do

        do ivar = 1, nvar2 ! if nvar is 0, then array a2 is not used
            do j = 1, ny
                do i = 1, np
                    i_loc = p(i)
                    b(i, j, ivar + nvar1, k) = a2(i_loc, j, k, ivar)
                end do
            end do
        end do

    end do

    return
end subroutine REDUCE_X_ALL

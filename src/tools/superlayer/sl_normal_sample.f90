!########################################################################
!#
!# The array b, with a number ofields set by nfields, is sampled along
!# direction (nx,ny,nz) and the obtained profiles are stored in array c.
!#
!# Assuming uniform spacing and periodic conditions. If parallel, the
!# points outside the local PE field are set to zero
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
subroutine SL_NORMAL_SAMPLE(imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield, &
                            scalex, scalez, factor, x, y, z, sl, b, c, nx, ny, nz)
    use TLab_Constants, only: wp, wi

#ifdef USE_MPI
    use TLabMPI_VARS
#endif

    implicit none

    integer(wi) imax, jmax, kmax, nmax, istep, kstep, nfield_loc, nfield
    real(wp) scalex, scalez, factor
    real(wp) x(*), y(*), z(*)
    real(wp) sl(imax, kmax)
    real(wp) b(imax, jmax, kmax, nfield_loc)
    real(wp) c(nfield, nmax, *)
    real(wp) nx(imax, jmax, kmax), ny(imax, jmax, kmax), nz(imax, jmax, kmax)

! -------------------------------------------------------------------
    integer(wi) i, k, n, ifield
    integer(wi) im, jm, km, ip, jp, kp
    integer(wi) iprofile
    real(wp) dx_u, dy_u, dz_u, dn_u
    real(wp) x_loc, y_loc, z_loc, nx_loc, ny_loc, nz_loc, norm, dy_loc, dn_loc
    real(wp) xr, yr, zr, xrc, yrc, zrc

! ###################################################################
! mean value of grid spacing
    dx_u = x(2) - x(1)
    dy_u = y(2) - y(1)
    dz_u = z(2) - z(1)
    dn_u = (dx_u + dy_u + dz_u)/C_3_R*factor

! ###################################################################
! Loop on the points
! ###################################################################
    iprofile = 0
    do k = kstep, kmax, kstep
        do i = istep, imax, istep
            iprofile = iprofile + 1

! -------------------------------------------------------------------
! Get normal direction at the point (pointing to the irrotational side)
! -------------------------------------------------------------------
            jm = int((sl(i, k) - y(1))/dy_u) + 1
            dy_loc = sl(i, k) - y(jm)
            nx_loc = nx(i, jm, k) + (nx(i, jm + 1, k) - nx(i, jm, k))/(y(jm + 1) - y(jm))*dy_loc
            ny_loc = ny(i, jm, k) + (ny(i, jm + 1, k) - ny(i, jm, k))/(y(jm + 1) - y(jm))*dy_loc
            nz_loc = nz(i, jm, k) + (nz(i, jm + 1, k) - nz(i, jm, k))/(y(jm + 1) - y(jm))*dy_loc
            norm = sqrt(nx_loc*nx_loc + ny_loc*ny_loc + nz_loc*nz_loc)
            nx_loc = -nx_loc/norm
            ny_loc = -ny_loc/norm
            nz_loc = -nz_loc/norm

! -------------------------------------------------------------------
! Sample fields along the normal
! -------------------------------------------------------------------
            do n = 1, nmax
                dn_loc = M_REAL(n - 1 - nmax/2)*dn_u

! sampling point in physical space
                x_loc = x(i) + dn_loc*nx_loc
                y_loc = sl(i, k) + dn_loc*ny_loc
#ifdef USE_MPI
                z_loc = z(k + ims_offset_k) + dn_loc*nz_loc
#else
                z_loc = z(k) + dn_loc*nz_loc
#endif

! periodic conditions in X
                if (x_loc > x(1) + scalex) x_loc = x_loc - scalex
                if (x_loc < x(1)) x_loc = x_loc + scalex
! periodic conditions in Z, if serial. If not, assign a given value
#ifdef USE_MPI
                if (z_loc > z(kmax + ims_offset_k) .or. z_loc < z(1 + ims_offset_k)) then
                    do ifield = 1, nfield_loc
                        c(ifield, n, iprofile) = C_0_R
                    end do

                else
#else
                    if (z_loc > z(1) + scalez) z_loc = z_loc - scalez
                    if (z_loc < z(1)) z_loc = z_loc + scalez
#endif

! get left corner
                    im = int((x_loc - x(1))/dx_u) + 1
                    jm = int((y_loc - y(1))/dy_u) + 1
#ifdef USE_MPI
                    km = int((z_loc - z(1 + ims_offset_k))/dz_u) + 1
#else
                    km = int((z_loc - z(1))/dz_u) + 1
#endif

! define relative distances for linear interpolation, and complemetaries
                    xr = (x_loc - x(im))/dx_u
                    yr = (y_loc - y(jm))/dy_u
#ifdef USE_MPI
                    zr = (z_loc - z(km + ims_offset_k))/dz_u
#else
                    zr = (z_loc - z(km))/dz_u
#endif
                    xrc = C_1_R - xr
                    yrc = C_1_R - yr
                    zrc = C_1_R - zr

! interpolate linearly all the fields
                    if (im == imax) then
                        ip = 1
                    else
                        ip = im + 1
                    end if
                    jp = jm + 1
                    if (km == kmax) then
                        kp = 1
                    else
                        kp = km + 1
                    end if
                    do ifield = 1, nfield_loc
                        c(ifield, n, iprofile) = &
                            b(im, jm, km, ifield)*xrc*yrc*zrc + &
                            b(im, jm, kp, ifield)*xrc*yrc*zr + &
                            b(im, jp, km, ifield)*xrc*yr*zrc + &
                            b(im, jp, kp, ifield)*xrc*yr*zr + &
                            b(ip, jm, km, ifield)*xr*yrc*zrc + &
                            b(ip, jm, kp, ifield)*xr*yrc*zr + &
                            b(ip, jp, km, ifield)*xr*yr*zrc + &
                            b(ip, jp, kp, ifield)*xr*yr*zr
                    end do

#ifdef USE_MPI
                end if
#endif

            end do

        end do
    end do

    return
end subroutine SL_NORMAL_SAMPLE

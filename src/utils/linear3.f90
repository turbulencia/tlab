#ifdef SINGLE_PREC
#define AXPY_LOC SAXPY
#define SCAL_LOC SSCAL
#else
#define AXPY_LOC DAXPY
#define SCAL_LOC DSCAL
#endif

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2009/09/01 - J.P. Mellado
!#              Reviewed, OpenMP
!# 2011/11/01 - C. Ansorge
!#              OpenMP Optimization
!#
!########################################################################
!# DESCRIPTION
!#
!# Vectorized tridiagonal solver based ob Thomas algorithm.
!#
!########################################################################

! #######################################################################
! LU factorization stage; L with diagonal unity
! #######################################################################
subroutine TRIDFS(nmax, a, b, c)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: nmax
    real(wp), dimension(nmax), intent(INOUT) :: a, b, c

! -----------------------------------------------------------------------
    integer(wi) n

! #######################################################################
    do n = 2, nmax
        a(n) = a(n)/b(n - 1)
        b(n) = b(n) - a(n)*c(n - 1)
    end do

! Final operations
    a = -a
    b = 1.0_wp/b
    c = -c

    return
end subroutine TRIDFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine TRIDSS(nmax, len, a, b, c, f)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

    integer(wi), intent(IN) :: nmax  ! dimension of tridiagonal systems
    integer(wi), intent(IN) :: len   ! number of systems to be solved
    real(wp), dimension(nmax), intent(IN) :: a, b, c ! factored LHS
    real(wp), dimension(len, nmax), intent(INOUT) :: f     ! RHS and solution

! -------------------------------------------------------------------
    integer(wi) :: n
    integer(wi) :: srt, end, siz
    real(wp) :: dummy1, dummy2

#ifdef USE_BLAS
    real(wp) alpha
    integer ilen
#else
    integer(wi) l
#endif

! ###################################################################
! -----------------------------------------------------------------------
! Forward sweep
! -----------------------------------------------------------------------
#ifdef USE_BLAS
    ilen = len
#endif

#ifdef USE_BLAS
!$omp parallel default(none) &
!$omp private(n,ilen,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,nmax,len)
#else
!$omp parallel default(none) &
!$omp private(n,l,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,nmax,len)
#endif

    call TLab_OMP_PARTITION(len, srt, end, siz)
    if (siz <= 0) then
        goto 999
    end if

#ifdef USE_BLAS
    ilen = siz
#endif

    do n = 2, nmax
        dummy1 = a(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n - 1), 1, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = f(l, n) + dummy1*f(l, n - 1)
        end do
#endif
    end do

! -----------------------------------------------------------------------
! Backward sweep
! -----------------------------------------------------------------------
    dummy1 = b(nmax)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, nmax), 1)
#else
    do l = srt, end
        f(l, nmax) = f(l, nmax)*dummy1
    end do
#endif

    do n = nmax - 1, 1, -1
        dummy1 = c(n)
        dummy2 = b(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n + 1), 1, f(srt, n), 1)
        call SCAL_LOC(ilen, dummy2, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = (f(l, n) + dummy1*f(l, n + 1))*dummy2
        end do

#endif
    end do
999 continue
!$omp end parallel

    return
end subroutine TRIDSS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine TRIDSS_ADD(nmax, len, a, b, c, f, g, h, d)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

    integer(wi), intent(IN) :: nmax  ! dimension of tridiagonal systems
    integer(wi), intent(IN) :: len   ! number of systems to be solved
    real(wp), dimension(nmax), intent(IN) :: a, b, c ! factored LHS
    real(wp), dimension(len, nmax), intent(INOUT) :: f     ! RHS and solution
    real(wp), dimension(len, nmax), intent(IN) :: g, h   ! Additive term
    real(wp), dimension(nmax), intent(IN) :: d     ! Scaled Jacobian

! -------------------------------------------------------------------
    integer(wi) :: n
    integer(wi) :: srt, end, siz
    real(wp) :: dummy1, dummy2

    integer(wi) l
#ifdef USE_BLAS
    real(wp) alpha
    integer ilen
#endif

! ###################################################################
! -----------------------------------------------------------------------
! Forward sweep
! -----------------------------------------------------------------------
#ifdef USE_BLAS
    ilen = len
#endif

#ifdef USE_BLAS
!$omp parallel default(none) &
!$omp private(n,l,ilen,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,d,g,h,nmax,len)
#else
!$omp parallel default(none) &
!$omp private(n,l,srt,end,siz,dummy1,dummy2) &
!$omp shared(f,a,b,c,d,g,h,nmax,len)
#endif

    call TLab_OMP_PARTITION(len, srt, end, siz)
    if (siz <= 0) then
        goto 999
    end if

#ifdef USE_BLAS
    ilen = siz
#endif

    do n = 2, nmax
        dummy1 = a(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n - 1), 1, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = f(l, n) + dummy1*f(l, n - 1)
        end do
#endif
    end do

! -----------------------------------------------------------------------
! Backward sweep
! -----------------------------------------------------------------------
    dummy1 = b(nmax)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, nmax), 1)
#else
    do l = srt, end
        f(l, nmax) = f(l, nmax)*dummy1
    end do
#endif

    do n = nmax - 1, 1, -1
        dummy1 = c(n)
        dummy2 = b(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n + 1), 1, f(srt, n), 1)
        call SCAL_LOC(ilen, dummy2, f(srt, n), 1)
        do l = srt, end
            f(l, n + 1) = f(l, n + 1) - (d(n + 1) + g(l, n + 1))*h(l, n + 1)  ! TO BE IMPLEMENTED USING BLAS
        end do
#else
        do l = srt, end
            f(l, n) = (f(l, n) + dummy1*f(l, n + 1))*dummy2
            f(l, n + 1) = f(l, n + 1) - (d(n + 1) + g(l, n + 1))*h(l, n + 1)
        end do
#endif
    end do
    do l = srt, end
        f(l, 1) = f(l, 1) - (d(1) + g(l, 1))*h(l, 1)
    end do

999 continue
!$omp end parallel

    return
end subroutine TRIDSS_ADD

!########################################################################
!# DESCRIPTION
!#
!# Circulant system (periodic bcs)
!#
!########################################################################

! #######################################################################
! LU factorization stage
! #######################################################################
subroutine TRIDPFS(nmax, a, b, c, d, e)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi) nmax
    real(wp), dimension(nmax) :: a, b, c, d, e

! -------------------------------------------------------------------
    integer(wi) n
    real(wp) sum

! ###################################################################
! Generate first elements of LU
    c(1) = c(1)/b(1)
    e(1) = a(1)/b(1)
    d(1) = c(nmax)

! Generate n=2 to n=n-2 elements of LU
    do n = 2, nmax - 2
        b(n) = b(n) - a(n)*c(n - 1)
        c(n) = c(n)/b(n)
        e(n) = -a(n)*e(n - 1)/b(n)
        d(n) = -d(n - 1)*c(n - 1)
    end do

! Generate n-1 elements
    b(nmax - 1) = b(nmax - 1) - a(nmax - 1)*c(nmax - 2)
    e(nmax - 1) = (c(nmax - 1) - a(nmax - 1)*e(nmax - 2))/b(nmax - 1)
    d(nmax - 1) = a(nmax) - d(nmax - 2)*c(nmax - 2)

! Generate the n-th element
    sum = 0.0_wp
    do n = 1, nmax - 1
        sum = sum + d(n)*e(n)
    end do
    b(nmax) = b(nmax) - sum

! Final operations
    do n = 1, nmax
        b(n) = 1.0_wp/b(n)
        a(n) = -a(n)*b(n)
        c(n) = -c(n)
        e(n) = -e(n)
    end do

    return
end subroutine TRIDPFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine TRIDPSS(nmax, len, a, b, c, d, e, f, wrk)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

    integer(wi), intent(IN) :: nmax, len
    real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
    real(wp), dimension(len, nmax), intent(INOUT) :: f
    real(wp), dimension(len) :: wrk

! -------------------------------------------------------------------
    real(wp) :: dummy1, dummy2
    integer(wi) :: srt, end, siz

    integer(wi) l, n
#ifdef USE_BLAS
    integer :: ilen
#endif

! -------------------------------------------------------------------
! Forward sweep
! -------------------------------------------------------------------

#ifdef USE_BLAS
!$omp parallel default( none ) &
!$omp private(n, l,ilen, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,wrk,nmax,a,b,c,d,e,len)
#else
!$omp parallel default( none ) &
!$omp private(n, l, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,wrk,nmax,a,b,c,d,e,len)
#endif

    call TLab_OMP_PARTITION(len, srt, end, siz)
    if (siz <= 0) then
        goto 999
    end if

#ifdef USE_BLAS
    ilen = siz
#endif

    dummy1 = b(1)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, 1), 1)
#else
    do l = srt, end
        f(l, 1) = f(l, 1)*dummy1
    end do
#endif

    do n = 2, nmax - 1
        dummy1 = a(n)
        dummy2 = b(n)
#ifdef USE_BLAS
        call SCAL_LOC(ilen, dummy2, f(srt, n), 1)
        call AXPY_LOC(ilen, dummy1, f(srt, n - 1), 1, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = f(l, n)*dummy2 + dummy1*f(l, n - 1)
        end do
#endif
    end do

    wrk(srt:end) = 0.0_wp

    do n = 1, nmax - 1
        dummy1 = d(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n), 1, wrk(srt), 1)
#else
        do l = srt, end
            wrk(l) = wrk(l) + dummy1*f(l, n)
        end do
#endif
    end do

    dummy1 = b(nmax)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, nmax), 1)
    call AXPY_LOC(ilen, -dummy1, wrk(srt), 1, f(srt, nmax), 1)
#else
    do l = srt, end
        f(l, nmax) = (f(l, nmax) - wrk(l))*dummy1
    end do
#endif

! -------------------------------------------------------------------
! Backward sweep
! -------------------------------------------------------------------
    dummy1 = e(nmax - 1)

#ifdef USE_BLAS
    call AXPY_LOC(ilen, dummy1, f(srt, nmax), 1, f(srt, nmax - 1), 1)
#else
    do l = srt, end
        f(l, nmax - 1) = dummy1*f(l, nmax) + f(l, nmax - 1)
    end do
#endif

    do n = nmax - 2, 1, -1
        dummy1 = c(n)
        dummy2 = e(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy2, f(srt, nmax), 1, f(srt, n), 1)
        call AXPY_LOC(ilen, dummy1, f(srt, n + 1), 1, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = f(l, n) + dummy1*f(l, n + 1) + dummy2*f(l, nmax)
        end do
#endif
    end do
999 continue
!$omp end parallel

    return
end subroutine TRIDPSS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine TRIDPSS_ADD(nmax, len, a, b, c, d, e, f, g, h, wrk)
    use TLab_Constants, only: wp, wi
    use TLab_OpenMP

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

    integer(wi), intent(IN) :: nmax, len
    real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
    real(wp), dimension(len, nmax), intent(INOUT) :: f
    real(wp), dimension(len, nmax), intent(IN) :: g, h
    real(wp), dimension(len) :: wrk

! -------------------------------------------------------------------
    real(wp) :: dummy1, dummy2
    integer(wi) :: srt, end, siz

    integer(wi) l, n
#ifdef USE_BLAS
    integer :: ilen
#endif

! -------------------------------------------------------------------
! Forward sweep
! -------------------------------------------------------------------

#ifdef USE_BLAS
!$omp parallel default( none ) &
!$omp private(n, l,ilen, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,g,h,wrk,nmax,a,b,c,d,e,len)
#else
!$omp parallel default( none ) &
!$omp private(n, l, dummy1, dummy2, srt, end,siz) &
!$omp shared(f,g,h,wrk,nmax,a,b,c,d,e,len)
#endif

    call TLab_OMP_PARTITION(len, srt, end, siz)
    if (siz <= 0) then
        goto 999
    end if

#ifdef USE_BLAS
    ilen = siz
#endif

    dummy1 = b(1)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, 1), 1)
#else
    do l = srt, end
        f(l, 1) = f(l, 1)*dummy1
    end do
#endif

    do n = 2, nmax - 1
        dummy1 = a(n)
        dummy2 = b(n)
#ifdef USE_BLAS
        call SCAL_LOC(ilen, dummy2, f(srt, n), 1)
        call AXPY_LOC(ilen, dummy1, f(srt, n - 1), 1, f(srt, n), 1)
#else
        do l = srt, end
            f(l, n) = f(l, n)*dummy2 + dummy1*f(l, n - 1)
        end do
#endif
    end do

    wrk(srt:end) = 0.0_wp

    do n = 1, nmax - 1
        dummy1 = d(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy1, f(srt, n), 1, wrk(srt), 1)
#else
        do l = srt, end
            wrk(l) = wrk(l) + dummy1*f(l, n)
        end do
#endif
    end do

    dummy1 = b(nmax)
#ifdef USE_BLAS
    call SCAL_LOC(ilen, dummy1, f(srt, nmax), 1)
    call AXPY_LOC(ilen, -dummy1, wrk(srt), 1, f(srt, nmax), 1)
#else
    do l = srt, end
        f(l, nmax) = (f(l, nmax) - wrk(l))*dummy1
    end do
#endif

! -------------------------------------------------------------------
! Backward sweep
! -------------------------------------------------------------------
    dummy1 = e(nmax - 1)

#ifdef USE_BLAS
    call AXPY_LOC(ilen, dummy1, f(srt, nmax), 1, f(srt, nmax - 1), 1)
#else
    do l = srt, end
        f(l, nmax - 1) = dummy1*f(l, nmax) + f(l, nmax - 1)
    end do
#endif

    do n = nmax - 2, 1, -1
        dummy1 = c(n)
        dummy2 = e(n)
#ifdef USE_BLAS
        call AXPY_LOC(ilen, dummy2, f(srt, nmax), 1, f(srt, n), 1)
        call AXPY_LOC(ilen, dummy1, f(srt, n + 1), 1, f(srt, n), 1)
        do l = srt, end
            f(l, n + 1) = f(l, n + 1) - g(l, n + 1)*h(l, n + 1)  ! TO BE IMPLEMENTED USING BLAS
        end do
#else
        do l = srt, end
            f(l, n) = f(l, n) + dummy1*f(l, n + 1) + dummy2*f(l, nmax)
            f(l, n + 1) = f(l, n + 1) - g(l, n + 1)*h(l, n + 1)
        end do
#endif
    end do
    do l = srt, end
        f(l, 1) = f(l, 1) - g(l, 1)*h(l, 1)
        f(l, nmax) = f(l, nmax) - g(l, nmax)*h(l, nmax)
    end do
999 continue
!$omp end parallel

    return
end subroutine TRIDPSS_ADD

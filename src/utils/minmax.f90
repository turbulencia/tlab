!########################################################################
!#
!# Determine Maximum and Minimum Values of a field
!#
!########################################################################
subroutine MINMAX(imax, jmax, kmax, a, amn, amx)
    use TLab_Constants, only: wp, wi
#ifdef USE_MPI
    use mpi_f08
#endif

    implicit none

    integer(wi) imax, jmax, kmax
    real(wp) a(imax*jmax*kmax), amn, amx

! -----------------------------------------------------------------------
#ifdef USE_MPI
    real(wp) pamn, pamx
    integer ims_err
#endif

! #######################################################################
    amn = minval(a)
    amx = maxval(a)

#ifdef USE_MPI
    pamn = amn
    pamx = amx
    call MPI_ALLREDUCE(pamn, amn, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(pamx, amx, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
#endif

    return
end subroutine MINMAX

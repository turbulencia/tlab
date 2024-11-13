module TLab_OpenMP
    use TLab_Constants, only: wi
    implicit none
    private

    integer, public :: TLab_OMP_numThreads
    integer, public :: TLab_OMP_error

    public :: TLab_OMP_PARTITION

contains
    subroutine TLab_OMP_PARTITION(len, omp_srt, omp_end, omp_siz)
#ifdef USE_OPENMP
        use OMP_LIB
#endif

        integer(wi), intent(IN) :: len
        integer(wi), intent(INOUT) :: omp_srt, omp_end, omp_siz

! -----------------------------------------------------------------------
#ifdef USE_OPENMP
        integer(wi) :: residual
        integer(wi) :: TLab_OMP_threadID
#endif

! #######################################################################
#ifdef USE_OPENMP
        TLab_OMP_threadID = omp_get_thread_num()

        if (len > TLab_OMP_numThreads) then
            residual = mod(len, TLab_OMP_numThreads)
            if (TLab_OMP_threadID >= residual) then
                omp_siz = len/TLab_OMP_numThreads
                omp_srt = TLab_OMP_threadID*(len/TLab_OMP_numThreads) + residual + 1
                omp_end = omp_srt + omp_siz - 1
            else
                omp_siz = len/TLab_OMP_numThreads + 1
                omp_srt = TLab_OMP_threadID*(len/TLab_OMP_numThreads + 1) + 1
                omp_end = omp_srt + omp_siz - 1
            end if
        else
            if (TLab_OMP_threadID == 0) then
                omp_siz = len
                omp_srt = 1
                omp_end = len
            else
                omp_srt = 1
                omp_end = -1
                omp_siz = -1
            end if
        end if

#else
        omp_siz = len
        omp_srt = 1
        omp_end = len

#endif

        return

    end subroutine TLab_OMP_PARTITION

end module TLab_OpenMP

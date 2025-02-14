module TLabMPI_VARS
    use TLab_Constants, only: wp, wi
    use mpi_f08
    implicit none

    type(MPI_Comm) :: ims_comm_xz                               ! Plane communicators
    type(MPI_Comm) :: ims_comm_x, ims_comm_z                    ! Directional communicators
    type(MPI_Comm) :: ims_comm_xz_aux
    type(MPI_Comm) :: ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_pro                                          ! Local PE number in global communicator
    integer :: ims_pro_i, ims_pro_j, ims_pro_k                  ! Local PE number in directional communicators
    integer :: ims_npro                                         ! Total number of PE
    integer :: ims_npro_i, ims_npro_j, ims_npro_k               ! Number of PEs in directional communicators
    integer :: ims_offset_i, ims_offset_j, ims_offset_k         ! Grid offsets

    integer :: ims_err, ims_tag

    real(wp) :: ims_time_min, ims_time_max, ims_time_trans      ! Profiling
    integer(wi) :: ims_bcs_imax, ims_bcs_jmax

    type(MPI_Datatype) :: TLAB_MPI_REAL_TYPE                    ! MPI Type control

#ifdef USE_PSFFT
    integer :: ims_nb_thrsupp_provided
    integer, dimension(2) :: ims_nb_proc_grid
    integer, dimension(3) :: ims_nb_msize
    integer, dimension(3) :: ims_nb_xsrt, ims_nb_xend, ims_nb_xsiz
    integer, dimension(3) :: ims_nb_ysrt, ims_nb_yend, ims_nb_ysiz
    integer, dimension(3) :: ims_nb_zsrt, ims_nb_zend, ims_nb_zsiz
#endif

end module TLabMPI_VARS
